# DMRseq Pipeline Refactor: Array Jobs + Threading Fix + Block Mode

## Context

The current `06_mch_dmrseq.R` processes all 21 chromosomes sequentially in a single SLURM job. Every chromosome crashes with the same BiocParallel error:

```
dmrseq error: error in evaluating the argument 'args' in selecting a method
for function 'do.call': wrong args for environment subassignment
```

This is not a timeout — it is a deterministic crash caused by the interaction between `MulticoreParam(workers = 40)` (which forks the R process) and `setDTthreads(24)` (which activates OpenMP threads inside each fork). Forking a process with active OpenMP state corrupts the child processes' environment handling.

The fix has three parts:
1. Eliminate the fork/OpenMP conflict by isolating data.table threading from MulticoreParam
2. Parallelize chromosomes as a SLURM array (one chromosome per task)
3. Add `block = TRUE` mode to detect megabase-scale methylation blocks alongside local DMRs

## Architecture

The monolithic script becomes three phases, each in its own R script and SLURM wrapper:

```
Phase 1 (split) ─── single job
    │
    ├── Phase 2a (local DMRs) ─── array 0-20, one chromosome each
    │       └── Phase 3a (gather local) ─── single job
    │
    └── Phase 2b (block DMRs) ─── array 0-20, one chromosome each
            └── Phase 3b (gather blocks) ─── single job
```

- **Phase 1** reads all 8 samples (from `combined_ch/`), extracts lambda rates, filters to canonical chromosomes, applies lambda correction, splits by chromosome, and saves per-chromosome RDS files. Identical to the current Phase 1. Saves lambda rates as a metadata file for the gather step.
- **Phase 2** loads one chromosome's RDS files (all 8 samples), merges to shared sites, builds BSseq, runs dmrseq. Two parallel arrays: one for local DMRs (`block = FALSE`), one for large-scale blocks (`block = TRUE`). Saves per-chromosome results as RDS.
- **Phase 3** loads all per-chromosome RDS results for one mode, combines them, writes the final TSV and summary. Two independent gather jobs, one per mode — local DMR results are available without waiting for blocks to finish.

## The threading fix

The root cause: Phase 1 sets `setDTthreads(24)`, then Phase 2 creates `MulticoreParam(workers = 40)`. MulticoreParam forks the R process, and each forked child inherits the parent's 24 OpenMP threads. The fork + OpenMP interaction corrupts environment state.

In the refactored version, Phase 1 and Phase 2 are separate R processes (separate SLURM jobs), so there is no inherited OpenMP state. The per-chromosome R script (`02_dmrseq_chr.R`) sets `setDTthreads(1)` at startup before creating MulticoreParam. data.table is only used for the merge step, which is fast at single-threaded — the bottleneck is dmrseq's smoothing and permutation, which MulticoreParam handles.

## Directory structure

```
scripts/dmrseq/
    01_split_samples.R          # Phase 1: read 8 samples, split by chromosome
    01_split_samples.sb         # SLURM wrapper (single job)
    02_dmrseq_chr.R             # Phase 2: per-chromosome dmrseq (takes --chr, --block)
    02_dmrseq_local.sb          # SLURM wrapper (array 0-20, local DMRs)
    02_dmrseq_blocks.sb         # SLURM wrapper (array 0-20, block = TRUE)
    03_gather_results.R         # Phase 3: combine per-chr results (takes --mode)
    03_gather_local.sb          # SLURM wrapper (single job, local results)
    03_gather_blocks.sb         # SLURM wrapper (single job, block results)
    submit_dmrseq.sh            # Dependency chaining across all phases

results/dmrseq/                 # Created by scripts at runtime
    .per_chr/local/             # Per-chromosome RDS results (local DMRs)
    .per_chr/blocks/            # Per-chromosome RDS results (blocks)
    mch_dmrs_local.tsv          # Combined local DMR table
    mch_dmrs_blocks.tsv         # Combined block table
    mch_dmr_summary.txt         # Summary covering both modes
    lambda_rates.tsv            # Sample lambda rates from Phase 1

logs/dmrseq/                    # SLURM output (add to create_log_dirs.sh)
```

Intermediate per-chromosome split RDS files (Phase 1 output, Phase 2 input) go in `DATA_DIR/.dmrseq_splits/` — these are large (8 samples x 21 chromosomes = 168 files, several GB total) and temporary. Phase 3 cleans them up after gathering results.

## R script details

### `01_split_samples.R`

Extracted from the current `06_mch_dmrseq.R` Phase 1. Reuses `split_sample_to_disk()` and the lambda extraction logic. CLI arguments:

- `--ch-dir` — path to combined_ch/ (DATA_DIR)
- `--out-dir` — path to results/dmrseq/ (BASEDIR/results/dmrseq)
- `--split-dir` — path for chromosome splits (DATA_DIR/.dmrseq_splits/)
- `--threads` — data.table threads (matches SLURM CPUs)

Outputs:
- `{split_dir}/{sample}_{chr}.rds` for each sample x chromosome (168 files)
- `{out_dir}/lambda_rates.tsv` with per-sample lambda rates

### `02_dmrseq_chr.R`

Extracted from the current `06_mch_dmrseq.R` Phase 2. Reuses `process_chromosome()` with modifications. CLI arguments:

- `--chr` — chromosome name (e.g., `chr1`)
- `--split-dir` — where Phase 1 saved the RDS splits
- `--out-dir` — where to save per-chromosome results
- `--mode` — `local` or `blocks` (controls `block` parameter to dmrseq)
- `--cutoff` — methylation difference cutoff (default 0.005)
- `--min-num-region` — minimum CH sites per candidate region (default 5)
- `--workers` — BiocParallel workers

Key changes from the current script:
- `setDTthreads(1)` at script startup — prevents the fork/OpenMP conflict
- Accepts `--chr` for a single chromosome instead of looping over all
- Accepts `--mode blocks` which passes `block = TRUE` to `dmrseq()`
- Saves results as RDS: `{out_dir}/{mode}/{chr}_dmrs.rds`
- Saves NULL (as an empty-result marker) when no DMRs are found, so the gather step knows the chromosome was processed vs. failed

### `03_gather_results.R`

New script. CLI arguments:

- `--results-dir` — path to results/dmrseq/
- `--mode` — `local` or `blocks`
- `--split-dir` — path to chromosome splits (for cleanup)
- `--clean` — whether to delete split files after gathering (default TRUE)

Reads all `{results_dir}/.per_chr/{mode}/chr*_dmrs.rds` files. Filters out NULLs (chromosomes with no DMRs). Combines with `do.call(c, ...)`. Writes:
- `{results_dir}/mch_dmrs_{mode}.tsv`
- `{results_dir}/mch_dmr_summary.txt` (appends if already exists, so both modes contribute)

Reads `{results_dir}/lambda_rates.tsv` for the summary. Reports per-chromosome DMR counts, width distribution, direction (hyper/hypo in mutant), q-value thresholds. If no DMRs found on any chromosome, writes the "no DMRs" summary (same language as the current script about CH spatial correlation).

Cleanup: if `--clean` is TRUE and this is the second mode to finish (both local/ and blocks/ have been gathered), deletes `{split_dir}/` (the Phase 1 splits) and `{results_dir}/.per_chr/`.

## SLURM wrappers

All wrappers follow the existing boilerplate pattern (header echo, TMPDIR setup, conda activation, input validation, timed execution, footer echo).

### `01_split_samples.sb`
- Single job (no array)
- 16 CPUs, 64G, 48h
- Log: `logs/dmrseq/split_%j.out`
- Calls `Rscript 01_split_samples.R --threads 16`

### `02_dmrseq_local.sb`
- Array 0-20 (21 chromosomes)
- 32 CPUs, 120G, 48h
- Log: `logs/dmrseq/local_%A_%a.out`
- Maps array index to chromosome via bash array:
  `CHRS=(chr1 chr2 ... chr19 chrX chrY)`
- Calls `Rscript 02_dmrseq_chr.R --chr $CHR --mode local --workers 28`

### `02_dmrseq_blocks.sb`
- Identical to `02_dmrseq_local.sb` except `--mode blocks`
- Log: `logs/dmrseq/blocks_%A_%a.out`

### `03_gather_local.sb` / `03_gather_blocks.sb`
- Single job (no array)
- 12 CPUs, 32G, 48h
- Log: `logs/dmrseq/gather_local_%j.out` / `logs/dmrseq/gather_blocks_%j.out`
- Calls `Rscript 03_gather_results.R --mode local` / `--mode blocks`

### `submit_dmrseq.sh`
- Validates that all 8 `combined_ch/*_CH.methylKit.gz` files exist
- Dependency chain:
  ```
  SPLIT_JID=$(sbatch --parsable 01_split_samples.sb)
  LOCAL_JID=$(sbatch --parsable --dependency=afterok:$SPLIT_JID 02_dmrseq_local.sb)
  BLOCK_JID=$(sbatch --parsable --dependency=afterok:$SPLIT_JID 02_dmrseq_blocks.sb)
  GATHER_L=$(sbatch --parsable --dependency=afterok:$LOCAL_JID 03_gather_local.sb)
  GATHER_B=$(sbatch --parsable --dependency=afterok:$BLOCK_JID 03_gather_blocks.sb)
  ```
- Prints dependency graph and job IDs

## Files to modify

- **`scripts/utils/create_log_dirs.sh`** — add `logs/dmrseq` to the `dirs` array
- **`scripts/submit_mch_pipeline.sh`** — replace the single `sbatch 06_mch_dmrseq.sb` with `bash "${SCRIPT_DIR}/dmrseq/submit_dmrseq.sh"`, or remove dmrseq from the main pipeline and note it as a standalone submission (since it has no dependencies on other pipeline steps)

## Files to archive

- `scripts/06_mch_dmrseq.R` → `archive/06_mch_dmrseq.R`
- `scripts/06_mch_dmrseq.sb` → `archive/06_mch_dmrseq.sb`

## Resource allocation rationale

| Job | CPUs | Memory | Time | Why |
|-----|------|--------|------|-----|
| Split | 16 | 64G | 48h | data.table reading; sequential samples; memory for one ~1B-row table at a time |
| Per-chr (local) | 32 | 120G | 48h | MulticoreParam(28 workers); BSseq for chr1 ~10GB + dmrseq working memory ~30-50GB |
| Per-chr (blocks) | 32 | 120G | 48h | Same as local; block mode changes smoothing bandwidth, not memory profile |
| Gather | 12 | 32G | 48h | Just loading small per-chr RDS files and combining |

Per the resource preference: all jobs get max walltime (48h), minimum 12 CPU / 32G, generous memory for data-heavy steps.

## Verification

**Smoke test first.** Before submitting the full 42-task array, test with a single chromosome:

1. Submit Phase 1 (split) normally — this has to run fully to produce the RDS files
2. Submit Phase 2 for one chromosome only: `sbatch --array=0 02_dmrseq_local.sb` (chr1, the largest chromosome at 72M sites — if it works here, it works everywhere)
3. Check the log: confirm `setDTthreads(1)` is set, MulticoreParam created successfully, no "wrong args for environment subassignment" error
4. If chr1 completes: submit the full arrays (local + blocks) and gather jobs
5. If chr1 still errors: the threading fix didn't resolve it; investigate further before burning 42 job slots

**Full run verification:**

1. Phase 1 log: all 8 samples read, 168 RDS files created in `.dmrseq_splits/`, `lambda_rates.tsv` written
2. Phase 2 array logs: each chromosome either finds DMRs or completes with "No DMRs" (both are valid outcomes for CH data). No errors.
3. Phase 3 logs: combined TSV(s) and summary written to `results/dmrseq/`
4. Inspect `results/dmrseq/mch_dmr_summary.txt` for per-chromosome counts and parameters
5. Verify split files cleaned up after both gather jobs complete
6. Finding zero DMRs across all chromosomes is a valid and expected outcome for CH data — the summary explains why (CH lacks the spatial correlation of CpG islands that dmrseq was designed for). The block mode results are the more scientifically interesting output for a genome-wide shift like BAP1-KO.
