# SLURM wrappers for the mCH sections

This directory holds one `.sb` wrapper per section and the submission script
that chains them.

| File | Purpose |
|------|---------|
| `<section>.sb` | One SLURM job. Activates `emseq_env`, sets `EMSEQ_ENV=expanse`, runs the section Rscript. |
| `submit_sections.sh` | Submits all 23 sections in dependency order with `--dependency=afterok`. |

Every wrapper writes its job log to `logs/sections/<section>_<jobid>.out`, a
path relative to the code directory
`/expanse/lustre/projects/csd940/zalibhai/emseq-repo`. Submit from that
directory. `submit_sections.sh` changes into it before it calls `sbatch`.

## Before the first submission

Complete steps 1 to 7 of [`docs/sections_setup.md`](../../../docs/sections_setup.md):
R packages installed, reference data copied into `data/`, feature BEDs built,
and `_shared_config.R` verified. Every section sources `_shared_config.R` at
startup, so one missing input file stops all 23 jobs the same way.

Section 50_01 also needs the step 02b output
(`results/02b_features/{sample}_feature_mch.tsv`). Step 02b is submitted
separately and is not part of this dependency graph.

## Submit the whole pipeline

```
cd /expanse/lustre/projects/csd940/zalibhai/emseq-repo
mkdir -p logs/sections
bash scripts/sections/slurm/submit_sections.sh 2>&1 | tee logs/sections/submit_sections.txt
```

The script prints one line per submission and a summary table of job id,
section, and what each job waits on. Keep the log: the table is the only record
of which job id belongs to which section.

Check the plan without submitting:

```
bash scripts/sections/slurm/submit_sections.sh --dry-run
```

`--dry-run` prints every `sbatch` command in submission order and substitutes
placeholder job ids (9000001 and up) in the dependency arguments. It runs
anywhere, including on a laptop, if you point it at a checkout:

```
EMSEQ_CODE_DIR=/Users/zakiralibhai/emseq bash scripts/sections/slurm/submit_sections.sh --dry-run
```

## Submit one section

Two equivalent ways. Direct `sbatch`:

```
cd /expanse/lustre/projects/csd940/zalibhai/emseq-repo
sbatch scripts/sections/slurm/10_01_chromatin_state.sb
```

Through the submission script, which resolves a unique name prefix and prints
the same summary table:

```
bash scripts/sections/slurm/submit_sections.sh 10_01_chromatin_state
bash scripts/sections/slurm/submit_sections.sh 10_01
```

Both forms skip every dependency. The named section is submitted alone, with no
`--dependency` argument. A section that reads another section's output stops
with a message naming the missing file and the section that writes it, for
example:

```
Gene-level mark table not found: .../gene_level_all_marks.tsv
Run section 20_02 (20_02_multi_mark_diffbind.R) first.
```

## Raise the permutation count

The four `40_*` sections take a permutation count. Defaults:

| Section | Default | Minimum |
|---------|---------|---------|
| `40_01_dmr_marks` | 5000 | none |
| `40_02_atac_loops` | 5000 | none |
| `40_03_domains` | 5000 | none |
| `40_04_gene_level` | 10000 | 100 |

The wrappers read `NTIMES` from the job environment and pass it to the R
script. `submit_sections.sh --ntimes N` exports it to the `40_*` jobs only:

```
bash scripts/sections/slurm/submit_sections.sh --ntimes 20000 --force-rerun 40_01_dmr_marks
bash scripts/sections/slurm/submit_sections.sh --ntimes 20000 --force-rerun
```

The same thing with plain `sbatch`:

```
sbatch --export=ALL,NTIMES=20000,FORCE_RERUN=1 scripts/sections/slurm/40_01_dmr_marks.sb
```

`--ntimes` and `--force-rerun` fail with an error when no `40_*` section is
being submitted, rather than being silently ignored.

Pair `--ntimes` with `--force-rerun` for 40_01, 40_03, and 40_04. Those three
cache their permutations under a fixed filename, so a higher count on its own
loads the old cache and reports the old count. 40_02 writes the count into the
cache filename (`40_02*_crosswise_n<ntimes>.rds`), so a new count builds a new
cache by itself.

Permutation runtime scales with the count. 20000 shuffles cost about four times
5000. Walltime is already at the 48 hour maximum.

## Force a permutation rerun

`FORCE_RERUN` set to any non-empty value makes the `40_*` sections ignore their
RDS caches and recompute. Only those four sections read it.

```
bash scripts/sections/slurm/submit_sections.sh --force-rerun 40_03_domains
sbatch --export=ALL,FORCE_RERUN=1 scripts/sections/slurm/40_03_domains.sb
```

Caches live in `results/sections/40_permutation/`:

| Section | Cache files |
|---------|-------------|
| `40_01` | `40_01_permutation_objects.rds` |
| `40_02` | `40_02*_crosswise_n<ntimes>.rds`, `40_02c_local_zscore_n<ntimes>.rds` |
| `40_03` | `40_03_perm_compartment.rds`, `40_03_perm_polycomb.rds`, `40_03_localz_compartment.rds`, `40_03_localz_polycomb.rds` |
| `40_04` | `40_04_gene_level_permutation.rds` |

Figures and tables are rebuilt on every run, cached or not.

40_04 compares its cache against `results/sections/fisher_test_registry.tsv`
and stops when the two cover different tests. So after any upstream section
adds, removes, or renames a Fisher test, rerun 40_04 with `--force-rerun`.

## Dependency graph

`submit_sections.sh` submits in four waves. Each wave passes
`--dependency=afterok:<job ids>` plus `--kill-on-invalid-dep=yes`, so a
dependent job is cancelled when its producer fails instead of waiting forever.

```
wave 1  (no dependency, 12 jobs)
    10_01_chromatin_state           10_02_ab_compartment
    10_04_subcompartment            20_01_mecp2_correlation
    30_01_loop_anchor_methylation   30_02_mecp2_loop_anchors
    40_01_dmr_marks                 50_01_feature_methylation
    60_01_methylation_scale         60_03_reconciliation
    60_04_aging                     70_01_k119ub_neuronal

wave 2  (6 jobs)
    10_01 -> 10_03_polycomb_enrichment    (gene_chromatin_state.tsv)
    10_01 -> 20_02_multi_mark_diffbind    (gene_chromatin_state.tsv)
    10_01 -> 40_02_atac_loops             (gene_chromatin_state.tsv)
    10_01 -> 40_03_domains                (gene_chromatin_state.tsv)
    60_01 -> 60_02_k119ub_unmethylated    (mecp2_no_meth_genes.tsv)
    70_01 -> 70_03_geneset_overlap        (neuronal_gene_set.tsv)

wave 3  (4 jobs)
    20_02        -> 20_03_quadrant_scatters     (gene_level_all_marks.tsv)
    20_02        -> 20_04_mch_mecp2_by_mark     (gene_level_all_marks.tsv)
    70_01, 20_02 -> 70_02_chromatin_remodeling  (both files)
    70_01, 20_02 -> 70_04_synapse_chromatin     (both files)

wave 4  (1 job)
    10_01, 10_02, 10_03, 20_01, 20_02, 30_01, 30_02,
    60_01, 60_02, 60_04, 70_01, 70_03 -> 40_04_gene_level
                                          (fisher_test_registry.tsv)
```

Sections that register a Fisher test append a row to
`results/sections/fisher_test_registry.tsv`. Section 40_04 reads that registry
and revalidates every recorded test by label shuffle, so it runs last. The list
of registering sections is the `FISHER_SECTIONS` array at the top of
`submit_sections.sh`.

## Resources

| Tier | Sections | CPUs | Memory | Walltime |
|------|----------|------|--------|----------|
| Standard | the other 19 | 16 | 64 GB | 48:00:00 |
| Permutation | 40_01, 40_02, 40_03, 40_04 | 32 | 128 GB | 48:00:00 |

Account `csd940`, partition `shared`. Each wrapper sets `TMPDIR` to
`/expanse/lustre/projects/csd940/zalibhai/emseq/.scratch/<jobid>`, on the same
lustre filesystem as the data.

## Watching the run

```
squeue -u $USER
squeue -u $USER --format="%.10i %.24j %.8T %.10M %R"
```

`State` of `Dependency` means the job is waiting for its producer.
`DependencyNeverSatisfied` cannot occur here: `--kill-on-invalid-dep=yes`
cancels those jobs instead.

Exit codes and elapsed time after the fact:

```
sacct -u $USER --starttime today --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS
```

Each job log ends with an `Exit code:` line and the elapsed seconds.

## When a job fails

1. Read `logs/sections/<section>_<jobid>.out`. Every section stops at the first
   missing or malformed input and names it.
2. Fix the cause, then resubmit that section alone
   (`submit_sections.sh <section>`).
3. Resubmit its consumers, in wave order, after it finishes. Cancelled
   dependents are not resubmitted automatically.
