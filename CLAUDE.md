# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

See [README.md](README.md) for project overview, pipeline architecture, samples, directory layout, and dependencies.

## Interpretation guidance

The BAP1-KO genome-wide mCH shift is expected biology, not a confounding artifact. High lambda_GC values in QQ plots are consistent with this. Do not treat lambda_GC > 1 as evidence of a statistical problem.

The upstream EM-seq Nextflow pipeline and MethylDackel extraction are complete and not rerun from here.

## Compute environment

All pipeline scripts run on SDSC Expanse (SLURM) via conda environment `emseq_env`. Account `csd940`.

On Expanse, the code and data live in separate directories:
- Code: `/expanse/lustre/projects/csd940/zalibhai/emseq-repo`
- Data: `/expanse/lustre/projects/csd940/zalibhai/emseq`
- Hi-C companion: `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c`
- BigWigs: `/expanse/lustre/projects/csd940/zalibhai/bigwigs`

`paths.yaml` at the repo root maps these roots for `local` and `expanse` environments. The section pipeline's shared config (`scripts/sections/_shared_config.R`) reads it via `EMSEQ_ENV`.

This local working directory mirrors the remote scripts and results. Scripts are developed locally and transferred to Expanse.

## Section pipeline

23 downstream analysis sections in `scripts/sections/`, organized in 7 thematic groups. Ported from the companion Biomodal CG 5mC/5hmC pipeline for single-modality mCH. Each section sources `_shared_config.R`, which pre-loads the mCH differential results, DiffBind tables, consensus peaks, and shared helpers.

Chromatin state is recorded as two columns per gene: `promoter_state` (over TSS +/- 2kb) and `body_state` (over the gene body). A single `chromatin_state` column does not exist.

Four sections produce handoff files that later sections consume:
- `10_01` writes `gene_chromatin_state.tsv` (19 columns)
- `20_02` writes `gene_level_all_marks.tsv` (27 columns)
- `60_01` writes `mecp2_no_meth_genes.tsv`
- `70_01` writes `neuronal_gene_set.tsv`

Gene-level Fisher tests go through `register_fisher_test()`, which writes per-test shards to `results/sections/fisher_registry/`. Section `40_04` reads every shard for permutation validation.

Setup guide: `docs/sections_setup.md`. Deviation log: `docs/biomodal-deviations.md`.

@README.md
@memories/MEMORY.md
@docs/abstract.md
@docs/aims.md
@docs/primary_handoff.md