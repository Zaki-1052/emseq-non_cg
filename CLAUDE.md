# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

See [README.md](README.md) for project overview, pipeline architecture, samples, directory layout, and dependencies.

## Interpretation guidance

The BAP1-KO genome-wide mCH shift is expected biology, not a confounding artifact. High lambda_GC values in QQ plots are consistent with this. Do not treat lambda_GC > 1 as evidence of a statistical problem.

The upstream EM-seq Nextflow pipeline and MethylDackel extraction are complete and not rerun from here.

## Compute environment

All pipeline scripts run on SDSC Expanse (SLURM) via conda environment `emseq_env`. Account `csd940`. Base directory on Expanse: `/expanse/lustre/projects/csd940/zalibhai/emseq`.

This local working directory mirrors the remote scripts and results. Scripts are developed locally and transferred to Expanse.

@README.md
@memories/MEMORY.md
@docs/abstract.md
@docs/aims.md
@docs/primary_handoff.md