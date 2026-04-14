# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics/systems biology research project building a **mechanistic immune digital twin of the salivary gland** in primary Sjögren's disease (pSS/SjD). It adapts the Boolean modelling framework from Zerrouk et al. (*npj Digital Medicine*, 2024) — originally developed for Rheumatoid Arthritis synovial tissue — to model SGECs and infiltrating immune cells in SjD.

## Environment Setup

```bash
# Create and activate the conda environment (file not yet committed)
conda env create -f envs/environment.yml
conda activate sjs-digitaltwin

# Python venv (local, gitignored)
source venv/bin/activate
```

The Python venv at `venv/` is gitignored. Primary dependency management will be via `envs/environment.yml` (conda). Key libraries: `numpy`, `scikit-learn`, `scipy` for binarisation/similarity scoring; R/Bioconductor (`limma`, `GEOquery`) for transcriptomic DEG analysis.

## Planned Repository Structure

The codebase is being built according to this pipeline (see README.md for full detail):

```
01_disease_map/       # CellDesigner XML maps (SjD Map + per-cell-type modules + assembled multicellular map)
02_boolean_model/     # CaSQ output, BMA input, raw attractors, filtered steady states
03_transcriptomics/   # Public GEO datasets (GSE23117, GSE40611, GSE84844, GSE157278) + processed outputs
04_validation/        # Literature validation notes, R scripts, calibrated model state
05_simulations/       # Drug KO, drug synergy, endotype-specific simulations
06_results/           # Figures and tables
envs/                 # Conda environment file
papers/               # Reference PDFs (disease_map.pdf, boolean_synovial.pdf)
```

## Scientific Pipeline

The modelling workflow proceeds in strict sequence:

1. **Disease map dissociation** — Parse SjD Map SBML from MINERVA; split into cell-type modules (SGECs, CD4+ Th1/Th17, B cells/plasmablasts, M1/M2 macrophages); add intercellular communication edges following CellDesigner conventions from Zerrouk et al., 2024.

2. **Boolean model generation** — Convert the assembled multi-cellular CellDesigner XML to an executable Boolean model using [CaSQ](https://github.com/sybila/casq):
   ```bash
   python casq.py --input SjD_multicellular_map.xml \
                  --output casq_output/SjS_boolean_model.json \
                  --threshold 0.5
   ```

3. **Attractor analysis** — Run [BMA](https://biomodelanalyzer.org/) on HPC to compute all fixed-point steady-state attractors.

4. **Transcriptomic validation** — DEG analysis on salivary gland GEO datasets using `limma` in R, then binarise expression vectors and compute similarity scores (1 − Hamming distance) against model steady states to identify the calibrated model state.

5. **In silico simulations** — Single and double gene KO experiments against the calibrated state; endotype-specific simulations using C1–C4 PRECISESADS signatures.

## Key Design Decisions

- **Intercellular edges**: secreted ligands use extracellular transport arrows to receptors on the receiving cell; cell-cell contact uses Heterodimer Complex Association — consistent with Zerrouk et al., 2024 CellDesigner conventions.
- **Validation metric**: 1 − Hamming distance between Boolean steady state and binarised DEG vector; nodes present in both the model and the expression data are used (intersection).
- **Primary validation datasets**: GSE23117 is the minimum benchmark (also used in the SjD Map paper). GSE40611, GSE84844, GSE157278 are secondary.
- **Phenotype readouts for simulations**: apoptosis, secretory function loss, lymphocytic infiltration markers, lymphomagenesis nodes.

## Reference Papers

- `papers/disease_map.pdf` — SjD Map paper (Silva-Saffar et al., *npj Syst Biol Appl*, 2026)
- `papers/boolean_synovial.pdf` — RA Boolean virtual twin (Zerrouk et al., *npj Digital Medicine*, 2024)

The RA model paper is the primary methodological template for this project.
