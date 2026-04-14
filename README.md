# SjS-DigitalTwin — Boolean Modelling of the Salivary Gland in Sjögren's Disease

> **Status:** Work in progress | **Started:** 2026  
> **Contact:** [nathan.foulquier@chu-brest.fr]

---

## Overview

This project builds a **mechanistic immune digital twin of the salivary gland** in primary Sjögren's disease (pSS/SjD). It adapts the modular Boolean modelling framework developed for Rheumatoid Arthritis synovial tissue (Zerrouk et al., *npj Digital Medicine*, 2024) to model the pathological interactions between salivary gland epithelial cells (SGECs) and infiltrating immune cells in SjD.

The starting point is the **SjD Map** (Silva-Saffar et al., *npj Syst Biol Appl*, 2026) — the first comprehensive Molecular Interaction Map for Sjögren's disease, available on MINERVA at [https://sjdmap.elixir-luxembourg.org/](https://sjdmap.elixir-luxembourg.org/).

The model is validated using:
1. **Public salivary gland transcriptomic datasets** (primary validation)

---

## Scientific Background

### Why a digital twin for SjD?

SjD is a prototypic systemic autoimmune disease with no approved treatment. A major obstacle to therapeutic development is **patient heterogeneity**: at least four molecular endotypes have been identified in blood (Sorret et al., *Nat Commun*, 2021), but the mechanisms driving salivary gland destruction across these endotypes remain poorly understood.

Static transcriptomic analyses describe *what* is dysregulated. A Boolean model describes *how* — it simulates the dynamic behaviour of signalling networks, identifies stable disease phenotypes (attractors), and predicts the effect of therapeutic perturbations in silico.

### Conceptual framework

```
SjD Map (SBML/CellDesigner)
        │
        ▼
Cell-type dissociation
  ├── Salivary gland epithelial cells (SGECs)
  ├── Infiltrating CD4+ T cells (Th1/Th17)
  ├── B cells / plasma cells
  └── Macrophages (M1/M2)
        │
        ▼
Multi-cellular map (CellDesigner XML)
  + Intercellular communication edges
        │
        ▼
CaSQ → Executable Boolean model
        │
        ▼
BMA (HPC) → Attractors / Steady states
        │
        ▼
Validation against transcriptomic data
  ├── Salivary gland datasets (GSE23117, GSE40611, ...)
        │
        ▼
In silico simulations
  ├── Drug KO experiments
  ├── Drug synergy prediction
  └── Endotype-specific phenotype modelling
```

---

## Repository Structure

```
SjS-DigitalTwin/
│
├── 01_disease_map/
│   ├── SjD_Map_original.xml          # CellDesigner XML from MINERVA (SjD Map)
│   ├── SjD_Map_celltype_modules/     # Cell-type-specific sub-maps
│   │   ├── SGEC_map.xml
│   │   ├── CD4_Th1_map.xml
│   │   ├── Bcell_map.xml
│   │   └── Macrophage_M1M2_map.xml
│   └── SjD_multicellular_map.xml     # Assembled multi-cellular map
│
├── 02_boolean_model/
│   ├── casq_output/                  # CaSQ-generated Boolean model files
│   ├── bma_input/                    # BMA-formatted model for attractor analysis
│   ├── attractors/                   # Raw attractor output from BMA (HPC)
│   └── steady_states_filtered.csv    # Filtered steady states used for validation
│
├── 03_transcriptomics/
│   ├── public_salivary_gland/        # Public SG datasets (see Data section)
│   │   ├── GSE23117/
│   │   ├── GSE40611/
│   │   └── GSE84844/
│   └── processed/
│       ├── DEG_tables/
│       ├── binary_vectors/           # Discretised expression for model validation
│       └── similarity_scores.csv
│
├── 04_validation/
│   ├── literature_validation.md      # Known biology reproduced by model
│   ├── transcriptomic_validation.R   # Similarity score computation
│   └── calibrated_state.csv         # Final calibrated model state
│
├── 05_simulations/
│   ├── drug_KO/                      # Single gene KO simulations
│   ├── drug_synergy/                 # Double KO combination experiments
│   └── endotype_simulations/         # Endotype-specific initial conditions
│
├── 06_results/
│   ├── figures/
│   └── tables/
│
├── envs/
│   └── environment.yml               # Conda environment
│
└── README.md
```

---

## Public Datasets Used

### Salivary Gland Transcriptomics (Primary Validation)

| Dataset | Platform | Design | N (SjS / Ctrl) | Access |
|---|---|---|---|---|
| **GSE23117** | Affymetrix microarray | Salivary gland biopsies, pSS vs controls | 10 / 5 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23117) |
| **GSE40611** | Affymetrix microarray | Minor salivary gland, pSS vs controls | 21 / 10 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40611) |
| **GSE84844** | Illumina RNA-seq | Parotid gland, anti-SSA+ pSS vs controls | 24 / 16 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84844) |
| **GSE157278** | Illumina RNA-seq | Labial salivary gland, pSS + scRNA context | 18 / 8 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157278) |

> **Note:** GSE23117 is already used in the SjD Map paper as a validation overlay. It is the minimum benchmark dataset for this project.

### Whole Blood RNA-seq (Secondary Validation — Systemic/Immune Components)

| Dataset | Cohort | Design | N | Access |
|---|---|---|---|---|
| **PRECISESADS** | Multi-centre EU | Whole blood RNA-seq, pSS + 5 other AID | 304 pSS / 341 ctrl | [ELIXIR Luxembourg](https://doi.org/10.17881/th9v-xt85) — controlled access |
| **GSE51092** | — | Microarray, whole blood pSS vs ctrl | 190 / 32 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51092) |
| **UKPSSR** | UK pSS registry | Microarray, whole blood | 151 / 29 | [EGA](https://www.ebi.ac.uk/ega/) — controlled access |

---

## Methods

### Step 1 — Disease Map Dissociation into Cell-Type Modules

The SjD Map integrates pathways from multiple cell types. We dissociate the integrated map into cell-type-specific modules following the approach of the RA Atlas:

- **SGECs** — ductal/acinar epithelial cell pathways: IFN signalling, apoptosis, aquaporin/secretory function, BAFF/APRIL production, MHC-I/II presentation
- **CD4+ Th1 / Th17 cells** — TCR signalling, cytokine production (IFN-γ, IL-17), co-stimulation
- **B cells / plasmablasts** — BCR signalling, antibody production, BTK pathway, lymphomagenesis nodes
- **M1/M2 macrophages** — TLR signalling, inflammasome, phagocytosis

Intercellular communication edges (cytokines, direct contact) are added manually following the CellDesigner convention used in Zerrouk et al., 2024:
- Secreted ligands → extracellular transport arrows → receptor on receiving cell
- Cell-cell contact → Heterodimer Complex Association

### Step 2 — Boolean Model Generation (CaSQ)

The multi-cellular CellDesigner XML is converted to an executable Boolean model using [CaSQ](https://github.com/sybila/casq) (Aghamiri et al., *Bioinformatics*, 2020):

```bash
# Example CaSQ conversion
python casq.py \
  --input SjD_multicellular_map.xml \
  --output casq_output/SjS_boolean_model.json \
  --threshold 0.5
```

### Step 3 — Attractor Analysis (BMA on HPC)

The Boolean model is analysed using [BMA (Bio Model Analyser)](https://biomodelanalyzer.org/) deployed on a high-performance computing cluster to identify all steady-state attractors. Fixed-point attractors are retained as candidate disease phenotypes.

### Step 4 — Transcriptomic Validation

**4a. DEG analysis on salivary gland datasets**

```r
# R — limma-based DEG analysis (consistent with SjD Map paper)
library(limma)
library(GEOquery)

# Load GSE23117
gse <- getGEO("GSE23117")
eset <- gse[[1]]

design <- model.matrix(~ 0 + group, data = pData(eset))
fit <- lmFit(eset, design)
contrast <- makeContrasts(pSS - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
DEGs <- topTable(fit2, adjust = "BH", p.value = 0.05, number = Inf)
```

**4b. Binarisation and similarity scoring**

Expression of each biomolecule present in the Boolean model is discretised (overexpressed = 1, underexpressed/unchanged = 0) and compared to model steady states using cosine similarity or Hamming distance, following the validation framework of Zerrouk et al., 2024.

```python
import numpy as np
from scipy.spatial.distance import hamming

def similarity_score(steady_state, binary_expression_vector):
    """Compute 1 - Hamming distance between model steady state and observed expression."""
    common_nodes = list(set(steady_state.keys()) & set(binary_expression_vector.keys()))
    ss_vec = np.array([steady_state[n] for n in common_nodes])
    exp_vec = np.array([binary_expression_vector[n] for n in common_nodes])
    return 1 - hamming(ss_vec, exp_vec)
```

The steady state with the highest similarity score defines the **calibrated model state**.

### Step 5 — In Silico Simulations

Using the calibrated model, we perform:

1. **Single KO simulations** — knock out individual nodes (drug targets) and observe phenotype readouts (apoptosis, secretory function loss, lymphocytic infiltration markers, lymphomagenesis)
2. **Double KO combinations** — screen for synergistic drug pairs
3. **Endotype-specific simulations** — initialise model with endotype-specific expression vectors (C1/C2/C3/C4 from PRECISESADS) to simulate endotype-specific disease phenotypes and predict differential drug responses

---

## Key Tools & Dependencies

| Tool | Purpose | Reference |
|---|---|---|
| [CellDesigner 4.4.2](http://www.celldesigner.org/) | Map building & visualisation | Funahashi et al., 2003 |
| [CaSQ](https://github.com/sybila/casq) | Map → Boolean model conversion | Aghamiri et al., *Bioinformatics*, 2020 |
| [BMA](https://biomodelanalyzer.org/) | Attractor analysis on HPC | Hall & Fisher, *Curr Protoc Bioinform*, 2020 |
| [MINERVA](https://minerva.pages.uni.lu/) | Map visualisation & overlays | — |
| R / Bioconductor | Transcriptomic analysis | — |
| Python (scikit-learn, numpy) | Binarisation, similarity scoring | — |

```bash
# Install conda environment
conda env create -f envs/environment.yml
conda activate sjs-digitaltwin
```

---

## Key References

1. **SjD Map:** Silva-Saffar E, Mariette X, [...] Niarakis A. *The SjD Map: an interactive pathway tour into Sjögren's disease signalling mechanisms.* npj Syst Biol Appl (2026). https://doi.org/10.1038/s41540-026-00670-x

2. **RA Boolean Virtual Twin:** Zerrouk N, Augé F, Niarakis A. *Building a modular and multi-cellular virtual twin of the synovial joint in Rheumatoid Arthritis.* npj Digital Medicine (2024). https://doi.org/10.1038/s41746-024-01396-y

3. **SjD Endotypes:** Cornec D, [...] Pers JO. *Identification and characterisation of endotypes in primary Sjögren's syndrome.* Nat Commun (2021). https://doi.org/10.1038/s41467-021-23472-7

4. **CaSQ:** Aghamiri SS, [...] Niarakis A. *Automated inference of Boolean models from molecular interaction maps using CaSQ.* Bioinformatics (2020). https://doi.org/10.1093/bioinformatics/btaa484

5. **RA macrophage Boolean models:** Zerrouk N, Alcraft R, Hall BA, Augé F, Niarakis A. *Large-scale computational modelling of the M1 and M2 synovial macrophages in rheumatoid arthritis.* npj Syst Biol Appl (2024).

6. **RA FLS Boolean model:** Singh V, Naldi A, Soliman S, Niarakis A. *A large-scale Boolean model of the rheumatoid arthritis fibroblast-like synoviocytes predicts drug synergies in the arthritic joint.* npj Syst Biol Appl (2023).

---

## Roadmap

- [ ] **Phase 1** — Multi-cellular map construction *(~6 months)*
  - [ ] Download and parse SjD Map SBML from MINERVA
  - [ ] Dissociate SjD Map into cell-type modules (SGECs, T, B, macrophages)
  - [ ] Add intercellular communication edges
  - [ ] Assemble and curate multi-cellular map in CellDesigner

- [ ] **Phase 2** — Boolean model generation & validation *(~6 months)*
  - [ ] Run CaSQ on multi-cellular map
  - [ ] Compute attractors with BMA on HPC
  - [ ] DEG analysis on salivary gland public datasets (GSE23117, GSE40611, GSE84844)
  - [ ] Compute similarity scores / calibrate model
  - [ ] Validate against known biology (IFN signatures, BAFF/APRIL, BTK)

- [ ] **Phase 3** — In silico simulations & manuscript *(~6 months)*
  - [ ] Drug KO simulations (JAK inhibitors, BTK inhibitors, anti-BAFF, rituximab)
  - [ ] Endotype-specific simulations with PRECISESADS signatures
  - [ ] Manuscript: "A Boolean immune digital twin of the salivary gland in Sjögren's disease"

---

## License

Code: MIT License. Maps and curated biological content: CC BY 4.0 (consistent with SjD Map licensing).
