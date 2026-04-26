# Dissociation Rules — SjD Map → modules cell-type

> **Phase 1.2 ROADMAP** — règles auditables d'assignation des nœuds de la SjD Map
> à un ou plusieurs types cellulaires (SGEC, TH1, TH17, TFH, TREG, BCELL, PLASMA, M1, M2, PDC).
>
> **Statut** : v2 — 2026-04-26 (ajout score plausibilité + sources).
> **Auteur** : SjS-DigitalTwin pipeline (`scripts/03_dissociate.py`).

---

## 1. Contexte & contrainte de départ

### 1.1 Topologie de la SjD Map (audit Phase 1.1)

La SjD Map MINERVA n'est **pas** organisée par type cellulaire mais par compartiment fonctionnel :

| compartmentId | Nom | n_species | Sémantique |
|---|---|---|---|
| 20513 | Cell | 624 | Intracellulaire générique (cytosol/membrane) |
| 21231 | nucleus | 308 | Noyau (Gene, RNA, TF) |
| 21555 | Extracellular_ligands | 73 | Ligands externes (cytokines reçues) |
| 21629 | Secreted_molecules | 40 | Molécules sécrétées (cytokines émises) |
| 21540 | Phenotypes | 14 | Nœuds de sortie phénotypiques |
| 20730 | Endoplasmic_reticulum | 5 | RE (Ca²⁺ signaling) |
| None | (no compartment) | 93 | Nœuds satellites (mix) |

**Conséquence** : les approches "compartiment-driven cell-type-specific" et "marker matching naïf" sont insuffisantes seules. On bascule sur une **approche hybride à 4 règles** avec multi-assignment autorisé.

### 1.2 Hypothèse fondamentale

Un nœud de pathway immunitaire générique (ex : `STAT1`, `JAK1`, `NFKB1`) est biologiquement actif dans **plusieurs cell-types**. Pour rester fidèle à la biologie, on autorise **multi-assignment** : un même nœud peut produire plusieurs clones cell-type-spécifiques (`STAT1_SGEC`, `STAT1_TH17`, etc.) avec des règles potentiellement différentes en aval.

C'est la convention adoptée par Zerrouk et al. 2024 (RA Atlas) pour les nœuds intracellulaires partagés.

---

## 2. Cell-types cibles

| Cell-type | Code | Inclusion |
|---|---|---|
| Salivary gland epithelial cells | `SGEC` | ✅ obligatoire (cible tissulaire) |
| CD4+ T helper 1 | `TH1` | ✅ |
| CD4+ T helper 17 | `TH17` | ✅ |
| T follicular helper | `TFH` | ✅ |
| T regulatory | `TREG` | ✅ |
| B cells (naïves/mémoire/GC) | `BCELL` | ✅ |
| Plasma cells | `PLASMA` | ✅ |
| Macrophage M1 | `M1` | ✅ |
| Macrophage M2 | `M2` | ✅ |
| Plasmacytoid DC | `PDC` | ✅ (audit : 9 marqueurs trouvés) |

**Décision** : 10 cell-types. SGEC obligatoire malgré 0/19 marqueurs canoniques natifs (cf. règle pathway-driven).

---

## 3. Règles d'assignation (ordre de priorité)

Les règles sont **appliquées dans l'ordre** ; la première qui produit un assignement non-vide donne le résultat. Une règle peut produire un set de cell-types (multi-assignment).

### Règle R1 — Compartiment extracellulaire → `EXTRA`

**Condition** : `compartmentId ∈ {21555 (Extracellular_ligands), 21629 (Secreted_molecules)}`

**Action** : assigner uniquement à `EXTRA` (pseudo-cell-type) → ce nœud sera **unique** dans la carte multi-cellulaire (pas de clonage), avec edges entrants (sécrétion) et sortants (réception) vers les cell-types appropriés.

**Confiance** : `HIGH`.

**Justification** : ces 113 nœuds sont par construction des intermédiaires intercellulaires (ligands sécrétés, antigènes, MHC sécrétés). Leur clonage par cell-type serait redondant.

### Règle R2 — Phénotype → `PHENOTYPE` (sortie globale)

**Condition** : `type == 'Phenotype'`

**Action** : assigner à un pseudo-cell-type `PHENOTYPE` (sortie globale, pas dupliquée).

**Confiance** : `HIGH`.

**Note** : pour les phénotypes cell-type-spécifiques (`Apoptosis` = SGEC + lymphocytes ; `B_Cell_Activation/Survival` = BCELL/PLASMA), un raffinement manuel sera fait en Phase 1.5 (assemblage). À Phase 1.2, on les garde globaux.

### Règle R3 — Marqueur cell-type exclusif → mono ou intra-lignée

Matching exact du nom (insensible à la casse, après strip) contre une liste de marqueurs canoniques par lignée :

| Lignée | Marqueurs exclusifs |
|---|---|
| **SGEC-only** | `AQP5`, `AQP3`, `MUC5B`, `MUC7`, `KRT7`, `KRT8`, `KRT18`, `KRT19`, `EPCAM`, `CFTR`, `CDH1`, `CLDN1`, `CLDN3`, `CLDN4`, `OCLN`, `AMY1A`, `PRB1`, `STATH` |
| **T-lineage (TH1+TH17+TFH+TREG)** | `CD3D`, `CD3E`, `CD3G`, `CD4`, `CD2`, `LCK`, `ZAP70`, `LAT`, `LCP2`, `CD28` (TCR signaling commun) |
| **TH1-only** | `TBX21`, `STAT4`, `IFNG` (mais cytokine → R1 si extracellulaire) |
| **TH17-only** | `IL17A`, `IL17F`, `RORC`, `RORA`, `IL23R`, `AHR` |
| **TFH-only** | `CXCR5`, `BCL6`, `ICOS`, `PDCD1`, `IL21` (ligand → R1) |
| **TREG-only** | `FOXP3`, `IL2RA`, `IKZF2` |
| **B-lineage (BCELL+PLASMA)** | `CD19`, `MS4A1`/`CD20`, `CD79A`, `CD79B`, `CR2`, `BLNK`, `BANK1` |
| **BCELL-only** | `BTK`, `BLK`, `TNFRSF13B`/`TACI`, `TNFRSF13C`/`BAFFR`, `FCRL3`, `MS4A1` |
| **PLASMA-only** | `XBP1`, `IRF4`, `PRDM1`/`BLIMP1`, `SDC1`/`CD138`, `MZB1`, `JCHAIN`, `TNFRSF17`/`BCMA` |
| **Macrophage-lineage** | `CD68`, `FCGR1A`, `FCGR2A`, `FCGR3A`, `CSF1R`, `MERTK`, `ITGAM`/`CD11B` |
| **M1-only** | `NOS2`, `CXCL10`, `CXCL11`, `CCR7` (en combinaison avec `CD68`) |
| **M2-only** | `CD163`, `MRC1`, `ARG1`, `MSR1`, `CCL22` (cytokine → R1) |
| **PDC-only** | `CLEC4C`/`BDCA2`, `IL3RA`/`CD123`, `LILRA4`, `IRF7`, `TLR7`, `TLR9` |

**Action** : assigner aux cell-types listés. **Confiance** : `HIGH`.

**Note** : si un nœud match un marqueur SGEC-only ET un marqueur d'un autre cell-type (cas rare), on flag pour revue manuelle.

### Règle R4 — Pathway-driven (depuis pathway hints des notes)

**Condition** : R1, R2, R3 non déclenchées **et** le nœud a un identifiant Reactome/KEGG/GO dans ses notes (extrait par `extract_pathway_hints_from_notes`).

**Action** : assigner aux cell-types associés au(x) pathway(s) trouvé(s) selon la table :

| Pathway / Reactome ID partiel | Cell-types assignés |
|---|---|
| `R-HSA-1280215` (Cytokine signaling) | TOUS |
| `R-HSA-877300` (Interferon signaling) | TOUS (avec emphasis SGEC, PDC, M1) |
| `R-HSA-877253` (Interferon alpha/beta signaling) | SGEC, PDC, BCELL, M1, M2 |
| `R-HSA-877312` (Interferon gamma signaling) | SGEC, M1, M2 |
| `R-HSA-983705` (Signaling by the BCR) | BCELL, PLASMA |
| `R-HSA-202403` (TCR signaling) | TH1, TH17, TFH, TREG |
| `R-HSA-168256` (Immune system) | TOUS |
| `R-HSA-168249` (Innate Immune System) | M1, M2, PDC, SGEC |
| `R-HSA-168256` (Adaptive Immune System) | TH1, TH17, TFH, TREG, BCELL, PLASMA |
| `R-HSA-449147` (TLR signaling) | M1, M2, PDC, SGEC |
| `R-HSA-983168` (Antigen processing/MHC-I) | TOUS sauf TREG |
| `R-HSA-2132295` (MHC-II) | M1, M2, BCELL, PDC, SGEC |
| `R-HSA-1474244` (Extracellular matrix) | SGEC |
| `R-HSA-372790` (Signaling by GPCR) | TOUS |
| `R-HSA-451927` (Interleukin-2 family signaling) | TH1, TH17, TFH, TREG |
| `R-HSA-512988` (Interleukin-6 signaling) | TOUS |
| `hsa04060` (KEGG cytokine-cytokine receptor) | TOUS |
| `hsa04062` (KEGG chemokine signaling) | TOUS |
| `hsa04630` (KEGG JAK-STAT) | TOUS |
| `hsa04064` (KEGG NF-κB) | TOUS |
| `hsa04668` (KEGG TNF signaling) | TOUS |
| `hsa04611` (KEGG platelet activation) | — (skip) |
| `hsa04920` (KEGG adipocytokine) | — (skip) |

**Confiance** : `MEDIUM`.

### Règle R5 — Propagation par voisinage (random walk)

**Condition** : R1–R4 non déclenchées **et** le nœud a ≥1 voisin (in ou out) déjà assigné par R1–R4.

**Action** : assigner aux cell-types majoritaires parmi les voisins assignés (vote pondéré par confiance R3 > R4).

**Confiance** : `LOW`.

**Justification** : un nœud connecté à un hub TCR sera probablement T-lineage. Heuristique imparfaite mais utile en l'absence d'autre signal.

### Règle R6 — Default fallback (pathway intracellulaire générique)

**Condition** : R1–R5 non déclenchées **et** `type ∈ {'Protein', 'Complex', 'Simple molecule'}` **et** nœud non isolé.

**Action** : assigner à `{SGEC, TH1, TH17, BCELL, M1, M2}` (cell-types principaux) — les nœuds générique de signaling sont présents partout.

**Confiance** : `LOW`.

**Justification** : pour ne pas perdre de nœuds, on suppose que le pathway est partagé. La curation manuelle Phase 1.3 affinera.

### Règle R6c — Raffinement expert (Phase 1.2bis)

**Condition** : nœud assigné par R6 **et** dont le nom (normalisé majuscules + alphanum) match une entrée de la table curée `scripts/lib/refinement.py`.

**Action** : remplacer l'assignation R6 par le set `cell-types` curé manuellement à partir de :

- **PanglaoDB** (Franzén 2019, PMID:30951143) — cell-type marker compendium
- **CellMarker 2.0** (Hu 2023, PMID:36300620) — manually curated single-cell markers
- **Human Protein Atlas** (Uhlén 2015, PMID:25613900) — tissue/cell expression
- **InnateDB** (Breuer 2013, PMID:23180781) — innate immunity curation
- **Reactome** (Jassal 2020, PMID:31691815) — pathway expression context
- **KEGG pathway maps** (Kanehisa 2023)
- **ImmGen** (Heng 2008, PMID:18800157) — immune cell transcriptional atlases
- **Littérature SjD primaire** : Mavragani 2017 (PMID:28604219) IFN signature ; Manganelli 2003 (PMID:12796328) Fas/FasL SGEC ; Lavie 2004 BAFF ; Verstappen 2021 review
- **Lineage TFs fondateurs** : Hori 2003 (FOXP3), Ivanov 2006 (RORγt), Szabo 2000 (T-bet), Johnston 2009 (BCL6), Steinfeld 2001 (AQP5)

**Score** : 60–88 selon spécificité du knowledge :
- 85–88 : restriction très forte avec consensus (ex. CD22 BCELL ; FLT3 PDC)
- 80–84 : restriction modérée validée par 2+ DBs
- 70–79 : axe biologique défini (ex. ISG → ALL ; Smad → ALL)
- 60–69 : housekeeping informé par HPA

**Confiance** : `MEDIUM` si score ≥ 75 sinon `LOW`.

**Source** : citation explicite dans la table (`DB1; DB2; PMID:…`).

**Cas spécial — DROP** : si le nœud est explicitement non-pertinent au modèle (ex. `IL5` Th2-only, `KDR` endothélial), la table renvoie un set vide, et le nœud bascule en R7 (drop expert).

**Sortie typique** : 589 nœuds R6 → R6c, 4 nœuds R6 → R7 (drop), 8 nœuds R6 sans match (conservés en R6 score 25).

### Règle R7 — Inassignable

**Condition** : R1–R6 non déclenchées (typiquement nœud isolé sans annotation, sans nom HGNC) **ou** drop explicite par R6c.

**Action** : `UNASSIGNED`.

**Décision Phase 1.3** : revue manuelle du fichier `unassigned_nodes.tsv` ; soit ajouter un marqueur, soit retirer le nœud.

---

## 4. Score de plausibilité et source

À chaque triplet (nœud, cell-type, règle) le dissociator attache :

- `plausibility_score` (0–100) — défendabilité biologique, indépendant de la confiance computationnelle
- `source` — référence primaire (PMID/doi) ou base de données (HGNC, UniProt, Reactome, KEGG, InnateDB, STRING, CellMarker2.0, PanglaoDB, HPA)

| Tranche | Sémantique | Règle typique |
|---|---|---|
| 95–100 | Marqueur lineage-defining, consensus immunologique | R1 (ECM), R3 mono-target |
| 80–94  | Marqueur fonctionnel à forte spécificité (DB curées) | R3 lignée, R4 pathway très spécifique (TCR, BCR) |
| 60–79  | Pathway-driven cohérent (Reactome/KEGG) | R4 ubiquitaire, R3 partagé large |
| 40–59  | Propagation par voisinage, validée par graphe | R5 |
| 20–39  | Fallback générique (signaling intracellulaire ubiquitaire) | R6 connecté |
| 1–19   | Fallback isolé, faible confiance | R6 isolé |
| 0      | Non assigné | R7 |

### Règles → score base + source par défaut

| Rule | Score base | Source par défaut |
|---|---|---|
| R1 | 95 | MINERVA SjD Map compartments 21555/21629 (Silva-Saffar 2026) |
| R2 | 90 | MINERVA SjD Map Phenotypes layer (Silva-Saffar 2026) |
| R3 | 95 (mono) / 82 (lignée ≤4) / 70 (large) | HGNC + PMID fondateur du marqueur (cf. `MARKER_SOURCE` dans `dissociator.py`) |
| R4 | 65 (override 70–82 pour TCR/BCR/Th17/IFN) | Reactome:R-HSA-… ou KEGG:hsaXXXXX |
| R5 | 45 + bonus vote (≤60) | propagation graphe + InnateDB Breuer 2013 PMID:23180781 + STRING v12 PMID:36370105 |
| R6 | 25 (connecté) / 15 (isolé) | InnateDB + STRING + HPA Uhlén 2015 PMID:25613900 |
| R6c | 60–88 selon spécificité | Table curée `refinement.py` (PanglaoDB, CellMarker 2.0, HPA, InnateDB, Reactome, KEGG, ImmGen, PMID primaires) |
| R7 | 0 | — |

## 5. Sortie attendue

### 5.1 `node_to_celltype.tsv`

Format (long) — une ligne par paire (nœud, cell-type) :

```tsv
node_id  node_name        node_type  compartment_id  celltype  confidence  rule  evidence            plausibility_score  source
20756    BTK              Protein    20513           BCELL     HIGH        R3    marker=BTK          95                  HGNC:1133; UniProt Q06187; KEGG:hsa04662 BCR
21240    STAT1 homodimer  Protein    21231           TH1       MEDIUM      R4    R-HSA-877300        65                  Reactome:R-HSA-877300 via MINERVA notes
21240    STAT1 homodimer  Protein    21231           BCELL     MEDIUM      R4    R-HSA-877300        65                  Reactome:R-HSA-877300 via MINERVA notes
...
```

### 5.2 `extracellular_nodes.tsv`

Sous-ensemble des nœuds R1 (extracellulaires) :

```tsv
node_id    node_name    node_type    compartment_id    canonical_form    g_r_p_variants
21580      CXCL13       Protein      21629             CXCL13            21490(Gene),21262(RNA)
...
```

### 5.3 `dissociation_summary.json`

Statistiques agrégées par règle, par cell-type, par compartiment + score moyen de plausibilité par cell-type, n_high (≥80), n_low (<40).

---

## 6. Gate 1.2 (passage Phase 1.3)

| Critère | Cible | Vérification |
|---|---|---|
| Couverture nœuds non-extracellulaires | ≥90% assignés (R1–R6) | `dissociation_summary.json` |
| Distribution cell-type | chaque cell-type a ≥80 nœuds | idem |
| Reviewers expert | revue de `node_to_celltype.tsv` signée | manuel |
| Cohérence biologique | aucun marqueur exclusif mal assigné (spot-check 20 nœuds aléatoires) | manuel |

---

## 7. Limites assumées

1. **R6 (default fallback)** sur-assigne par défaut à 6 cell-types — corrigé en Phase 1.2bis par R6c (table `refinement.py`).
2. **R5 (random walk)** propage le bruit si un voisin a été mal assigné. Mitigation : ne déclencher R5 que si ≥2 voisins concordants.
3. **Phénotypes** restent globaux à Phase 1.2 ; cell-type-specificité reportée à Phase 1.5.
4. **Triplets Gene+RNA+Protein** sont assignés indépendamment — si discordance, alerter (rare attendu).

---

*Dernière mise à jour : 2026-04-26*
