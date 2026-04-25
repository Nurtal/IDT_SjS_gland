# ROADMAP — SjS-DigitalTwin

> **Projet** : Jumeau numérique mécaniste de la glande salivaire dans la maladie de Sjögren primaire (pSS/SjD), basé sur un modèle Booléen multi-cellulaire.
> **Méthodologie de référence** : Zerrouk et al., *npj Digital Medicine*, 2024 (Atlas RA synovial).
> **Point de départ** : SjD Map (Silva-Saffar et al., *npj Syst Biol Appl*, 2026) — MINERVA.
> **Durée estimée** : 12 à 18 mois.
> **Statut** : Phase 0 démarrée (`scripts/01_download_sjd_map.py` opérationnel).

---

## Évaluation stratégique de l'étude

### Originalité scientifique

| Dimension | Évaluation | Justification |
|---|---|---|
| **Nouveauté du modèle** | Forte | Premier jumeau numérique Booléen *multi-cellulaire* dédié à la glande salivaire SjD. |
| **Nouveauté méthodologique** | Modérée | Méthode de Zerrouk (2024) déjà publiée et validée sur la RA — la valeur ajoutée méthodologique tient à l'adaptation tissue-spécifique (épithélium sécrétoire vs synoviocytes), pas à un cadre fondamentalement nouveau. |
| **Niche thérapeutique** | Forte | SjD n'a aucun traitement approuvé. Plusieurs essais récents (ianalumab, iscalimab, dazodalibep, remibrutinib) → besoin pressant de stratification mécaniste et de prédiction de combinaisons. |
| **Pont avec données patients** | Modérée à forte | Couplage attracteurs ↔ DEGs salivaire (4 datasets) + endotypes PRECISESADS = lien rare à grande échelle pour SjD. |

### Forces

1. **Réutilisation directe** d'une méthodologie validée par publication majeure (npj Digital Medicine, IF ≈ 15) — réduit le risque méthodologique.
2. **Carte de référence pré-existante** (SjD Map) : économie de 12–18 mois de curation manuelle.
3. **Datasets de validation accessibles** : 4 cohortes salivaire publiques + PRECISESADS (accès contrôlé) couvrant à la fois biopsies et profils sang.
4. **Application clinique tangible** : prédictions de KO ↔ médicaments en cours d'essai = pertinence translationnelle directe.
5. **Pipeline auditable** : `01_download_sjd_map.py` produit déjà des gates de validation chiffrés (≥800 espèces, 14 phénotypes, validation libsbml).

### Risques majeurs

| Risque | Probabilité | Impact | Mitigation prévue |
|---|---|---|---|
| **Tractabilité computationnelle** : ≈1000–1500 nœuds dépassent la capacité usuelle de BMA. | Élevée | Bloquant Phase 2 | POC multi-outils dès Phase 2.2 (mpbn, MaBoSS, pyboolnet, bioLQM-réduction). |
| **Dissociation cell-type ambiguë** dans la SjD Map (carte intégrée, pas par compartiment cellulaire). | Élevée | Curation très lourde | Stratégie hybride markers + pathway-enrichment + curation littérature (Phase 1.2). |
| **Validation signal/bruit** : score de Hamming sans baseline aléatoire = risque de p-hacking. | Moyenne | Critique reviewer | Distribution null obligatoire (1000 vecteurs aléatoires), p-value empirique. |
| **Endotypes blood→tissue** : signatures PRECISESADS issues du sang, modèle de la glande. | Moyenne | Limite interprétative | Discussion explicite + restriction aux nœuds présents simultanément dans l'expression sang et le modèle salivaire. |
| **Comparaison RA vs SjD** : Zerrouk 2024 est le benchmark — un reviewer demandera ce que cette étude apporte de plus. | Élevée | Acceptation | Cadrer l'apport sur (i) tissu épithélial sécrétoire absent de la RA, (ii) endotypes, (iii) prédiction de combinaisons d'essais en cours. |
| **Cibles thérapeutiques mal mappées** sur les nœuds (ex : iscalimab cible CD40L sur T, pas sur la cellule présentatrice). | Moyenne | Crédibilité clinique | Mapping drogue→nœud documenté DrugBank/ChEMBL avec PMID (Phase 3.1). |

### Chances de publication

**Estimation réaliste** :
- *npj Systems Biology and Applications* (cohérence avec Silva-Saffar 2026) — **probable** si Phase 2 + Phase 3 sont menées proprement, même sans validation expérimentale. IF ≈ 4. Audience appropriée.
- *npj Digital Medicine* (cohérence avec Zerrouk 2024) — **possible mais exigeant** : nécessitera (i) au moins une prédiction validée *a posteriori* contre un essai clinique récent (ex : succès/échec ianalumab Phase 2), (ii) endotypes interprétés avec des prédictions différentielles solides. IF ≈ 15.
- *Bioinformatics* / *Journal of Autoimmunity* — **alternatives sûres** si la dimension digital twin reste descriptive.

**Conditions critiques pour acceptation top-tier** :
1. Au moins une **prédiction nouvelle** non triviale (non déjà documentée par Zerrouk ou Silva-Saffar).
2. **Reproduction qualitative** d'au moins un essai clinique récent (rituximab inefficace en monothérapie SjD ; BAFF-inhibition prometteuse) → preuve de concept rétrospective.
3. **Différenciation des endotypes** mécaniquement crédible — pas seulement clustering.
4. **Robustesse** documentée (sensibilité aux règles, conditions initiales, bruit).

### Verdict global

L'étude est **pertinente, faisable et publiable**, à condition de traiter les risques tractabilité et validation statistique en amont. Elle se positionne dans la lignée d'un programme de recherche productif et reconnu (groupe Niarakis : RA Atlas, SjD Map, modèles macrophages M1/M2, FLS) plutôt que comme rupture conceptuelle. La fenêtre d'opportunité est ouverte tant qu'aucun autre groupe n'a publié l'équivalent SjD — surveiller régulièrement bioRxiv et la littérature `Niarakis A.[Author]`.

---

## Vue d'ensemble

```
Phase 0   Infrastructure & repro                                 [2 sem]
Phase 1   Construction de la carte multi-cellulaire              [4 mois]
  1.1     Audit topologique & sémantique de la SjD Map           [2 sem]
  1.2     Stratégie de dissociation cell-type                    [3 sem]
  1.3     Modules cell-type (4 sous-cartes)                      [6 sem]
  1.4     Edges intercellulaires                                 [3 sem]
  1.5     Assemblage & QC                                        [2 sem]
Phase 2   Modèle Booléen & analyse d'attracteurs                 [4 mois]
  2.1     Conversion CaSQ + curation des règles                  [4 sem]
  2.2     POC multi-outils (BMA / mpbn / MaBoSS / bioLQM)        [2 sem]
  2.3     Énumération & filtrage des attracteurs                 [3 sem]
  2.4     Validation transcriptomique                            [4 sem]
  2.5     Validation littérature                                 [2 sem]
Phase 3   Simulations in silico & manuscrit                      [4-6 mois]
  3.1     Mapping drogue → nœud                                  [2 sem]
  3.2     KO simples                                             [3 sem]
  3.3     Synergies (KO doubles)                                 [3 sem]
  3.4     Robustesse & sensibilité                               [2 sem]
  3.5     Endotypes (PRECISESADS)                                [3 sem]
  3.6     Rédaction & soumission                                 [6 sem]
```

---

## Phase 0 — Infrastructure & reproductibilité *(2 semaines)* — ✅ COMPLÈTE (2026-04-25)

### Objectif

Verrouiller l'environnement et les conventions **avant** de produire de la science, pour garantir la reproductibilité bit-à-bit.

### Statut

| Sous-tâche | Statut | Livrable |
|---|---|---|
| 0.1 Environnement scellé | ✅ | `envs/environment.yml` complété (pyboolnet, mpbn, edgeR, DESeq2, sva, MaBoSS) |
| 0.2 Conventions de nommage | ✅ | `CONVENTIONS.md` |
| 0.3 Dockerfile + compose | ✅ | `Dockerfile`, `docker-compose.yml`, `.dockerignore` |
| 0.4 CI/CD | ✅ | `.github/workflows/ci.yml`, `pyproject.toml`, 5 smoke tests OK |
| 0.5 Journal | ✅ | `journal.md` |

**Gate 0** : `[ ]` reproduction Docker à valider sur 2 machines distinctes (test reporté — pas de Docker accessible localement à cette session).

### Tâches

#### 0.1 Environnement scellé
- Compléter `envs/environment.yml` avec versions épinglées :
  - Ajouter : `pyboolnet`, `mpbn`, `colomoto-jupyter`, `bioconductor-edger`, `bioconductor-deseq2`, `bioconductor-sva`, `r-readr`, `r-survival`
  - Pin exact des versions critiques : `casq==1.4.3`, `python-libsbml==5.20.2`, `lxml==4.9.x`
- `Dockerfile` basé sur `colomoto/colomoto-docker:latest` (inclut CaSQ, MaBoSS, pyboolnet, bioLQM, GINsim) + couche projet avec `environment.yml`
- `docker-compose.yml` avec montage `./scripts:/work/scripts` pour développement live

#### 0.2 Conventions de nommage
- Document `CONVENTIONS.md` :
  - Identifiants espèces : préfixe HGNC (gènes), CHEBI (métabolites), GO (phénotypes)
  - Suffixe cell-type : `_SGEC`, `_TH1`, `_TH17`, `_BCELL`, `_PLASMA`, `_M1`, `_M2`, `_PDC` (pour pDC infiltrantes)
  - Compartiments : `extracellular`, `cytosol_<celltype>`, `nucleus_<celltype>`, `membrane_<celltype>`
  - Phénotypes (sortie) : `phen_<name>` (ex : `phen_apoptosis_SGEC`, `phen_secretory_loss`)

#### 0.3 CI/CD minimaliste
- GitHub Actions : `.github/workflows/ci.yml`
  - Job 1 : lint (ruff, mypy --strict sur `scripts/lib/`)
  - Job 2 : exécution `01_download_sjd_map.py` avec mocked HTTP (`responses` ou cache JSON commité)
  - Job 3 : pytest sur fonctions pures (`scripts/lib/celldesigner_xml.py` notamment)

#### 0.4 Versionning des artéfacts
- `01_disease_map/SjD_Map_original.xml` versionné via Git LFS (taille >1 Mo)
- `data/` (transcriptomes bruts) **gitignored** — récupération via script
- Tags Git : `v0.x` après chaque gate franchi

### Gate 0 (passage Phase 1)
- [ ] `docker run` reproduit `01_download_sjd_map.py` à l'identique sur 2 machines distinctes
- [ ] CI verte sur `main`
- [ ] `CONVENTIONS.md` validé et commité

### Livrables
- `Dockerfile`, `docker-compose.yml`
- `envs/environment.yml` (complété)
- `CONVENTIONS.md`
- `.github/workflows/ci.yml`

---

## Phase 1 — Construction de la carte multi-cellulaire *(4 mois)*

### Objectif

Transformer la SjD Map intégrée (mono-bloc) en une carte multi-cellulaire structurée par type cellulaire avec communications intercellulaires explicites, suivant les conventions CellDesigner de Zerrouk et al. 2024.

### Sous-phase 1.1 — Audit topologique & sémantique *(2 semaines)* — ✅ COMPLÈTE (2026-04-25)

**Statut** :

| Tâche | Résultat |
|---|---|
| Téléchargement SjD Map | ✅ 1157 espèces, 598 réactions, 14 phénotypes |
| Bug pagination MINERVA corrigé | ✅ déduplication par id |
| Bug `aliasId` dans `make_reaction` | ✅ corrigé (resolve fallback) |
| Audit topologique | ✅ 1157 nœuds, 1151 edges, 316 composantes faibles |
| Couverture annotations éléments | ✅ 94.21% (gate ≥80%) |
| Couverture PMID réactions | ⚠️ 16.72% (gate FAIL ≥50%) |
| Détection marqueurs cell-type | ⚠️ SGEC 0/19 — module à construire from-scratch |
| Détection cytokines clés | ✅ 25/36 cytokines clés présentes (92 nœuds avec triplets Gene/RNA/Protein) |
| Pathway hints depuis notes | ✅ 235 éléments avec ≥1 ID Reactome/KEGG/GO |

**Findings critiques pour Phase 1.2** :
1. **27% nœuds isolés** (315/1157, dont 290 Protein) → carte très fragmentée
2. **SGEC quasi-absent** (0 marqueurs canoniques AQP5/MUC5B/CFTR) → la SjD Map se concentre sur l'immunité, pas la fonction épithéliale → **module SGEC à construire manuellement** depuis littérature SjD spécifique
3. **125 triplets Gene+RNA+Protein** (sur 526 noms uniques) → **réduction formelle obligatoire** avant cell-type cloning, sinon explosion combinatoire (1157 × N_celltypes)
4. **6 compartiments** présents → compartiment-driven viable mais 93 nodes en `None`

**Livrables produits** :
- `01_disease_map/SjD_Map_original.xml` (2.0 Mo)
- `01_disease_map/cache/elements_raw.json` (3.4 Mo)
- `01_disease_map/cache/reactions_raw.json` (1.3 Mo)
- `01_disease_map/audit_summary.json`
- `01_disease_map/audit_report.md`
- `01_disease_map/audit_celltype_hits.tsv`
- `01_disease_map/audit_intercellular_pivots.tsv`

**Justification** : avant toute dissociation, il faut connaître précisément le matériau de départ. Aucun audit n'a été publié au-delà du papier original.

**Tâches** :
- Notebook `scripts/02_audit_map.ipynb` (sortie HTML) :
  - Topologie : composantes connexes, nœuds isolés, in/out-degree distribution, cycles fondamentaux, hubs (degree > 20)
  - Couverture sémantique :
    - % espèces avec annotation HGNC / UniProt / ChEBI / Ensembl
    - % réactions avec ≥1 PMID
    - Compartiments présents et nombre d'espèces par compartiment
    - Types CellDesigner (PROTEIN, RNA, GENE, COMPLEX, SIMPLE_MOLECULE, PHENOTYPE, DEGRADED…)
  - Cartographie pathway : enrichissement Reactome / KEGG via API REST sur les HGNC symbols
  - Vérification croisée : intersection avec les 14 nœuds phénotypes attendus (déjà implémentée)
- Identification des **nœuds-pivots** entre cell-types : ligands sécrétés (cytokines), récepteurs membranaires, complexes MHC

**Gate 1.1** :
- [ ] ≥80 % des espèces ont au moins une annotation externe
- [ ] Toutes les composantes connexes >5 nœuds documentées
- [ ] Mapping pathway Reactome/KEGG produit pour ≥70 % des espèces

**Livrable** : `01_disease_map/audit_report.html` + `01_disease_map/audit_summary.json`

### Sous-phase 1.2 — Stratégie de dissociation *(3 semaines)*

**Justification** : la SjD Map n'a pas (ou peu) d'annotation cell-type explicite — c'est l'étape la plus risquée de Phase 1. Approche hybride à valider.

**Trois approches combinées** :

1. **Compartiment-driven** (à vérifier en 1.1)
   - Si la SjD Map a déjà un compartiment par cell-type → assignation directe
   - Sinon → fallback sur 2/3

2. **Marker-driven** (forte spécificité)
   - Sentinelles cellulaires :
     - SGEC : `AQP5`, `MUC5B`, `KRT7`, `EPCAM`, `CFTR`
     - CD4+ Th1 : `CD3D`, `CD3E`, `IFNG`, `TBX21`, `STAT1`
     - CD4+ Th17 : `CD3D`, `IL17A`, `IL17F`, `RORC`, `STAT3`
     - B cell : `CD19`, `CD20/MS4A1`, `CD79A/B`, `BTK`
     - Plasma : `XBP1`, `IRF4`, `PRDM1`, `CD138/SDC1`
     - M1 : `CD68`, `IL6`, `TNF`, `IL1B`, `NOS2`
     - M2 : `CD163`, `MRC1`, `IL10`, `ARG1`
     - pDC (à inclure ?) : `CLEC4C/BDCA2`, `IL3RA`, `IRF7`
   - Propagation de l'assignation aux voisins par diffusion sur le graphe (random walk pondéré)

3. **Pathway-driven** (forte sensibilité)
   - Enrichissement par cell-type via signatures Reactome/MSigDB :
     - SGEC : `Salivary_secretion`, `Tight_junction`, `IFN_alpha_response`
     - T cell : `TCR_signaling`, `TH17_differentiation`, `JAK_STAT`
     - B cell : `BCR_signaling`, `Antibody_production`, `BAFF_signaling`
     - Macrophage : `TLR_signaling`, `Inflammasome`, `Phagocytosis`

**Combinaison** : score de confiance par nœud par cell-type (markers > pathway > compartment). Multi-assignation autorisée (un nœud peut appartenir à plusieurs cell-types — il sera **cloné** lors de l'assemblage).

**Livrables** :
- `01_disease_map/dissociation_rules.md` — règles auditables
- `01_disease_map/node_to_celltype.tsv` — `node_id, celltype, confidence, rule_triggered, evidence`
- `01_disease_map/extracellular_nodes.tsv` — ligands/cytokines partagés (cas spécial à ne pas dupliquer)

**Gate 1.2** :
- [ ] ≥90 % des nœuds non-extracellulaires assignés à ≥1 cell-type
- [ ] Les 10 % restants documentés et justifiés (orphelins / stroma / matrice / inconnu)
- [ ] Revue par expert (immunologiste / SjD specialist) du fichier `node_to_celltype.tsv`

### Sous-phase 1.3 — Modules cell-type *(6 semaines)*

**Tâches** :
- Pour chaque cell-type (4–5 modules : SGEC, Th1/Th17 fusionnés, B/Plasma fusionnés, M1/M2 fusionnés ; pDC en option) :
  1. Extraction du sous-graphe correspondant (`node_to_celltype.tsv`)
  2. Export CellDesigner XML : `01_disease_map/SjD_Map_celltype_modules/<celltype>_map.xml`
  3. Revue littérature ciblée : ≥10 PMID par cell-type pour identifier les pathways manquants
     - SGEC : autophagie, sénescence cellulaire (Manoussakis 2020), pyroptose induite par IFN
     - Tfh (sous-type T) : axe IL-21/CXCR5 (à ajouter si absent)
     - Plasma : différenciation BLIMP1/IRF4/XBP1
     - Macrophages : polarisation M1/M2 et plasticité
  4. Compléments manuels dans CellDesigner 4.4.2 (interface graphique)
  5. Vérification de la **fonctionnalité** : ≥1 boucle de feedback identifiable par module, ≥1 phénotype de sortie attaché

- Outils :
  - Script `scripts/03_split_celltype.py` : extraction automatique sous-graphes en SBML
  - Script `scripts/04_validate_module.py` : checks topologiques par module

**Gate 1.3** :
- [ ] 4 modules ouverts proprement dans CellDesigner 4.4.2
- [ ] Chaque module : ≥1 feedback loop, ≥1 phénotype de sortie
- [ ] Chaque module : ≥80 nœuds, ≥50 réactions (échelle Zerrouk)
- [ ] Revue par expert SjD signée

**Livrable** : `01_disease_map/SjD_Map_celltype_modules/*.xml` + rapport de complétude

### Sous-phase 1.4 — Edges intercellulaires *(3 semaines)*

**Justification** : c'est ici que la valeur ajoutée multi-cellulaire se construit. Doit être systématique, pas ad hoc.

**Sources** :
1. **CellPhoneDB v4** (Efremova et al., *Nat Protoc*, 2020) — base ligand-récepteur curée
2. **OmniPath** (Türei et al., *Mol Syst Biol*, 2016) — interactions inter-cellulaires intégrées
3. **Curation littérature SjD-spécifique** :
   - IFN-α (pDC → SGEC) — signature IFN type I
   - BAFF / APRIL (SGEC, Mφ → B cell) — survie des lymphocytes B
   - CXCL13 (SGEC → B/T cells) — formation des structures lymphoïdes ectopiques (TLS)
   - IL-21 (Tfh → B cell) — différenciation plasmocytaire
   - IL-6, IL-7, IL-15, IL-22 (SGEC → T cell) — survie/recrutement
   - LIGHT / LTα1β2 (T → SGEC) — destruction épithéliale
   - Anti-Ro/La (Plasma → SGEC, autocrine) — réaction immune locale

**Conventions CellDesigner (cf. Zerrouk 2024)** :
- Ligand sécrété : `extracellular_transport` arrow → récepteur
- Contact cell-cell : `Heterodimer Complex Association`
- Auto-cycle (ligand sécrété + récepteur sur la même cellule, ex : autocrine BAFF) : transport via `extracellular`

**Tâches** :
- Script `scripts/05_intercellular_edges.py` : intersection CellPhoneDB ∩ OmniPath filtrée sur les paires cell-types présentes
- Curation manuelle SjD : table `01_disease_map/intercellular_edges.tsv` avec colonnes :
  - `source_celltype, target_celltype, ligand, receptor, mechanism, evidence_PMID, in_cellphonedb, in_omnipath, sjs_specific`

**Gate 1.4** :
- [ ] ≥30 edges intercellulaires (cible : 40–60, échelle Zerrouk RA)
- [ ] Chaque edge a ≥1 PMID ou est référencé dans CellPhoneDB/OmniPath
- [ ] Couverture des 4 cytokines/voies clés en SjD : IFN-α, BAFF, CXCL13, IL-21

### Sous-phase 1.5 — Assemblage & QC *(2 semaines)*

**Tâches** :
- Script `scripts/06_assemble_multicellular.py` :
  - Merge des 4 modules + edges intercellulaires
  - Clonage des nœuds multi-assignés avec préfixes `_SGEC`, `_TH17`, etc.
  - Création des compartiments cell-type-spécifiques
- QC topologique automatique :
  - Pas de nœud orphelin (degree 0)
  - Pas de réaction sans réactif/produit valide
  - Cohérence de signe sur les chemins critiques (ex : IFN-α → STAT1 → IRF7 → IFN-α doit être cohérent)
- Visualisation : ouverture dans CellDesigner 4.4.2 + export PNG layouté

**Gate 1.5 (passage Phase 2)** :
- [ ] Carte ouvre dans CellDesigner 4.4.2 sans erreur
- [ ] ≥1200 espèces, ≥800 réactions, ≥30 edges intercellulaires
- [ ] Validation libsbml : 0 erreur fatale
- [ ] Toutes les références d'espèces résolues
- [ ] 14 phénotypes présents
- [ ] Revue scientifique par 2 experts (immunologue + modélisateur)

**Livrable** : `01_disease_map/SjD_multicellular_map.xml` + `01_disease_map/multicellular_map.png`

---

## Phase 2 — Modèle Booléen & analyse d'attracteurs *(4 mois)*

### Objectif

Convertir la carte CellDesigner en modèle Booléen exécutable, énumérer les attracteurs, et identifier l'attracteur correspondant à l'état SjD via validation transcriptomique.

### Sous-phase 2.1 — Conversion CaSQ & curation des règles *(4 semaines)*

**Justification** : la sortie brute de CaSQ n'est jamais directement utilisable. C'est l'étape la plus chronophage de Phase 2.

**Tâches** :
- Run CaSQ avec paramètres canoniques :
  ```bash
  casq --input 01_disease_map/SjD_multicellular_map.xml \
       --output 02_boolean_model/casq_output/SjS_boolean.sbml \
       --csv 02_boolean_model/casq_output/SjS_boolean.csv \
       --bnet 02_boolean_model/casq_output/SjS_boolean.bnet
  ```
- Audit automatique des règles générées (`scripts/07_audit_rules.py`) :
  - Nœuds **sources** sans entrée (input nodes) → forcer état initial ou supprimer
  - Nœuds **puits** sans sortie (terminal nodes hors phénotypes) → vérifier
  - Modifieurs CellDesigner sans signe explicite → curation manuelle
  - Détection des règles triviales (`x = x`) ou inconsistantes
- Curation manuelle ciblée :
  - Pour chaque hub (degree > 10) : revue de la règle générée vs littérature
  - Conversion d'opérateurs OR par défaut → AND quand justifié biologiquement (ex : activation T cell = TCR AND co-stimulation)
- Réduction formelle via `bioLQM` : élimination des nœuds tampons / simplification de constantes

**Outils** :
- `bioLQM` (Naldi 2018) pour réduction formelle
- `pyboolnet` pour parsing et manipulation
- `colomoto-jupyter` pour notebooks interactifs

**Gate 2.1** :
- [ ] Modèle Booléen valide, importable dans `pyboolnet` et `MaBoSS`
- [ ] ≥95 % des règles auditées (manuellement ou par script)
- [ ] Log de curation `02_boolean_model/curation_log.md` : pour chaque règle modifiée → ancienne règle, nouvelle règle, justification, PMID

**Livrable** : `02_boolean_model/casq_output/SjS_boolean.{sbml,bnet,csv}` + `curation_log.md`

### Sous-phase 2.2 — POC multi-outils *(2 semaines)*

**Justification** : BMA est le choix par défaut Zerrouk, mais sa scalabilité est limitée (~500 nœuds confortables, dégradation forte au-delà). Modèle multi-cellulaire ≈ 1000–1500 nœuds → **BMA seul est risqué**.

**POC parallèle sur le modèle Phase 2.1** :

| Outil | Méthode | Scaling | Avantage | Limite |
|---|---|---|---|---|
| **BMA** | Stabilization Theorem | ≤500 nœuds confortable | Méthode du papier RA (Zerrouk) | Sync seulement, scaling limité |
| **mpbn** | Trap-spaces minimaux (≃ attracteurs async) | 10⁴+ nœuds | Très scalable, async | Pas de dynamique temporelle |
| **MaBoSS** | Continuous-time Markov, async | 10³+ nœuds | Stochastique, populations, dynamique | Configuration .bnd/.cfg lourde |
| **pyboolnet** | SAT-based attractor detection | ≈10³ nœuds | Async, point fixes + cycles | Plus lent que mpbn pour énumération |
| **bioLQM** | Réduction formelle | n/a | Préprocessing avant énumération | Outil de réduction, pas de simulation |

**Stratégie** :
1. Run **mpbn** en premier → trap-spaces minimaux comme proxy d'attracteurs async (rapide, scalable)
2. Run **MaBoSS** pour la dynamique temporelle et la robustesse stochastique
3. Run **BMA** sur HPC en parallèle pour cohérence avec papier de référence (si tractable)
4. Run **pyboolnet** comme cross-check sur point fixes

**Gate 2.2** :
- [ ] ≥1 outil produit un set complet d'attracteurs en <24h sur le modèle complet
- [ ] Cohérence des point fixes entre ≥2 outils (les point fixes async = sync = trap-spaces ponctuels)
- [ ] Document `02_boolean_model/tool_choice.md` justifiant l'outil retenu pour la suite

**Livrable** : `02_boolean_model/poc_results/` (3–4 sous-dossiers, un par outil) + `tool_choice.md`

### Sous-phase 2.3 — Énumération & filtrage *(3 semaines)*

**Tâches** :
- Énumération exhaustive des attracteurs (point fixes + trap-spaces minimaux + cycles) via l'outil retenu
- Annotation de chaque attracteur par état des 14 phénotypes (vecteur ON/OFF)
- Filtrage en 3 catégories :
  1. **Disease attractors** : ≥1 phénotype "disease" actif (Inflammation, Apoptosis, Fibrosis, B_Cell_Activation, T_Cell_Activation/Differentiation)
  2. **Healthy-like attractors** : tous phénotypes "disease" inactifs
  3. **Mixed / ambigus** : à curer cas par cas

**Gate 2.3** :
- [ ] ≥1 attracteur "disease" identifié (sinon le modèle ne capture pas SjD → revoir Phase 1.4 / 2.1)
- [ ] ≥1 attracteur "healthy-like" identifié (sinon pas de baseline pour les KO)
- [ ] Stabilité numérique : 2 runs indépendants (random seed) produisent les mêmes attracteurs

**Livrable** : `02_boolean_model/attractors_disease.csv`, `attractors_healthy.csv`, `phenotype_signature_per_attractor.csv`

### Sous-phase 2.4 — Validation transcriptomique *(4 semaines)*

**Justification** : c'est ici que le modèle est calibré sur la réalité. Étape critique pour la crédibilité du papier.

**4 datasets** : GSE23117 (microarray, n=15), GSE40611 (microarray, n=31), GSE84844 (RNA-seq, n=40), GSE157278 (RNA-seq, n=26).

**Pipeline R harmonisé** (`04_validation/scripts/deg_pipeline.R`) :
- Chargement via `GEOquery` ; normalisation : `RMA` (microarray) ou `vst`/`DESeq2` (RNA-seq)
- DEG :
  - Microarray : `limma` avec `voom` si raw counts disponibles
  - RNA-seq : `edgeR` ou `DESeq2`
- Batch correction : `sva::ComBat` si effets batch détectés
- Méta-analyse : Fisher combined p-values ou `metafor` ; FDR < 0.05
- Output : table DEG par dataset + table consensus

**Binarisation** :
- Seuil **par dataset** : top 25 % overexpressed (FC > Q3 des FC positifs significatifs) → 1 ; bottom 25 % → 0 ; reste → masqué
- Justification : seuil global = sensibilité aux différences de plateforme
- Mapping HGNC symbol ↔ nom de nœud (table de correspondance dans `04_validation/node_gene_mapping.tsv`)

**Score de similarité** :
- 1 − Hamming distance sur intersection {nœuds modèle} ∩ {gènes binarisés non masqués}
- **Baseline aléatoire obligatoire** : 1000 vecteurs binaires aléatoires de même densité → distribution null par attracteur
- p-value empirique : P(score_random ≥ score_observed)
- Correction multiple sur le nombre d'attracteurs testés (Bonferroni ou BH)

**Critère de calibration** :
- Attracteur retenu = celui avec p < 0.01 sur GSE23117 **et** p < 0.05 sur ≥1 autre dataset
- Si aucun attracteur ne passe → diagnostic : revoir Phase 1.4 (intercellulaires manquants ?) ou 2.1 (règles trop OR ?)

**Gate 2.4** :
- [ ] DEG tables produites pour les 4 datasets, FDR < 0.05
- [ ] ≥1 attracteur calibré (p<0.01 sur GSE23117 + p<0.05 sur ≥1 autre)
- [ ] Reproductibilité : permutation des labels patients/contrôles → score aléatoire, p > 0.5 (sanity check)

**Livrables** : `03_transcriptomics/processed/DEG_tables/*.tsv`, `binary_vectors/*.tsv`, `04_validation/similarity_scores.csv`, `04_validation/calibrated_state.csv`

### Sous-phase 2.5 — Validation littérature *(2 semaines)*

**Justification** : un score statistique seul ne suffit pas — l'attracteur calibré doit reproduire des faits biologiques connus.

**Checklist de 15–20 faits attendus pour SjD** (à valider sur l'attracteur calibré) :

| # | Fait biologique attendu | Référence type |
|---|---|---|
| 1 | Signature IFN type I active (STAT1, IRF7, ISGs ON) | Mavragani 2017 |
| 2 | Signature IFN type II cohérente (IFNG-driven) | Nezos 2015 |
| 3 | BAFF/TNFSF13B surexprimé | Mariette 2003 |
| 4 | BTK actif dans B cells | Verstappen 2019 |
| 5 | AQP5 down/délocalisé dans SGEC | Steinfeld 2001 |
| 6 | Axe Th17/IL-17 actif | Lin 2014 |
| 7 | FoxP3/Tregs réduits ou dysfonctionnels | Liu 2014 |
| 8 | CXCL13 élevé (TLS markers) | Salomonsson 2002 |
| 9 | Apoptose SGEC active (Fas/FasL, caspases) | Manganelli 2003 |
| 10 | MHC-II up sur SGEC | Tsunawaki 2002 |
| 11 | Plasma cells / IRF4 / XBP1 actif | Hansen 2003 |
| 12 | M1 polarisation prédominante | Greenwell-Wild 2011 |
| 13 | Lymphomagenesis nodes (NF-κB, MYC) actifs ou activables | Nezos 2014 |
| 14 | JAK-STAT axis actif | Sjöstrand 2024 |
| 15 | Type III IFN (IFNL) cohérent | Apostolou 2016 |

**Score qualitatif** : nombre de faits reproduits / 15. Cible : ≥12/15 (80 %).

**Gate 2.5 (passage Phase 3)** :
- [ ] Score littérature ≥ 12/15
- [ ] Discussion explicite des faits non reproduits + cause probable (manquant Phase 1, règle CaSQ, donnée trop spécifique)
- [ ] Validation par expert SjD signée

**Livrable** : `04_validation/literature_validation.md`

---

## Phase 3 — Simulations in silico & manuscrit *(4–6 mois)*

### Objectif

Utiliser le modèle calibré pour prédire les effets de perturbations thérapeutiques (KO simples, synergies), simuler l'hétérogénéité endotypique, démontrer la robustesse, et publier.

### Sous-phase 3.1 — Mapping drogue → nœud *(2 semaines)*

**Justification** : les "KO drogues" doivent reposer sur un mapping documenté, pas sur des suppositions.

**Tâches** :
- Table de mapping `05_simulations/drug_target_mapping.tsv` :
  - Colonnes : `drug_name, status_in_SjD (approved/Phase3/Phase2/preclinical), target_gene, target_node, action (inhibit/activate), evidence_PMID, drugbank_id, chembl_id`
- Drogues prioritaires (essais SjD récents) :
  - **Ianalumab (anti-BAFF)** → KO `TNFSF13B` ou `TNFRSF13C/BAFF-R`
  - **Iscalimab (anti-CD40L)** → KO `CD40LG` (sur T) ou `CD40` (sur B/SGEC)
  - **Dazodalibep (anti-CD40L)** → idem iscalimab
  - **Rituximab (anti-CD20)** → KO `MS4A1` puis cascade B cell
  - **Tofacitinib (JAK1/3i)** → KO `JAK1` + `JAK3`
  - **Filgotinib (JAK1i sélectif)** → KO `JAK1`
  - **Remibrutinib (BTKi)** → KO `BTK`
  - **Nipocalimab / Efgartigimod (anti-FcRn)** → KO `FCGRT`
  - **Abatacept (CTLA4-Ig)** → activation `CTLA4` ou KO `CD80`/`CD86`
  - **Hydroxychloroquine** → modulation `TLR7`/`TLR9` (pour pDC)
  - **Anifrolumab (anti-IFNAR)** → KO `IFNAR1`/`IFNAR2`
  - **Belimumab (anti-BAFF)** → similaire ianalumab

**Gate 3.1** :
- [ ] ≥10 drogues mappées avec PMID + DrugBank ID
- [ ] Pour chaque drogue : justification du choix de nœud cible (un seul vs plusieurs)

**Livrable** : `05_simulations/drug_target_mapping.tsv`

### Sous-phase 3.2 — KO simples *(3 semaines)*

**Tâches** :
- Pour chaque drogue : forcer le nœud cible à 0 (ou 1 pour activateurs) à partir de l'attracteur calibré
- Re-énumération des attracteurs accessibles (subset basin)
- Quantification des **changements de phénotypes** :
  - Vecteur phénotype `pre = [phen_apoptosis_SGEC, phen_secretory_loss, phen_inflammation, ...]`
  - Vecteur phénotype `post`
  - Distance L1 normalisée + breakdown par phénotype
- Définition d'un "score thérapeutique" agrégé : somme pondérée des phénotypes "disease" éteints (avec poids cliniques : secretory_loss > inflammation > infiltration)

**Visualisations** :
- Heatmap drogues × phénotypes (ON→OFF en vert, OFF→ON en rouge, inchangé en gris)
- Scatter : score thérapeutique prédit vs efficacité clinique connue (validation rétrospective : rituximab faible, ianalumab fort, etc.)

**Gate 3.2** :
- [ ] Heatmap complète pour ≥10 drogues
- [ ] Reproduction qualitative d'au moins 2 résultats d'essais cliniques connus
- [ ] Identification de ≥1 cible non-évidente (drogue avec score élevé non encore en essai)

**Livrables** : `05_simulations/drug_KO/results.csv`, `06_results/figures/single_KO_heatmap.svg`

### Sous-phase 3.3 — Synergies (KO doubles) *(3 semaines)*

**Tâches** :
- Screening exhaustif sur les n cibles cliniquement pertinentes (n ≈ 10 → C(10,2) = 45 paires)
- Pour chaque paire (A, B) :
  - Score combiné `s_AB`
  - Score additif attendu `s_A + s_B - (s_A * s_B)` (modèle Bliss)
  - Synergie = `s_AB - s_Bliss`
- Top-10 paires synergiques avec interprétation mécanistique :
  - Quel chemin de la carte explique la synergie ?
  - Cohérence avec données précliniques connues (ex : anti-BAFF + anti-CD20 = synergie attendue)

**Gate 3.3** :
- [ ] ≥45 paires testées
- [ ] Top-10 paires synergiques identifiées et interprétées mécaniquement
- [ ] ≥1 synergie nouvelle non publiée (apport original du papier)

**Livrables** : `05_simulations/drug_synergy/synergy_matrix.csv`, `06_results/figures/synergy_heatmap.svg`

### Sous-phase 3.4 — Robustesse & sensibilité *(2 semaines)*

**Justification** : sans cette section, un reviewer top-tier rejettera. C'est le test de stress du modèle.

**Tests** :
1. **Sensibilité aux conditions initiales** :
   - 1000 inits aléatoires (densité 0.5)
   - Distribution des bassins d'attraction
   - % d'inits convergeant vers l'attracteur calibré
   - Cible : ≥10 % (pas trop niche, pas trop dominant)

2. **Sensibilité aux règles** :
   - Flip aléatoire de 1 %, 2 %, 5 % des règles → re-énumération attracteurs
   - L'attracteur calibré reste-t-il accessible ?
   - Cible : ≥80 % à 1 %, ≥60 % à 5 %

3. **Sensibilité au bruit MaBoSS** :
   - Run MaBoSS avec différents taux de bruit (`time_tick`, `max_time`)
   - Probabilité stationnaire de l'attracteur calibré
   - Cible : probabilité stable sur 3 ordres de grandeur de paramètres

4. **Sensibilité au choix de l'outil** :
   - Re-faire les KO simples top-5 avec un second outil (ex : MaBoSS si mpbn principal)
   - Cohérence des prédictions

**Gate 3.4** :
- [ ] Tous les tests passent les seuils
- [ ] Document `04_validation/robustness.md` avec figures de sensibilité

**Livrable** : `06_results/figures/robustness_panels.svg`

### Sous-phase 3.5 — Endotypes *(3 semaines)*

**Justification** : c'est l'apport translationnel le plus original — différencier mécaniquement les endotypes C1–C4 PRECISESADS.

**⚠ Limite à discuter explicitement** : signatures PRECISESADS issues du **sang**, modèle de la **glande**. Restreindre l'analyse aux nœuds présents dans :
- L'expression sang (gènes mesurés en RNA-seq)
- ET le modèle salivaire
- ET avec une expression non-tissue-spécifique documentée

**Tâches** :
- Initialiser le modèle avec le vecteur d'expression binarisé de chaque endotype C1–C4 (sur l'intersection)
- Énumérer attracteurs accessibles depuis chaque init endotype
- Comparer signatures phénotypiques par endotype
- Tester les top-10 drogues par endotype → matrice drogue × endotype × score thérapeutique
- Identifier les drogues **endotype-spécifiques** (efficaces sur un endotype, pas sur un autre)

**Gate 3.5** :
- [ ] Différenciation mécanistique des 4 endotypes (signatures phénotypiques distinctes)
- [ ] ≥2 prédictions d'endotype-specific drug response interprétables
- [ ] Discussion claire des limites blood-vs-tissue

**Livrables** :
- `05_simulations/endotype_simulations/endotype_attractors.csv`
- `06_results/figures/endotype_drug_response.svg`

### Sous-phase 3.6 — Rédaction & soumission *(6 semaines, en parallèle)*

**Cibles de revues** (par ordre d'ambition décroissante) :
1. *npj Digital Medicine* (IF ≈ 15) — si ≥1 prédiction validée rétrospectivement contre essai clinique
2. *npj Systems Biology and Applications* (IF ≈ 4) — cible réaliste, cohérent avec Silva-Saffar 2026
3. *Bioinformatics* / *Briefings in Bioinformatics* — alternative outil-centré
4. *Journal of Autoimmunity* — alternative clinical-centré

**Structure du manuscrit** :
- **Fig 1** : Workflow général (carte → modèle → validation → simulations)
- **Fig 2** : Carte multi-cellulaire (statistiques + visualisation)
- **Fig 3** : Attracteur calibré + score similarité + validation littérature (panels A-D)
- **Fig 4** : Heatmap KO simples + reproduction rétrospective d'essais cliniques
- **Fig 5** : Top synergies avec interprétation mécanistique
- **Fig 6** : Endotypes — signatures phénotypiques + drug response différentielle
- **Fig S1-Sn** : Robustesse, sensibilité, validation détaillée par dataset

**Code & données** :
- Tag GitHub release `v1.0`
- Dépôt Zenodo avec DOI (carte XML, modèle Booléen, attracteurs, scripts)
- BioModels deposit (modèle SBML-qual)
- MINERVA upload (carte multi-cellulaire) en complément de la SjD Map originale

**Gate 3.6** :
- [ ] Manuscrit soumis
- [ ] Code/données archivés (Zenodo + BioModels)

---

## Récapitulatif des gates scientifiques

| Phase | Gate critique | Critère bloquant si échec |
|---|---|---|
| 0 | Repro Docker | Ne pas démarrer Phase 1 |
| 1.2 | ≥90 % nœuds assignés | Revoir stratégie de dissociation |
| 1.5 | Carte assemblée + ≥1200 espèces | Revoir Phase 1.3 ou 1.4 |
| 2.1 | Modèle Booléen valide | Revoir conversion / curation règles |
| 2.2 | ≥1 outil tractable | Réduction de modèle obligatoire avant Phase 2.3 |
| 2.3 | ≥1 attracteur disease + ≥1 healthy | Revoir règles ou edges intercellulaires |
| 2.4 | ≥1 attracteur calibré (p<0.01) | Revoir Phase 1.4 ou 2.1 — projet en risque |
| 2.5 | ≥12/15 faits littérature | Revoir Phase 1.3 — projet en risque |
| 3.4 | Robustesse ≥80 % | Revoir Phase 2.1 |
| 3.6 | Manuscrit soumis | — |

---

## Risques transverses & mitigations

| Risque | Détection | Mitigation |
|---|---|---|
| Dérive temporelle (>18 mois) | Suivi mensuel des gates | Réduction du périmètre : pDC retirées si Phase 1 dérape |
| Concurrence (autre groupe publie SjD digital twin) | Veille bioRxiv mensuelle (`Niarakis`, `Sjogren`, `digital twin`, `Boolean model`) | Soumission preprint à Gate 2.5 pour priorité |
| Indisponibilité PRECISESADS (accès contrôlé) | Demande accès dès Phase 0 | Backup : GSE51092 (whole blood, public) + signatures publiées de Tarn et al. 2019 |
| Bug critique dans la SjD Map MINERVA | Audit Phase 1.1 | Communication directe avec auteurs Silva-Saffar / Niarakis |
| Échec scaling Booléen | POC Phase 2.2 | Réduction formelle bioLQM + restriction aux 5–6 cell-types core |

---

## Veille active à maintenir tout au long du projet

- **bioRxiv** : alertes hebdomadaires sur "Sjogren digital twin", "Sjogren Boolean model", "salivary gland computational model"
- **PubMed** : alertes mensuelles sur Niarakis A., Mariette X., Pers JO., Cornec D.
- **ClinicalTrials.gov** : suivi des essais Phase 2/3 SjD pour valider rétrospectivement les prédictions
- **MINERVA** : nouvelles versions de la SjD Map
- **CellPhoneDB / OmniPath** : nouvelles releases pour edges intercellulaires

---

## Décisions ouvertes (à trancher avant Phase 1.3)

1. **Inclure pDC ?** — pertinent biologiquement (signature IFN-α majoritaire) mais ajoute un cell-type → coût Phase 1.3 + Phase 1.4. **Recommandation** : oui, modulo gate 1.1 montre couverture suffisante des marqueurs pDC.
2. **Th1 et Th17 fusionnés ou séparés ?** — fusion = simplification ; séparation = différenciation mécanistique. **Recommandation** : séparés, car endotypes peuvent différer sur cet axe.
3. **B cells et Plasma cells fusionnés ou séparés ?** — différenciation BLIMP1/IRF4 importante. **Recommandation** : séparés.
4. **Validation expérimentale (wet lab) ?** — hors périmètre actuel mais boost majeur publication. **Recommandation** : si collaboration accessible (CHU Brest), prioriser validation d'1 prédiction non-évidente in vitro sur lignée SGEC (HSY ou A253).

---

*Dernière mise à jour : 2026-04-25*
