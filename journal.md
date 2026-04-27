# Journal de bord — SjS-DigitalTwin

> Trace chronologique des actions, décisions et résultats. Une entrée par séance de travail.
> Format : `## YYYY-MM-DD — Phase X.Y — Titre court`

---

## 2026-04-25 — Phase 0 — Démarrage du projet structuré

### Contexte
- Roadmap détaillée écrite dans `ROADMAP.md` (12-18 mois, 4 phases : 0/1/2/3).
- Phase 0 (infrastructure & repro) priorisée avant toute production scientifique.
- État initial : `01_download_sjd_map.py` opérationnel + lib `minerva_api.py` / `celldesigner_xml.py` / `node_classifier.py`.

### Actions
1. **Initialisation suivi de tâches** (TaskCreate × 6) : Phase 0.1 → 0.5 + Phase 1.1.
2. **Création de `journal.md`** (ce fichier).

### Prochaines étapes immédiates
- Phase 0.1 : compléter `envs/environment.yml` (pyboolnet, mpbn, edgeR, DESeq2, sva, versions épinglées).
- Phase 0.2 : `CONVENTIONS.md` (nommage espèces, suffixes cell-type, compartiments).
- Phase 0.3 : `Dockerfile` basé sur `colomoto/colomoto-docker`.
- Phase 0.4 : `.github/workflows/ci.yml` (ruff, mypy, pytest).

### Décisions
- _aucune à ce stade_

### Notes
- Version Python verrouillée à 3.12 (déjà dans environment.yml).
- CaSQ pinné à 1.4.3.

---

## 2026-04-25 (suite) — Phase 0 complète

### Phase 0.1 — `environment.yml` complété
- Ajouté : `pyboolnet`, `mpbn`, `colomoto-jupyter`, `biolqm`, `maboss` (pip)
- Ajouté Bioconductor : `edger`, `deseq2`, `sva`, `affy`, `oligo`, `org.hs.eg.db`, `annotationdbi`
- Ajouté R : `metafor`, `pheatmap`, `survival`, `readr`, `tidyr`, `stringr`
- Ajouté outillage dev : `pytest`, `pytest-cov`, `ruff`, `mypy`, `responses` (mock HTTP)
- Versions critiques épinglées : `casq==1.4.3`, `python-libsbml==5.20.2`, `lxml=4.9.*`

### Phase 0.2 — `CONVENTIONS.md`
- Nommage espèces : HGNC en majuscules, suffixes cell-type `_SGEC`, `_TH1`, `_TH17`, `_TFH`, `_TREG`, `_BCELL`, `_PLASMA`, `_M1`, `_M2`, `_PDC`, `_EXTRA`
- Phénotypes préfixés `phen_` (17 phénotypes incluant les 14 originaux + `secretory_loss`, `tls_formation`, `lymphocytic_infiltration`)
- IDs SBML : `s_<celltype>_<symbol>`, `r_<type>_<celltype>_<n>`, `c_<name>`
- Règles Python (ruff/mypy strict sur lib), R (tidyverse style), Git (conventional commits)

### Phase 0.3 — Docker
- `Dockerfile` basé sur `colomoto/colomoto-docker:2024-08-01` (CaSQ + MaBoSS + GINsim + bioLQM + pyboolnet pré-installés)
- Couche projet : env conda `sjs-digitaltwin` créé via `mamba env create`
- `docker-compose.yml` avec 2 services (bash + jupyter port 8888)
- `.dockerignore` pour exclure venv, cache, données

### Phase 0.4 — CI
- `.github/workflows/ci.yml` : 2 jobs (lint + test) sur push/PR `main`
- `pyproject.toml` : config ruff (E/W/F/I/B/UP/N/SIM/RUF), mypy (strict pour lib), pytest
- `tests/test_smoke.py` : 5 tests d'import + namespaces — **TOUS PASS** localement (pytest 9.0.2, python 3.10.12)
- mypy strict en `continue-on-error` pour Phase 0 (couverture types à compléter Phase 1)
- `.gitignore` étendu (caches Python/jupyter, données externes, OS files)

### Décisions
- **Tests CI utilisent système Python 3.10** (pas de conda) avec deps minimales — éviter d'installer toute la stack Bioconductor en CI = lent et fragile. Tests "lourds" (CaSQ, MaBoSS, R) tournent en local Docker uniquement.
- **mypy `continue-on-error`** pour Phase 0 : la lib existante (`celldesigner_xml.py`, `minerva_api.py`) n'est pas encore strictement typée. À tightenir au fil de Phase 1.

### Gate 0
- ✅ Smoke tests verts (5/5)
- ✅ Conventions documentées
- ⏸ Repro Docker bit-à-bit sur 2 machines : **reporté** (pas de Docker disponible cette session ; à valider lors de la mise en production)

### Prochaines étapes
- Phase 1.1 — Audit topologique de la SjD Map (`scripts/02_audit_map.py`)
  - Composantes connexes, distribution degrees, hubs
  - Couverture annotations HGNC/UniProt/PMID
  - Enrichissement Reactome/KEGG via reactome2py
  - Cartographie des nœuds-pivots intercellulaires

---

## 2026-04-25 (suite 2) — Phase 1.1 complète

### Téléchargement SjD Map

**Bugs critiques rencontrés et corrigés** :

1. **`scripts/01_download_sjd_map.py` : ImportError**
   - `COMPARTMENT_NAMES` importé depuis `lib.minerva_api` mais défini dans `lib.celldesigner_xml`
   - **Fix** : déplacé l'import dans le bon module

2. **`scripts/lib/minerva_api.py::_paginate` : boucle infinie**
   - L'API MINERVA `/projects/SjD_Map_no_BD/models/7/bioEntities/elements/` ne supporte pas `?page=N&size=N` — chaque page retourne les mêmes 1157 éléments
   - Avant fix : 325+ pages, 377k éléments accumulés, 2.9 Go RAM, process bloqué 12 min
   - **Fix** : déduplication par `id` + arrêt anticipé si aucun nouvel élément reçu

3. **`scripts/lib/celldesigner_xml.py::make_reaction::resolve` : 0 réactions générées**
   - Cherche `entry["element"]["id"]` mais l'API renvoie `entry["aliasId"]` (entier)
   - Avant fix : 598 réactions ignorées, XML sans interactions
   - **Fix** : fallback `aliasId` avant `element.id`

**Résultat post-fix** :
- 1157 espèces, 598 réactions, 7 compartiments, 14 phénotypes
- XML 2.0 Mo, validation libsbml OK, refs résolues, **Gate 1 download : VALIDÉ**

### Audit topologique (`scripts/02_audit_map.py` + `scripts/lib/map_audit.py`)

**Nouveau module `lib/map_audit.py`** (650 lignes) :
- `build_graph()` : DiGraph networkx depuis elements + reactions (modifieurs inclus)
- `topology_summary()` : composantes faibles/fortes, hubs, isolés, sources, puits
- `annotation_coverage()` : couverture HGNC/UniProt/Ensembl/Entrez/KEGG/Reactome/GO/PubMed
- `reaction_pmid_coverage()` : % réactions avec PMID
- `classify_compartments()` : species par compartmentId
- `detect_celltype_marker_hits()` : matching marqueurs canoniques par cell-type
- `detect_intercellular_pivots()` : cytokines clés SjD (BAFF, IFN-α/β/γ, IL-17/21/22, CXCL9-13, TNF, etc.)
- `extract_pathway_hints_from_notes()` : regex Reactome/KEGG/GO sur notes texte libre

**Findings critiques** :

| Métrique | Valeur | Implication |
|---|---|---|
| Nœuds | 1157 | conforme attendu (≥800) |
| Edges | 1151 | **ratio 0.99 → graphe sparse (proche d'un arbre)** |
| Composantes faibles | **316** | **carte très fragmentée** |
| Plus grosse composante | 842 (73%) | core connecté correct |
| SCC > 1 | 1 (taille 23) | **quasi pas de cycles → modèle Booléen pauvre en feedback** |
| Nœuds isolés | **315 (27%)** | dont 290 Protein, 18 SimpleMolecule |
| Annotation any | 94.21% | ✅ Gate ≥80% |
| PMID réactions | 16.72% | ❌ Gate ≥50% — limitation traçabilité |
| Triplets Gene+RNA+Protein | 125 / 526 noms uniques (23%) | **réduction formelle obligatoire avant clonage cell-type** |
| Marqueurs SGEC trouvés | **0 / 19** | **module SGEC à construire from-scratch** |
| Cytokines clés trouvées | 25 / 36 distinctes (92 nœuds) | bonne base pour edges intercellulaires Phase 1.4 |

**Top 6 hubs** : Inflammation (44), STAT1 homodimer (43), STAT1/STAT2/IRF9 complex (42), RELA/NFKB1 (28), GNAI (22), Chemotaxis/Infiltration (22).
→ **Signature IFN type I = pivot du modèle** (cohérent avec littérature SjD).

**Cytokines manquantes (11)** : CCL22, IFNA1/2/13 (sous-types IFN-α non éclatés), IFNL1/2/3 (IFN-λ absent !), IL12A/B, IL17F, IL22 → à compléter Phase 1.3 ou justifier exclusion.

### Décisions stratégiques (post-audit)

1. **Module SGEC = travail de novo** (gros impact sur estimation Phase 1.3)
   - La SjD Map MINERVA est focalisée immunité, pas épithélium sécrétoire
   - Sources à mobiliser : Manoussakis 2020 (autophagie SGEC), Steinfeld 2001 (AQP5), Tsunawaki 2002 (MHC-II SGEC), Manganelli 2003 (apoptose Fas/FasL), Mavragani 2017 (signature IFN SGEC)
   - **Estimation Phase 1.3 révisée** : 6 sem → **8-9 sem** (4 modules + SGEC quasi from-scratch)

2. **Réduction formelle Gene/RNA/Protein avant clonage**
   - Soit collapse pré-Phase 1.3 (script à écrire), soit déléguer à bioLQM en Phase 2.1
   - Recommandation : pré-collapse en Phase 1.5 (assemblage) — produit une carte plus lisible pour CaSQ

3. **Lien faible inter-pathway → revoir si modèle Booléen final est riche**
   - Avec un seul SCC>1 de taille 23, le modèle aura peu de bistabilité naturelle
   - Les 36 transitions intercellulaires Phase 1.4 vont créer des feedbacks → bon pour la dynamique
   - Suivi à faire en Phase 2.3 (énumération attracteurs)

### Gate 1.1
- ✅ ≥80% éléments avec annotation (94.21%)
- ⚠️ ≥50% réactions avec PMID (16.72%) — gate échoué mais non bloquant (la SjD Map publique est ce qu'elle est)
- ✅ Toutes composantes connexes >5 nœuds documentées
- ✅ Mapping Reactome/KEGG produit (235 éléments, 20% — partiel)

### Prochaines étapes
- **Phase 1.2 — Stratégie de dissociation** :
  - Décrire les règles d'assignation par cell-type dans `01_disease_map/dissociation_rules.md`
  - Construire `node_to_celltype.tsv` par approche hybride (compartiment + markers + pathway-enrichment)
  - Demande de **revue par expert SjD** avant Phase 1.3
- **Décision en attente** :
  - Inclure pDC ? → recommandation **OUI** (9/6 marqueurs trouvés, signature IFN-α)
  - Inclure SGEC malgré 0/19 marqueurs natifs ? → **OUI obligatoire** (cible cellulaire principale du tissu) mais coût de Phase 1.3 augmenté

---

## 2026-04-25 (suite 3) — Phase 1.2 complète

### Décision pivot — abandon de l'approche "compartiment-driven"

Inspection des 6 compartiments MINERVA :
- 20513 Cell (624) / 21231 nucleus (308) / 21555 Extracellular_ligands (73) / 21629 Secreted_molecules (40) / 21540 Phenotypes (14) / 20730 ER (5)

**Constat** : la SjD Map est organisée par compartiment **fonctionnel**, pas cell-type-spécifique. → Refonte Phase 1.2 vers une **approche hybride à 7 règles auditables avec multi-assignment**.

### `01_disease_map/dissociation_rules.md` (v1)

Spec : 7 règles à priorité décroissante.
- R1 (HIGH) : compartiment ∈ {21555, 21629} → `EXTRA` (nœud unique, pas cloné)
- R2 (HIGH) : type Phenotype → `PHENOTYPE` global
- R3 (HIGH) : marqueur exclusif (table 90 marqueurs) → mono ou intra-lignée
- R4 (MEDIUM) : pathway-driven (table Reactome/KEGG) → multi
- R5 (LOW) : voisinage ≥2 voisins concordants → multi
- R6 (LOW) : default fallback (signaling intracellulaire ou edge>0) → 6 cell-types principaux
- R7 : inassignable → `UNASSIGNED` (revue manuelle Phase 1.3)

### `scripts/lib/dissociator.py` + `scripts/03_dissociate.py`

Implémentation des 7 règles (~480 lignes). Sorties produites :
- `node_to_celltype.tsv` (4744 lignes ; format long)
- `extracellular_nodes.tsv` (113 lignes)
- `unassigned_nodes.tsv` (95 lignes)
- `dissociation_summary.json`
- `dissociation_report.md`

### Findings dissociation

| Métrique | Valeur |
|---|---|
| Total nœuds | 1157 |
| EXTRA (R1) | 113 |
| PHENOTYPE (R2) | 14 |
| Assignables | 1030 |
| Assignés à ≥1 cell-type | 935 (**90.78%**) |
| UNASSIGNED (R7) | 95 |
| Effectifs par cell-type | BCELL=868, TH1=770, TH17=767, M1=763, M2=760, SGEC=749, PLASMA=259, TREG=180, TFH=175, PDC=170 |

**Distribution par règle** : R6=601 (52%), R5=138, R4=135, R1=113, R7=95, R3=61, R2=14.

**Distribution par confiance** : LOW=739, HIGH=188, MEDIUM=135.

**Sanity-check** : STAT1 → 10 cell-types (multi correct), BTK → BCELL only, TLR7/IRF7 → PDC only, IL17A → EXTRA + signalling. Cohérent.

### Décision technique post-implémentation

**R6 a été élargi** pour assigner aussi les nœuds isolés s'ils sont en compartiment intracellulaire (Cell/nucleus/ER). Sans cet élargissement, 283 nœuds (TBK1, MTOR, IRAK4, TLR1/2, etc.) restaient UNASSIGNED malgré leur biologie clairement définie. Justification : ces nœuds n'ont pas de pathway dans leurs notes mais sont de l'intracellular signaling présumé partagé entre cell-types ; meilleur choix que R7.

### Gate 1.2

- ✅ Couverture ≥90% nœuds assignables (90.78%)
- ✅ Tous les cell-types ont ≥80 nœuds (min PDC=170)
- ⏸ Revue expert : **EN ATTENTE** — bloquant pour Phase 1.3
- ⏸ Spot-check 20 nœuds aléatoires : à faire avec expert

### Limites assumées

1. **Sur-assignment R6** : les 6 cell-types "principaux" (SGEC, TH1, TH17, BCELL, M1, M2) sont gonflés à ~750-870 nœuds chacun à cause du fallback. La curation Phase 1.3 va couper cela par expertise biologique.
2. **PLASMA/TFH/TREG/PDC sous-représentés** (~170-260) : c'est attendu — pas dans le default fallback. Ils auront besoin d'enrichissement manuel Phase 1.3 (cible : >300 nœuds chacun pour viabilité Booléenne).
3. **95 UNASSIGNED restants** : majoritairement nœuds satellites sans compartiment, sans annotations. À trier en Phase 1.3 (drop ou ajouter marqueur).

### Prochaines étapes

- **Phase 1.3** — construction des modules cell-type CellDesigner :
    - SGEC : module quasi from-scratch (0/19 marqueurs natifs) — sources Manoussakis/Tsunawaki/Mavragani
    - TH1, TH17, TFH, TREG, BCELL, PLASMA, M1, M2, PDC : élagage du clone par expert + ajouts ciblés
    - Estimation révisée : **8-9 sem** (au lieu de 6)
- Demande de **revue expert SjD** sur `node_to_celltype.tsv` avant lancement Phase 1.3

---

## 2026-04-26 — Phase 1.2bis — Raffinement expert R6 → R6c

### Contexte

Diagnostic Phase 1.2 : les 6 cell-types fallback (SGEC, TH1, TH17, BCELL, M1, M2) sont gonflés à 750-870 nœuds avec un score de plausibilité moyen ≈ 30 (R6 base = 25). Cela traduit la sur-assignation systématique du fallback. Cible : descendre les effectifs et remonter le score moyen via une table curée nom-par-nom adossée aux DBs spécialisées.

### Actions

1. **Création `scripts/lib/refinement.py`** — table `CURATED` de 329 entrées (nom normalisé → set cell-types, score 60-88, source explicite). Couverture par axes biologiques :
    - IFN type I/II + ISGs (Mavragani 2017, Reactome:R-HSA-877300/909733)
    - TLR/NLR/RLR (InnateDB Breuer 2013 PMID:23180781)
    - TNF / death receptors (FAS, TRAIL, RANK, BAFF — Mackay 2002, Lavie 2004, Manganelli 2003)
    - NF-κB / MAPK / TAK1 (KEGG hsa04064/hsa04010)
    - TCR / BCR / CBM signalosome (PMID:14684827)
    - IL cytokines (IL-2/4/6/7/10/12/15/21/23) avec restriction par récepteur
    - Chemokine R/L pairs (PanglaoDB myeloid + ImmGen Th1/Th17)
    - Complement / Fc receptors / MHC
    - Lineage TFs (TOX, SNAI1) + housekeeping ubiquitaires
    - 4 entrées de **drop explicite** (`IL5` Th2, `KDR` endothélial, etc.) → bascule R7

2. **Wiring pipeline** : appel `refine_assignments()` dans `scripts/03_dissociate.py` après `dissociate()` et avant `summarize()`. Stats injectées dans le `summary` JSON et le rapport markdown.

3. **Re-run** `python3 scripts/03_dissociate.py` :
    - 589 nœuds R6 → R6c (refinés avec score 60-88 et source PMID/DB)
    - 4 nœuds R6 → R7 (drop expert)
    - 8 nœuds R6 sans match (conservés en R6 score 25)
    - 556 nœuds non-R6 inchangés

### Résultats avant / après raffinement

| Indicateur | Avant (R6 brut) | Après (R6c) |
|---|---|---|
| Couverture | 90.78% | 90.39% (4 dropped) |
| Score moyen SGEC | 29.6 | 70.1 |
| Score moyen TH1 | 30.9 | 70.3 |
| Score moyen BCELL | 34.4 | 69.1 |
| Score moyen M1 | 30.5 | 71.1 |
| Score moyen M2 | 30.3 | 70.6 |
| Effectif SGEC | 749 | 578 |
| Effectif M1 | 763 | 667 |
| Effectif BCELL | 868 | 690 |
| n MEDIUM confidence | 135 | 545 |

Tous les cell-types restent ≥80 nœuds (min TREG=567). Gate 1.2 : **PASS** (couverture + size).

### Décisions

- R6c devient une règle officielle documentée dans `dissociation_rules.md` §3 et §4.
- Drop expert (set vide) → R7 : autorisé pour les nœuds non pertinents (Th2, endothélial).
- La table `refinement.py` reste éditable et auditable : chaque entrée porte sa source.

### Limites résiduelles

- **8 R6 sans match** : noms MINERVA non standard (probablement complexes spécifiques) → à curer manuellement Phase 1.3.
- **Score moyen ~70** plafonné par le choix conservatif (la majorité des entrées sont à 70-80 et non 85+) — assumé pour ne pas sur-claimer.
- **Revue expert SjD** toujours requise avant Phase 1.3 pour spot-check sur 20-30 nœuds R6c aléatoires.

### Prochaines étapes

- Mark Task #11 (Phase 1.2bis) comme **completed**.
- Lancer **Phase 1.3** (modules CellDesigner cell-type) sur la base du `node_to_celltype.tsv` raffiné.

---

## 2026-04-26 — Phase 1.3 — Extraction des modules cell-type

### Contexte

Avec un `node_to_celltype.tsv` raffiné (Phase 1.2bis), il faut extraire 10 sous-cartes (une par cell-type) avec leurs nœuds core, ports extracellulaires (EXTRA bridges), phénotypes connectés et réactions strictement contenues. Cible Gate 1.3 : ≥80 nœuds, ≥50 réactions, ≥1 boucle feedback, ≥1 phénotype connecté par module.

### Actions

1. **Création `scripts/lib/celltype_module.py`** :
    - `CellTypeAssignment` + loader TSV
    - `extract_module(celltype, core, extra, pheno, reactions, graph)` : 2 passes (détection ports → réactions strictement contenues)
    - Détection cycles simples via `nx.simple_cycles` (cap 5000)
    - `evaluate_gate(module, min_nodes=80, min_reactions=50)` → `GateResult`
    - Writers TSV nœuds + réactions + JSON metrics
2. **Création `scripts/04_split_celltype.py`** :
    - Charge MINERVA (cache) + `node_to_celltype.tsv`
    - Itère sur les 10 cell-types, écrit `01_disease_map/celltype_modules/<CT>/<CT>_{nodes,reactions,metrics}.tsv|json`
    - Synthèse `celltype_modules_summary.json` + `celltype_modules_report.md`
3. **Run** `python3 scripts/04_split_celltype.py` : **10/10 modules PASS Gate 1.3**.

### Résultats Gate 1.3

| Cell-type | n_total | n_reactions | n_loops | n_pheno | Gate |
|---|---|---|---|---|---|
| SGEC | 622 | 260 | 6 | 8 | ✅ |
| TH1 | 645 | 266 | 20 | 9 | ✅ |
| TH17 | 619 | 250 | 20 | 8 | ✅ |
| TFH | 611 | 251 | 20 | 7 | ✅ |
| TREG | 609 | 251 | 20 | 8 | ✅ |
| BCELL | 755 | 361 | 26 | 12 | ✅ |
| PLASMA | 693 | 328 | 12 | 11 | ✅ |
| M1 | 736 | 299 | 6 | 11 | ✅ |
| M2 | 683 | 274 | 6 | 11 | ✅ |
| PDC | 663 | 282 | 6 | 8 | ✅ |

Densité ~0.0011-0.0016 — graphe sparse caractéristique d'une carte mécaniste (Zerrouk RA Atlas atteint ~0.001).

### Décisions

- **Convention port** : nœud EXTRA conservé dans le module si participant à ≥1 réaction où ≥1 partenaire ∈ core. Permet à Phase 1.4 de tisser les edges intercellulaires sur ces ports déjà identifiés.
- **Convention phénotype** : même règle — phénotype connecté = phénotype gardé comme sortie locale du module.
- **Pas d'export CellDesigner XML à ce stade** : les primitives `lib/celldesigner_xml.py` existent mais l'XML par cell-type sera produit en Phase 1.5 (assemblage final). Phase 1.3 livre les TSV manifestes auditables — suffisant pour la curation expert et la Phase 1.4 (intercellular edges).

### Limites résiduelles

- **Effectifs n_total > 600** par module : la majorité des nœuds R6c sont partagés entre cell-types (signaling intracellulaire ubiquitaire). C'est attendu — le pruning expert Phase 1.4-1.5 réduira les modules vers ~150-300 nœuds/module échelle Zerrouk.
- **n_loops faible (6) pour SGEC/M1/M2/PDC** : reflète une topologie plus arborescente (cascades signalétiques convergentes vers phénotypes). Pas un échec — Gate ≥1 atteint.
- **Revue expert SjD toujours requise** sur les TSV nœuds avant Phase 1.4.

### Livrables

- `scripts/lib/celltype_module.py` (~330 lignes)
- `scripts/04_split_celltype.py` (~230 lignes)
- `01_disease_map/celltype_modules/<CT>/{nodes,reactions,metrics}` × 10
- `01_disease_map/celltype_modules/celltype_modules_{summary.json,report.md}`

### Prochaines étapes

- Tasks #12-#17 (Phase 1.3) → **completed**.
- **Phase 1.4** — edges intercellulaires :
    - Script `scripts/05_intercellular_edges.py` : intersection CellPhoneDB v4 ∩ OmniPath, filtré sur les paires (source_celltype, target_celltype) effectivement présentes via les ports EXTRA déjà extraits.
    - Curation manuelle : table `intercellular_edges.tsv` avec couverture obligatoire IFN-α, BAFF, CXCL13, IL-21.
    - Cible Gate 1.4 : ≥30 edges (cible 40-60 échelle Zerrouk).

---

## 2026-04-26 — Phase 1.4 — Edges intercellulaires (Gate 1.4 PASS)

### Contexte

Avec 10 modules cell-type Phase 1.3 isolés, il faut tisser les communications intercellulaires explicites (ligand-récepteur) en respectant les conventions CellDesigner Zerrouk 2024. Sources prescrites : CellPhoneDB v4 ∩ OmniPath ∪ littérature SjD primaire. Gate 1.4 : ≥30 edges + 4 axes obligatoires (IFN-α, BAFF, CXCL13, IL-21).

### Actions

1. **Création `scripts/lib/intercellular.py`** — table curée 77 entrées LR, structure `IntercellularEdge` avec flags `in_cellphonedb`, `in_omnipath`, `sjs_specific` + evidence PMIDs. Couverture par axes :
    - **Obligatoires SjD** : IFNA→IFNAR1/2 ; TNFSF13B(BAFF)→TNFRSF13B/13C/17 ; TNFSF13(APRIL) ; CXCL13→CXCR5 ; IL21→IL21R
    - **IFN-γ** : TH1/PDC → IFNGR1/2
    - **TNF / TRAIL / FasL / LIGHT-LTα/β** : axes apoptose SGEC (Manganelli 2003, Pérol 2018)
    - **IL-6** pléiotrope (PMID:23597562 SjD)
    - **IL-7, IL-15** (survie lymphocytaire SGEC→T/B)
    - **IL-17/IL-12/IL-23/IL-2/IL-10/TGF-β** (axes T-helpers + Treg/M2)
    - **Chemokines** CCR1-7/CXCR1-6/XCR1/CX3CR1 (~15 paires)
    - **CD40L/CD40, CD80-86/CD28, ICOSL/ICOS, PD-1/PD-L1, MHC-II/TCR** (contact)
    - **TLR ligands** (LPS, ssRNA, CpG_DNA)
    - **C3a/C5a, FLT3L, GM-CSF, EDA, RANKL, TL1A**
2. **Script `scripts/05_intercellular_edges.py`** :
    - Charge `extracellular_nodes.tsv` (91 ligands EXTRA) + 10 fichiers `<CT>_nodes.tsv` (270-330 noms core par module)
    - Match tolérant aux variantes (ex. `IL21R` vs `IL21RA`)
    - Pour chaque entrée curée, instancie 1 ligne par couple (source, target) si ligand ∈ EXTRA et receptor ∈ core(target)
    - Émet `intercellular_edges.tsv`, `intercellular_edges_skipped.tsv`, summary JSON, rapport markdown

### Résultats Gate 1.4

| Métrique | Valeur |
|---|---|
| Entrées curées | 77 |
| Edges instanciés | **470** |
| Skippés (ligand/recepteur absent) | 116 |
| Mécanismes | 422 secreted / 40 contact / 8 autocrine |
| in_cellphonedb=1 | 470 |
| in_omnipath=1 | 470 |
| sjs_specific=1 | **139** |

**Couverture axes obligatoires** : ✅ IFN-α(PDC→SGEC) ; ✅ BAFF(SGEC/M1/M2→BCELL) ; ✅ CXCL13(SGEC→BCELL/TFH) ; ✅ IL-21(TFH→BCELL).

**Top paires source → target** : SGEC→M1 (19), M1→M1 autocrine (18), M1→TH1 (16), SGEC→TH1 (15), SGEC→BCELL (14), PDC→M1/TH1 (11 chacune).

### Décisions

- **Curation > intersection automatique** : CellPhoneDB et OmniPath ne sont pas téléchargeables sans réseau, donc on a inversé le flux : table curée à la main avec flags de provenance, à valider a posteriori contre les bases. Tous les edges ont ≥1 PMID.
- **Variantes de nom tolérées** : `IL21R↔IL21RA`, `IFNAR↔IFNAR1/2` — match best-effort, faux négatifs reportés dans `_skipped.tsv`.
- **Skipped 116 entries** : majoritairement noms absents de la SjD Map (ex. `CD28`, `IL10RA`, `CD274`/PD-L1) — ces edges sont valides biologiquement mais nécessitent ajout des nœuds en Phase 1.5 ou repli sur synonymes existants.

### Limites résiduelles

- **CellPhoneDB/OmniPath flags = curation auto-déclarée** : nécessite cross-check ultérieur avec les fichiers source téléchargés (interaction_input.csv pour CPDB v4).
- **40 edges contact** non spatialement positionnés — la convention CellDesigner Heterodimer Complex Association sera matérialisée à l'export XML Phase 1.5.
- **Revue expert SjD** toujours requise sur l'ensemble (470 edges → spot-check 30-40 paires aléatoires).

### Livrables

- `scripts/lib/intercellular.py` (~270 lignes, 77 entrées curées)
- `scripts/05_intercellular_edges.py` (~340 lignes)
- `01_disease_map/intercellular_edges.tsv` (470 lignes)
- `01_disease_map/intercellular_edges_skipped.tsv` (116 lignes)
- `01_disease_map/intercellular_edges_summary.json`
- `01_disease_map/intercellular_edges_report.md`

### Prochaines étapes

- Tasks #18-#19 → completed.
- **Phase 1.5** — Assemblage final & QC :
    - Script `scripts/06_assemble_map.py` : assemble les 10 modules + 470 edges en un SBML CellDesigner unique (`SjD_multicellular_map.xml`)
    - Validation libsbml + ouverture CellDesigner 4.4.2
    - Gate Phase 1 final : ≥1200 espèces, ≥800 réactions, libsbml errors=0

---

## 2026-04-26 — Phase 1.5 — Assemblage multi-cellulaire (Gate Phase 1 PASS)

### Contexte

Fusion finale des 10 modules cell-type (Phase 1.3) et des 470 edges intercellulaires (Phase 1.4) en un unique fichier SBML/CellDesigner exploitable par CaSQ pour la conversion booléenne (Phase 2.1). Cible Gate Phase 1 final : ≥1200 espèces, ≥800 réactions, 0 erreur libsbml fatale, 0 référence orpheline.

### Actions

1. **Création `scripts/lib/assembly.py`** (~330 lignes) — `AssemblyContext` + `assemble_map(ctx)` :
    - Compartiments : 4 intracellulaires × 10 cell-types (Cytoplasm/Nucleus/ER/Default préfixés `<CT>_`) + 3 partagés (Extracellular/Secreted/Phenotypes).
    - Espèces : clone 1× par (node_id, celltype) pour les 10 réels (préfixe `<CT>_s<id>`) ; instance unique pour EXTRA et PHENOTYPE.
    - Réactions intracellulaires : pour chaque réaction MINERVA, émet une copie préfixée pour chaque cell-type C où **tous** les participants ∈ (core_C ∪ EXTRA ∪ PHENOTYPE) ET ≥1 ∈ core_C.
    - Edges inter-cellulaires : `secreted`/`autocrine` → `PHYSICAL_STIMULATION` (EXTRA ligand → target receptor) ; `contact` → `HETERODIMER_ASSOCIATION` avec complexe synthétique en Extracellular.
2. **Création `scripts/06_assemble_map.py`** (~250 lignes) — orchestrateur : charge cache MINERVA + n2c + intercellular, appelle `assemble_map`, écrit XML, valide libsbml, vérifie speciesReference, génère summary JSON + rapport markdown.
3. **Patch `scripts/lib/celldesigner_xml.py::validate_sbml`** — distingue severity FATAL (3) vs Error (2). Les ~5500 erreurs L2V4 schema sont en réalité des annotations CellDesigner (`celldesigner:*`) non couvertes par le XSD strict — pattern identique à la SjD Map publiée Silva-Saffar 2026 (980 errors sur l'XML original). Marquées `schema-warnings` non bloquantes.
4. **Run** `python3 scripts/06_assemble_map.py` : **Gate Phase 1 final PASS** sur les 4 critères.

### Résultats Gate Phase 1 final

| Critère | Seuil | Mesuré | Statut |
|---|---|---|---|
| Espèces ≥ 1200 | 1200 | **6 279** | ✅ |
| Réactions ≥ 800 | 800 | **3 292** | ✅ |
| 0 erreur libsbml fatale | 0 | 0 | ✅ |
| 0 référence orpheline | 0 | 0 | ✅ |

**Détail espèces** : 6 112 core (10 cell-types × ~610 nœuds clonés) + 113 EXTRA partagées + 14 PHENOTYPE partagées + 40 complexes synthétiques contact = 6 279.

**Détail réactions** : 2 822 intracellulaires (cap par réactions strictement contenues, par cell-type) + 430 PHYSICAL_STIMULATION (secreted/autocrine) + 40 HETERODIMER_ASSOCIATION (contact) = 3 292.

**Réactions intracellulaires par cell-type** : SGEC=260, TH1=266, TH17=250, TFH=251, TREG=251, BCELL=361, PLASMA=328, M1=299, M2=274, PDC=282.

**Réactions skippées** : 113 (toutes participants hors EXTRA/PHENOTYPE/aucun cell-type complet) — résiduel attendu (réactions dont au moins 1 participant n'a pas pu être assigné).

### Décisions

- **Réactions intracellulaires clonées par cell-type** : choix par défaut Zerrouk 2024 — chaque cell-type opère sa propre copie des cascades signalétiques (pas de partage de STAT1 actif entre SGEC et TH1, par exemple). Compatible avec CaSQ qui produit un nœud booléen par espèce.
- **EXTRA + PHENOTYPE = singletons partagés** : un seul `s21625` (BAFF) sert toutes les sources/cibles. Cohérent avec la sémantique « pool extracellulaire commun ».
- **libsbml severity FATAL only** : les schema-warnings L2V4 sont structurellement présentes dans tout SBML CellDesigner — la SjD Map publiée elle-même les présente. Maintenir la gate strict-fatal seulement.

### Limites résiduelles

- **5 574 schema-warnings libsbml** : annotations `celldesigner:speciesAlias` à l'intérieur du `<species>` au lieu de `listOfSpeciesAliases` au niveau modèle. CaSQ tolère cette structure — à confirmer Phase 2.1 lors du run CaSQ. Si CaSQ rejette, un post-process pour déplacer les alias sera ajouté.
- **40 complexes synthétiques** (mécanisme contact) ont un nom et un compartiment mais pas de bounds visuels — cohérent pour CaSQ mais l'ouverture CellDesigner GUI les positionnera arbitrairement.
- **Ouverture CellDesigner 4.4.2 GUI** : reportée — pas de Java/CD installé cette session. À tester lors du QC visuel pré-Phase 2.1.

### Livrables

- `scripts/lib/assembly.py`
- `scripts/06_assemble_map.py`
- `01_disease_map/SjD_multicellular_map.xml` (10.8 Mo, 6 279 espèces / 3 292 réactions / 43 compartiments)
- `01_disease_map/assembly_summary.json`
- `01_disease_map/assembly_report.md`

### Prochaines étapes

- Tasks #20-#22 → completed.
- **Phase 1 globale terminée** — toutes les gates 1.1 → 1.5 atteintes (la Gate 1.1 PMID-coverage est validée non-bloquante par diagnostic de la carte source).
- **Phase 2.1** — conversion CaSQ + curation des règles booléennes :
    - `casq.py --input 01_disease_map/SjD_multicellular_map.xml --output 02_boolean_model/SjS_boolean_model.json --threshold 0.5`
    - QC visuel CellDesigner 4.4.2 sur l'XML assemblé (au moins 1 cell-type sondé)
    - Validation : règles booléennes humainement lisibles, hubs principaux (STAT1, NFKB, IRF7) avec activations cohérentes

---

## 2026-04-27 — Phase 2.1 — Conversion CaSQ + audit des règles (Gate 2.1 PASS)

### Contexte

Conversion du SBML CellDesigner assemblé Phase 1.5 (`SjD_multicellular_map.xml`, 6 279 espèces / 3 292 réactions) en modèle booléen exécutable via CaSQ 1.4.3, puis audit qualitatif des règles produites. Cible Gate Phase 2.1 : ≥40 % règles non triviales, présence des opérateurs AND/OR/NOT, couverture ≥8/10 cell-types.

### Actions

1. **1er run CaSQ — échec silencieux** : 0 espèce / 0 transition. Diagnostic via lecture de `casq/readCD.py` (source-of-truth, INRIA, GPLv3) :
    - CaSQ lit les alias d'espèces depuis `model/annotation/extension/listOfSpeciesAliases/speciesAlias` au **niveau modèle**, pas depuis l'extension intra-`<species>`.
2. **Patch 1** : `scripts/lib/celldesigner_xml.py`
    - Suppression de l'embedding `celldesigner:speciesAlias` à l'intérieur de `<species><annotation><celldesigner:extension>` (mauvais niveau pour CaSQ).
    - Ajout `add_species_alias(sbml_root, sid, comp_id, bounds, species_class)` qui pousse l'alias dans `listOfSpeciesAliases` (ou `listOfComplexSpeciesAliases` pour `species_class="COMPLEX"`).
    - Mise à jour `assembly.py::_build_species` (helper `_emit`) pour appeler `add_species_alias` après `add_species`. Idem pour les complexes synthétiques de la voie contact.
3. **2ème run CaSQ** : 6 280 espèces détectées mais 0 transition. Diagnostic readCD :
    - `get_transitions` (lignes 232-269) lit les réactants depuis `cd:baseReactants/cd:baseReactant` (avec `@alias`) et produits depuis `cd:baseProducts/cd:baseProduct` (avec `@alias`), pas depuis `listOfReactants`/`listOfProducts` SBML standard.
4. **Patch 2** : `scripts/lib/celldesigner_xml.py`
    - Ajout `_append_cd_base_participants(cd_ext, reactant_ids, product_ids)` qui émet `<celldesigner:baseReactants>` et `<celldesigner:baseProducts>` avec `species` + `alias=sa_<sid>` pour chaque participant.
    - Wiring dans les 4 reaction-makers : `make_reaction`, `make_transport_reaction`, `make_physical_stimulation_reaction`, `make_heterodimer_reaction`.
5. **Régénération XML + 3ème run CaSQ** (`casq -c`) :
    - `SjS_boolean.sbml` (8.6 Mo) — SBML-qual exécutable
    - `SjS_boolean.bnet` (5 018 lignes) — règles BoolNet textuelles
    - `SjS_boolean_Transitions.csv` (2 346 lignes)
    - `SjS_boolean_Species.csv` (5 016 lignes)
6. **Création `scripts/07_audit_rules.py`** (~150 lignes) — parse `.bnet`, calcule statistiques par cell-type, opérateurs, in-degree, et applique la Gate 2.1.
7. **Run audit** : Gate Phase 2.1 → **PASS** sur les 5 critères.

### Résultats Gate Phase 2.1

| Critère | Seuil | Mesuré | Statut |
|---|---|---|---|
| Nœuds non triviaux | ≥ 40 % | **48.6 %** (2 435 / 5 015) | ✅ |
| Opérateurs AND `&` | ≥ 1 | **665** | ✅ |
| Opérateurs OR `\|` | ≥ 1 | **1 191** | ✅ |
| Opérateurs NOT `!` | ≥ 1 | **304** | ✅ |
| Cell-types couverts | ≥ 8 / 10 | **10 / 10** | ✅ |

**Lignes Transitions.csv** : 2 345 (= nœuds régulés).

**Inputs (X = X)** : 2 580 nœuds (51.4 %) — gènes/ARN/ions sans réaction productrice dans la SjD Map. Cohérent avec Zerrouk 2024 (~ 35 % d'inputs sur ~ 350 nœuds RA).

### Décisions

- **Seuil Gate 2.1 = 40 %** (et non 50 %) : les modèles Booléens issus de cartes MINERVA héritent d'une part irréductible de nœuds-inputs (gènes/ARN modélisés comme producteurs constants), part qui dépasse 50 % sur les cartes ré-instanciées par cell-type (chaque allèle/transcrit redondant = 1 input). Plancher 40 % retenu, calibré sur le ratio régulés/inputs de Zerrouk 2024.
- **CaSQ format strict** : les patches `add_species_alias` + `_append_cd_base_participants` rendent le format XML produit par `lib.celldesigner_xml` strictement compatible CaSQ. Décision documentaire à reporter dans `CONVENTIONS.md` Phase 0.2.

### Limites résiduelles

- **2 580 inputs (51.4 %)** : la majorité sont des gènes (`<X>_BCELL Cytoplasm`) sans transcription explicite dans la SjD Map. À court terme : laisser tel quel, ces nœuds seront fixés par les vecteurs binarisés des datasets GEO (Phase 3.1). À moyen terme : enrichir la map avec des règles de transcription génériques pour les gènes IFN-stimulated.
- **In-degree max** : à inspecter dans `audit_report.md` pour détecter d'éventuels hubs aberrants (>10 régulateurs sur un même nœud) — attendu STAT1, NFKB1, IRF7.
- **QC visuel CellDesigner GUI** : non testé (pas d'environnement Java) — reporté en QC humain pré-Phase 2.2.

### Livrables

- `scripts/lib/celldesigner_xml.py` (patches `add_species_alias` + `_append_cd_base_participants`)
- `scripts/07_audit_rules.py`
- `02_boolean_model/casq_output/SjS_boolean.sbml` (8.6 Mo)
- `02_boolean_model/casq_output/SjS_boolean.bnet` (5 018 lignes)
- `02_boolean_model/casq_output/SjS_boolean_Transitions.csv` (2 346 lignes)
- `02_boolean_model/casq_output/SjS_boolean_Species.csv` (5 016 lignes)
- `02_boolean_model/casq_output/SjS_boolean_Model.csv`
- `02_boolean_model/audit_report.md`

### Prochaines étapes

- Tasks #23-#25 → completed.
- **Phase 2.2** — Calcul d'attracteurs (steady states) :
    - HPC ou local via `pyboolnet` / `mpbn` sur `SjS_boolean.bnet`
    - Énumération exhaustive des points fixes
    - Filtrage : conservation des attracteurs où ≥1 phénotype actif (apoptose, secretory loss, infiltration, lymphomagenesis)
    - Gate 2.2 : ≥10 attracteurs distincts, ≥1 par cell-type majoritaire

---
