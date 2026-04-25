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
