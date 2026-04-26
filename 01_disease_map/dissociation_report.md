# Dissociation Report — SjD Map → cell-types

> **Phase 1.2 ROADMAP** — sortie de `scripts/03_dissociate.py`.
> Règles définies dans `dissociation_rules.md`.

## 1. Couverture globale

- Nœuds totaux : **1157**
- Extracellulaires (R1) : **113**
- Phénotypes (R2) : **14**
- Inassignables (R7) : **99**
- Assignables (hors EXTRA/PHENOTYPE) : **1030**
- Assignés à ≥1 cell-type : **931** (90.39%)

## 2. Distribution par règle

| Règle | n_nœuds |
|---|---|
| R1 | 113 |
| R2 | 14 |
| R3 | 61 |
| R4 | 135 |
| R5 | 138 |
| R6 | 8 |
| R6c | 589 |
| R7 | 99 |

## 3. Distribution par confiance

| Confidence | n_nœuds |
|---|---|
| MEDIUM | 545 |
| LOW | 325 |
| HIGH | 188 |
| (R7) | 99 |

## 4. Effectifs et plausibilité par cell-type

| Cell-type | n_nœuds | seuil 80 | score moyen | n ≥80 | n <40 |
|---|---|---|---|---|---|
| SGEC | 578 | ✅ | 70.1 | 35 | 8 |
| TH1 | 596 | ✅ | 70.3 | 59 | 8 |
| TH17 | 577 | ✅ | 70.1 | 54 | 8 |
| TFH | 571 | ✅ | 70.8 | 54 | 0 |
| TREG | 567 | ✅ | 70.5 | 48 | 0 |
| BCELL | 690 | ✅ | 69.1 | 88 | 8 |
| PLASMA | 634 | ✅ | 68.8 | 62 | 0 |
| M1 | 667 | ✅ | 71.1 | 78 | 8 |
| M2 | 625 | ✅ | 70.6 | 58 | 8 |
| PDC | 607 | ✅ | 71.0 | 59 | 0 |
| EXTRA | 113 | — | 95.0 | 113 | 0 |
| PHENOTYPE | 14 | — | 90.0 | 14 | 0 |

## 5. Gate 1.2

- Couverture ≥90% (assignables) : **PASS** (90.39%)
- Tous les cell-types ont ≥80 nœuds : **PASS**
- Revue expert : **EN ATTENTE** (manuel)

## 6. Raffinement expert R6 → R6c (Phase 1.2bis)

Application post-dissociation de la table curée `scripts/lib/refinement.py` (PanglaoDB, CellMarker 2.0, HPA, InnateDB, Reactome, KEGG, ImmGen + littérature SjD primaire).

| Indicateur | Valeur |
|---|---|
| Entrées curées | 329 |
| Nœuds R6 → R6c (raffinés) | 589 |
| Nœuds R6 → R7 (drop expert) | 4 |
| Nœuds R6 sans match (conservés R6) | 8 |

## 7. Recommandations Phase 1.3

- Réviser les **99 nœuds inassignables** (`unassigned_nodes.tsv`) — soit ajouter un marqueur, soit retirer.

---

_Généré par `scripts/03_dissociate.py`_