# Dissociation Report — SjD Map → cell-types

> **Phase 1.2 ROADMAP** — sortie de `scripts/03_dissociate.py`.
> Règles définies dans `dissociation_rules.md`.

## 1. Couverture globale

- Nœuds totaux : **1157**
- Extracellulaires (R1) : **113**
- Phénotypes (R2) : **14**
- Inassignables (R7) : **95**
- Assignables (hors EXTRA/PHENOTYPE) : **1030**
- Assignés à ≥1 cell-type : **935** (90.78%)

## 2. Distribution par règle

| Règle | n_nœuds |
|---|---|
| R1 | 113 |
| R2 | 14 |
| R3 | 61 |
| R4 | 135 |
| R5 | 138 |
| R6 | 601 |
| R7 | 95 |

## 3. Distribution par confiance

| Confidence | n_nœuds |
|---|---|
| LOW | 739 |
| HIGH | 188 |
| MEDIUM | 135 |
| (R7) | 95 |

## 4. Effectifs par cell-type (clones potentiels)

| Cell-type | n_nœuds | seuil 80 |
|---|---|---|
| SGEC | 749 | ✅ |
| TH1 | 770 | ✅ |
| TH17 | 767 | ✅ |
| TFH | 175 | ✅ |
| TREG | 180 | ✅ |
| BCELL | 868 | ✅ |
| PLASMA | 259 | ✅ |
| M1 | 763 | ✅ |
| M2 | 760 | ✅ |
| PDC | 170 | ✅ |
| EXTRA | 113 | — |
| PHENOTYPE | 14 | — |

## 5. Gate 1.2

- Couverture ≥90% (assignables) : **PASS** (90.78%)
- Tous les cell-types ont ≥80 nœuds : **PASS**
- Revue expert : **EN ATTENTE** (manuel)

## 6. Recommandations Phase 1.3

- Réviser les **95 nœuds inassignables** (`unassigned_nodes.tsv`) — soit ajouter un marqueur, soit retirer.

---

_Généré par `scripts/03_dissociate.py`_