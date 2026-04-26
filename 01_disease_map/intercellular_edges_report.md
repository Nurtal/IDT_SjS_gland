# Intercellular Edges Report — Phase 1.4

> Sortie de `scripts/05_intercellular_edges.py` — table curée `scripts/lib/intercellular.py` croisée avec les modules Phase 1.3.

## 1. Synthèse globale

- Entrées curées : **77**
- Edges instanciés : **470**
- Entrées skippées : **116**

## 2. Couverture des axes obligatoires SjD (Gate 1.4)

| Axe | Statut |
|---|---|
| IFN-α (PDC→SGEC) | ✅ |
| BAFF (SGEC/M1/M2→BCELL) | ✅ |
| CXCL13 (SGEC→BCELL/TFH) | ✅ |
| IL-21 (TFH→BCELL) | ✅ |

## 3. Distribution par mécanisme

| Mécanisme | n_edges |
|---|---|
| secreted | 422 |
| contact | 40 |
| autocrine | 8 |

## 4. Distribution par paire (source → target)

| Source | Target | n_edges |
|---|---|---|
| SGEC | M1 | 19 |
| M1 | M1 | 18 |
| M1 | TH1 | 16 |
| SGEC | TH1 | 15 |
| SGEC | BCELL | 14 |
| PDC | M1 | 11 |
| PDC | TH1 | 11 |
| M1 | BCELL | 11 |
| M1 | M2 | 11 |
| SGEC | M2 | 11 |
| SGEC | PDC | 11 |
| SGEC | TFH | 9 |
| TH1 | M1 | 9 |
| TH1 | SGEC | 9 |
| M1 | TH17 | 9 |
| SGEC | SGEC | 9 |
| SGEC | TH17 | 9 |
| SGEC | TREG | 9 |
| TH1 | BCELL | 8 |
| TH1 | M2 | 8 |
| M1 | TFH | 8 |
| M1 | TREG | 8 |
| PDC | BCELL | 6 |
| PDC | M2 | 6 |
| SGEC | PLASMA | 6 |
| M1 | SGEC | 6 |
| M1 | PDC | 6 |
| TH1 | PDC | 6 |
| TH1 | TH1 | 6 |
| M2 | M1 | 6 |

_… et 58 autres paires._

## 5. Gate 1.4

- ≥30 edges intercellulaires : **PASS** (470)
- 4 axes SjD obligatoires couverts : **PASS**
- Edges sjs_specific : 139
- Revue expert SjD : **EN ATTENTE** (manuel)

## 6. Sources / provenance

| Provenance | n_edges |
|---|---|
| in_cellphonedb=1 | 470 |
| in_omnipath=1 | 470 |
| sjs_specific=1 | 139 |

---

_Généré par `scripts/05_intercellular_edges.py`_