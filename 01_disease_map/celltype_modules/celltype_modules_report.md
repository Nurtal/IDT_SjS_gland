# Cell-type Modules Report — Phase 1.3

> Sortie de `scripts/04_split_celltype.py`. Extraction des sous-cartes à partir de `node_to_celltype.tsv` (Phase 1.2bis raffiné R6c).

## 1. Synthèse globale

- Cell-types extraits : **10**
- Réactions sources    : **598**
- Nœuds sources        : **1157**
- Modules valides Gate 1.3 : **10/10**

## 2. Métriques par module

| Cell-type | n_core | n_ports | n_pheno | n_total | n_reactions | n_edges | n_loops | density |
|---|---|---|---|---|---|---|---|---|
| SGEC | 578 | 36 | 8 | 622 | 260 | 522 | 6 | 0.0014 |
| TH1 | 596 | 40 | 9 | 645 | 266 | 545 | 20 | 0.0013 |
| TH17 | 577 | 34 | 8 | 619 | 250 | 519 | 20 | 0.0014 |
| TFH | 571 | 33 | 7 | 611 | 251 | 520 | 20 | 0.0014 |
| TREG | 567 | 34 | 8 | 609 | 251 | 527 | 20 | 0.0014 |
| BCELL | 690 | 53 | 12 | 755 | 361 | 819 | 26 | 0.0014 |
| PLASMA | 634 | 48 | 11 | 693 | 328 | 756 | 12 | 0.0016 |
| M1 | 667 | 58 | 11 | 736 | 299 | 615 | 6 | 0.0011 |
| M2 | 625 | 47 | 11 | 683 | 274 | 566 | 6 | 0.0012 |
| PDC | 607 | 48 | 8 | 663 | 282 | 590 | 6 | 0.0013 |

## 3. Gate 1.3

Critères : ≥80 nœuds, ≥50 réactions, ≥1 boucle de feedback, ≥1 phénotype connecté.

| Cell-type | ≥80 nœuds | ≥50 réactions | ≥1 loop | ≥1 phénotype | Statut |
|---|---|---|---|---|---|
| SGEC | ✅ | ✅ | ✅ | ✅ | **PASS** |
| TH1 | ✅ | ✅ | ✅ | ✅ | **PASS** |
| TH17 | ✅ | ✅ | ✅ | ✅ | **PASS** |
| TFH | ✅ | ✅ | ✅ | ✅ | **PASS** |
| TREG | ✅ | ✅ | ✅ | ✅ | **PASS** |
| BCELL | ✅ | ✅ | ✅ | ✅ | **PASS** |
| PLASMA | ✅ | ✅ | ✅ | ✅ | **PASS** |
| M1 | ✅ | ✅ | ✅ | ✅ | **PASS** |
| M2 | ✅ | ✅ | ✅ | ✅ | **PASS** |
| PDC | ✅ | ✅ | ✅ | ✅ | **PASS** |

## 4. Recommandations

- ✅ Tous les modules passent Gate 1.3 — passage Phase 1.4 autorisé.

- Revue expert SjD signée requise sur `<celltype>_nodes.tsv` avant Phase 1.4.

---

_Généré par `scripts/04_split_celltype.py`_