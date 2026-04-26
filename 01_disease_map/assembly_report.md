# Assembly Report — Phase 1.5

> Sortie de `scripts/06_assemble_map.py` — fusion des 10 modules cell-type + 470 edges inter-cellulaires en un unique CellDesigner SBML.

## 1. Synthèse globale

- Compartiments : **43**
- Espèces       : **6279**
- Réactions     : **3292**

## 2. Détail des espèces

| Catégorie | n |
|---|---|
| core (clonées par cell-type) | 6112 |
| EXTRA (partagées) | 113 |
| PHENOTYPE (partagées) | 14 |
| complexes inter-cellulaires (contact) | 40 |

## 3. Réactions intracellulaires par cell-type

| Cell-type | n_reactions |
|---|---|
| SGEC | 260 |
| TH1 | 266 |
| TH17 | 250 |
| TFH | 251 |
| TREG | 251 |
| BCELL | 361 |
| PLASMA | 328 |
| M1 | 299 |
| M2 | 274 |
| PDC | 282 |
| **Total intra** | **2822** |
| Réactions skippées (participants non résolus) | 113 |

## 4. Edges inter-cellulaires matérialisés

| Mécanisme | n |
|---|---|
| secreted/autocrine (PHYSICAL_STIMULATION) | 430 |
| contact (HETERODIMER_ASSOCIATION) | 40 |
| skipped (ligand/recepteur non résolus) | 0 |

## 5. Validation SBML

- Erreurs FATAL libsbml : **0**
- Schema-warnings (annotations CellDesigner non couvertes par L2V4) : 5574 _(non bloquantes — pattern identique à la SjD Map publiée Silva-Saffar 2026)_
- speciesReference non résolues : **0**

Premiers messages libsbml :

```
[Error] L190:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L218:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L246:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L372:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L452:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L480:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L509:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L538:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L566:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

[Error] L595:C6 — An SBML XML document must conform to the XML Schema for the corresponding SBML Level, Version and Release. The XML Schema for SBML defines the basic SBML object structure, the data types used by those objects, and the order in which the objects may appear in an SBML document.
Reference: L2V4 Section 4.1
 Incorrect ordering of <annotation> and <notes> elements -- <notes> must come before <annotation> due to the way that the XML Schema for SBML is defined.

```

## 6. Gate Phase 1 final

| Critère | Seuil | Mesuré | Statut |
|---|---|---|---|
| Espèces ≥ 1200 | 1200 | 6279 | ✅ |
| Réactions ≥ 800 | 800 | 3292 | ✅ |
| 0 erreur libsbml fatale | 0 | 0 | ✅ |
| 0 référence orpheline | 0 | 0 | ✅ |

**Statut global Phase 1 : PASS**

---

_Généré par `scripts/06_assemble_map.py`_