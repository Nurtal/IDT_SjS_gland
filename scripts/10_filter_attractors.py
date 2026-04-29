"""
10_filter_attractors.py — Phase 2.3 — scénarios + filtrage par phénotype.

Pourquoi pas un simple filtre du `trap_spaces.tsv` Phase 2.2 ?
    Sur 5 015 nœuds dont ~ 2 580 inputs, le cap=5 000 trap-spaces de mpbn
    n'échantillonne qu'**un coin** de l'espace : tous les attracteurs
    énumérés partagent la même configuration d'inputs (Fibrosis=
    Matrix_degradation=Regulated_Necrosis=ON, autres OFF). Pour discriminer
    `disease` vs `healthy_like`, on doit **fixer les inputs** correspondant
    aux scénarios biologiques attendus avant de re-énumérer.

Stratégie scénarios :
    SCENARIO `healthy`     : tous les triggers SjD à 0 (pas d'inputs viraux,
                             pas d'auto-Ag, pas de cytokines pro-inflam).
    SCENARIO `disease_ifn` : dsRNA, CpG_DNA, ssRNA → IFN-α (axe IFN type I,
                             Mavragani 2017) — endotype IFN.
    SCENARIO `disease_ag`  : auto-antigènes + cytokines (Ag/IgG, IFNG, TNF,
                             IL6) — endotype inflammatoire chronique.
    SCENARIO `disease_full`: tous triggers ON (worst case).

Pour chaque scénario : run mpbn capé à 100 trap-spaces (suffit pour les
diversifier sur les phénotypes downstream), annote, classe.

Sortie :
    02_boolean_model/attractors_<scenario>.tsv (4 fichiers)
    02_boolean_model/attractors_disease.tsv  (concat scénarios disease_*)
    02_boolean_model/attractors_healthy.tsv  (= scénario healthy)
    02_boolean_model/phenotype_signature_per_attractor.tsv (tous scénarios)
    02_boolean_model/attractor_summary.json
"""
from __future__ import annotations

import csv
import json
import logging
import time
from collections import Counter
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
BNET_NORM = ROOT / "02_boolean_model" / "poc_results" / "mpbn" / "SjS_boolean_normalized.bnet"
OUT = ROOT / "02_boolean_model"
DISEASE_TSV = OUT / "attractors_disease.tsv"
HEALTHY_TSV = OUT / "attractors_healthy.tsv"
PHENO_TSV = OUT / "phenotype_signature_per_attractor.tsv"
SUMMARY_JSON = OUT / "attractor_summary.json"

# Inputs SjD-pertinents (existant dans le bnet normalisé)
# `healthy` = état de référence : aucun trigger ET les 9 disease phenotypes
# forcés OFF. Sans cette contrainte, les self-loops du bnet (ex. RIPK3 = RIPK3)
# laissent mpbn libre de fixer les effecteurs ON, ce qui pollue tous les
# attracteurs avec Regulated_Necrosis/Inflammation systématiquement ON.
# Forcer les phénotypes OFF dans healthy revient à demander : « quelle config
# upstream est compatible avec une glande sans aucun signe pathologique ».
_TRIGGERS_OFF = {
    "dsRNA_polyIC": 0,
    "CpG_DNA": 0,
    "Allergen": 0,
    "TNF_Extracellular": 0,
    "IFNG_Extracellular": 0,
    "IFNB1_Extracellular": 0,
    "IL6_Extracellular": 0,
}
_DISEASE_PHENOS_OFF = {
    "Apoptosis_phenotype": 0,
    "Inflammation_phenotype": 0,
    "Fibrosis_phenotype": 0,
    "B_Cell_Activation_Survival_phenotype": 0,
    "T_Cell_Activation_Differentiation_phenotype": 0,
    "Chemotaxis_Infiltration_phenotype": 0,
    "Lymphoid_organ_development_phenotype": 0,
    "Matrix_degradation_phenotype": 0,
    "Regulated_Necrosis_phenotype": 0,
}
SCENARIOS: dict[str, dict[str, int]] = {
    "healthy": {**_TRIGGERS_OFF, **_DISEASE_PHENOS_OFF},
    "disease_ifn": {
        "dsRNA_polyIC": 1,
        "CpG_DNA": 1,
        "IFNB1_Extracellular": 1,
        "Allergen": 0,
        "TNF_Extracellular": 0,
    },
    "disease_ag": {
        "Allergen": 1,
        "TNF_Extracellular": 1,
        "IFNG_Extracellular": 1,
        "IL6_Extracellular": 1,
        "dsRNA_polyIC": 0,
    },
    "disease_full": {
        "dsRNA_polyIC": 1,
        "CpG_DNA": 1,
        "Allergen": 1,
        "TNF_Extracellular": 1,
        "IFNG_Extracellular": 1,
        "IFNB1_Extracellular": 1,
        "IL6_Extracellular": 1,
    },
}

CAP_PER_SCENARIO = 100
TIME_BUDGET_S = 300

# Phénotypes considérés comme "disease" (suit Silva-Saffar 2026 + Zerrouk 2024)
DISEASE_PHENOTYPES = (
    "Apoptosis_phenotype",
    "Inflammation_phenotype",
    "Fibrosis_phenotype",
    "B_Cell_Activation_Survival_phenotype",
    "T_Cell_Activation_Differentiation_phenotype",
    "Chemotaxis_Infiltration_phenotype",
    "Lymphoid_organ_development_phenotype",
    "Matrix_degradation_phenotype",
    "Regulated_Necrosis_phenotype",
)

# Phénotypes "homéostatiques" — informationnels uniquement, pas utilisés dans la classification
HOMEOSTATIC_PHENOTYPES = (
    "Angiogenesis_phenotype",
    "Cell_Proliferation_Survival_phenotype",
    "MHC_Class_1_Activation_phenotype",
    "MHC_Class_2_Activation_phenotype",
    "Phagocytosis_phenotype",
)

ALL_PHENOTYPES = DISEASE_PHENOTYPES + HOMEOSTATIC_PHENOTYPES

CELL_TYPES = ("SGEC", "TH1", "TH17", "TFH", "TREG", "BCELL", "PLASMA", "M1", "M2", "PDC")


def cell_type_of(node: str) -> str | None:
    for ct in CELL_TYPES:
        # Préfixes safe-id : `<gène>_<CT>_<comp>` ou `<gène>_<CT>`
        if f"_{ct}_" in node or node.endswith(f"_{ct}"):
            return ct
    return None


def write_tsv(path: Path, records: list[dict]) -> None:
    if not records:
        path.write_text("", encoding="utf-8")
        return
    cols = list(records[0].keys())
    with path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            fh.write("\t".join(
                str(r[c]) if r[c] is not None else "*" for c in cols
            ) + "\n")


def enumerate_scenario(mbn_load_path: Path, name: str, fixings: dict[str, int],
                       cap: int, budget_s: int) -> list[dict]:
    """Charge le bnet, fixe les inputs spécifiés, énumère les trap-spaces."""
    import mpbn

    mbn = mpbn.MPBooleanNetwork.load(str(mbn_load_path))
    applied = []
    skipped = []
    for k, v in fixings.items():
        if k in mbn:
            mbn[k] = int(v)
            applied.append((k, v))
        else:
            skipped.append(k)
    logger.info("[%s] fixings appliqués : %d, skipped : %d (%s)",
                name, len(applied), len(skipped),
                skipped[:5] if skipped else "")

    t0 = time.time()
    traps: list[dict] = []
    for ts in mbn.attractors():
        traps.append(dict(ts))
        if len(traps) >= cap:
            break
        if time.time() - t0 > budget_s:
            logger.warning("[%s] budget %ds dépassé — arrêt", name, budget_s)
            break
    logger.info("[%s] → %d trap-spaces (%.1fs)", name, len(traps), time.time() - t0)
    return traps


def main() -> None:
    # Index nœuds (lecture du bnet pour extraire la liste des nœuds)
    nodes: list[str] = []
    with BNET_NORM.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#") or line.startswith("targets"):
                continue
            if "," in line:
                nodes.append(line.split(",", 1)[0].strip())
    logger.info("Nœuds bnet : %d", len(nodes))

    missing = [p for p in ALL_PHENOTYPES if p not in nodes]
    if missing:
        raise SystemExit(f"Phénotypes manquants dans bnet : {missing}")

    # Pré-compute cell-type membership pour annotation
    ct_members: dict[str, list[str]] = {ct: [] for ct in CELL_TYPES}
    for n in nodes:
        ct = cell_type_of(n)
        if ct:
            ct_members[ct].append(n)

    pheno_signatures: list[dict] = []
    disease_attractors: list[dict] = []
    healthy_attractors: list[dict] = []
    per_scenario_files: dict[str, Path] = {}

    for scenario_name, fixings in SCENARIOS.items():
        traps = enumerate_scenario(
            BNET_NORM, scenario_name, fixings,
            cap=CAP_PER_SCENARIO, budget_s=TIME_BUDGET_S,
        )
        scenario_records: list[dict] = []
        for i, ts in enumerate(traps):
            trap_id = f"{scenario_name}_t{i}"
            pheno_vec = {p: ts.get(p) for p in ALL_PHENOTYPES}
            pheno_int = {
                p: (int(v) if v in (0, 1) else None) for p, v in pheno_vec.items()
            }
            n_disease_on = sum(1 for p in DISEASE_PHENOTYPES if pheno_int[p] == 1)
            n_disease_off = sum(1 for p in DISEASE_PHENOTYPES if pheno_int[p] == 0)
            n_disease_unk = sum(1 for p in DISEASE_PHENOTYPES if pheno_int[p] is None)
            if n_disease_on >= 1:
                klass = "disease"
            elif n_disease_off == len(DISEASE_PHENOTYPES):
                klass = "healthy_like"
            else:
                klass = "mixed"

            ct_activity = {}
            for ct in CELL_TYPES:
                ct_activity[ct] = sum(1 for n in ct_members[ct] if ts.get(n) == 1)

            record = {
                "trap_id": trap_id,
                "scenario": scenario_name,
                "class": klass,
                "n_disease_on": n_disease_on,
                "n_disease_off": n_disease_off,
                "n_disease_unk": n_disease_unk,
                **{f"pheno_{p}": pheno_int[p] for p in ALL_PHENOTYPES},
                **{f"ct_{ct}_on": ct_activity[ct] for ct in CELL_TYPES},
            }
            scenario_records.append(record)
            pheno_signatures.append(record)
            if klass == "disease":
                disease_attractors.append(record)
            elif klass == "healthy_like":
                healthy_attractors.append(record)

        per_path = OUT / f"attractors_{scenario_name}.tsv"
        write_tsv(per_path, scenario_records)
        per_scenario_files[scenario_name] = per_path
        logger.info("[%s] %d records → %s", scenario_name, len(scenario_records), per_path)

    write_tsv(PHENO_TSV, pheno_signatures)
    write_tsv(DISEASE_TSV, disease_attractors)
    write_tsv(HEALTHY_TSV, healthy_attractors)

    # --- Stats ---
    classes = Counter(r["class"] for r in pheno_signatures)
    by_scenario = Counter(r["scenario"] for r in pheno_signatures)
    pheno_on = {p: sum(1 for r in pheno_signatures if r[f"pheno_{p}"] == 1) for p in ALL_PHENOTYPES}

    summary = {
        "n_attractors_total": len(pheno_signatures),
        "n_disease": classes.get("disease", 0),
        "n_healthy_like": classes.get("healthy_like", 0),
        "n_mixed": classes.get("mixed", 0),
        "by_scenario": dict(by_scenario),
        "phenotype_on_counts": pheno_on,
        "disease_phenotypes_used": list(DISEASE_PHENOTYPES),
        "homeostatic_phenotypes_used": list(HOMEOSTATIC_PHENOTYPES),
        "scenarios_definitions": SCENARIOS,
        "cap_per_scenario": CAP_PER_SCENARIO,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print()
    print("=" * 70)
    print("  PHASE 2.3 — Filtrage par phénotype")
    print("=" * 70)
    print(f"  Total attracteurs        : {summary['n_attractors_total']}")
    print(f"  Disease (≥1 disease ON)  : {summary['n_disease']}")
    print(f"  Healthy-like (tous OFF)  : {summary['n_healthy_like']}")
    print(f"  Mixed                    : {summary['n_mixed']}")
    print("-" * 70)
    print("  Top phénotypes ON :")
    for p, c in sorted(pheno_on.items(), key=lambda x: -x[1])[:8]:
        print(f"    {c:5d} | {p}")
    print("-" * 70)
    gate = (summary["n_disease"] >= 1) and (summary["n_healthy_like"] >= 1)
    print(f"  Gate Phase 2.3 (filtrage) : {'PASS' if gate else 'FAIL'}")
    print(f"  Sorties : {OUT}/attractors_*.tsv + phenotype_signature_per_attractor.tsv")
    print("=" * 70)


if __name__ == "__main__":
    main()
