"""
08_mpbn_attractors.py — Phase 2.2 mpbn — énumération des trap-spaces minimaux.

Les trap-spaces minimaux d'un réseau booléen async = sur-approximation des
attracteurs asynchrones (cf. Pauleve 2020, Klarner 2015). mpbn calcule cela
via solving ASP, scalable jusqu'à 10⁴+ nœuds.

Le `.bnet` brut produit par CaSQ contient des espaces dans les noms (cause :
suffixe ` Cytoplasm` collé au nom MINERVA). On commence par produire un
`.bnet` normalisé (espaces → `__`) avant de le passer à mpbn.

Sortie :
    02_boolean_model/poc_results/mpbn/
        SjS_boolean_normalized.bnet
        trap_spaces.tsv (1 ligne / trap-space, 1 col / nœud, valeurs 0/1/*)
        summary.json
"""
from __future__ import annotations

import csv
import json
import logging
import re
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
BNET_RAW = ROOT / "02_boolean_model" / "casq_output" / "SjS_boolean.bnet"
SPECIES_CSV = ROOT / "02_boolean_model" / "casq_output" / "SjS_boolean_Species.csv"
OUT = ROOT / "02_boolean_model" / "poc_results" / "mpbn"
OUT.mkdir(parents=True, exist_ok=True)
BNET_NORM = OUT / "SjS_boolean_normalized.bnet"
NAME_MAP = OUT / "name_map.tsv"
TRAPS_TSV = OUT / "trap_spaces.tsv"
SUMMARY = OUT / "summary.json"


def safe(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_]", "_", name.strip())
    s = re.sub(r"_+", "_", s).strip("_")
    if s and s[0].isdigit():
        s = "_" + s
    return s or "_x"


def collect_names(species_csv: Path, bnet_raw: Path) -> dict[str, str]:
    """Construit le mapping `nom canonique` → `nom safe bnet`.

    Les noms sont collectés à la fois depuis le CSV des espèces CaSQ
    (canon) et depuis les targets/littéraux du `.bnet` brut, pour couvrir
    les littéraux apparaissant uniquement à droite (constantes d'entrée).
    """
    raw_names: set[str] = set()

    # Source 1 : Species CSV (champ 0 — CSV quoté car les noms peuvent
    # contenir des virgules, ex. `PI(3,4,5)P3_BCELL Cytoplasm`).
    with species_csv.open() as fh:
        reader = csv.reader(fh)
        next(reader)  # header
        for row in reader:
            if row and row[0].strip():
                raw_names.add(row[0].strip())

    # Source 2 : targets du bnet
    with bnet_raw.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("targets"):
                continue
            if "," not in line:
                continue
            target = line.split(",", 1)[0].strip()
            if target:
                raw_names.add(target)

    # Résolution des collisions de safe-name : suffixe numérique stable.
    seen: dict[str, str] = {}
    used: set[str] = set()
    # Tri pour mapping reproductible
    for raw in sorted(raw_names):
        base = safe(raw)
        candidate = base
        i = 1
        while candidate in used:
            i += 1
            candidate = f"{base}__{i}"
        seen[raw] = candidate
        used.add(candidate)
    return seen


def normalize_bnet(src: Path, dst: Path, name_map: dict[str, str]) -> None:
    """Réécrit le `.bnet` en remplaçant chaque nom canonique par son alias safe.

    Algorithme : longest-match-first, on parcourt la règle de gauche à droite
    et à chaque position on tente de matcher le plus long nom du mapping.
    Les opérateurs `& | ! ( )` qui ne font pas partie d'un nom canonique
    sont conservés tels quels.
    """
    # Liste triée par longueur décroissante pour matching glouton
    sorted_names = sorted(name_map.keys(), key=len, reverse=True)

    def normalize_rule(rule: str) -> str:
        out: list[str] = []
        i = 0
        n = len(rule)
        while i < n:
            ch = rule[i]
            # Opérateurs Booléens passent tels quels
            if ch in "&|!()":
                out.append(ch)
                i += 1
                continue
            if ch.isspace():
                i += 1
                continue
            # Tentative de match d'un nom canonique
            matched = None
            for nm in sorted_names:
                if rule.startswith(nm, i):
                    matched = nm
                    break
            if matched is not None:
                out.append(name_map[matched])
                i += len(matched)
            else:
                # Fallback : prendre une rafale alphanumérique brute
                m = re.match(r"[A-Za-z0-9_\-./,:+]+", rule[i:])
                if m:
                    out.append(safe(m.group(0)))
                    i += len(m.group(0))
                else:
                    i += 1
        return "".join(out)

    def split_target_rule(line: str) -> tuple[str, str] | None:
        """Sépare target et rule en tenant compte des virgules internes aux
        noms (ex. `PI(4,5)P2_BCELL Cytoplasm`). Stratégie : longest-match
        depuis la position 0, puis `, ` séparateur."""
        for nm in sorted_names:
            if line.startswith(nm):
                rest = line[len(nm):]
                rest_stripped = rest.lstrip()
                if rest_stripped.startswith(","):
                    return nm, rest_stripped[1:]
        # Fallback : split sur la première virgule
        if "," in line:
            t, r = line.split(",", 1)
            return t, r
        return None

    out_lines: list[str] = []
    with src.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                out_lines.append(line)
                continue
            if line.startswith("targets"):
                out_lines.append(line)
                continue
            split = split_target_rule(line)
            if split is None:
                out_lines.append(line)
                continue
            target_raw, rule = split
            target_raw = target_raw.strip()
            target = name_map.get(target_raw, safe(target_raw))
            new_rule = normalize_rule(rule).strip()
            out_lines.append(f"{target}, {new_rule}")
    dst.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    logger.info("bnet normalisé écrit : %s", dst)


def main() -> None:
    name_map = collect_names(SPECIES_CSV, BNET_RAW)
    logger.info("Mapping nom→safe collecté : %d entrées", len(name_map))
    with NAME_MAP.open("w", encoding="utf-8") as fh:
        fh.write("raw\tsafe\n")
        for raw, sf in sorted(name_map.items()):
            fh.write(f"{raw}\t{sf}\n")
    normalize_bnet(BNET_RAW, BNET_NORM, name_map)

    import mpbn

    logger.info("Chargement modèle mpbn…")
    t0 = time.time()
    mbn = mpbn.MPBooleanNetwork.load(str(BNET_NORM))
    logger.info("Modèle chargé : %d nœuds (%.1fs)", len(mbn), time.time() - t0)

    # Énumération bornée : sur un modèle 5k+ nœuds, le nombre de trap-spaces
    # peut atteindre plusieurs dizaines de milliers. On échantillonne les
    # `MAX_TRAPS` premiers (clingo les énumère en ordre stable) — suffisant
    # pour caractériser la structure et calibrer la Phase 2.4. L'énumération
    # exhaustive sera reportée à 2.1bis après réduction bioLQM.
    MAX_TRAPS = 5000
    TIME_BUDGET_S = 600
    logger.info(
        "Énumération des trap-spaces minimaux (cap %d, budget %ds)…",
        MAX_TRAPS, TIME_BUDGET_S,
    )
    t0 = time.time()
    trap_spaces: list[dict] = []
    truncated = False
    for ts in mbn.attractors():
        trap_spaces.append(dict(ts))
        elapsed_now = time.time() - t0
        if len(trap_spaces) >= MAX_TRAPS:
            truncated = True
            logger.info("Cap %d atteint — arrêt", MAX_TRAPS)
            break
        if elapsed_now > TIME_BUDGET_S:
            truncated = True
            logger.info("Budget temps %ds dépassé — arrêt", TIME_BUDGET_S)
            break
    elapsed = time.time() - t0
    logger.info("→ %d trap-spaces (%.1fs, truncated=%s)",
                len(trap_spaces), elapsed, truncated)

    if trap_spaces:
        nodes = sorted(trap_spaces[0].keys())
        with TRAPS_TSV.open("w") as fh:
            fh.write("trap_id\t" + "\t".join(nodes) + "\n")
            for i, ts in enumerate(trap_spaces):
                vals = [str(ts.get(n, "*")) for n in nodes]
                fh.write(f"trap_{i}\t" + "\t".join(vals) + "\n")
        logger.info("Trap-spaces écrits : %s", TRAPS_TSV)

    fixed_count = 0
    free_distrib: list[int] = []
    for ts in trap_spaces:
        free = sum(1 for v in ts.values() if v not in (0, 1))
        free_distrib.append(free)
        if free == 0:
            fixed_count += 1

    summary = {
        "n_nodes": len(mbn),
        "n_trap_spaces": len(trap_spaces),
        "n_full_fixed_points": fixed_count,
        "free_dim_min": min(free_distrib) if free_distrib else 0,
        "free_dim_max": max(free_distrib) if free_distrib else 0,
        "free_dim_mean": (sum(free_distrib) / len(free_distrib)) if free_distrib else 0,
        "elapsed_seconds": round(elapsed, 2),
        "truncated": truncated,
        "cap_max_traps": MAX_TRAPS,
        "cap_time_budget_s": TIME_BUDGET_S,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print()
    print("=" * 70)
    print("  PHASE 2.2 — mpbn (trap-spaces minimaux)")
    print("=" * 70)
    print(f"  Nœuds                   : {summary['n_nodes']}")
    print(f"  Trap-spaces minimaux    : {summary['n_trap_spaces']}")
    print(f"  Point fixes complets    : {summary['n_full_fixed_points']}")
    print(f"  Free-dim (min/mean/max) : {summary['free_dim_min']} / "
          f"{summary['free_dim_mean']:.1f} / {summary['free_dim_max']}")
    print(f"  Temps                   : {summary['elapsed_seconds']} s")
    print(f"  Sortie                  : {OUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()
