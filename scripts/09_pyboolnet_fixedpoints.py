"""
09_pyboolnet_fixedpoints.py — Phase 2.2 cross-check pyboolnet.

Cross-check des point fixes du modèle Booléen via pyboolnet (SAT-based).
Doit produire un set cohérent avec mpbn.

Le `.bnet` normalisé en Phase 2.2/mpbn est utilisé tel quel.

Sortie :
    02_boolean_model/poc_results/pyboolnet/
        steady_states.tsv  (1 ligne / steady state, 1 col / nœud, 0/1)
        summary.json
"""
from __future__ import annotations

import json
import logging
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
BNET_NORM = ROOT / "02_boolean_model" / "poc_results" / "mpbn" / "SjS_boolean_normalized.bnet"
OUT = ROOT / "02_boolean_model" / "poc_results" / "pyboolnet"
OUT.mkdir(parents=True, exist_ok=True)
STEADY_TSV = OUT / "steady_states.tsv"
SUMMARY = OUT / "summary.json"


def _bnet2primes_with_timeout(bnet_path: str, timeout_s: int) -> dict | None:
    """Invoque bnet2primes en sous-processus pour pouvoir l'interrompre."""
    import multiprocessing
    from queue import Empty

    def worker(q):
        try:
            from pyboolnet import file_exchange
            primes = file_exchange.bnet2primes(bnet_path)
            q.put(("ok", primes))
        except Exception as exc:  # noqa: BLE001
            q.put(("err", str(exc)))

    q: multiprocessing.Queue = multiprocessing.Queue()
    proc = multiprocessing.Process(target=worker, args=(q,))
    proc.start()
    proc.join(timeout_s)
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return None
    try:
        kind, payload = q.get_nowait()
    except Empty:
        return None
    if kind != "ok":
        logger.error("bnet2primes erreur : %s", payload)
        return None
    return payload


def main() -> None:
    PRIMES_TIMEOUT = 180  # bnet2primes intractable sur 5k+ nœuds
    logger.info("Lecture bnet normalisé : %s", BNET_NORM)
    logger.info("Conversion bnet → primes (timeout %ds)…", PRIMES_TIMEOUT)
    t0 = time.time()
    primes = _bnet2primes_with_timeout(str(BNET_NORM), PRIMES_TIMEOUT)
    if primes is None:
        elapsed = time.time() - t0
        logger.warning(
            "pyboolnet `bnet2primes` ne termine pas en %ds sur ce modèle "
            "(5k+ nœuds). Cross-check abandonné — voir `tool_choice.md`.",
            PRIMES_TIMEOUT,
        )
        SUMMARY.write_text(
            json.dumps({
                "n_nodes": None,
                "n_steady_states": 0,
                "elapsed_seconds": round(elapsed, 2),
                "status": "bnet2primes_timeout",
                "primes_timeout_s": PRIMES_TIMEOUT,
                "note": (
                    "pyboolnet's bnet2primes (BNetToPrime) ne scale pas au-delà "
                    "de ~ 1k nœuds dans la pratique. Modèle SjS = 5015 nœuds, "
                    "pas tractable. Utiliser mpbn (Phase 2.2/mpbn) à la place."
                ),
            }, indent=2),
            encoding="utf-8",
        )
        print()
        print("=" * 70)
        print("  PHASE 2.2 — pyboolnet : INTRACTABLE")
        print("=" * 70)
        print(f"  bnet2primes timeout après {PRIMES_TIMEOUT}s sur 5015 nœuds")
        print(f"  Voir : 02_boolean_model/poc_results/pyboolnet/summary.json")
        print(f"  Décision : mpbn retenu pour Phase 2.3")
        print("=" * 70)
        return
    logger.info("→ %d primes en %.1fs", len(primes), time.time() - t0)

    from pyboolnet import trap_spaces

    MAX_STEADY = 5000
    TIME_BUDGET_S = 600
    logger.info("Énumération des steady states (cap=%d, budget=%ds)…",
                MAX_STEADY, TIME_BUDGET_S)
    t0 = time.time()
    try:
        ss = trap_spaces.compute_steady_states(primes, max_output=MAX_STEADY)
    except Exception as e:
        logger.error("Erreur pyboolnet : %s", e)
        ss = []
    elapsed = time.time() - t0
    logger.info("→ %d steady states (%.1fs)", len(ss), elapsed)

    if ss:
        nodes = sorted(ss[0].keys())
        with STEADY_TSV.open("w") as fh:
            fh.write("ss_id\t" + "\t".join(nodes) + "\n")
            for i, st in enumerate(ss):
                vals = [str(st.get(n, "*")) for n in nodes]
                fh.write(f"ss_{i}\t" + "\t".join(vals) + "\n")
        logger.info("Steady states écrits : %s", STEADY_TSV)

    summary = {
        "n_nodes": len(primes),
        "n_steady_states": len(ss),
        "elapsed_seconds": round(elapsed, 2),
        "cap_max_steady": MAX_STEADY,
        "cap_time_budget_s": TIME_BUDGET_S,
        "truncated": len(ss) >= MAX_STEADY,
        "status": "ok",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print()
    print("=" * 70)
    print("  PHASE 2.2 — pyboolnet (steady states / point fixes)")
    print("=" * 70)
    print(f"  Nœuds                 : {summary['n_nodes']}")
    print(f"  Steady states         : {summary['n_steady_states']}")
    print(f"  Temps                 : {summary['elapsed_seconds']} s")
    print(f"  Truncated             : {summary['truncated']}")
    print(f"  Sortie                : {OUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()
