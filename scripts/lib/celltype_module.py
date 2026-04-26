"""
celltype_module.py — Extraction des sous-cartes par cell-type (Phase 1.3).

Pour chaque cell-type ∈ {SGEC, TH1, TH17, TFH, TREG, BCELL, PLASMA, M1, M2, PDC},
extrait :

    1. La liste des nœuds assignés (depuis `node_to_celltype.tsv`).
    2. Les nœuds extracellulaires (EXTRA, R1) connectés à ces nœuds via une
       réaction MINERVA → conservés comme "ports" intercellulaires (in/out).
    3. Les nœuds phénotype (PHENOTYPE, R2) connectés → conservés comme sorties.
    4. Le sous-graphe induit (réactions dont tous les participants —
       réactants/produits/modifiers — sont dans le set étendu).
    5. Les boucles de feedback détectées (cycles simples du DiGraph induit).
    6. Métriques : n_nodes, n_reactions, n_edges, n_feedback_loops,
       n_phenotype_links, n_extra_ports, density.

Le format de sortie est volontairement TSV/JSON (pas de XML CellDesigner ici) :
l'export CellDesigner est délégué à `scripts/lib/celldesigner_xml.py` pour
qu'on puisse driver avec ou sans CellDesigner GUI selon le besoin expert.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import networkx as nx

logger = logging.getLogger(__name__)


CELL_TYPES: tuple[str, ...] = (
    "SGEC", "TH1", "TH17", "TFH", "TREG",
    "BCELL", "PLASMA", "M1", "M2", "PDC",
)


# ---------------------------------------------------------------------------
# Chargement node_to_celltype.tsv
# ---------------------------------------------------------------------------


@dataclass
class CellTypeAssignment:
    node_id: int
    name: str
    type: str
    compartment_id: int | None
    celltype: str
    confidence: str
    rule: str
    plausibility_score: int
    source: str


def load_node_to_celltype(path: Path) -> list[CellTypeAssignment]:
    rows: list[CellTypeAssignment] = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            try:
                comp = int(r["compartment_id"]) if r.get("compartment_id") else None
            except ValueError:
                comp = None
            try:
                score = int(r.get("plausibility_score") or 0)
            except ValueError:
                score = 0
            rows.append(CellTypeAssignment(
                node_id=int(r["node_id"]),
                name=r.get("node_name") or "",
                type=r.get("node_type") or "",
                compartment_id=comp,
                celltype=r["celltype"],
                confidence=r.get("confidence") or "",
                rule=r.get("rule") or "",
                plausibility_score=score,
                source=r.get("source") or "",
            ))
    return rows


def index_by_celltype(
    rows: list[CellTypeAssignment],
) -> dict[str, set[int]]:
    """Retourne {celltype → set(node_id)} pour les 10 vrais cell-types."""
    out: dict[str, set[int]] = {ct: set() for ct in CELL_TYPES}
    for r in rows:
        if r.celltype in out:
            out[r.celltype].add(r.node_id)
    return out


def collect_extra_and_phenotype(
    rows: list[CellTypeAssignment],
) -> tuple[set[int], set[int]]:
    """Set des node_id EXTRA (R1) et PHENOTYPE (R2)."""
    extra = {r.node_id for r in rows if r.celltype == "EXTRA"}
    phen = {r.node_id for r in rows if r.celltype == "PHENOTYPE"}
    return extra, phen


# ---------------------------------------------------------------------------
# Extraction du sous-module
# ---------------------------------------------------------------------------


@dataclass
class CellTypeModule:
    celltype: str
    core_nodes: set[int]                    # nœuds assignés au cell-type
    extra_ports: set[int] = field(default_factory=set)  # ligands EXTRA connectés
    phenotype_outputs: set[int] = field(default_factory=set)
    reactions: list[dict[str, Any]] = field(default_factory=list)
    induced_graph: nx.DiGraph | None = None
    metrics: dict[str, Any] = field(default_factory=dict)

    @property
    def all_nodes(self) -> set[int]:
        return self.core_nodes | self.extra_ports | self.phenotype_outputs


def _resolve_id(entry: dict[str, Any]) -> int | None:
    eid = entry.get("aliasId")
    if eid is None:
        elem = entry.get("element") or {}
        eid = elem.get("id")
    return int(eid) if eid is not None else None


def _reaction_participants(rxn: dict[str, Any]) -> tuple[set[int], set[int], set[int]]:
    """(reactants, products, modifiers) en sets d'IDs résolus."""
    reactants = {
        rid for rid in (_resolve_id(r) for r in rxn.get("reactants") or [])
        if rid is not None
    }
    products = {
        rid for rid in (_resolve_id(p) for p in rxn.get("products") or [])
        if rid is not None
    }
    modifiers = {
        rid for rid in (_resolve_id(m) for m in rxn.get("modifiers") or [])
        if rid is not None
    }
    return reactants, products, modifiers


def extract_module(
    celltype: str,
    core_nodes: set[int],
    extra_set: set[int],
    phenotype_set: set[int],
    reactions: list[dict[str, Any]],
    full_graph: nx.DiGraph,
) -> CellTypeModule:
    """
    Construit le module pour `celltype` :

    1. Identifie les EXTRA ports = ligands extracellulaires participant à une
       réaction où ≥1 partenaire est dans `core_nodes`.
    2. Identifie les phénotypes connectés = nœuds PHENOTYPE participant à une
       réaction où ≥1 partenaire est dans `core_nodes`.
    3. Garde les réactions où **tous** les participants ∈ core_nodes ∪ ports.
    4. Construit l'induced subgraph DiGraph.
    5. Calcule métriques + feedback loops.
    """
    extra_ports: set[int] = set()
    phen_outputs: set[int] = set()
    kept_reactions: list[dict[str, Any]] = []

    # Pass 1 — détection des ports
    for rxn in reactions:
        r_set, p_set, m_set = _reaction_participants(rxn)
        all_p = r_set | p_set | m_set
        if not (all_p & core_nodes):
            continue
        # ports = participants extracellulaires/phénotype
        for nid in all_p:
            if nid in extra_set:
                extra_ports.add(nid)
            elif nid in phenotype_set:
                phen_outputs.add(nid)

    accepted_set = core_nodes | extra_ports | phen_outputs

    # Pass 2 — réactions strictement contenues dans le set étendu
    for rxn in reactions:
        r_set, p_set, m_set = _reaction_participants(rxn)
        all_p = r_set | p_set | m_set
        if not all_p:
            continue
        if not all_p.issubset(accepted_set):
            continue
        if not (all_p & core_nodes):
            # ne doit pas arriver vu pass 1, mais sécurité
            continue
        kept_reactions.append(rxn)

    # Sous-graphe induit
    sub_nodes = [n for n in accepted_set if full_graph.has_node(n)]
    induced = full_graph.subgraph(sub_nodes).copy()

    # Feedback loops (cycles simples) — limité pour éviter explosion combinatoire
    n_loops = 0
    sample_loops: list[list[int]] = []
    try:
        for i, cyc in enumerate(nx.simple_cycles(induced)):
            n_loops += 1
            if i < 10:
                sample_loops.append(cyc)
            if i >= 5000:
                break
    except Exception as exc:  # noqa: BLE001
        logger.warning("simple_cycles failed for %s: %s", celltype, exc)

    metrics = {
        "celltype": celltype,
        "n_core_nodes": len(core_nodes),
        "n_extra_ports": len(extra_ports),
        "n_phenotype_outputs": len(phen_outputs),
        "n_total_nodes": len(accepted_set),
        "n_reactions": len(kept_reactions),
        "n_edges": induced.number_of_edges(),
        "n_feedback_loops": n_loops,
        "density": round(nx.density(induced), 4) if induced.number_of_nodes() else 0.0,
        "sample_loops": sample_loops[:5],
    }

    return CellTypeModule(
        celltype=celltype,
        core_nodes=set(core_nodes),
        extra_ports=extra_ports,
        phenotype_outputs=phen_outputs,
        reactions=kept_reactions,
        induced_graph=induced,
        metrics=metrics,
    )


# ---------------------------------------------------------------------------
# Sorties TSV
# ---------------------------------------------------------------------------


def write_module_nodes_tsv(
    module: CellTypeModule,
    elements_by_id: dict[int, dict[str, Any]],
    assignments_by_node: dict[int, list[CellTypeAssignment]],
    path: Path,
) -> None:
    fields = [
        "node_id", "node_name", "node_type", "compartment_id",
        "role",  # core / extra_port / phenotype_output
        "rule", "confidence", "plausibility_score", "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for nid in sorted(module.all_nodes):
            el = elements_by_id.get(nid, {})
            if nid in module.core_nodes:
                role = "core"
            elif nid in module.extra_ports:
                role = "extra_port"
            else:
                role = "phenotype_output"
            # trouver la ligne d'assignation correspondant à ce module
            rule = confidence = source = ""
            score = 0
            for a in assignments_by_node.get(nid, []):
                if (role == "core" and a.celltype == module.celltype) or \
                   (role == "extra_port" and a.celltype == "EXTRA") or \
                   (role == "phenotype_output" and a.celltype == "PHENOTYPE"):
                    rule, confidence = a.rule, a.confidence
                    score, source = a.plausibility_score, a.source
                    break
            w.writerow({
                "node_id": nid,
                "node_name": el.get("name") or "",
                "node_type": el.get("type") or "",
                "compartment_id": el.get("compartmentId") or "",
                "role": role,
                "rule": rule,
                "confidence": confidence,
                "plausibility_score": score,
                "source": source,
            })


def write_module_reactions_tsv(
    module: CellTypeModule,
    elements_by_id: dict[int, dict[str, Any]],
    path: Path,
) -> None:
    fields = [
        "reaction_id", "reaction_type",
        "reactants", "products", "modifiers",
    ]

    def fmt(ids: set[int]) -> str:
        parts = []
        for nid in sorted(ids):
            nm = (elements_by_id.get(nid) or {}).get("name") or str(nid)
            parts.append(f"{nid}:{nm}")
        return ";".join(parts)

    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for rxn in module.reactions:
            r_set, p_set, m_set = _reaction_participants(rxn)
            w.writerow({
                "reaction_id": rxn.get("reactionId") or rxn.get("id") or "",
                "reaction_type": rxn.get("type") or "",
                "reactants": fmt(r_set),
                "products": fmt(p_set),
                "modifiers": fmt(m_set),
            })


# ---------------------------------------------------------------------------
# Validation Gate 1.3
# ---------------------------------------------------------------------------


@dataclass
class GateResult:
    celltype: str
    n_nodes: int
    n_reactions: int
    n_feedback_loops: int
    n_phenotype_outputs: int
    pass_min_nodes: bool       # ≥80
    pass_min_reactions: bool   # ≥50 (cible Zerrouk)
    pass_feedback: bool        # ≥1 loop
    pass_phenotype: bool       # ≥1 output

    @property
    def pass_all(self) -> bool:
        return all([
            self.pass_min_nodes, self.pass_min_reactions,
            self.pass_feedback, self.pass_phenotype,
        ])


def evaluate_gate(
    module: CellTypeModule,
    min_nodes: int = 80,
    min_reactions: int = 50,
) -> GateResult:
    n = module.metrics["n_total_nodes"]
    r = module.metrics["n_reactions"]
    fb = module.metrics["n_feedback_loops"]
    ph = module.metrics["n_phenotype_outputs"]
    return GateResult(
        celltype=module.celltype,
        n_nodes=n,
        n_reactions=r,
        n_feedback_loops=fb,
        n_phenotype_outputs=ph,
        pass_min_nodes=n >= min_nodes,
        pass_min_reactions=r >= min_reactions,
        pass_feedback=fb >= 1,
        pass_phenotype=ph >= 1,
    )


# ---------------------------------------------------------------------------
# Construction id_map pour export CellDesigner
# ---------------------------------------------------------------------------


def build_id_map_for_module(
    module: CellTypeModule,
    elements_by_id: dict[int, dict[str, Any]],
) -> dict[int, str]:
    """Construit {minerva_id → species_sbml_id} préfixé par celltype."""
    from lib.celldesigner_xml import _species_id  # noqa: PLC0415
    out: dict[int, str] = {}
    for nid in module.all_nodes:
        el = elements_by_id.get(nid)
        if el:
            out[nid] = _species_id(el, module.celltype)
    return out
