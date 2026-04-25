"""
map_audit.py — Fonctions d'audit topologique et sémantique de la SjD Map.

Phase 1.1 de la ROADMAP : produit les statistiques nécessaires pour informer
la stratégie de dissociation cell-type (Phase 1.2).

Fonctions principales :
    build_graph(elements, reactions)      → networkx.DiGraph
    topology_summary(graph)               → dict
    annotation_coverage(elements)         → dict
    detect_intercellular_pivots(elements) → list[dict]
    pathway_enrichment_via_notes(elements) → dict[str, list[str]]
    classify_compartments(elements)       → dict[int|None, dict]
"""

from __future__ import annotations

import logging
import re
from collections import Counter, defaultdict
from typing import Any

import networkx as nx

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Familles de marqueurs cell-type — cf. CONVENTIONS.md §1.2
# ---------------------------------------------------------------------------

CELLTYPE_MARKERS: dict[str, set[str]] = {
    "SGEC": {
        "AQP5", "AQP3", "MUC5B", "MUC7", "KRT7", "KRT8", "KRT18", "KRT19",
        "EPCAM", "CFTR", "CDH1", "CLDN1", "CLDN3", "CLDN4", "OCLN", "ZO1",
        "AMY1A", "PRB1", "STATH",
    },
    "TH1": {
        "CD3D", "CD3E", "CD3G", "CD4", "IFNG", "TBX21", "STAT4",
    },
    "TH17": {
        "IL17A", "IL17F", "IL21", "RORC", "RORA", "STAT3",
        "AHR", "IL23R",
    },
    "TFH": {
        "CXCR5", "ICOS", "BCL6", "PDCD1",
    },
    "TREG": {
        "FOXP3", "IL2RA", "CTLA4", "IKZF2",
    },
    "BCELL": {
        "CD19", "MS4A1", "CD79A", "CD79B", "CR2", "BTK", "BLK", "BANK1",
        "TNFRSF13B", "TNFRSF13C", "FCRL3",
    },
    "PLASMA": {
        "XBP1", "IRF4", "PRDM1", "SDC1", "MZB1", "JCHAIN", "TNFRSF17",
    },
    "M1": {
        "CD68", "TNF", "IL1B", "IL6", "NOS2", "CXCL10", "CCR7", "FCGR1A",
    },
    "M2": {
        "CD163", "MRC1", "IL10", "ARG1", "MSR1", "CCL22", "TGFB1",
    },
    "PDC": {
        "CLEC4C", "IL3RA", "IRF7", "TLR7", "TLR9", "LILRA4",
    },
}

# Cytokines et chemokines clés en SjD (signature partagée extracellulaire)
SJD_CYTOKINES_OF_INTEREST: set[str] = {
    "IFNA1", "IFNA2", "IFNA13", "IFNB1", "IFNG", "IFNL1", "IFNL2", "IFNL3",
    "TNF", "IL1B", "IL6", "IL10", "IL12A", "IL12B", "IL15", "IL17A", "IL17F",
    "IL21", "IL22", "IL23A", "TNFSF13", "TNFSF13B",  # APRIL, BAFF
    "CXCL9", "CXCL10", "CXCL11", "CXCL13",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL19", "CCL21", "CCL22",
    "LTA", "LTB", "TNFSF14",  # LTα, LTβ, LIGHT
}


# ---------------------------------------------------------------------------
# Construction du graphe à partir des données API MINERVA
# ---------------------------------------------------------------------------


def build_graph(
    elements: list[dict[str, Any]],
    reactions: list[dict[str, Any]],
) -> nx.DiGraph:
    """
    Construit un graphe orienté à partir des éléments et réactions.

    Convention :
        - Un nœud = un élément MINERVA (clé : id entier)
        - Un edge réactant→produit pour chaque réaction
        - Les modifieurs créent des edges supplémentaires modifier→produit
        - Attributs nœud : name, type, compartmentId, references (list)
        - Attributs edge : reaction_id, reaction_type, role (substrate/modifier)

    Args:
        elements   : liste de dicts (sortie de get_all_elements)
        reactions  : liste de dicts (sortie de get_all_reactions)

    Returns:
        networkx.DiGraph
    """
    graph: nx.DiGraph = nx.DiGraph()

    for el in elements:
        eid = el.get("id")
        if eid is None:
            continue
        graph.add_node(
            int(eid),
            name=el.get("name") or "",
            type=el.get("type") or "Unknown",
            compartment_id=el.get("compartmentId"),
            references=el.get("references") or [],
            notes=el.get("notes") or "",
        )

    def _resolve_id(entry: dict[str, Any]) -> int | None:
        """MINERVA API renvoie soit aliasId (int), soit element:{id:...}."""
        eid = entry.get("aliasId")
        if eid is None:
            elem = entry.get("element") or {}
            eid = elem.get("id")
        return int(eid) if eid is not None else None

    for rxn in reactions:
        rxn_id = rxn.get("reactionId") or rxn.get("id")
        rxn_type = rxn.get("type") or "Unknown"

        reactants = [_resolve_id(r) for r in rxn.get("reactants") or []]
        products = [_resolve_id(p) for p in rxn.get("products") or []]
        modifiers = [
            (_resolve_id(m), m.get("type"))
            for m in rxn.get("modifiers") or []
        ]

        for src in reactants:
            for dst in products:
                if src is None or dst is None:
                    continue
                if graph.has_node(src) and graph.has_node(dst):
                    graph.add_edge(
                        src, dst,
                        reaction_id=rxn_id,
                        reaction_type=rxn_type,
                        role="substrate",
                    )

        for mod_id, mod_type in modifiers:
            if mod_id is None:
                continue
            for dst in products:
                if dst is None:
                    continue
                if graph.has_node(mod_id) and graph.has_node(dst):
                    graph.add_edge(
                        mod_id, dst,
                        reaction_id=rxn_id,
                        reaction_type=rxn_type,
                        role=f"modifier:{mod_type}",
                    )

    logger.info(
        "Graphe construit : %d nœuds, %d edges",
        graph.number_of_nodes(), graph.number_of_edges(),
    )
    return graph


# ---------------------------------------------------------------------------
# Topologie
# ---------------------------------------------------------------------------


def topology_summary(graph: nx.DiGraph, hub_threshold: int = 20) -> dict[str, Any]:
    """
    Résumé topologique : composantes, degrees, hubs, isolés, cycles.

    Args:
        graph         : graphe orienté
        hub_threshold : seuil de degree pour qu'un nœud soit "hub"

    Returns:
        dict structuré pour rapport JSON
    """
    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()

    weak_components = list(nx.weakly_connected_components(graph))
    strong_components = list(nx.strongly_connected_components(graph))

    in_degrees = dict(graph.in_degree())
    out_degrees = dict(graph.out_degree())
    total_degrees = {n: in_degrees[n] + out_degrees[n] for n in graph.nodes}

    isolated = [n for n, d in total_degrees.items() if d == 0]
    sources = [n for n, d in in_degrees.items() if d == 0 and out_degrees[n] > 0]
    sinks = [n for n, d in out_degrees.items() if d == 0 and in_degrees[n] > 0]
    hubs = sorted(
        [(n, d) for n, d in total_degrees.items() if d >= hub_threshold],
        key=lambda x: -x[1],
    )

    deg_dist = Counter(total_degrees.values())

    return {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "n_weakly_connected_components": len(weak_components),
        "largest_weak_component_size": max((len(c) for c in weak_components), default=0),
        "n_strongly_connected_components": len(strong_components),
        "largest_strong_component_size": max((len(c) for c in strong_components), default=0),
        "n_isolated_nodes": len(isolated),
        "n_source_nodes": len(sources),
        "n_sink_nodes": len(sinks),
        "n_hubs": len(hubs),
        "hub_threshold": hub_threshold,
        "top_20_hubs": [
            {
                "id": int(n),
                "name": graph.nodes[n].get("name"),
                "type": graph.nodes[n].get("type"),
                "in_degree": in_degrees[n],
                "out_degree": out_degrees[n],
                "total_degree": d,
            }
            for n, d in hubs[:20]
        ],
        "degree_distribution": dict(sorted(deg_dist.items())),
        "isolated_node_names": [
            graph.nodes[n].get("name") for n in isolated[:50]
        ],
    }


# ---------------------------------------------------------------------------
# Annotations & couverture sémantique
# ---------------------------------------------------------------------------

# Préfixes MIRIAM courants dans MINERVA (resource → identifier scheme)
ANNOTATION_PREFIXES: dict[str, str] = {
    "hgnc": "HGNC",
    "hgnc.symbol": "HGNC",
    "uniprot": "UniProt",
    "ensembl": "Ensembl",
    "entrez": "Entrez",
    "ncbigene": "Entrez",
    "chebi": "ChEBI",
    "kegg": "KEGG",
    "reactome": "Reactome",
    "go": "GO",
    "pubmed": "PubMed",
    "doi": "DOI",
}


def _classify_reference(ref: dict[str, Any]) -> str | None:
    """Mappe une référence MINERVA à un préfixe normalisé."""
    rtype = (ref.get("type") or "").lower()
    link = (ref.get("link") or "").lower()
    resource = (ref.get("resource") or "").lower()

    for key, label in ANNOTATION_PREFIXES.items():
        if key in rtype or key in link or key in resource:
            return label
    return None


def annotation_coverage(elements: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Compte la couverture des annotations externes.

    Returns:
        {
            "total_elements": int,
            "by_annotation": {HGNC: count, UniProt: count, ...},
            "with_any_annotation": int,
            "without_annotation": int,
            "with_pmid": int,
        }
    """
    total = len(elements)
    by_label: Counter = Counter()
    with_any = 0
    with_pmid = 0
    sample_no_annot: list[str] = []

    for el in elements:
        refs = el.get("references") or []
        labels: set[str] = set()
        has_pmid = False
        for r in refs:
            label = _classify_reference(r)
            if label:
                labels.add(label)
            if label == "PubMed" or "pubmed" in (r.get("type") or "").lower():
                has_pmid = True

        if labels:
            with_any += 1
            for label in labels:
                by_label[label] += 1
        elif len(sample_no_annot) < 30:
            sample_no_annot.append(el.get("name") or f"id={el.get('id')}")

        if has_pmid:
            with_pmid += 1

    return {
        "total_elements": total,
        "by_annotation": dict(by_label.most_common()),
        "with_any_annotation": with_any,
        "without_annotation": total - with_any,
        "with_pmid": with_pmid,
        "coverage_pct_any": round(100 * with_any / total, 2) if total else 0.0,
        "sample_without_annotation": sample_no_annot,
    }


def reaction_pmid_coverage(reactions: list[dict[str, Any]]) -> dict[str, Any]:
    """Couverture PMID sur les réactions."""
    total = len(reactions)
    with_ref = 0
    with_pmid = 0
    for r in reactions:
        refs = r.get("references") or []
        if refs:
            with_ref += 1
        for ref in refs:
            if "pubmed" in (ref.get("type") or "").lower():
                with_pmid += 1
                break

    return {
        "total_reactions": total,
        "with_any_reference": with_ref,
        "with_pmid": with_pmid,
        "pmid_coverage_pct": round(100 * with_pmid / total, 2) if total else 0.0,
    }


# ---------------------------------------------------------------------------
# Compartiments
# ---------------------------------------------------------------------------


def classify_compartments(elements: list[dict[str, Any]]) -> dict[Any, dict[str, Any]]:
    """
    Pour chaque compartmentId présent : nombre d'espèces et top types.

    Returns:
        {compartment_id: {"n_species": int, "by_type": {...}, "sample_names": [...]}}
    """
    out: dict[Any, dict[str, Any]] = defaultdict(
        lambda: {"n_species": 0, "by_type": Counter(), "sample_names": []}
    )
    for el in elements:
        cid = el.get("compartmentId")
        out[cid]["n_species"] += 1
        out[cid]["by_type"][el.get("type")] += 1
        if len(out[cid]["sample_names"]) < 10:
            name = el.get("name")
            if name:
                out[cid]["sample_names"].append(name)

    # Convertir Counter → dict pour sérialisation JSON
    for cid in out:
        out[cid]["by_type"] = dict(out[cid]["by_type"].most_common())
    return dict(out)


# ---------------------------------------------------------------------------
# Pivots intercellulaires
# ---------------------------------------------------------------------------


def detect_celltype_marker_hits(
    elements: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    """
    Pour chaque cell-type, liste les nœuds dont le nom matche un marqueur.

    Returns:
        {celltype: [{"id": ..., "name": ..., "type": ...}, ...]}
    """
    hits: dict[str, list[dict]] = {ct: [] for ct in CELLTYPE_MARKERS}
    for el in elements:
        name = (el.get("name") or "").upper().strip()
        # Strip common suffixes/prefixes
        name_clean = re.sub(r"[_\-\s]+", "", name)
        for ct, markers in CELLTYPE_MARKERS.items():
            for marker in markers:
                if marker == name or marker == name_clean:
                    hits[ct].append({
                        "id": el.get("id"),
                        "name": el.get("name"),
                        "type": el.get("type"),
                        "matched_marker": marker,
                    })
                    break
    return hits


def detect_intercellular_pivots(
    elements: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """
    Identifie les nœuds correspondant à des cytokines / chemokines clés en SjD.
    Ce sont les candidats naturels pour les edges intercellulaires (Phase 1.4).

    Returns:
        liste de dicts {id, name, type, matched_cytokine}
    """
    out = []
    for el in elements:
        name = (el.get("name") or "").upper().strip()
        if name in SJD_CYTOKINES_OF_INTEREST:
            out.append({
                "id": el.get("id"),
                "name": el.get("name"),
                "type": el.get("type"),
                "matched_cytokine": name,
            })
    return out


# ---------------------------------------------------------------------------
# Analyse des notes (pathways encodés en texte libre)
# ---------------------------------------------------------------------------

PATHWAY_PATTERNS = [
    re.compile(r"REACT_\w+", re.I),
    re.compile(r"R-HSA-\d+", re.I),
    re.compile(r"hsa\d{5}", re.I),  # KEGG
    re.compile(r"GO:\d{7}"),
]


def extract_pathway_hints_from_notes(
    elements: list[dict[str, Any]],
) -> dict[str, list[str]]:
    """
    Extraction d'identifiants Reactome / KEGG / GO depuis le champ `notes`
    des éléments MINERVA (texte libre).

    Returns:
        {element_name: [pathway_id, ...]}
    """
    out: dict[str, list[str]] = {}
    for el in elements:
        notes = el.get("notes") or ""
        if not notes:
            continue
        ids: list[str] = []
        for pat in PATHWAY_PATTERNS:
            ids.extend(pat.findall(notes))
        if ids:
            name = el.get("name") or f"id={el.get('id')}"
            out[name] = sorted(set(ids))
    return out


# ---------------------------------------------------------------------------
# Synthèse globale
# ---------------------------------------------------------------------------


def full_audit(
    elements: list[dict[str, Any]],
    reactions: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Lance toutes les analyses d'audit et retourne un dict consolidé.
    """
    graph = build_graph(elements, reactions)
    return {
        "topology": topology_summary(graph),
        "element_annotations": annotation_coverage(elements),
        "reaction_annotations": reaction_pmid_coverage(reactions),
        "compartments": classify_compartments(elements),
        "celltype_marker_hits": detect_celltype_marker_hits(elements),
        "intercellular_pivots": detect_intercellular_pivots(elements),
        "pathway_hints_from_notes": extract_pathway_hints_from_notes(elements),
    }
