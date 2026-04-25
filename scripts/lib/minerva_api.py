"""
minerva_api.py — Wrapper pour l'API REST MINERVA (sjdmap.elixir-luxembourg.org).

Fonctions publiques :
    get_all_elements(session)  → list[dict]
    get_all_reactions(session) → list[dict]
    get_project_info(session)  → dict
    get_model_info(session)    → dict
    make_session()             → requests.Session  (retry + headers)
    load_or_fetch_elements(cache_dir, session) → list[dict]
    load_or_fetch_reactions(cache_dir, session) → list[dict]

Constantes :
    BASE_URL, PROJECT_ID, MODEL_ID
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration MINERVA
# ---------------------------------------------------------------------------

BASE_URL = "https://sjdmap.elixir-luxembourg.org/minerva/api"
PROJECT_ID = "SjD_Map_no_BD"
MODEL_ID = 7

# Taille de page pour la pagination (l'API renvoie un tableau plat sans totalCount)
PAGE_SIZE = 500

# Délai poli entre les pages (secondes)
REQUEST_DELAY = 0.2

# ---------------------------------------------------------------------------
# Session HTTP avec retry automatique
# ---------------------------------------------------------------------------

def make_session() -> requests.Session:
    """
    Crée une Session requests configurée avec :
    - 5 tentatives avec backoff exponentiel (1, 2, 4, 8, 16 s)
    - Retry sur les codes HTTP 429, 500, 502, 503, 504
    - User-Agent identifiant le projet
    """
    session = requests.Session()

    retry = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    session.headers.update({
        "User-Agent": "SjS-DigitalTwin/0.1 (research; contact: nathan.foulquier@chu-brest.fr)",
        "Accept": "application/json",
    })
    return session


# ---------------------------------------------------------------------------
# Fonctions d'interrogation de l'API
# ---------------------------------------------------------------------------

def get_project_info(session: requests.Session) -> dict[str, Any]:
    """Retourne les métadonnées du projet MINERVA (version, nombre de modèles, etc.)."""
    url = f"{BASE_URL}/projects/{PROJECT_ID}/"
    resp = session.get(url, timeout=30)
    resp.raise_for_status()
    return resp.json()


def get_model_info(session: requests.Session) -> dict[str, Any]:
    """Retourne les métadonnées du modèle SjD Map (id=7)."""
    url = f"{BASE_URL}/projects/{PROJECT_ID}/models/{MODEL_ID}/"
    resp = session.get(url, timeout=30)
    resp.raise_for_status()
    return resp.json()


def _paginate(
    session: requests.Session,
    endpoint: str,
    page_size: int = PAGE_SIZE,
    delay: float = REQUEST_DELAY,
) -> list[dict[str, Any]]:
    """
    Itère sur un endpoint MINERVA. L'API peut :
      - soit renvoyer tout en une seule réponse (ignore les query params)
      - soit paginer correctement avec ?page=N&size=N

    La fin est détectée par :
      1. La réponse est plus petite que page_size → dernière page
      2. La page suivante contient les mêmes IDs que la page précédente → l'API
         ignore la pagination, on a déjà tout

    Args:
        session   : Session requests configurée
        endpoint  : URL complète (sans paramètres de pagination)
        page_size : Nombre d'éléments par page
        delay     : Pause entre les pages (secondes)

    Returns:
        Liste dédupliquée par champ `id` de tous les objets JSON.

    Raises:
        requests.HTTPError si une page renvoie un code >= 400.
        ValueError si la réponse n'est pas un tableau JSON.
    """
    results: list[dict] = []
    seen_ids: set = set()
    page = 0

    while True:
        params = {"size": page_size, "page": page}
        logger.debug("GET %s page=%d size=%d", endpoint, page, page_size)

        resp = session.get(endpoint, params=params, timeout=60)
        resp.raise_for_status()

        batch = resp.json()
        if not isinstance(batch, list):
            raise ValueError(
                f"Réponse inattendue (attendu list, reçu {type(batch).__name__}) "
                f"sur {endpoint} page {page}: {str(batch)[:200]}"
            )

        # Filtrer les doublons (cas où l'API ignore la pagination)
        new_items = []
        for item in batch:
            iid = item.get("id")
            if iid is None or iid not in seen_ids:
                if iid is not None:
                    seen_ids.add(iid)
                new_items.append(item)

        results.extend(new_items)
        logger.info(
            "Page %d : %d éléments reçus, %d nouveaux (total unique : %d)",
            page, len(batch), len(new_items), len(results),
        )

        # Conditions d'arrêt
        if len(batch) < page_size:
            logger.debug("Dernière page (batch < size)")
            break
        if len(new_items) == 0:
            logger.info(
                "Aucun nouvel élément sur la page %d → l'API ne pagine pas, arrêt.",
                page,
            )
            break

        page += 1
        if delay > 0:
            time.sleep(delay)

    return results


def get_all_elements(session: requests.Session) -> list[dict[str, Any]]:
    """
    Récupère tous les bioEntities/elements du modèle SjD Map.

    Un élément représente une espèce moléculaire (Protein, Gene, RNA, Complex,
    SimpleMolecule, Phenotype, Drug, Ion, Degraded).

    Champs clés retournés par l'API :
        id, elementId, name, type, compartmentId,
        bounds (x, y, w, h), references (list of {type, link}),
        notes (texte libre avec annotations Reactome/KEGG)

    Returns:
        Liste de dicts, un par élément. Taille attendue : >= 800.
    """
    endpoint = (
        f"{BASE_URL}/projects/{PROJECT_ID}/models/{MODEL_ID}"
        f"/bioEntities/elements/"
    )
    elements = _paginate(session, endpoint)

    assert len(elements) >= 800, (
        f"Nombre d'éléments suspicieusement bas : {len(elements)} "
        f"(attendu >= 800 d'après la publication)"
    )
    logger.info("Total éléments : %d", len(elements))
    return elements


def get_all_reactions(session: requests.Session) -> list[dict[str, Any]]:
    """
    Récupère toutes les réactions du modèle SjD Map.

    Une réaction représente une interaction moléculaire (State transition,
    Transport, Heterodimer association, Physical stimulation, Catalysis,
    Inhibition, Negative influence, Transcription).

    Champs clés retournés par l'API :
        id, reactionId, type,
        reactants (list of {element, stoichiometry}),
        products  (list of {element, stoichiometry}),
        modifiers (list of {element, type}),
        notes, references

    Returns:
        Liste de dicts, une par réaction. Taille attendue : >= 590.
    """
    endpoint = (
        f"{BASE_URL}/projects/{PROJECT_ID}/models/{MODEL_ID}"
        f"/bioEntities/reactions/"
    )
    reactions = _paginate(session, endpoint)

    assert len(reactions) >= 590, (
        f"Nombre de réactions suspicieusement bas : {len(reactions)} "
        f"(attendu >= 590 d'après la publication)"
    )
    logger.info("Total réactions : %d", len(reactions))
    return reactions


# ---------------------------------------------------------------------------
# Cache local (évite de ré-interroger l'API à chaque run)
# ---------------------------------------------------------------------------

def _cache_path(cache_dir: Path, name: str) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{name}.json"


def _save_cache(data: list[dict], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)
    logger.info("Cache sauvegardé : %s (%d entrées)", path, len(data))


def _load_cache(path: Path) -> list[dict] | None:
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    logger.info("Cache chargé : %s (%d entrées)", path, len(data))
    return data


def load_or_fetch_elements(
    cache_dir: Path,
    session: requests.Session | None = None,
    force_refresh: bool = False,
) -> list[dict[str, Any]]:
    """
    Charge les éléments depuis le cache local si disponible,
    sinon les récupère depuis MINERVA et les met en cache.

    Args:
        cache_dir     : Répertoire de cache (ex. Path("01_disease_map/cache"))
        session       : Session requests (créée automatiquement si None)
        force_refresh : Si True, ignore le cache et ré-interroge l'API

    Returns:
        Liste complète des éléments.
    """
    path = _cache_path(cache_dir, "elements_raw")
    if not force_refresh:
        cached = _load_cache(path)
        if cached is not None:
            return cached

    if session is None:
        session = make_session()

    elements = get_all_elements(session)
    _save_cache(elements, path)
    return elements


def load_or_fetch_reactions(
    cache_dir: Path,
    session: requests.Session | None = None,
    force_refresh: bool = False,
) -> list[dict[str, Any]]:
    """
    Charge les réactions depuis le cache local si disponible,
    sinon les récupère depuis MINERVA et les met en cache.

    Args:
        cache_dir     : Répertoire de cache (ex. Path("01_disease_map/cache"))
        session       : Session requests (créée automatiquement si None)
        force_refresh : Si True, ignore le cache et ré-interroge l'API

    Returns:
        Liste complète des réactions.
    """
    path = _cache_path(cache_dir, "reactions_raw")
    if not force_refresh:
        cached = _load_cache(path)
        if cached is not None:
            return cached

    if session is None:
        session = make_session()

    reactions = get_all_reactions(session)
    _save_cache(reactions, path)
    return reactions


# ---------------------------------------------------------------------------
# Utilitaires d'inspection
# ---------------------------------------------------------------------------

def summarize_elements(elements: list[dict]) -> dict[str, Any]:
    """
    Retourne un résumé statistique des éléments pour vérification rapide.

    Returns:
        Dict avec :
            total         : nombre total d'éléments
            by_type       : {type: count}
            by_compartment: {compartmentId: count}
            missing_name  : nombre d'éléments sans nom
            missing_bounds: nombre d'éléments sans coordonnées spatiales
    """
    from collections import Counter

    by_type = Counter(e.get("type") for e in elements)
    by_compartment = Counter(e.get("compartmentId") for e in elements)
    missing_name = sum(1 for e in elements if not e.get("name"))
    missing_bounds = sum(1 for e in elements if not e.get("bounds"))

    return {
        "total": len(elements),
        "by_type": dict(by_type.most_common()),
        "by_compartment": dict(by_compartment.most_common()),
        "missing_name": missing_name,
        "missing_bounds": missing_bounds,
    }


def summarize_reactions(reactions: list[dict]) -> dict[str, Any]:
    """
    Retourne un résumé statistique des réactions.

    Returns:
        Dict avec :
            total           : nombre total de réactions
            by_type         : {type: count}
            with_modifiers  : nombre de réactions avec au moins un modifier
            modifier_types  : {modifier.type: count}
    """
    from collections import Counter

    by_type = Counter(r.get("type") for r in reactions)
    with_modifiers = sum(1 for r in reactions if r.get("modifiers"))
    modifier_types: Counter = Counter()
    for r in reactions:
        for m in r.get("modifiers", []):
            modifier_types[m.get("type")] += 1

    return {
        "total": len(reactions),
        "by_type": dict(by_type.most_common()),
        "with_modifiers": with_modifiers,
        "modifier_types": dict(modifier_types.most_common()),
    }


def get_phenotype_nodes(elements: list[dict]) -> list[dict]:
    """Retourne uniquement les éléments de type Phenotype."""
    return [e for e in elements if e.get("type") == "Phenotype"]


def find_element_by_name(
    elements: list[dict], name: str, exact: bool = False
) -> list[dict]:
    """
    Recherche des éléments par nom.

    Args:
        name  : Terme de recherche (insensible à la casse)
        exact : Si True, correspondance exacte ; sinon, sous-chaîne

    Returns:
        Liste des éléments correspondants.
    """
    name_lower = name.lower()
    if exact:
        return [e for e in elements if (e.get("name") or "").lower() == name_lower]
    return [e for e in elements if name_lower in (e.get("name") or "").lower()]


# ---------------------------------------------------------------------------
# Point d'entrée pour test rapide en CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import pprint
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stdout,
    )

    cache_dir = Path(__file__).parents[2] / "01_disease_map" / "cache"
    session = make_session()

    print("=== Infos projet ===")
    pprint.pprint(get_project_info(session))

    print("\n=== Infos modèle ===")
    pprint.pprint(get_model_info(session))

    print("\n=== Chargement des éléments ===")
    elements = load_or_fetch_elements(cache_dir, session)
    pprint.pprint(summarize_elements(elements))

    print("\n=== Chargement des réactions ===")
    reactions = load_or_fetch_reactions(cache_dir, session)
    pprint.pprint(summarize_reactions(reactions))

    print("\n=== Nœuds Phenotype ===")
    for p in get_phenotype_nodes(elements):
        print(f"  [{p.get('id')}] {p.get('name')}")
