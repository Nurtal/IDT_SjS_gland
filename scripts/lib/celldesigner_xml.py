"""
celldesigner_xml.py — Primitives lxml pour construire un XML CellDesigner 4.4.2.

Ce module fournit toutes les fonctions bas-niveau nécessaires pour reconstruire
un fichier SBML/CellDesigner valide depuis les données JSON de l'API MINERVA.

Le format doit être strictement conforme au schéma CellDesigner 4.4.2 car CaSQ
est très sensible aux noms d'éléments, aux namespaces et aux valeurs des attributs.

Fonctions publiques principales :
    make_sbml_root(model_id, model_name, width, height)  → Element
    add_compartment(model, comp_id, comp_name, ...)      → Element
    make_species(elem, compartment_id, module_prefix)    → Element
    make_reaction(rxn, id_map, reaction_id_prefix)       → Element
    make_transport_reaction(...)                          → Element
    make_heterodimer_reaction(...)                        → Element
    write_celldesigner_xml(root, path)
    validate_sbml(path)                                   → (bool, list[str])

Constantes de mapping :
    MINERVA_TYPE_TO_CD, MINERVA_REACTION_TO_CD,
    MODIFIER_TYPE_TO_CD, COMPARTMENT_NAMES, MIRIAM_PREFIX
"""

from __future__ import annotations

import re
import logging
from pathlib import Path
from typing import Any

from lxml import etree

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Namespaces — DOIVENT être exacts pour la compatibilité CaSQ
# ---------------------------------------------------------------------------

NS_SBML = "http://www.sbml.org/sbml/level2/version4"
NS_CD = "http://www.sbml.org/2001/ns/celldesigner"
NS_RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
NS_BQBIOL = "http://biomodels.net/biology-qualifiers/"
NS_DC = "http://purl.org/dc/elements/1.1/"

NSMAP: dict[str | None, str] = {
    None: NS_SBML,
    "celldesigner": NS_CD,
    "rdf": NS_RDF,
    "bqbiol": NS_BQBIOL,
    "dc": NS_DC,
}

def _cd(tag: str) -> str:
    return f"{{{NS_CD}}}{tag}"

def _rdf(tag: str) -> str:
    return f"{{{NS_RDF}}}{tag}"

def _bq(tag: str) -> str:
    return f"{{{NS_BQBIOL}}}{tag}"

def _sbml(tag: str) -> str:
    return f"{{{NS_SBML}}}{tag}"

# ---------------------------------------------------------------------------
# Mappings MINERVA → CellDesigner
# ---------------------------------------------------------------------------

MINERVA_TYPE_TO_CD: dict[str, str] = {
    "Protein":        "PROTEIN",
    "Gene":           "GENE",
    "RNA":            "RNA",
    "Complex":        "COMPLEX",
    "SimpleMolecule": "SIMPLE_MOLECULE",
    "Phenotype":      "PHENOTYPE",
    "Degraded":       "DEGRADED",
    "Drug":           "DRUG",
    "Ion":            "ION",
    "Unknown":        "UNKNOWN",
}

# Réactions : le type est posé à la fois sur l'attribut SBML <reaction type=...>
# ET dans <celldesigner:reactionType> — les deux doivent correspondre.
MINERVA_REACTION_TO_CD: dict[str, str] = {
    "State transition":       "STATE_TRANSITION",
    "Transcription":          "TRANSCRIPTION",
    "Translation":            "TRANSLATION",
    "Transport":              "TRANSPORT",
    "Heterodimer association":"HETERODIMER_ASSOCIATION",
    "Physical stimulation":   "PHYSICAL_STIMULATION",
    "Catalysis":              "CATALYSIS",
    "Inhibition":             "INHIBITION",
    "Negative influence":     "UNKNOWN_NEGATIVE_INFLUENCE",
    "Positive influence":     "UNKNOWN_POSITIVE_INFLUENCE",
    "Boolean logic gate":     "BOOLEAN_LOGIC_GATE_AND",
    "Reduced physical stimulation": "REDUCED_PHYSICAL_STIMULATION",
    "Reduced modulation":     "REDUCED_MODULATION",
    "Trigger":                "TRIGGER",
}

# Types de modifiers dans les réactions
MODIFIER_TYPE_TO_CD: dict[str, str] = {
    "Catalysis":              "CATALYSIS",
    "Inhibition":             "INHIBITION",
    "Physical stimulation":   "PHYSICAL_STIMULATION",
    "Negative influence":     "UNKNOWN_NEGATIVE_INFLUENCE",
    "Positive influence":     "UNKNOWN_POSITIVE_INFLUENCE",
    "Modulation":             "MODULATION",
    "Trigger":                "TRIGGER",
    "Reduced physical stimulation": "REDUCED_PHYSICAL_STIMULATION",
}

# Compartiments MINERVA → noms CellDesigner
COMPARTMENT_NAMES: dict[int | None, str] = {
    20513: "Cytoplasm",
    21231: "Nucleus",
    21555: "Extracellular",
    21629: "Secreted",
    20730: "EndoplasmicReticulum",
    21540: "Phenotypes",
    None:  "Default",
}

# Préfixes MIRIAM pour les annotations de références
MIRIAM_PREFIX: dict[str, str] = {
    "UNIPROT":       "urn:miriam:uniprot:",
    "HGNC_SYMBOL":   "urn:miriam:hgnc.symbol:",
    "ENSEMBL":       "urn:miriam:ensembl:",
    "CHEBI":         "urn:miriam:chebi:",
    "PUBMED":        "urn:miriam:pubmed:",
    "KEGG_PATHWAY":  "urn:miriam:kegg.pathway:",
    "KEGG_COMPOUND": "urn:miriam:kegg.compound:",
    "GO":            "urn:miriam:go:",
    "INTERPRO":      "urn:miriam:interpro:",
    "PFAM":          "urn:miriam:pfam:",
    "REACTOME":      "urn:miriam:reactome:",
}

# ---------------------------------------------------------------------------
# Utilitaires internes
# ---------------------------------------------------------------------------

def _safe_id(raw: str) -> str:
    """
    Convertit une chaîne en identifiant SBML valide (commence par lettre/_,
    ne contient que lettres/chiffres/_ et -).
    """
    s = re.sub(r"[^A-Za-z0-9_\-]", "_", str(raw))
    if s and s[0].isdigit():
        s = "_" + s
    return s or "_unknown"


def _species_id(elem: dict, prefix: str = "") -> str:
    """Génère un ID SBML unique pour une espèce."""
    raw = f"{prefix}_s{elem['id']}" if prefix else f"s{elem['id']}"
    return _safe_id(raw)


def _reaction_id(rxn: dict, prefix: str = "") -> str:
    """Génère un ID SBML unique pour une réaction."""
    raw = f"{prefix}_re{rxn['id']}" if prefix else f"re{rxn['id']}"
    return _safe_id(raw)


def _compartment_id(comp_id: int | None, prefix: str = "") -> str:
    name = COMPARTMENT_NAMES.get(comp_id, f"comp_{comp_id}")
    return _safe_id(f"{prefix}_{name}" if prefix else name)


# ---------------------------------------------------------------------------
# Construction de la racine SBML
# ---------------------------------------------------------------------------

def make_sbml_root(
    model_id: str,
    model_name: str,
    width: float = 10815.0,
    height: float = 6000.0,
) -> etree._Element:
    """
    Crée l'élément racine <sbml> avec le namespace map complet.

    Structure :
        <sbml level="2" version="4">
          <model id="..." name="...">
            <annotation>
              <celldesigner:extension>
                <celldesigner:modelVersion>4.0</celldesigner:modelVersion>
                <celldesigner:modelDisplay .../>
                <celldesigner:listOfCompartmentAliases/>
                <celldesigner:listOfComplexSpeciesAliases/>
                <celldesigner:listOfSpeciesAliases/>
              </celldesigner:extension>
            </annotation>
            <listOfCompartments/>
            <listOfSpecies/>
            <listOfReactions/>
          </model>
        </sbml>

    Returns:
        Tuple (sbml_root, model_element) pour faciliter les ajouts ultérieurs.
    """
    sbml = etree.Element(
        _sbml("sbml"),
        attrib={"level": "2", "version": "4"},
        nsmap=NSMAP,
    )

    model = etree.SubElement(
        sbml, _sbml("model"),
        attrib={"id": _safe_id(model_id), "name": model_name},
    )

    # --- annotation CellDesigner (métadonnées du modèle) ---
    annotation = etree.SubElement(model, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))

    mv = etree.SubElement(cd_ext, _cd("modelVersion"))
    mv.text = "4.0"

    etree.SubElement(
        cd_ext, _cd("modelDisplay"),
        attrib={
            "sizeX": str(int(width)),
            "sizeY": str(int(height)),
            "color": "ffffff",
            "baseFont": "Default",
            "baseFontSize": "12",
        },
    )

    # Listes d'alias CellDesigner (remplies plus tard)
    for tag in (
        "listOfCompartmentAliases",
        "listOfComplexSpeciesAliases",
        "listOfSpeciesAliases",
        "listOfReactionAliases",
    ):
        etree.SubElement(cd_ext, _cd(tag))

    # Listes SBML (remplies plus tard)
    for tag in ("listOfCompartments", "listOfSpecies", "listOfReactions"):
        etree.SubElement(model, _sbml(tag))

    return sbml


def _get_list(sbml_root: etree._Element, list_tag: str) -> etree._Element:
    """Retourne la première liste SBML portant le tag donné (ex. 'listOfSpecies')."""
    result = sbml_root.find(f".//{_sbml(list_tag)}")
    if result is None:
        raise ValueError(f"<{list_tag}> introuvable dans le document SBML")
    return result


# ---------------------------------------------------------------------------
# Compartiments
# ---------------------------------------------------------------------------

def add_compartment(
    sbml_root: etree._Element,
    comp_minerva_id: int | None,
    module_prefix: str = "",
    outside: str | None = None,
) -> etree._Element:
    """
    Ajoute un <compartment> SBML + son alias CellDesigner.

    Args:
        sbml_root       : Racine <sbml> du document
        comp_minerva_id : ID numérique MINERVA du compartiment (ou None)
        module_prefix   : Préfixe du module (ex. "SGEC") pour noms uniques
        outside         : ID SBML du compartiment parent (optionnel)

    Returns:
        L'élément <compartment> créé.
    """
    comp_name = COMPARTMENT_NAMES.get(comp_minerva_id, f"comp_{comp_minerva_id}")
    comp_id = _compartment_id(comp_minerva_id, module_prefix)

    list_of_comp = _get_list(sbml_root, "listOfCompartments")

    # Éviter les doublons
    existing = list_of_comp.find(f"{_sbml('compartment')}[@id='{comp_id}']")
    if existing is not None:
        return existing

    attrib: dict[str, str] = {
        "id": comp_id,
        "name": f"{module_prefix} {comp_name}".strip() if module_prefix else comp_name,
    }
    if outside:
        attrib["outside"] = outside

    compartment = etree.SubElement(list_of_comp, _sbml("compartment"), attrib=attrib)

    # Alias CellDesigner
    cd_aliases = sbml_root.find(f".//{_cd('listOfCompartmentAliases')}")
    if cd_aliases is not None:
        alias = etree.SubElement(
            cd_aliases, _cd("compartmentAlias"),
            attrib={"id": f"ca_{comp_id}", "compartment": comp_id},
        )
        etree.SubElement(
            alias, _cd("bounds"),
            attrib={"x": "0", "y": "0", "w": "800", "h": "600"},
        )

    return compartment


# ---------------------------------------------------------------------------
# Espèces (species)
# ---------------------------------------------------------------------------

def _make_miriam_annotation(species_sbml_id: str, references: list[dict]) -> etree._Element | None:
    """
    Construit le bloc RDF MIRIAM pour les annotations biologiques.

    Structure :
        <rdf:RDF>
          <rdf:Description rdf:about="#species_id">
            <bqbiol:isVersionOf>
              <rdf:Bag>
                <rdf:li rdf:resource="urn:miriam:uniprot:P12345"/>
              </rdf:Bag>
            </bqbiol:isVersionOf>
          </rdf:Description>
        </rdf:RDF>
    """
    if not references:
        return None

    miriam_items = []
    for ref in references:
        ref_type = ref.get("type", "")
        ref_link = ref.get("link", "") or ref.get("resource", "")
        prefix = MIRIAM_PREFIX.get(ref_type)
        if prefix and ref_link:
            # Extraire l'ID depuis le lien (ex. "https://www.uniprot.org/uniprot/P12345")
            ref_id = ref_link.rstrip("/").split("/")[-1]
            miriam_items.append(f"{prefix}{ref_id}")

    if not miriam_items:
        return None

    rdf_root = etree.Element(_rdf("RDF"), nsmap={
        "rdf": NS_RDF,
        "bqbiol": NS_BQBIOL,
    })
    description = etree.SubElement(
        rdf_root, _rdf("Description"),
        attrib={_rdf("about"): f"#{species_sbml_id}"},
    )
    is_version_of = etree.SubElement(description, _bq("isVersionOf"))
    bag = etree.SubElement(is_version_of, _rdf("Bag"))
    for uri in miriam_items:
        etree.SubElement(bag, _rdf("li"), attrib={_rdf("resource"): uri})

    return rdf_root


def make_species(
    elem: dict[str, Any],
    compartment_sbml_id: str,
    module_prefix: str = "",
) -> etree._Element:
    """
    Construit un élément <species> SBML complet depuis un dict MINERVA.

    Structure produite :
        <species id="SGEC_s42" name="STAT1" compartment="SGEC_Cytoplasm"
                 initialAmount="0" boundaryCondition="false" hasOnlySubstanceUnits="false">
          <annotation>
            <celldesigner:extension>
              <celldesigner:positionToCompartment>inside</celldesigner:positionToCompartment>
              <celldesigner:speciesIdentity>
                <celldesigner:class>PROTEIN</celldesigner:class>
                <celldesigner:name>STAT1</celldesigner:name>
              </celldesigner:speciesIdentity>
            </celldesigner:extension>
            <rdf:RDF> ... annotations MIRIAM ... </rdf:RDF>
          </annotation>
        </species>

    Args:
        elem               : Dict élément MINERVA
        compartment_sbml_id: ID SBML du compartiment cible
        module_prefix      : Préfixe du module (ex. "SGEC")

    Returns:
        Élément <species> prêt à insérer dans <listOfSpecies>.
    """
    sid = _species_id(elem, module_prefix)
    name = elem.get("name") or f"unknown_{elem.get('id', 'x')}"
    cd_class = MINERVA_TYPE_TO_CD.get(elem.get("type", ""), "PROTEIN")

    species = etree.Element(
        _sbml("species"),
        attrib={
            "id": sid,
            "name": name,
            "compartment": compartment_sbml_id,
            "initialAmount": "0",
            "boundaryCondition": "false",
            "hasOnlySubstanceUnits": "false",
        },
    )

    # --- annotation CellDesigner ---
    annotation = etree.SubElement(species, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))

    pos = etree.SubElement(cd_ext, _cd("positionToCompartment"))
    pos.text = "inside"

    identity = etree.SubElement(cd_ext, _cd("speciesIdentity"))
    cls_el = etree.SubElement(identity, _cd("class"))
    cls_el.text = cd_class

    name_el = etree.SubElement(identity, _cd("name"))
    name_el.text = name

    # NB : l'alias de position spatiale n'est PAS ajouté ici. CaSQ lit
    # `speciesAlias` depuis `model/annotation/extension/listOfSpeciesAliases/`,
    # pas depuis l'extension intra-species. Utiliser `add_species_alias` après
    # `add_species` pour exposer l'alias au niveau modèle.

    # Notes textuelles (annotations libres MINERVA)
    notes_text = elem.get("notes", "")
    if notes_text:
        notes_el = etree.SubElement(species, _sbml("notes"))
        p_el = etree.SubElement(
            notes_el,
            "{http://www.w3.org/1999/xhtml}p",
        )
        p_el.text = notes_text

    # Annotations MIRIAM (UniProt, HGNC, Ensembl, ChEBI…)
    references = elem.get("references", []) or []
    miriam = _make_miriam_annotation(sid, references)
    if miriam is not None:
        annotation.append(miriam)

    return species


def add_species(sbml_root: etree._Element, species_el: etree._Element) -> None:
    """Insère un élément <species> dans <listOfSpecies>."""
    _get_list(sbml_root, "listOfSpecies").append(species_el)


def add_species_alias(
    sbml_root: etree._Element,
    species_sbml_id: str,
    compartment_sbml_id: str,
    bounds: dict[str, Any] | None = None,
    species_class: str = "PROTEIN",
) -> etree._Element | None:
    """
    Ajoute un `celldesigner:speciesAlias` (ou `complexSpeciesAlias` pour les
    complexes) au niveau modèle, dans `listOfSpeciesAliases` —
    **format attendu par CaSQ** (`./sbml:annotation/cd:extension/
    cd:listOfSpeciesAliases/cd:speciesAlias`).

    Args:
        sbml_root           : racine <sbml> du document
        species_sbml_id     : id de l'espèce (attribut `species`)
        compartment_sbml_id : id du compartiment associé
        bounds              : dict {x, y, w, h} (sinon défaut)
        species_class       : "PROTEIN" / "COMPLEX" / etc. — détermine
                              speciesAlias vs complexSpeciesAlias.

    Returns:
        L'élément `<celldesigner:speciesAlias>` (ou complexSpeciesAlias) créé.
    """
    list_tag = (
        "listOfComplexSpeciesAliases"
        if species_class == "COMPLEX"
        else "listOfSpeciesAliases"
    )
    alias_tag = (
        "complexSpeciesAlias"
        if species_class == "COMPLEX"
        else "speciesAlias"
    )

    list_alias = sbml_root.find(f".//{_cd(list_tag)}")
    if list_alias is None:
        logger.warning("Liste d'alias %s introuvable", list_tag)
        return None

    alias = etree.SubElement(
        list_alias, _cd(alias_tag),
        attrib={
            "id": f"sa_{species_sbml_id}",
            "species": species_sbml_id,
            "compartmentAlias": f"ca_{compartment_sbml_id}",
        },
    )
    b = bounds or {}
    etree.SubElement(
        alias, _cd("bounds"),
        attrib={
            "x": str(b.get("x", 0.0)),
            "y": str(b.get("y", 0.0)),
            "w": str(b.get("w", 60.0)),
            "h": str(b.get("h", 25.0)),
        },
    )
    return alias


# ---------------------------------------------------------------------------
# Réactions
# ---------------------------------------------------------------------------

def _make_species_ref(species_sbml_id: str, stoich: float = 1.0) -> etree._Element:
    return etree.Element(
        _sbml("speciesReference"),
        attrib={"species": species_sbml_id, "stoichiometry": str(stoich)},
    )


def _make_modifier_ref(species_sbml_id: str) -> etree._Element:
    return etree.Element(
        _sbml("modifierSpeciesReference"),
        attrib={"species": species_sbml_id},
    )


def _append_cd_base_participants(
    cd_ext: etree._Element,
    reactant_ids: list[str],
    product_ids: list[str],
) -> None:
    """
    Ajoute `celldesigner:baseReactants` et `celldesigner:baseProducts` à
    l'extension de la réaction — **format requis par CaSQ** pour parser les
    transitions (cf. `casq.readCD.get_transitions`).

    Chaque baseReactant/baseProduct référence à la fois `species` (ID SBML)
    et `alias` (ID de speciesAlias `sa_<species_id>`).
    """
    if reactant_ids:
        lor_cd = etree.SubElement(cd_ext, _cd("baseReactants"))
        for sid in reactant_ids:
            etree.SubElement(
                lor_cd, _cd("baseReactant"),
                attrib={"species": sid, "alias": f"sa_{sid}"},
            )
    if product_ids:
        lop_cd = etree.SubElement(cd_ext, _cd("baseProducts"))
        for sid in product_ids:
            etree.SubElement(
                lop_cd, _cd("baseProduct"),
                attrib={"species": sid, "alias": f"sa_{sid}"},
            )


def make_reaction(
    rxn: dict[str, Any],
    id_map: dict[int, str],
    module_prefix: str = "",
) -> etree._Element | None:
    """
    Construit un élément <reaction> SBML complet depuis un dict MINERVA.

    Structure produite :
        <reaction id="re123" name="..." type="STATE_TRANSITION" reversible="false">
          <listOfReactants>
            <speciesReference species="SGEC_s42"/>
          </listOfReactants>
          <listOfProducts>
            <speciesReference species="SGEC_s43"/>
          </listOfProducts>
          <listOfModifiers>
            <modifierSpeciesReference species="SGEC_s10"/>
          </listOfModifiers>
          <annotation>
            <celldesigner:extension>
              <celldesigner:reactionType>STATE_TRANSITION</celldesigner:reactionType>
              <celldesigner:listOfModification>
                <celldesigner:modification type="CATALYSIS" modifiers="SGEC_s10"/>
              </celldesigner:listOfModification>
            </celldesigner:extension>
          </annotation>
        </reaction>

    Args:
        rxn           : Dict réaction MINERVA
        id_map        : {minerva_element_id (int) → species_sbml_id (str)}
        module_prefix : Préfixe du module pour l'ID de réaction

    Returns:
        Élément <reaction>, ou None si les espèces sont introuvables dans id_map.
    """
    rid = _reaction_id(rxn, module_prefix)
    cd_type = MINERVA_REACTION_TO_CD.get(rxn.get("type", ""), "STATE_TRANSITION")
    name = rxn.get("name") or rid

    # --- Résolution des IDs espèces ---
    reactants_raw = rxn.get("reactants") or []
    products_raw = rxn.get("products") or []
    modifiers_raw = rxn.get("modifiers") or []

    def resolve(entry: dict) -> str | None:
        # MINERVA API renvoie soit "aliasId" (entier), soit
        # "element": {"id": ...} selon les endpoints/versions.
        eid = entry.get("aliasId")
        if eid is None:
            elem = entry.get("element") or {}
            eid = elem.get("id")
        if eid is None:
            return None
        return id_map.get(int(eid))

    reactant_ids = [r for r in (resolve(e) for e in reactants_raw) if r]
    product_ids = [r for r in (resolve(e) for e in products_raw) if r]
    modifier_data = [
        (resolve(m), m.get("type", ""))
        for m in modifiers_raw
    ]
    modifier_data = [(sid, mtype) for sid, mtype in modifier_data if sid]

    if not reactant_ids and not product_ids:
        logger.warning("Réaction %s ignorée : aucun réactant/produit résolu", rid)
        return None

    reaction = etree.Element(
        _sbml("reaction"),
        attrib={
            "id": rid,
            "name": name,
            "type": cd_type,
            "reversible": "false",
        },
    )

    # Réactants
    if reactant_ids:
        lor = etree.SubElement(reaction, _sbml("listOfReactants"))
        for sid in reactant_ids:
            lor.append(_make_species_ref(sid))

    # Produits
    if product_ids:
        lop = etree.SubElement(reaction, _sbml("listOfProducts"))
        for sid in product_ids:
            lop.append(_make_species_ref(sid))

    # Modifiers SBML
    if modifier_data:
        lom = etree.SubElement(reaction, _sbml("listOfModifiers"))
        for sid, _ in modifier_data:
            lom.append(_make_modifier_ref(sid))

    # Annotation CellDesigner
    annotation = etree.SubElement(reaction, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))

    rt_el = etree.SubElement(cd_ext, _cd("reactionType"))
    rt_el.text = cd_type

    _append_cd_base_participants(cd_ext, reactant_ids, product_ids)

    if modifier_data:
        lom_cd = etree.SubElement(cd_ext, _cd("listOfModification"))
        modifier_ids_str = ",".join(sid for sid, _ in modifier_data)
        for sid, mtype in modifier_data:
            cd_mtype = MODIFIER_TYPE_TO_CD.get(mtype, "MODULATION")
            etree.SubElement(
                lom_cd, _cd("modification"),
                attrib={
                    "type": cd_mtype,
                    "modifiers": sid,
                    "aliases": f"sa_{sid}",
                    "targetLineIndex": "-1,0",
                },
            )
        # Ligne visuelle générique
        _append_reaction_line(cd_ext, reactant_ids, product_ids)

    return reaction


def _append_reaction_line(
    cd_ext: etree._Element,
    reactant_ids: list[str],
    product_ids: list[str],
) -> None:
    """Ajoute un élément <celldesigner:line> basique pour le rendu visuel."""
    line = etree.SubElement(
        cd_ext, _cd("line"),
        attrib={"color": "ff000000", "width": "1.0"},
    )
    connect_scheme = etree.SubElement(
        line, _cd("connectScheme"),
        attrib={"connectPolicy": "direct", "rectangleIndex": "0"},
    )
    if reactant_ids:
        etree.SubElement(
            connect_scheme, _cd("listOfLineDirection"),
        )


def add_reaction(sbml_root: etree._Element, reaction_el: etree._Element) -> None:
    """Insère un élément <reaction> dans <listOfReactions>."""
    _get_list(sbml_root, "listOfReactions").append(reaction_el)


# ---------------------------------------------------------------------------
# Réactions inter-cellulaires spécialisées
# ---------------------------------------------------------------------------

def make_transport_reaction(
    reaction_id: str,
    source_species_id: str,
    target_species_id: str,
    reaction_name: str = "",
) -> etree._Element:
    """
    Crée une réaction TRANSPORT entre deux espèces dans des compartiments différents.

    Utilisé pour les arêtes inter-cellulaires de type cytokine sécrétée :
        [Ligand_Cell] → Transport → [Extracellular_Ligand]

    Args:
        reaction_id       : ID SBML de la réaction (ex. "inter_re_001")
        source_species_id : ID SBML de l'espèce source (compartiment Secreted)
        target_species_id : ID SBML de l'espèce cible (compartiment Extracellular)
        reaction_name     : Nom lisible (optionnel)

    Returns:
        Élément <reaction> de type TRANSPORT.
    """
    reaction = etree.Element(
        _sbml("reaction"),
        attrib={
            "id": _safe_id(reaction_id),
            "name": reaction_name or reaction_id,
            "type": "TRANSPORT",
            "reversible": "false",
        },
    )

    lor = etree.SubElement(reaction, _sbml("listOfReactants"))
    lor.append(_make_species_ref(source_species_id))

    lop = etree.SubElement(reaction, _sbml("listOfProducts"))
    lop.append(_make_species_ref(target_species_id))

    annotation = etree.SubElement(reaction, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))
    rt = etree.SubElement(cd_ext, _cd("reactionType"))
    rt.text = "TRANSPORT"

    _append_cd_base_participants(
        cd_ext, [source_species_id], [target_species_id],
    )

    return reaction


def make_physical_stimulation_reaction(
    reaction_id: str,
    source_species_id: str,
    target_species_id: str,
    reaction_name: str = "",
) -> etree._Element:
    """
    Crée une réaction PHYSICAL_STIMULATION.

    Utilisé pour la 2ème étape des arêtes inter-cellulaires :
        [Extracellular_Ligand] → PhysicalStimulation → [Receptor_TargetCell]

    Returns:
        Élément <reaction> de type PHYSICAL_STIMULATION.
    """
    reaction = etree.Element(
        _sbml("reaction"),
        attrib={
            "id": _safe_id(reaction_id),
            "name": reaction_name or reaction_id,
            "type": "PHYSICAL_STIMULATION",
            "reversible": "false",
        },
    )

    lor = etree.SubElement(reaction, _sbml("listOfReactants"))
    lor.append(_make_species_ref(source_species_id))

    lop = etree.SubElement(reaction, _sbml("listOfProducts"))
    lop.append(_make_species_ref(target_species_id))

    annotation = etree.SubElement(reaction, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))
    rt = etree.SubElement(cd_ext, _cd("reactionType"))
    rt.text = "PHYSICAL_STIMULATION"

    _append_cd_base_participants(
        cd_ext, [source_species_id], [target_species_id],
    )

    return reaction


def make_heterodimer_reaction(
    reaction_id: str,
    species1_id: str,
    species2_id: str,
    complex_id: str,
    reaction_name: str = "",
) -> etree._Element:
    """
    Crée une réaction HETERODIMER_ASSOCIATION pour les contacts cellule-cellule.

    Convention CellDesigner : exactement 2 réactants, 1 produit (le complexe).

    Exemple :
        CD4::CD40LG + BCELL::CD40 → Intercellular_CD40LG_CD40

    Args:
        reaction_id : ID SBML de la réaction
        species1_id : ID SBML de la première espèce (ex. "CD4_CD40LG")
        species2_id : ID SBML de la seconde espèce (ex. "BCELL_CD40")
        complex_id  : ID SBML du complexe produit

    Returns:
        Élément <reaction> de type HETERODIMER_ASSOCIATION.
    """
    reaction = etree.Element(
        _sbml("reaction"),
        attrib={
            "id": _safe_id(reaction_id),
            "name": reaction_name or reaction_id,
            "type": "HETERODIMER_ASSOCIATION",
            "reversible": "false",
        },
    )

    lor = etree.SubElement(reaction, _sbml("listOfReactants"))
    lor.append(_make_species_ref(species1_id))
    lor.append(_make_species_ref(species2_id))

    lop = etree.SubElement(reaction, _sbml("listOfProducts"))
    lop.append(_make_species_ref(complex_id))

    annotation = etree.SubElement(reaction, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))
    rt = etree.SubElement(cd_ext, _cd("reactionType"))
    rt.text = "HETERODIMER_ASSOCIATION"

    _append_cd_base_participants(
        cd_ext, [species1_id, species2_id], [complex_id],
    )

    return reaction


def make_extracellular_species(
    gene_name: str,
    compartment_sbml_id: str,
    species_id: str | None = None,
) -> etree._Element:
    """
    Crée une espèce intermédiaire dans le compartiment extracellulaire
    pour les arêtes inter-cellulaires de type Transport.

    Args:
        gene_name          : Nom du gène/protéine (ex. "TNFSF13B")
        compartment_sbml_id: ID SBML du compartiment Extracellular
        species_id         : ID SBML forcé (sinon généré depuis gene_name)

    Returns:
        Élément <species> pour le ligand extracellulaire.
    """
    sid = species_id or _safe_id(f"extracellular_{gene_name}")
    species = etree.Element(
        _sbml("species"),
        attrib={
            "id": sid,
            "name": gene_name,
            "compartment": compartment_sbml_id,
            "initialAmount": "0",
            "boundaryCondition": "false",
            "hasOnlySubstanceUnits": "false",
        },
    )
    annotation = etree.SubElement(species, _sbml("annotation"))
    cd_ext = etree.SubElement(annotation, _cd("extension"))
    pos = etree.SubElement(cd_ext, _cd("positionToCompartment"))
    pos.text = "inside"
    identity = etree.SubElement(cd_ext, _cd("speciesIdentity"))
    cls_el = etree.SubElement(identity, _cd("class"))
    cls_el.text = "PROTEIN"
    name_el = etree.SubElement(identity, _cd("name"))
    name_el.text = gene_name
    return species


# ---------------------------------------------------------------------------
# Écriture et validation
# ---------------------------------------------------------------------------

def write_celldesigner_xml(root: etree._Element, path: Path | str) -> None:
    """
    Sérialise l'arbre XML vers un fichier avec déclaration XML correcte.

    Args:
        root : Racine <sbml> de l'arbre lxml
        path : Chemin de sortie
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tree = etree.ElementTree(root)
    tree.write(
        str(path),
        xml_declaration=True,
        encoding="UTF-8",
        pretty_print=True,
    )
    logger.info("XML écrit : %s (%.1f Ko)", path, path.stat().st_size / 1024)


def validate_sbml(path: Path | str) -> tuple[bool, list[str]]:
    """
    Valide un fichier SBML avec libsbml.

    Note : les schémas L2V4 stricts considèrent les annotations CellDesigner
    (`celldesigner:*`) comme non conformes — la SjD Map originale présente
    ~1000 « Error » de catégorie *General SBML conformance* qui sont en fait
    des artefacts du namespace `celldesigner` non couvert par le XSD. On
    distingue donc :
        - Severity FATAL (3) : erreurs structurelles fatales → bloquantes.
        - Severity Error sur catégorie *General SBML conformance* :
          marquées comme "schema_warnings" et non bloquantes.

    Returns:
        (valid: bool, errors: list[str]) — valid=True si 0 erreur FATAL.
    """
    try:
        import libsbml
    except ImportError:
        logger.warning("libsbml non disponible — validation SBML ignorée")
        return True, []

    reader = libsbml.SBMLReader()
    doc = reader.readSBMLFromFile(str(path))

    errors = []
    fatal = 0
    for i in range(doc.getNumErrors()):
        err = doc.getError(i)
        msg = (
            f"[{err.getSeverityAsString()}] "
            f"L{err.getLine()}:C{err.getColumn()} — {err.getMessage()}"
        )
        errors.append(msg)
        if err.getSeverity() >= libsbml.LIBSBML_SEV_FATAL:
            fatal += 1

    if fatal > 0:
        logger.error("%d erreur(s) fatale(s) SBML dans %s", fatal, path)
    else:
        logger.info(
            "Validation SBML OK : %s (%d schema-warnings non bloquantes)",
            path, len(errors),
        )

    return fatal == 0, errors


def count_elements(sbml_root: etree._Element) -> dict[str, int]:
    """
    Compte les espèces et réactions dans l'arbre XML.

    Returns:
        {"species": N, "reactions": M, "compartments": K}
    """
    return {
        "compartments": len(_get_list(sbml_root, "listOfCompartments")),
        "species":    len(_get_list(sbml_root, "listOfSpecies")),
        "reactions":  len(_get_list(sbml_root, "listOfReactions")),
    }


def build_species_id_map(sbml_root: etree._Element) -> dict[str, str]:
    """
    Retourne {species_id → name} depuis <listOfSpecies>.
    Utile pour les vérifications post-assemblage.
    """
    result = {}
    for sp in _get_list(sbml_root, "listOfSpecies"):
        sid = sp.get("id", "")
        name = sp.get("name", "")
        if sid:
            result[sid] = name
    return result


def check_species_references(sbml_root: etree._Element) -> list[str]:
    """
    Vérifie que toutes les <speciesReference> et <modifierSpeciesReference>
    pointent vers des espèces existantes dans <listOfSpecies>.

    Returns:
        Liste des IDs non résolus (vide si tout est OK).
    """
    known_ids = set(build_species_id_map(sbml_root).keys())
    unresolved = []

    for rxn in _get_list(sbml_root, "listOfReactions"):
        for ref in rxn.iter():
            tag_local = etree.QName(ref.tag).localname if ref.tag else ""
            if tag_local in ("speciesReference", "modifierSpeciesReference"):
                sid = ref.get("species", "")
                if sid and sid not in known_ids:
                    unresolved.append(f"{rxn.get('id', '?')} → {sid}")

    return unresolved
