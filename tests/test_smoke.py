"""Smoke tests — vérifie que les modules lib s'importent."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def test_import_celldesigner_xml() -> None:
    from lib import celldesigner_xml  # noqa: F401


def test_import_minerva_api() -> None:
    from lib import minerva_api  # noqa: F401


def test_import_node_classifier() -> None:
    from lib import node_classifier  # noqa: F401


def test_celldesigner_namespaces() -> None:
    from lib.celldesigner_xml import NS_CD, NS_SBML

    assert NS_SBML == "http://www.sbml.org/sbml/level2/version4"
    assert NS_CD == "http://www.sbml.org/2001/ns/celldesigner"


def test_minerva_session_factory() -> None:
    from lib.minerva_api import make_session

    session = make_session()
    assert session.headers["Accept"] == "application/json"
    assert "SjS-DigitalTwin" in session.headers["User-Agent"]
