"""Check documented and committed observable words without computing moments.

Canonical formatting must preserve the frozen AST, including transform
sandwiches and sea wrappers. Pitfall: arbitrary static identifiers are valid
operator names; the parser cannot decide whether a model supplies them.
"""

import json
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.words.catalogue import CATALOGUE, match_catalogue
from lsquant.words.grammar import parse


def _words(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "word" and isinstance(item, str):
                yield item
            else:
                yield from _words(item)
    elif isinstance(value, list):
        for item in value:
            yield from _words(item)


CASES = [(entry.shape_id, entry.canonical_example) for entry in CATALOGUE]
for path in sorted(set(REPO.glob("runs/goldens/*.json")) |
                   set(REPO.glob("examples/**/run*.json"))):
    CASES.extend((f"{path.relative_to(REPO)}:{i}", word)
                 for i, word in enumerate(_words(json.loads(path.read_text()))))


@pytest.mark.parametrize("label, text", CASES, ids=[c[0] for c in CASES])
def test_documented_word_roundtrip(label, text):
    """Preserve slots, rank and sea status after canonical formatting."""
    word = parse(text)
    assert parse(word.canonical()) == word
    assert parse(str(word)).to_dict() == word.to_dict()


@pytest.mark.parametrize("entry", CATALOGUE, ids=lambda e: e.shape_id)
def test_catalogue_examples_match(entry):
    """Resolve every documented shape back to its own catalogue entry."""
    assert match_catalogue(parse(entry.canonical_example)).entry == entry


@pytest.mark.parametrize("text", [
    "Tr[delta(E-H)", "Tr[delta(E-H)]]", "sea[Tr[delta(E-H)]",
    "Tr[]", "Tr[sea[]]", "Tr[Unknown(E)]", "Tr[$ delta(E-H)]",
    "Tr[[Q,U(t)] delta(E-H)]", "Tr[VX@unknown delta(E-H)]",
])
def test_malformed_word_rejected(text):
    """Reject broken brackets, empty traces and unknown token syntax."""
    with pytest.raises(ValueError):
        parse(text)


