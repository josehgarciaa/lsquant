"""Small public-contract regressions; no build or external model input is needed."""

from __future__ import annotations

import ast
from pathlib import Path
import re
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))


def test_package_readme_is_local_and_shippable():
    metadata = (REPO / "python/pyproject.toml").read_text()
    name = re.search(r'^readme = "([^"]+)"$', metadata, re.MULTILINE).group(1)
    package_root = REPO / "python"
    readme = (package_root / name).resolve()
    assert readme.is_relative_to(package_root)
    assert readme.is_file() and readme.read_text().strip()


def test_engine_instruction_wording_matches_plan():
    stub = ast.parse((REPO / "python/lsquant/_core.pyi").read_text())
    instruction = next(node.value for node in stub.body
                       if isinstance(node, ast.Assign)
                       and any(isinstance(t, ast.Name) and t.id == "Instruction"
                               for t in node.targets))
    assert len(instruction.slice.elts) == 7
    for name in ("python/lsquant/moments.py", "lsquant/core/chebyshev/moment_engine.hpp"):
        text = (REPO / name).read_text()
        assert "seven-instruction" in text
        assert "five-instruction" not in text


def test_legacy_dipole_cache_read_location_is_documented():
    text = (REPO / "docs/how-to/use-siesta.md").read_text()
    assert "model.dir/.lsquant-cache/dipoles-<key>.npz" in text
    assert "same cache-key and certificate validation" in text


def test_unsupported_generator_is_not_catalogued():
    pytest.importorskip("lsquant")
    from lsquant.runfile import GENERATOR_CATALOGUE, resolve_generator_path

    assert "gen_models" not in GENERATOR_CATALOGUE
    with pytest.raises(ValueError, match="Unknown generator name"):
        resolve_generator_path("gen_models", str(REPO))


def test_printed_defaults_match_minimal_resolution():
    pytest.importorskip("lsquant")
    from lsquant.runfile import RunFile, defaults_for_basis

    resolved = RunFile.from_dict({
        "model": {"kind": "csr", "dir": "goldens/stage0/chain1d", "label": "chain1d"},
        "word": "Tr[delta(E-H)]",
        "expansion": {"moments": 64},
    }, repo_root=str(REPO))

    def compare(shared_defaults, actual, path):
        if isinstance(shared_defaults, dict) and isinstance(actual, dict):
            for key in shared_defaults.keys() & actual.keys():
                compare(shared_defaults[key], actual[key], f"{path}.{key}")
        else:
            assert shared_defaults == actual, path

    for section, defaults in defaults_for_basis().items():
        compare(defaults, getattr(resolved, section), section)


def test_loaded_core_public_names_are_declared_in_stub():
    pytest.importorskip("lsquant")
    from lsquant import _core

    tree = ast.parse((REPO / "python/lsquant/_core.pyi").read_text())
    declared = {node.name for node in tree.body
                if isinstance(node, (ast.ClassDef, ast.FunctionDef))}
    declared.update(node.target.id for node in tree.body
                    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name))
    declared.update(target.id for node in tree.body if isinstance(node, ast.Assign)
                    for target in node.targets if isinstance(target, ast.Name))
    assert {name for name in dir(_core) if not name.startswith("_")} <= declared


def test_stub_signatures_and_readonly_profile():
    tree = ast.parse((REPO / "python/lsquant/_core.pyi").read_text())
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    ranged = classes["RangedRandomPhaseStates"]
    constructor = next(node for node in ranged.body if isinstance(node, ast.FunctionDef))
    assert [arg.arg for arg in constructor.args.args] == [
        "self", "dim", "num_states", "seed", "row_begin", "row_end"]
    block_methods = [node for node in classes["BlockLocalStates"].body
                     if isinstance(node, ast.FunctionDef)]
    assert ["self", "n1", "n2", "n3", "block_size"] in [
        [arg.arg for arg in node.args.args] for node in block_methods if node.name == "__init__"]
    index = next(node for node in block_methods if node.name == "block_index")
    assert [arg.arg for arg in index.args.args] == ["m1", "m2", "m3", "n1", "n2"]
    assert any(isinstance(d, ast.Name) and d.id == "staticmethod" for d in index.decorator_list)
    profile = next(node for node in classes["ChebyshevMoments"].body
                   if isinstance(node, ast.FunctionDef) and node.name == "build_profile")
    assert [d.id for d in profile.decorator_list] == ["property"]
    assert not any(isinstance(node, ast.FunctionDef) and node.name == "build_profile"
                   and any(isinstance(d, ast.Attribute) and d.attr == "setter"
                           for d in node.decorator_list)
                   for node in classes["ChebyshevMoments"].body)
    sync = next(node for node in tree.body
                if isinstance(node, ast.FunctionDef) and node.name == "synchronize_cuda")
    assert not sync.args.args
