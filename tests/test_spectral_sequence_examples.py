from __future__ import annotations

import builtins
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TESTS = ROOT / "tests"
for path in (SRC, ROOT, TESTS):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from matrices import IV  # noqa: E402
from element import HomoElem  # noqa: E402
from spectral_sequence_examples import all_examples, build_all_examples  # noqa: E402


@pytest.mark.parametrize("record", all_examples(), ids=lambda r: r.name)
def test_examples_build_and_metadata(record):
    data = record.build()
    assert data["name"]
    assert data["ss"] is not None
    assert data["field"]
    assert data["pages"]
    assert data["expected_stable_page"] in data["pages"]


@pytest.mark.parametrize(
    "example_key,page_key,src_name,tgt_name",
    [
        ("example_2_hopf_fibration", "p_2", "u", "a"),
        ("example_3_universal_circle_bundle", "p_2", "u", "c"),
        ("example_4_quaternionic_hopf", "p_4", "u", "a"),
        ("example_5_path_space_kz3", "p_3", "t", "a"),
    ],
)
def test_declared_nonzero_differentials_are_present(example_key, page_key, src_name, tgt_name):
    data = build_all_examples()[example_key]
    page = data["pages"][page_key]
    gens = data["generators"]
    src = HomoElem(page, gens[src_name])
    tgt = HomoElem(page, gens[tgt_name])
    assert src in page.d.info
    assert page.d.info[src] == tgt


@pytest.mark.parametrize(
    "record,should_succeed",
    [
        (all_examples()[0], True),
        (all_examples()[1], True),
        (all_examples()[2], True),
        (all_examples()[3], True),
        (all_examples()[4], True),
    ],
    ids=lambda x: x.name if hasattr(x, "name") else str(x),
)
def test_expected_stable_page_query_at_unit(record, should_succeed):
    """
    Regression check at (0,0):
    non-first-quadrant targets are treated as zero modules, so these queries should succeed.
    """
    data = record.build()
    stable_page = data["pages"][data["expected_stable_page"]]
    if should_succeed:
        _ = stable_page[IV([0, 0])]
    else:
        with pytest.raises(AssertionError):
            _ = stable_page[IV([0, 0])]


@pytest.mark.parametrize("record", all_examples(), ids=lambda r: r.name)
def test_higher_page_queries_with_auto_zero_completion(record, monkeypatch):
    """
    Query several first-quadrant bidegrees on the expected stable page.
    Missing differential values are auto-filled by mocking input() to "0".
    """
    monkeypatch.setattr(builtins, "input", lambda _prompt="": "")

    data = record.build()
    stable_page = data["pages"][data["expected_stable_page"]]

    queries = [
        IV([0, 1]),
        IV([1, 1]),
        IV([2, 0]),
    ]
    for bideg in queries:
        module = stable_page[bideg]
        assert module is not None


def test_higher_page_query_invokes_completion_input_when_needed(monkeypatch):
    """
    Example 5 at p_4(0,4) needs extra differential info.
    Verify pytest can handle this by mocking input() and completing the query.
    """
    calls = {"n": 0}

    def fake_input(_prompt=""):
        calls["n"] += 1
        return ""

    monkeypatch.setattr(builtins, "input", fake_input)

    data = build_all_examples()["example_5_path_space_kz3"]
    page = data["pages"]["p_4"]
    module = page[IV([0, 4])]

    assert module is not None
    assert calls["n"] >= 1


@pytest.mark.parametrize("record", all_examples(), ids=lambda r: r.name)
def test_expected_stable_page_grid_queries(record, monkeypatch):
    """
    Stronger regression: query a 10x10 first-quadrant grid on each example's expected stable page.
    """
    monkeypatch.setattr(builtins, "input", lambda _prompt="": "")

    data = record.build()
    stable_page = data["pages"][data["expected_stable_page"]]

    for p in range(10):
        for q in range(10):
            module = stable_page[IV([p, q])]
            assert module is not None
