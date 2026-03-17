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


def _rank_mod(rows, mod):
    if not rows:
        return 0

    m = [r[:] for r in rows]
    nrows = len(m)
    ncols = len(m[0]) if nrows else 0
    rank = 0

    for c in range(ncols):
        pivot = None
        for i in range(rank, nrows):
            if m[i][c] % mod != 0:
                pivot = i
                break
        if pivot is None:
            continue

        m[rank], m[pivot] = m[pivot], m[rank]
        inv = pow(m[rank][c], -1, mod)
        m[rank] = [(v * inv) % mod for v in m[rank]]

        for i in range(nrows):
            if i == rank:
                continue
            factor = m[i][c] % mod
            if factor == 0:
                continue
            m[i] = [(m[i][j] - factor * m[rank][j]) % mod for j in range(ncols)]

        rank += 1
        if rank == nrows:
            break

    return rank


def _module_dim_over_field(module):
    if module.S is None:
        return 0

    mod = int(module.domain.mod)
    s_rows = [[int(x) % mod for x in row] for row in module.S.to_list()]
    rank_s = _rank_mod(s_rows, mod)

    if module.R is None:
        rank_r = 0
    else:
        r_rows = [[int(x) % mod for x in row] for row in module.R.to_list()]
        rank_r = _rank_mod(r_rows, mod)

    return rank_s - rank_r


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


def test_leibniz_derives_t_squared_on_example_5():
    data = build_all_examples()["example_5_path_space_kz3"]
    p3 = data["pages"]["p_3"]
    a = data["generators"]["a"]
    t = data["generators"]["t"]

    t2 = HomoElem(p3, t**2)
    assert t2 not in p3.d.info

    p3.d.extend_by_forward_leibniz(IV([0, 4]))

    assert t2 in p3.d.info
    assert p3.d.info[t2] == HomoElem(p3, 2 * a * t)


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


def test_higher_page_query_uses_leibniz_before_prompt(monkeypatch):
    """
    Example 5 at p_4(0,4): d(t^2) should be derived from d(t) by Leibniz,
    so no interactive prompt is needed for this query.
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
    assert calls["n"] == 0


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


@pytest.mark.parametrize(
    "example_key,stable_page_key,expected_nonzero",
    [
        (
            "example_1_trivial_product",
            "p_2",
            {(0, 0): 1, (0, 1): 1, (2, 0): 1, (2, 1): 1},
        ),
        (
            "example_2_hopf_fibration",
            "p_3",
            {(0, 0): 1, (2, 1): 1},
        ),
        (
            "example_3_universal_circle_bundle",
            "p_3",
            {(0, 0): 1},
        ),
        (
            "example_4_quaternionic_hopf",
            "p_5",
            {(0, 0): 1, (4, 3): 1},
        ),
        (
            "example_5_path_space_kz3",
            "p_4",
            {(0, 0): 1, (0, 6): 1, (3, 4): 1},
        ),
    ],
)
def test_guide_ground_truth_dimensions_10x10(example_key, stable_page_key, expected_nonzero, monkeypatch):
    """
    Ground-truth checks from the guide's "what the stable page should look like" sections,
    compared on a 10x10 first-quadrant window.
    """
    monkeypatch.setattr(builtins, "input", lambda _prompt="": "")

    data = build_all_examples()[example_key]
    stable_page = data["pages"][stable_page_key]

    for p in range(10):
        for q in range(10):
            module = stable_page[IV([p, q])]
            got = _module_dim_over_field(module)
            expected = expected_nonzero.get((p, q), 0)
            assert got == expected, (
                f"{example_key} {stable_page_key} wrong dim at ({p},{q}): "
                f"expected {expected}, got {got}"
            )
