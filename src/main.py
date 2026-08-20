#!/usr/bin/env python3
"""
Build the sample algebraic spectral sequence and export an HTML chart.

The sample is the spectral sequence over Z[a, t]/(a^2) used in the project
write-up. The script takes explicit chart bounds and scans the corresponding
rectangle of bidegrees before rendering through the bundled SeqSee frontend.
"""

import os
import sys

# Support running as a script: `python src/main.py ...`
if __package__ is None or __package__ == "":
    repo_root = os.path.dirname(os.path.dirname(__file__))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

from sympy import ZZ
from sympy.abc import symbols

from src.spectral_sequence import SpectralSequence
from src.element import HomoElem, Bidegree
from seqsee_main import process_data


def build_data_dict(min_x, max_x, min_y, max_y):
    a, t = symbols("a t")
    ss = SpectralSequence(
        ZZ,
        [a, t],
        [[3, 0], [0, 2]],
        [[1, 0],
         [-1, 1]],
    )

    ss.kill(a**2)
    ss.add_page({a: 0, t: 0})
    ss.add_page({a: 0, t: 0})
    # d3(a)=0 is mathematically redundant because d3(t)=a and d3^2=0,
    # but the current implementation does not infer d^2=0 automatically.
    ss.add_page({t: a, a: 0})
    p4 = ss.add_page()   # page 4

    data = {
        "$schema": "https://raw.githubusercontent.com/JoeyBF/SeqSee/refs/heads/master/seqsee/input_schema.json",
        "header": {
            "metadata": {
                "htmltitle": "Sample algebraic spectral sequence over Z",
            },
            "aliases": {
                "attributes": {
                    "defaultNode": [{"color": "gray"}],
                    "defaultEdge": [{"color": "gray"}],
                },
                "colors": {
                    "gray": "#666",
                },
            },
        },
        "nodes": {},
        "edges": [],
    }

    # Walk the given rectangle of bidegrees.
    for x in range(min_x, max_x + 1):
        for y in range(min_y, max_y + 1):
            try:
                module = p4[x, y]
            except KeyError:
                # No module at this bidegree -> skip.
                continue

            structural = module.get_structural_information()
            if structural is None:
                continue
            gens, torsion = structural
            if callable(gens):
                gens = gens()

            # Keep only genuinely nontrivial quotient classes.
            # `classify(g) == 1` excludes generators already trivial by relations.
            nontrivial_coords = [
                g for g, t_info in zip(gens, torsion)
                if (not ss.domain.is_unit(t_info)) and module.classify(g) == 1
            ]
            if len(nontrivial_coords) == 0:
                continue

            for col, abs_coord in enumerate(nontrivial_coords):
                elem = HomoElem(
                    p4,
                    abs_coordinate=abs_coord,
                    abs_bideg=Bidegree([x, y]),
                )
                node_id = f"n_{x}_{y}_{col}"

                from sympy import latex
                data["nodes"][node_id] = {
                    "x": x,
                    "y": y,
                    "label": f"${latex(elem.poly.as_expr())}$",
                }

    print(f"Debug: found {len(data['nodes'])} total nodes")
    return data


def main():
    if len(sys.argv) != 6:
        print("Usage: python src/main.py <min_x> <max_x> <min_y> <max_y> <out.html>")
        sys.exit(1)

    min_x, max_x, min_y, max_y = map(int, sys.argv[1:5])
    out_path = sys.argv[5]

    data = build_data_dict(min_x, max_x, min_y, max_y)
    process_data(data, out_path)

    print(f"Wrote {out_path} covering x in [{min_x},{max_x}], y in [{min_y},{max_y}]")


if __name__ == "__main__":
    main()
