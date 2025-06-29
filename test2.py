#!/usr/bin/env python3
"""
convert_to_seqsee_html.py

Builds your spectral sequence exactly as in Test.py, then calls
main.generate_html(data) directly to emit the SVG+HTML output,
without ever serializing to intermediate JSON.

Now takes explicit chart bounds (min_x, max_x, min_y, max_y)
and iterates over that grid, pulling out each module’s basis.
"""

import sys
from sympy.abc import symbols
from sympy import GF
from spectral_sequence import SpectralSequence
from element import HomoElem, Bidegree

from seqsee_main import process_data

SCHEMA = {
  "$schema": "https://raw.githubusercontent.com/JoeyBF/SeqSee/refs/heads/master/seqsee/input_schema.json",
  "header": {
    "metadata": {
      "htmltitle": "First few C-motivic stable stems",
      "title": "First few $\\mathbb{C}$-motivic stable stems",
      "authors": ["Joey Beauvais-Feisthauer", "Daniel C. Isaksen"]
    },
    "aliases": {
      "attributes": {
        "defaultNode": [ {"color": "gray"         } ],
        "defaultEdge": [ {"color": "gray"         } ],
        "tau1"       : [ {"color": "tau1color"    } ],
        "tau1extn"   : [ {"color": "tau1extncolor"} ]
      },
      "colors": {
        "gray"         : "#666"   ,
        "tau1color"    : "#DD0000",
        "tau1extncolor": "magenta"
      }
    }
  },
  "nodes": {
    "1"   : {"x": 0, "y": 0, "label": "$1$ (0)"}                              ,
    "h0"  : {"x": 0, "y": 1, "label": "$h_0$ (0)"}                            ,
    "h0^2": {"x": 0, "y": 2, "label": "$h_0^2$ (0)"}                          ,
    "h0^3": {"x": 0, "y": 3, "label": "$h_0^3$ (0)"}                          ,
    "h1"  : {"x": 1, "y": 1, "label": "$h_1$ (1)"}                            ,
    "h1^2": {"x": 2, "y": 2, "label": "$h_1^2$ (2)"}                          ,
    "h2"  : {"x": 3, "y": 1, "label": "$h_2$ (2)"}                            ,
    "h0h2": {"x": 3, "y": 2, "label": "$h_0 h_2$ (2)"}                        ,
    "h1^3": {"x": 3, "y": 3, "label": "$h_1^3$ (3)"}                          ,
    "h1^4": { "x": 4, "y": 4, "label": "$h_1^4$ (4)", "attributes": ["tau1"] }
  },
  "edges": [
  ]
}


def build_data_dict(min_x, max_x, min_y, max_y):
    # --- replicate Test.py setup ---
    a, t = symbols("a t")
    ss = SpectralSequence(
        GF(3),
        [a, t],
        [[3, 0], [0, 2]],
        [[1, 0],
         [-1, 1]]
    )
    ss.kill(a**2)
    ss.add_page({a: 0, t: 0})
    ss.add_page({a: 0, t: 0})
    ss.add_page({t: a})
    p4 = ss.add_page()   # page 4

    data = {
      "$schema": "https://raw.githubusercontent.com/JoeyBF/SeqSee/refs/heads/master/seqsee/input_schema.json",
      "header": {
        "metadata": {
          "htmltitle": "First few C-motivic stable stems",
          "title": "First few $\\mathbb{C}$-motivic stable stems",
          "authors": ["Joey Beauvais-Feisthauer", "Daniel C. Isaksen"]
        },
        "aliases": {
          "attributes": {
            "defaultNode": [ {"color": "gray"         } ],
            "defaultEdge": [ {"color": "gray"         } ],
            "tau1"       : [ {"color": "tau1color"    } ],
            "tau1extn"   : [ {"color": "tau1extncolor"} ]
          },
          "colors": {
            "gray"         : "#666"   ,
            "tau1color"    : "#DD0000",
            "tau1extncolor": "magenta"
          }
        }
      },
      "nodes": {
      },
      "edges": [
      ]
    }

    # Walk the given rectangle of bidegrees
    for x in range(min_x, max_x + 1):
        for y in range(min_y, max_y + 1):
            try:
                module = p4[x, y]
            except KeyError:
                # no module at this bidegree → skip
                continue

            for col in range(module.sp_basis.cols):
                elem = HomoElem(
                    p4,
                    abs_coordinate=module.sp_basis.col(col),
                    abs_bideg= Bidegree([x, y])
                )
                node_id = f"n_{x}_{y}_{col}"

                from sympy import latex
                data["nodes"][node_id] = {
                    "x": x,
                    "y": y,
                    "label": f"${latex(elem.poly.as_expr())}$"
                }
    print(f"Debug: found {len(data['nodes'])} total nodes")
    return data

def main():
    global global_css

    if len(sys.argv) != 6:
        print("Usage: convert_to_seqsee_html.py <min_x> <max_x> <min_y> <max_y> <out.html>")
        sys.exit(1)

    min_x, max_x, min_y, max_y = map(int, sys.argv[1:5])
    out_path = sys.argv[5]

    data = build_data_dict(min_x, max_x, min_y, max_y)
    process_data(data, out_path)

    print(f"Wrote {out_path} covering x ∈ [{min_x},{max_x}], y ∈ [{min_y},{max_y}]")


if __name__ == "__main__":
    main()
