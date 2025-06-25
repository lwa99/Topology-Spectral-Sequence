"""
convert_to_seqsee_html.py

Builds your spectral sequence exactly as in Test.py, then calls
main.generate_html(data) directly to emit the SVG+HTML output,
without ever serializing to intermediate JSON.
"""

import sys
from sympy.abc import symbols
from sympy import GF
from spectral_sequence import SpectralSequence
from element import HomoElem
from utilities import Vector

# Import the rendering functions directly from main.py
# (assumes convert_to_seqsee_html.py sits alongside main.py)
from seqsee_main import generate_html, compute_chart_dimensions, generate_css_styles


def build_data_dict():
    # --- replicate Test.py setup ---
    a, t = symbols("a t")
    ss = SpectralSequence(
        GF(3),
        [a, t],
        [[3, 0], [0, 2]],
        [[1, 0],
         [-1, 1]]
    )

    # impose relation a^2 = 0 and build through page 4
    ss.kill(a**2)
    ss.add_page({a: 0, t: 0})
    ss.add_page({a: 0, t: 0})
    ss.add_page({t: a})
    p4 = ss.add_page()   # page 4

    # now assemble the JSON-style dict that main.generate_html expects
    data = {
        "header": {
            "metadata": {
                "title": f"Page {p4.page_num} of your spectral sequence"
            },
            # leave chart bounds null so compute_chart_dimensions fills them in
            "chart": {
                "width":  {"min": None, "max": None},
                "height": {"min": None, "max": None},
                # scale, nodeSize, etc. will all be pulled from schema defaults
                "scale": None
            },
            # you can also set aliases/color overrides here if you like
            "aliases": {}
        },
        "nodes": {},
        "edges": []
    }

    # collect all surviving basis elements as nodes
    for (x_deg, y_deg), module in p4.items():
        for col in range(module.sp_basis.cols):
            elem = HomoElem(
                p4,
                abs_coordinate=module.sp_basis.col(col),
                abs_bideg= Vector(x_deg, y_deg)
            )
            node_id = f"n_{x_deg}_{y_deg}_{col}"
            data["nodes"][node_id] = {
                "x": int(x_deg),
                "y": int(y_deg),
                "label": f"${elem.poly.as_expr()}$",
                # you can add styling attributes here, e.g.
                # "attributes": [{"color": "h0", "size": 1.2}]
            }

    # (Optional) populate data["edges"] by inspecting ss.differentials or extensions

    return data

def main():
    if len(sys.argv) != 2:
        print(sys.argv)
        print("Usage: --- <output.html>")
        sys.exit(1)
    out_path = sys.argv[1]

    # Build the data dict
    data = build_data_dict()

    # Prepare CSS defaults & compute chart dimensions
    generate_css_styles(data)            # fills in global_css
    compute_chart_dimensions(data)       # sets width/height bounds

    # Render the full HTML
    html = generate_html(data)

    # Write it out
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()