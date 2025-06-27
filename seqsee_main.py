"""
https://github.com/JoeyBF/SeqSee/tree/master

MIT License

Copyright (c) 2024 Joey Beauvais-Feisthauer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
import copy
import math
from collections import defaultdict
from jinja2 import Environment, FileSystemLoader

# Scale (pixels per unit)
scale = None

class CssStyle:
    """A list of CSS styles stored as a dict."""
    def __init__(self, *args, **kwargs):
        self._styles = {}
        for a in args:
            self.append(a)
        for name, value in kwargs.items():
            self._styles[name] = value

    def __getitem__(self, key):
        return self._styles[key]

    def keys(self):
        return self._styles.keys()

    def items(self):
        return self._styles.items()

    def append(self, other):
        if isinstance(other, dict):
            self._styles.update(other)
        elif isinstance(other, CssStyle):
            self._styles.update(other._styles)
        else:
            raise TypeError("Unsupported style type")

    def generate(self, parent="", indent=4):
        subnodes, stylenodes = [], []
        lines = []
        for name, value in self._styles.items():
            if isinstance(value, dict) or isinstance(value, CssStyle):
                subnodes.append((name, CssStyle(value) if isinstance(value, dict) else value))
            else:
                stylenodes.append((name, value))
        if stylenodes:
            lines.append(f"{parent} {{")
            for k, v in stylenodes:
                val = v if isinstance(v, str) else f"{v}px"
                lines.append(" "*indent + f"{k}: {val};")
            lines.append("}\n")
        for subname, substyle in subnodes:
            nested = substyle.generate(parent=(parent + " " + subname).strip(), indent=indent)
            if isinstance(nested, str):
                lines.extend(nested.split("\n"))
            else:
                lines.extend(nested)
        return "\n".join(lines)

# Global CSS accumulator
global_css = CssStyle()

# Load template
def load_template():
    env = Environment(loader=FileSystemLoader(searchpath="."))
    return env.get_template("template.html.jinja")

# Style helpers
def style_and_aliases_from_attributes(attributes):
    new_style = CssStyle()
    aliases = []
    for attr in attributes:
        if isinstance(attr, dict):
            for k, v in attr.items():
                if k == "color":
                    new_style.append({"fill": v, "stroke": v})
                elif k == "size":
                    new_style.append({"r": scale * float(v)})
                elif k == "thickness":
                    new_style.append({"stroke-width": scale * float(v)})
                else:
                    new_style.append({k: v})
        elif isinstance(attr, str):
            aliases.append(attr)
    return new_style, aliases

# Generate CSS classes
def generate_css_styles(data):
    global global_css
    global_css = CssStyle()
    # default node/edge
    global_css.append({".defaultNode": {"stroke-width": 0}})
    global_css.append({".defaultEdge": {"stroke": "#000", "stroke-width": "1px"}})

# Chart dimensions autodetect
def compute_chart_dimensions(data):
    nodes = data.get("nodes", {})
    xs = [n.get("x", 0) for n in nodes.values()]
    ys = [n.get("y", 0) for n in nodes.values()]
    w_min, w_max = (min(xs), max(xs)) if xs else (0, 2)
    h_min, h_max = (min(ys), max(ys)) if ys else (0, 2)
    hdr = data.setdefault("header", {}).setdefault("chart", {})
    hdr["width"] = {"min": w_min, "max": w_max}
    hdr["height"] = {"min": h_min, "max": h_max}

# Position nodes (simple grid)
def calculate_absolute_positions(data):
    for nid, node in data.get("nodes", {}).items():
        node["absoluteX"] = node.get("x", 0)
        node["absoluteY"] = node.get("y", 0)

# SVG generation
def generate_nodes_svg(data):
    out = '<g id="nodes">\n'
    for nid, node in data.get("nodes", {}).items():
        cx = node["absoluteX"] * scale
        cy = node["absoluteY"] * scale
        style, aliases = style_and_aliases_from_attributes(node.get("attributes", []))
        s = style.generate(indent=0).replace("\n", " ").strip()
        cls = " ".join(aliases)
        out += f'<circle id="{nid}" class="defaultNode {cls}" cx="{cx}" cy="{cy}" {s}></circle>\n'
    out += '</g>\n'
    return out


def generate_edges_svg(data):
    out = '<g id="edges">\n'
    for edge in data.get("edges", []):
        src = data["nodes"][edge.get("source")]
        if "target" in edge:
            dst = data["nodes"][edge.get("target")]
            x1, y1 = src["absoluteX"] * scale, src["absoluteY"] * scale
            x2, y2 = dst["absoluteX"] * scale, dst["absoluteY"] * scale
            out += f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="defaultEdge"></line>\n'
    out += '</g>\n'
    return out


def generate_svg(data):
    calculate_absolute_positions(data)
    return generate_edges_svg(data) + generate_nodes_svg(data)

# Main HTML renderer
def generate_html(data):
    global scale
    scale = data.get("header", {}).get("chart", {}).get("scale") or 50
    generate_css_styles(data)
    compute_chart_dimensions(data)
    svg_content = generate_svg(data)
    template = load_template()
    return template.render(data=data, css_styles=global_css.generate(), svg_content=svg_content)

# # CLI entry
if __name__ == "__main__":
    # Expect test2.py to write data.py automatically
    if len(sys.argv) != 2:
        print("Usage: python seqsee_main.py <output.html>")
        sys.exit(1)
    output_path = sys.argv[1]
    # Dynamically import test2 module to get data
    import test2
    data = test2.build_data_dict(0, 5, 0, 5)  # use default window or modify as needed
    html = generate_html(data)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Wrote {output_path}")
