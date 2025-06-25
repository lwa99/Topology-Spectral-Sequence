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
from jinja2 import Environment, FileSystemLoader

# Scale (pixels per unit)
scale = None


class CssStyle:
    """A simple CSS style container."""
    def __init__(self, *args, **kwargs):
        self._styles = {}
        for a in args:
            self.append(a)
        for name, value in kwargs.items():
            self._styles[name] = value

    def append(self, other):
        if isinstance(other, dict):
            self._styles.update(other)
        elif isinstance(other, CssStyle):
            self._styles.update(other._styles)

    def generate(self):
        lines = []
        for selector, props in self._styles.items():
            lines.append(f"{selector} {{")
            for k, v in props.items():
                lines.append(f"  {k}: {v};")
            lines.append("}\n")
        return "\n".join(lines)

# Global CSS accumulator
global_css = CssStyle()

def load_template():
    env = Environment(loader=FileSystemLoader(searchpath=sys.path[0]))
    return env.get_template("template.html.jinja")


def generate_css_styles(data):
    """Populate CSS classes from data["aliases"]."""
    global global_css
    global_css = CssStyle()
    for name, attrs in data.get("aliases", {}).items():
        selector = f".{name}"
        global_css._styles[selector] = attrs


def compute_chart_dimensions(data):
    """Automatically fill width/height in header.chart."""
    nodes = data.get("nodes", {})
    xs = [node["x"] for node in nodes.values()]
    ys = [node["y"] for node in nodes.values()]
    w_min, w_max = (min(xs), max(xs)) if xs else (0, 2)
    h_min, h_max = (min(ys), max(ys)) if ys else (0, 2)
    data["header"]["chart"]["width"] = {"min": w_min, "max": w_max}
    data["header"]["chart"]["height"] = {"min": h_min, "max": h_max}


def calculate_absolute_positions(data):
    """Compute absoluteX/absoluteY for each node."""
    nodes = data.get("nodes", {})
    for node_id, node in nodes.items():
        node["absoluteX"] = node["x"]
        node["absoluteY"] = node["y"]


def generate_nodes_svg(data):
    svg = ['<g id="nodes">']
    for nid, node in data.get("nodes", {}).items():
        cx = node["absoluteX"] * scale
        cy = node["absoluteY"] * scale
        style = ''
        svg.append(f'<circle id="{nid}" cx="{cx}" cy="{cy}" r="{scale}" {style}/>' )
    svg.append('</g>')
    return "\n".join(svg)


def generate_edges_svg(data):
    svg = ['<g id="edges">']
    # No edges by default
    svg.append('</g>')
    return "\n".join(svg)


def generate_svg(data):
    calculate_absolute_positions(data)
    return generate_edges_svg(data) + generate_nodes_svg(data)


def generate_html(data):
    global scale
    scale = data.get("header", {}).get("chart", {}).get("scale", 50)
    generate_css_styles(data)
    compute_chart_dimensions(data)
    svg_content = generate_svg(data)
    template = load_template()
    return template.render(data=data, css_styles=global_css.generate(), svg_content=svg_content)
