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
import copy
import json
import jsonschema
import math
import sys
import os
from collections import defaultdict
from jsonschema.exceptions import ValidationError
from jinja2 import Environment, FileSystemLoader

# The distance between successive x or y coordinates. Units are in pixels.
# This will be fixed throughout the html file, but zooming is implemented
# through a transformation matrix applied to the <g> element.
scale = None

# Lifted/adapted from MIT-licensed https://github.com/slacy/pyssed/
class CssStyle:
    """A list of CSS styles, but stored as a dict.
    Can contain nested styles."""

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
        self._styles = (self + other)._styles

    def __add__(self, other):
        summed = copy.deepcopy(self)
        if isinstance(other, str):
            sel, val = other.split(":", 1)
            summed._styles[sel] = val
        elif isinstance(other, dict):
            summed._styles.update(other)
        elif isinstance(other, CssStyle):
            summed._styles.update(other._styles)
        else:
            raise TypeError("Unsupported type for style addition")
        return summed

    def generate(self, parent="", indent=4):
        subnodes, stylenodes, result = [], [], []
        for name, value in self.items():
            if isinstance(value, dict):
                subnodes.append((name, CssStyle(value)))
            elif isinstance(value, CssStyle):
                subnodes.append((name, value))
            elif isinstance(value, (str, int, float)):
                stylenodes.append((name, value))
            else:
                raise TypeError("Bad style node type")

        if stylenodes:
            result.append(parent.strip() + " {")
            for k, v in stylenodes:
                attr = k.strip(" ;:")
                val = v if isinstance(v, str) else f"{v}px"
                result.append(" " * indent + f"{attr}: {val};")
            result.append("}")
            result.append("")

        for subname, substyle in subnodes:
            result += substyle.generate(parent=(parent + " " + subname).strip(), indent=indent)

        return "\n".join(result) if parent == "" else result

# Global CSS accumulator
global_css = CssStyle()

# Load JSON schema for input validation
def load_schema():
    schema_path = os.path.join(os.path.dirname(__file__), "input_schema.json")
    with open(schema_path, "r", encoding="utf-8") as f:
        return json.load(f)

schema = load_schema()

# Load Jinja2 template
def load_template():
    env = Environment(loader=FileSystemLoader(searchpath=os.path.dirname(__file__)))
    return env.get_template("template.html.jinja")

# Schema-default helpers
def get_schema_default(data, path):
    default_node = schema
    for key in path:
        default_node = default_node["properties"][key]
    return default_node.get("default")

def get_value_or_schema_default(data, path):
    try:
        cur = data
        for key in path:
            cur = cur[key]
        return cur
    except (KeyError, TypeError):
        return get_schema_default(data, path)

# CSS-safe name
def cssify_name(name):
    if str(name).isnumeric():
        name = "n" + name
    return "." + name

# Attributes → style & aliases
def style_and_aliases_from_attributes(attributes):
    new_style = CssStyle()
    aliases = []
    for attr in attributes:
        if isinstance(attr, dict):
            for k, v in attr.items():
                if k == "color":
                    key = cssify_name(v)
                    if key in global_css.keys():
                        new_style += global_css[key]
                    else:
                        new_style += {"fill": v, "stroke": v}
                elif k == "size":
                    new_style += {"r": scale * float(v)}
                elif k == "thickness":
                    new_style += {"stroke-width": scale * float(v)}
                elif k == "arrowTip":
                    if v == "none":
                        new_style += {"marker-end": "none"}
                    else:
                        new_style += {"marker-end": f"url(#arrow-{v})"}
                elif k == "pattern":
                    if v == "solid": new_style += {"stroke-dasharray": "none"}
                    elif v == "dashed": new_style += {"stroke-dasharray": "5, 5"}
                    elif v == "dotted": new_style += {"stroke-dasharray": "0, 2", "stroke-linecap": "round"}
                else:
                    new_style += {k: v}
        elif isinstance(attr, str):
            aliases.append(cssify_name(attr).lstrip("."))
    return new_style, aliases

# Combine style + alias classes
def generate_style(style, aliases):
    st = copy.deepcopy(style)
    for alias in aliases:
        st.append(global_css[cssify_name(alias)])
    return st

# Ensure nested JSON path exists
def ensure_json_path_is_defined(data, path):
    cur = data
    for key in path:
        cur = cur.setdefault(key, {})

# Chart dimensions autodetection
def compute_chart_dimensions(data):
    nodes = data.get("nodes", {})
    def bounds(dim, coord, default):
        ensure_json_path_is_defined(data, ["header","chart",dim])
        mn = data["header"]["chart"][dim].get("min")
        if mn is None:
            vmin = min((n[coord] for n in nodes.values()), default=default)
            data["header"]["chart"][dim]["min"] = 2*((vmin//2)-1)
        mx = data["header"]["chart"][dim].get("max")
        if mx is None:
            vmax = max((n[coord] for n in nodes.values()), default=default)
            data["header"]["chart"][dim]["max"] = 2*((vmax//2)+1)
    bounds("width","x",0)
    bounds("height","y",0)

# Positioning with spacing, slope, size
def calculate_absolute_positions(data):
    nodes_by = defaultdict(list)
    for nid,node in data.get("nodes",{}).items():
        nodes_by[(node["x"],node["y"])].append(nid)
    default_pos = schema["properties"]["nodes"]["additionalProperties"]["properties"]["position"]["default"]
    for bd, lst in nodes_by.items():
        nodes_by[bd] = sorted(lst, key=lambda id: data["nodes"][id].get("position", default_pos))
    size = get_value_or_schema_default(data, ["header","chart","nodeSize"])
    spacing = get_value_or_schema_default(data, ["header","chart","nodeSpacing"])
    slope = get_value_or_schema_default(data, ["header","chart","nodeSlope"])
    dist = spacing + 2*size
    theta = math.atan(slope) if slope is not None else math.pi/2
    for (x,y), lst in nodes_by.items():
        n=len(lst)
        total = (n-1)*dist
        for i,nid in enumerate(lst):
            off = -total/2 + i*dist
            data["nodes"][nid]["absoluteX"] = x + off*math.cos(theta)
            data["nodes"][nid]["absoluteY"] = y + off*math.sin(theta)

# SVG generation
def generate_nodes_svg(data):
    out = '<g id="nodes-group">\n'
    for nid,node in data.get("nodes",{}).items():
        cx,cy = node["absoluteX"]*scale, node["absoluteY"]*scale
        style,aliases = style_and_aliases_from_attributes(node.get("attributes",[]))
        s = style.generate(indent=0).replace("\n"," ").strip(" {}")
        if s: s=f'style="{s}"'
        cls=" ".join(aliases)
        lbl=node.get("label","")
        out+=f'<circle id="{nid}" class="defaultNode {cls}" cx="{cx}" cy="{cy}" {s} data-label="{lbl}"></circle>\n'
    out+='</g>\n'
    return out

def generate_edges_svg(data):
    out = '<g id="edges-group">\n'
    for edge in data.get("edges",[]):
        src=data["nodes"][edge["source"]]
        if "target" in edge:
            tgt=data["nodes"][edge["target"]]
            tx,ty=tgt["absoluteX"]*scale,tgt["absoluteY"]*scale
        elif "offset" in edge:
            tx=(src["absoluteX"]+edge["offset"]["x"])*scale
            ty=(src["absoluteY"]+edge["offset"]["y"])*scale
        else: raise NotImplementedError
        x1,y1=src["absoluteX"]*scale,src["absoluteY"]*scale
        style,aliases=style_and_aliases_from_attributes(edge.get("attributes",[]))
        s=style.generate(indent=0).replace("\n"," ").strip(" {}")
        cls=" ".join(aliases)
        if edge.get("bezier"):
            cp=edge["bezier"]
            if len(cp)==1:
                c0x,c0y=cp[0]["x"]*scale+x1,cp[0]["y"]*scale+y1
                d=f"Q {c0x} {c0y} {tx} {ty}"
            else:
                c0x,c0y=cp[0]["x"]*scale+x1,cp[0]["y"]*scale+y1
                c1x,c1y=cp[1]["x"]*scale+tx,cp[1]["y"]*scale+ty
                d=f"C {c0x} {c0y} {c1x} {c1y} {tx} {ty}"
            line=f'<path d="M {x1} {y1} {d}" class="{cls}" style="fill: none;{s}"></path>\n'
        else:
            line=f'<line x1="{x1}" y1="{y1}" x2="{tx}" y2="{ty}" class="defaultEdge {cls}" style="{s}"></line>\n'
        out+=line.replace(' style=""','')
    out+='</g>\n'
    return out

# Wrapper
def generate_svg(data):
    calculate_absolute_positions(data)
    return generate_edges_svg(data) + generate_nodes_svg(data)

# Populate CSS classes
def generate_css_styles(data):
    global global_css
    base = ["header","aliases"]
    # colors
    cp=base+["colors"]
    colors={
        "backgroundColor":get_schema_default(data,cp+["backgroundColor"]),
        "textColor":get_schema_default(data,cp+["textColor"]),
        "borderColor":get_schema_default(data,cp+["borderColor"]),
    }
    colors.update(get_value_or_schema_default(data,cp))
    # attributes
    ap=base+["attributes"]
    attr_alias={
        "grid":get_schema_default(data,ap+["grid"]),
        "defaultNode":get_schema_default(data,ap+["defaultNode"]),
        "defaultEdge":get_schema_default(data,ap+["defaultEdge"]),
    }
    user_attrs=get_value_or_schema_default(data,ap)
    for k,v in user_attrs.items(): attr_alias[k]=attr_alias.get(k,[])+v
    # generate color styles
    for nm,val in colors.items(): global_css += {cssify_name(nm):{"fill":val,"stroke":val}}
    global_css += {".backgroundStyle":{"background-color":colors["backgroundColor"],"fill":colors["backgroundColor"]},
                   "#tooltip":{"color":colors["textColor"],"border-color":colors["borderColor"]},
                   ".axis":{"stroke":colors["borderColor"],"stroke-width":"2px"},
                   ".tick text, .katex":{"color":colors["textColor"],"fill":"currentColor"}}
    size=get_value_or_schema_default(data,["header","chart","nodeSize"])
    global_css += {"circle":{"stroke-width":0,"r":scale*size}}
    for alias,attrs in attr_alias.items(): st,als=style_and_aliases_from_attributes(attrs); global_css += {cssify_name(alias):generate_style(st,als)}

# Main HTML generation
def generate_html(data):
    generate_css_styles(data)
    compute_chart_dimensions(data)
    svg = generate_svg(data)
    tpl = load_template()
    return tpl.render(data=data, spacing=scale, css_styles=global_css.generate(), static_svg_content=svg)

# CLI entry
def process_json(input_file, output_file):
    global global_css, scale
    with open(input_file,'r',encoding='utf-8') as f: data=json.load(f)
    try: jsonschema.validate(instance=data, schema=schema)
    except ValidationError as e:
        print("JSON validation error:", e)
        sys.exit(1)
    scale = get_value_or_schema_default(data,["header","chart","scale"])
    html = generate_html(data)
    with open(output_file,'w',encoding='utf-8') as f: f.write(html)
    print(f"Generated {output_file} successfully.")
    global_css = CssStyle()

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Usage: seqsee <input.json> <output.html>")
        sys.exit(1)
    process_json(sys.argv[1], sys.argv[2])
