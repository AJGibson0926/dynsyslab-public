"""Architecture-diagram tool: spec contract, nested layout and SVG output."""
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from architecture_diagrams import ICONS, SCHEMA, SpecError, load, parse, render  # noqa: E402


def spec(**root):
    base = {"id": "pp", "kind": "module", "label": "Post-Processing", "icon": "monitor", "layout": "column",
            "children": [{"id": "hdmd", "label": "H-DMD", "icon": "matrix", "notes": {"rank": "0"}},
                         {"id": "sindy", "label": "SINDy", "icon": "library"}]}
    base.update(root)
    return {"schema": SCHEMA, "title": "t", "root": base,
            "edges": [{"from": "hdmd", "to": "sindy", "from_side": "bottom", "to_side": "top"}]}


class ArchitectureDiagramTests(unittest.TestCase):
    def test_children_are_inside_their_container(self):
        diagram = parse(spec(), set(ICONS))
        ET.fromstring(render(diagram))
        outer = diagram.blocks["pp"]
        for name in ("hdmd", "sindy"):
            inner = diagram.blocks[name]
            self.assertGreaterEqual(inner.x, outer.x)
            self.assertGreaterEqual(inner.y, outer.y)
            self.assertLessEqual(inner.x + inner.w, outer.x + outer.w)
            self.assertLessEqual(inner.y + inner.h, outer.y + outer.h)
        self.assertLess(diagram.blocks["hdmd"].y, diagram.blocks["sindy"].y)  # column stacks vertically

    def test_contract_rejects_bad_specs(self):
        with self.assertRaises(SpecError):
            parse(spec(icon="no-such-icon"), set(ICONS))
        with self.assertRaises(SpecError):
            parse({**spec(), "schema": "other/v1"}, set(ICONS))
        bad = spec()
        bad["edges"].append({"from": "hdmd", "to": "missing"})
        with self.assertRaises(SpecError):
            parse(bad, set(ICONS))
        bad = spec()
        bad["edges"][0]["route"] = [{"x": "nowhere.left"}]
        with self.assertRaises(SpecError):
            render(parse(bad, set(ICONS)))

    def test_every_icon_is_valid_svg(self):
        for name, draw in ICONS.items():
            ET.fromstring(f'<svg xmlns="http://www.w3.org/2000/svg">{draw()}</svg>')

    def test_templates_scaffold_valid_specs(self):
        from architecture_diagrams.templates import TEMPLATES, scaffold
        for name in TEMPLATES:
            ET.fromstring(render(parse(scaffold(name, name), set(ICONS))))

    def test_every_committed_spec_is_valid(self):
        for path in (ROOT / "docs" / "diagrams").glob("*.diagram.json"):
            ET.fromstring(render(load(path, set(ICONS))))

    def test_pipeline_overview_spec_renders(self):
        diagram = load(ROOT / "docs" / "diagrams" / "pipeline-overview.diagram.json", set(ICONS))
        ET.fromstring(render(diagram))
        self.assertIn("hdmd", [c.id for c in diagram.blocks["hdmd-row"].children])


if __name__ == "__main__":
    unittest.main()
