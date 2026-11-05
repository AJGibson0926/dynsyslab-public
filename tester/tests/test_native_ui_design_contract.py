from __future__ import annotations

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[2]


def _luminance(rgb: tuple[int, int, int]) -> float:
    channels = []
    for value in rgb:
        normalized = value / 255.0
        channels.append(normalized / 12.92 if normalized <= 0.04045
                        else ((normalized + 0.055) / 1.055) ** 2.4)
    return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2]


def _contrast(first: tuple[int, int, int], second: tuple[int, int, int]) -> float:
    light, dark = sorted((_luminance(first), _luminance(second)), reverse=True)
    return (light + 0.05) / (dark + 0.05)


class NativeUiDesignContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.theme = (ROOT / "viewer/native/include/dynsyslab_viewer/Theme.hpp").read_text(
            encoding="utf-8")
        cls.colors = {
            match.group(1): tuple(int(match.group(index)) for index in range(2, 5))
            for match in re.finditer(
                r"constexpr COLORREF (\w+) = RGB\((\d+),\s*(\d+),\s*(\d+)\);",
                cls.theme,
            )
        }

    def test_light_theme_text_and_focus_pairs_meet_wcag_aa(self):
        text_pairs = (
            ("text", "background"),
            ("text", "panelRaised"),
            ("muted", "panelRaised"),
            ("accentStrong", "selection"),
            ("accentStrong", "panel"),
            ("blue", "panelRaised"),
            ("danger", "panelRaised"),
            ("warning", "panelRaised"),
        )
        for foreground, background in text_pairs:
            with self.subTest(pair=(foreground, background)):
                self.assertGreaterEqual(
                    _contrast(self.colors[foreground], self.colors[background]), 4.5)
        self.assertGreaterEqual(
            _contrast(self.colors["blue"], self.colors["panelRaised"]), 3.0)

    def test_shared_controls_are_keyboard_focusable_and_dark_theme_is_retired(self):
        panel = (ROOT / "compute_workbench/native/src/PanelCommon.cpp").read_text(
            encoding="utf-8")
        viewer = (ROOT / "viewer/native/src/Workbench.cpp").read_text(encoding="utf-8")
        native_sources = "\n".join(
            path.read_text(encoding="utf-8", errors="ignore")
            for root in (ROOT / "viewer/native", ROOT / "compute_workbench/native")
            for path in root.rglob("*.cpp")
        )
        self.assertGreaterEqual(panel.count("WS_TABSTOP"), 4)
        self.assertIn("WS_TABSTOP | BS_OWNERDRAW", viewer)
        self.assertIn("ODS_FOCUS", viewer)
        self.assertNotIn("DarkMode_Explorer", native_sources)

    def test_solver_and_viewer_layouts_preserve_the_approved_information_architecture(self):
        solver = (ROOT / "compute_workbench/native/src/WorkbenchLayout.cpp").read_text(
            encoding="utf-8")
        solver_shell = (ROOT / "compute_workbench/native/src/Workbench.cpp").read_text(
            encoding="utf-8")
        viewer = (ROOT / "viewer/native/src/WorkbenchLayout.cpp").read_text(
            encoding="utf-8")
        for stage in ("1  System", "2  ICs", "3  Control", "4  Integrator", "5  Review"):
            self.assertIn(stage, solver_shell)
        self.assertIn("Review, output, and provenance", solver_shell)
        self.assertIn("collapsedReviewH", solver)
        self.assertIn("contentHeight * 28 / 100", viewer)
        self.assertLess(viewer.index("{sortDate_"), viewer.index("{navigation_.handle()"))
        self.assertIn("selectionContext_", viewer)


if __name__ == "__main__":
    unittest.main()
