"""Render the low-amplitude echo-jet success case through Viewer PlotCanvas."""

from pathlib import Path

import plot_rpe_echo_jet_success_detail as viewer_plot


ROOT = Path(__file__).resolve().parents[1]
viewer_plot.TARGET_PATH = ROOT / (
    "archives/solutions/rpe/"
    "rpe-periodic-100khz-20kpa-phase-portrait__20260911T214426Z.npz"
)
viewer_plot.CONTROL_PATH = ROOT / (
    "archives/solutions/rpe/"
    "rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300__"
    "20260914T171839Z.npz"
)
viewer_plot.VIEWER_BUNDLE = (
    ROOT / "viewer_outputs/rpe-echo-jet2-20kpa-success-viewer"
)
viewer_plot.PAYLOAD_PATH = viewer_plot.VIEWER_BUNDLE / "payload.json"
viewer_plot.CATALOG_PATH = viewer_plot.VIEWER_BUNDLE / "catalog.json"
viewer_plot.RENDER_REPORT_PATH = (
    ROOT / "viewer_outputs/rpe-echo-jet2-20kpa-success-detail.render.json"
)
viewer_plot.OUTPUT_PATH = (
    ROOT / "viewer_outputs/rpe-echo-jet2-20kpa-success-detail.png"
)
viewer_plot.CATALOG_ID = "rpe-echo-jet2-20kpa-success-detail"
viewer_plot.TILE_NAME = "RPE 20 kPa echo-jet KLQR — final six cycles"
viewer_plot.CATALOG_LABEL = "20 kPa echo-jet KLQR"
viewer_plot.PAGE_LABEL = "Successful low-amplitude control"


if __name__ == "__main__":
    viewer_plot.main()
