"""Small, serializable contracts shared with the C++ Viewer application."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


CATALOG_SCHEMA = "dynsyslab.viewer.catalog/v1"
PAYLOAD_SCHEMA = "dynsyslab.viewer.payload/v1"


@dataclass(frozen=True)
class Relation:
    target: str
    kind: str
    label: str

    def as_dict(self) -> dict[str, str]:
        return {"target": self.target, "kind": self.kind, "label": self.label}


@dataclass
class CatalogItem:
    identifier: str
    kind: str
    category: str
    label: str
    source: Path
    system: str = ""
    subtitle: str = ""
    modified: float = 0.0
    byte_size: int = 0
    views: list[str] = field(default_factory=list)
    tags: list[str] = field(default_factory=list)
    properties: dict[str, Any] = field(default_factory=dict)
    relation_paths: list[tuple[str, str]] = field(default_factory=list)
    relations: list[Relation] = field(default_factory=list)
    role: str = "primary"
    lineage_parent_path: str = ""
    parent_identifier: str = ""
    lineage_identifier: str = ""
    primary: bool = True
    variant_count: int = 0
    # Python-only scientific identity used while assembling the catalog.  The
    # public subset is serialized below, but native clients do not need to
    # understand it in order to render the hierarchy.
    grouping: dict[str, Any] = field(default_factory=dict)

    def as_dict(self, root: Path) -> dict[str, Any]:
        try:
            source = self.source.relative_to(root).as_posix()
        except ValueError:
            source = str(self.source)
        return {
            "id": self.identifier,
            "kind": self.kind,
            "category": self.category,
            "label": self.label,
            "source": source,
            "system": self.system,
            "subtitle": self.subtitle,
            "modified": self.modified,
            "bytes": self.byte_size,
            "views": self.views,
            "tags": self.tags,
            "properties": self.properties,
            "relations": [item.as_dict() for item in self.relations],
            "role": self.role,
            "parent": self.parent_identifier,
            "lineage": self.lineage_identifier or self.identifier,
            "primary": self.primary,
            "variants": self.variant_count,
            "grouping": self.grouping,
        }


VIEW_LABELS = {
    "overview": "Overview",
    "time": "Time series",
    "phase": "Phase space",
    "three_d": "3D state",
    "fourier": "Spectrum",
    "control": "Control",
    "observable_time": "Observable",
    "observable_heatmap": "Observable map",
    "koopman_spectrum": "Eigenvalues",
    "observable_spectrum_comparison": "Compare spectra",
    "krom": "KROM modes",
    "pod_mode": "POD mode",
    "control_error": "B-psi error",
    "bifurcation": "Bifurcation",
    "observability": "Observability",
    "klqe": "KLQE",
    "data_overview": "Dataset",
    "open": "Open report",
}
