"""Independent durable archive package for dynsyslab numerical solutions.

``solutions/`` stores completed numerical outputs. ``history/`` stores optional
in-progress checkpoints. Future exploration and visualization tools belong here,
not in the solver or system library.
"""
from .store import ARCHIVE_ROOT, ArchiveSession, compact_hierarchy, deduplicate_solution_time_grids, extend_solution, extend_solution_from_archive, load_control_target_time_series, load_control_time_series, load_hierarchy_time_series, load_time_series, recompress_archive, rewrite_archive_metadata, save_hierarchy, update_archive_metadata
from .index import INDEX_PATH, build_index, find_efficiency_candidate, index_solution, index_archive_object, register_campaign_run, record_campaign_job, record_campaign_solutions, campaign_record, campaign_solution_paths, list_campaign_runs
from time_series import TimeSeries
from .schema import ArchiveSchemaError, require_format, require_keys

__all__ = ["ARCHIVE_ROOT", "ArchiveSession", "TimeSeries", "ArchiveSchemaError", "require_format", "require_keys", "load_time_series", "load_control_time_series", "load_control_target_time_series", "load_hierarchy_time_series", "extend_solution", "extend_solution_from_archive", "save_hierarchy", "compact_hierarchy", "deduplicate_solution_time_grids", "recompress_archive", "update_archive_metadata", "rewrite_archive_metadata", "INDEX_PATH", "build_index", "index_solution", "index_archive_object", "find_efficiency_candidate", "register_campaign_run", "record_campaign_job", "record_campaign_solutions", "campaign_record", "campaign_solution_paths", "list_campaign_runs"]
