"""SQLite metadata index for fast, non-destructive archive discovery."""
from __future__ import annotations

import json
import os
import hashlib
from contextlib import contextmanager
from pathlib import Path
import sqlite3
from typing import Any
import numpy as np

from .store import ARCHIVE_ROOT

INDEX_PATH = ARCHIVE_ROOT / "management" / "archive_index.sqlite"
ARCHIVE_SUMMARY_VERSION = 2
_SUMMARY_FIELDS = (
    "label", "system_label", "system", "method", "analysis", "format",
    "system_preset", "equations_version", "equations_fingerprint",
    "integration_form", "krom_mode_count", "eigenfunction_representation",
    "control_map_time_domain", "control_coordinate_system", "control_maps",
    "selected_modes", "independent_archive", "delta_t", "dataset_role",
    "attempt_count", "stop_reason", "controller_kind", "source_archive",
    "source_archives", "source_spectrum", "source_krom", "result_archives",
)


def summarize_archive_metadata(metadata_json: str) -> str:
    """Return the compact, versioned metadata projection consumed by UIs."""
    try:
        source=json.loads(metadata_json) if metadata_json else {}
        summary={key:source[key] for key in _SUMMARY_FIELDS if key in source}
        if source.get("format")=="dynsyslab.hierarchy/v1":
            summary.update({key:source[key] for key in ("job","systems","tier") if key in source})
        configuration=source.get("system_configuration",{})
        if "identifier" in configuration:
            summary["system_configuration"]={"identifier":configuration["identifier"]}
        if source.get("analysis_metadata",{}).get("krom") is not None:
            summary["analysis_metadata"]={"krom":True}
        krom=source.get("controller_configuration",{}).get("settings",{}).get("krom")
        if krom is not None:
            summary["controller_configuration"]={"settings":{"krom":krom}}
        source_krom=source.get("training_preflight",{}).get("source_krom")
        if source_krom is not None:
            summary["training_preflight"]={"source_krom":source_krom}
        return json.dumps(summary,separators=(",",":"),sort_keys=True)
    except (TypeError,ValueError):
        return "{}"


@contextmanager
def _connect(archive_root: str | Path | None=None):
    index_path=(Path(archive_root).resolve()/"management"/"archive_index.sqlite") if archive_root else INDEX_PATH
    index_path.parent.mkdir(parents=True, exist_ok=True)
    con=sqlite3.connect(index_path,timeout=30); con.execute("PRAGMA journal_mode=WAL"); con.execute("PRAGMA synchronous=NORMAL")
    con.execute("""CREATE TABLE IF NOT EXISTS solutions (
        path TEXT PRIMARY KEY, mtime_ns INTEGER NOT NULL, size_bytes INTEGER NOT NULL,
        label TEXT, system TEXT, driver_amplitude REAL, driver_frequency REAL,
        time_start REAL, time_end REAL, dimensionless INTEGER, status TEXT,
        metadata_json TEXT NOT NULL)""")
    con.execute("CREATE INDEX IF NOT EXISTS solutions_lookup ON solutions(label, driver_amplitude, driver_frequency, dimensionless, time_end)")
    con.execute("""CREATE TABLE IF NOT EXISTS archive_objects (
        relative_path TEXT PRIMARY KEY, modified_ticks INTEGER NOT NULL,
        size_bytes INTEGER NOT NULL, sha256 TEXT NOT NULL,
        arrays_json TEXT NOT NULL,
        summary_version INTEGER NOT NULL DEFAULT 0,
        summary_json TEXT NOT NULL DEFAULT '')""")
    columns={row[1] for row in con.execute("PRAGMA table_info(archive_objects)")}
    if "summary_version" not in columns:
        con.execute("ALTER TABLE archive_objects ADD COLUMN summary_version INTEGER NOT NULL DEFAULT 0")
    if "summary_json" not in columns:
        con.execute("ALTER TABLE archive_objects ADD COLUMN summary_json TEXT NOT NULL DEFAULT ''")
    if "metadata_json" in columns:
        con.executescript("""BEGIN IMMEDIATE;
            DROP TABLE IF EXISTS archive_objects_compact;
            CREATE TABLE archive_objects_compact (
                relative_path TEXT PRIMARY KEY, modified_ticks INTEGER NOT NULL,
                size_bytes INTEGER NOT NULL, sha256 TEXT NOT NULL,
                arrays_json TEXT NOT NULL, summary_version INTEGER NOT NULL DEFAULT 0,
                summary_json TEXT NOT NULL DEFAULT '');
            INSERT INTO archive_objects_compact SELECT relative_path,modified_ticks,
                size_bytes,sha256,arrays_json,summary_version,summary_json
                FROM archive_objects;
            DROP TABLE archive_objects;
            ALTER TABLE archive_objects_compact RENAME TO archive_objects;
            COMMIT;""")
        try:
            con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            con.execute("VACUUM")
        except sqlite3.OperationalError:
            pass  # The compact schema is committed; a later maintenance pass can vacuum.
    con.execute("""CREATE TABLE IF NOT EXISTS campaign_runs (
        campaign_id TEXT PRIMARY KEY, manifest_path TEXT NOT NULL, label TEXT,
        goal TEXT, campaign_type TEXT, system TEXT, preset TEXT,
        created_utc TEXT, submitted_utc TEXT, last_updated_utc TEXT,
        status TEXT NOT NULL, metadata_json TEXT NOT NULL)""")
    con.execute("""CREATE TABLE IF NOT EXISTS campaign_jobs (
        campaign_id TEXT NOT NULL, job_id TEXT NOT NULL, job_label TEXT,
        status TEXT NOT NULL, manifest_path TEXT, updated_utc TEXT NOT NULL,
        PRIMARY KEY (campaign_id, job_id))""")
    con.execute("""CREATE TABLE IF NOT EXISTS campaign_solutions (
        campaign_id TEXT NOT NULL, solution_path TEXT NOT NULL, job_id TEXT,
        recorded_utc TEXT NOT NULL, PRIMARY KEY (campaign_id, solution_path))""")
    con.execute("CREATE INDEX IF NOT EXISTS campaign_solutions_lookup ON campaign_solutions(campaign_id, solution_path)")
    con.execute("""CREATE TABLE IF NOT EXISTS archive_index_state (
        singleton INTEGER PRIMARY KEY CHECK(singleton=1),
        generation INTEGER NOT NULL)""")
    con.execute("INSERT OR IGNORE INTO archive_index_state(singleton,generation) VALUES(1,0)")
    con.execute("""CREATE TABLE IF NOT EXISTS archive_catalog_snapshots (
        consumer TEXT PRIMARY KEY, generation INTEGER NOT NULL,
        schema_version INTEGER NOT NULL, payload_json TEXT NOT NULL)""")
    try:
        yield con
        con.commit()
    except Exception:
        con.rollback()
        raise
    finally:
        con.close()

def canonical_label(value: str) -> str:
    import re
    return re.sub(r"-c\d+-t\d+$", "", value)


def _invalidate_catalog(con: sqlite3.Connection) -> None:
    con.execute("UPDATE archive_index_state SET generation=generation+1 WHERE singleton=1")
    con.execute("DELETE FROM archive_catalog_snapshots")

def _row(path: Path) -> tuple[Any,...] | None:
    try:
        with np.load(path,allow_pickle=False) as data: metadata=json.loads(str(data["metadata_json"]))
        source=((metadata.get("forcing",{}).get("sources") or [{}])[0].get("parameters",{})); span=metadata.get("time_span",[None,None]); stat=path.stat()
        return (str(path.resolve()),stat.st_mtime_ns,stat.st_size,canonical_label(str(metadata.get("label",""))),str(metadata.get("system","")),source.get("amplitude"),source.get("frequency_hz"),span[0],span[-1],int(bool(metadata.get("options",{}).get("use_dimensionless_form",False))),metadata.get("status","complete"),json.dumps(metadata,sort_keys=True))
    except Exception: return None

def index_solution(path: str | Path, *, archive_root: str | Path | None=None) -> bool:
    row=_row(Path(path));
    if row is None: return False
    with _connect(archive_root) as con:
        con.execute("INSERT INTO solutions VALUES (?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(path) DO UPDATE SET mtime_ns=excluded.mtime_ns,size_bytes=excluded.size_bytes,label=excluded.label,system=excluded.system,driver_amplitude=excluded.driver_amplitude,driver_frequency=excluded.driver_frequency,time_start=excluded.time_start,time_end=excluded.time_end,dimensionless=excluded.dimensionless,status=excluded.status,metadata_json=excluded.metadata_json",row)
    return index_archive_object(path,archive_root=archive_root)

def index_archive_object(path: str | Path, *, archive_root: str | Path | None=None) -> bool:
    """Publish one immutable NPZ pointer into the shared Archiver index.

    The exact digest is paid once at publication. Later reconciliations compare
    path, byte size, and nanosecond mtime before reopening or rehashing a file.
    """
    root=Path(archive_root or ARCHIVE_ROOT).resolve(); target=Path(path).resolve()
    try:
        relative=target.relative_to(root).as_posix()
        stat=target.stat()
        with np.load(target,allow_pickle=False) as data:
            arrays=list(data.files)
            metadata=str(data["metadata_json"]) if "metadata_json" in data.files else ""
        digest=hashlib.sha256()
        with target.open("rb") as stream:
            for block in iter(lambda:stream.read(1024*1024),b""): digest.update(block)
        row=(relative,stat.st_mtime_ns,stat.st_size,digest.hexdigest(),
             json.dumps(arrays,separators=(",",":")),
             ARCHIVE_SUMMARY_VERSION,summarize_archive_metadata(metadata))
    except (OSError,ValueError,KeyError):
        return False
    with _connect(root) as con:
        previous=con.execute("""SELECT relative_path,modified_ticks,size_bytes,
            sha256,arrays_json,summary_version,summary_json FROM archive_objects
            WHERE relative_path=?""",(relative,)).fetchone()
        if previous==row:
            return True
        con.execute("""INSERT INTO archive_objects(
            relative_path,modified_ticks,size_bytes,sha256,arrays_json,
            summary_version,summary_json) VALUES (?,?,?,?,?,?,?)
            ON CONFLICT(relative_path) DO UPDATE SET
            modified_ticks=excluded.modified_ticks,size_bytes=excluded.size_bytes,
            sha256=excluded.sha256,arrays_json=excluded.arrays_json,
            summary_version=excluded.summary_version,summary_json=excluded.summary_json""",row)
        _invalidate_catalog(con)
    return True

def build_index(*, reconcile: bool=True, archive_root: str | Path | None=None) -> dict[str,int]:
    root=Path(archive_root or ARCHIVE_ROOT).resolve(); index_path=root/"management"/"archive_index.sqlite"
    scanned=indexed=unchanged=removed=failed=0
    with _connect(root) as con:
        known={row[0]:(int(row[1]),int(row[2]),int(row[3]),bool(row[4])) for row in
               con.execute("SELECT relative_path,modified_ticks,size_bytes,summary_version,length(summary_json) FROM archive_objects")}
        known_solutions={row[0]:(int(row[1]),int(row[2])) for row in
                         con.execute("SELECT path,mtime_ns,size_bytes FROM solutions")}
    present:set[str]=set()
    for path in root.rglob("*.npz"):
        if path.name.endswith(".partial.npz"): continue
        scanned+=1
        try:
            stat=path.stat(); absolute=str(path.resolve())
            relative_path=path.resolve().relative_to(root); relative=relative_path.as_posix(); present.add(relative)
            signature=(stat.st_mtime_ns,stat.st_size)
            solution=bool(relative_path.parts and relative_path.parts[0].lower()=="solutions")
            indexed_signature=known.get(relative)
            summary_current=bool(indexed_signature and
                indexed_signature[2]==ARCHIVE_SUMMARY_VERSION and indexed_signature[3])
            if indexed_signature and indexed_signature[:2]==signature and summary_current and (not solution or known_solutions.get(absolute)==signature):
                unchanged+=1; continue
            success=index_solution(path,archive_root=root) if solution else index_archive_object(path,archive_root=root)
            if not success and solution: success=index_archive_object(path,archive_root=root)
            if success: indexed+=1
            else: failed+=1
        except (OSError,ValueError):
            failed+=1
    if reconcile:
        with _connect(root) as con:
            removed_archives=0
            for relative in known:
                if relative not in present:
                    con.execute("DELETE FROM archive_objects WHERE relative_path=?",(relative,)); removed+=1; removed_archives+=1
            for absolute in known_solutions:
                if not Path(absolute).is_file():
                    con.execute("DELETE FROM solutions WHERE path=?",(absolute,)); removed+=1
            if removed_archives:
                _invalidate_catalog(con)
    return {"scanned":scanned,"indexed":indexed,"unchanged":unchanged,"removed":removed,"failed":failed,
            "database_bytes":index_path.stat().st_size if index_path.exists() else 0}

def _compatible(metadata: dict[str,Any], spec: dict[str,Any], *, allow_shorter: bool) -> bool:
    driver=spec.get("driver",{}); source=((metadata.get("forcing",{}).get("sources") or [{}])[0].get("parameters",{}))
    kind=driver.get("kind","zero")
    expected_class={"sinusoid":"SinusoidalDriver","zero":None}.get(kind)
    if kind not in {"sinusoid","zero"}: return False  # no unsafe continuation of arbitrary/tabulated forcing
    if expected_class is not None and source.get("_dynsyslab_class") != expected_class: return False
    if kind=="zero" and metadata.get("forcing_type","free") != "free": return False
    for key in ("amplitude","frequency_hz"):
        if key in driver and not np.isclose(float(driver[key]),float(source.get(key,float("nan"))),rtol=1e-12,atol=1e-12): return False
    if kind=="sinusoid" and not np.isclose(float(driver.get("phase_rad",0.)),float(source.get("phase_rad",0.)),rtol=1e-12,atol=1e-12): return False
    if bool(metadata.get("options",{}).get("use_dimensionless_form",False)) != bool(spec.get("solver",{}).get("options",{}).get("use_dimensionless_form",False)): return False
    have=metadata.get("time_span",[]); need=spec.get("time_span",[])
    if len(have)!=2 or len(need)!=2 or not np.isclose(float(have[0]),float(need[0]),rtol=1e-12,atol=1e-12): return False
    return float(have[-1]) < float(need[-1]) if allow_shorter else np.isclose(float(have[-1]),float(need[-1]),rtol=1e-12,atol=1e-12)

def find_efficiency_candidate(spec: dict[str,Any]) -> dict[str,Any] | None:
    label=str(spec.get("label","")); driver=spec.get("driver",{}); dim=int(bool(spec.get("solver",{}).get("options",{}).get("use_dimensionless_form",False)))
    with _connect() as con:
        rows=con.execute("SELECT path,metadata_json FROM solutions WHERE label=? AND dimensionless=? AND status IN ('complete','completed')",(label,dim)).fetchall()
    exact=[]; extend=[]
    for path,raw in rows:
        metadata=json.loads(raw)
        if _compatible(metadata,spec,allow_shorter=False): exact.append((path,metadata))
        elif _compatible(metadata,spec,allow_shorter=True): extend.append((path,metadata))
    if exact: return {"kind":"reuse","archive_path":exact[-1][0],"metadata":exact[-1][1]}
    if extend:
        path,metadata=max(extend,key=lambda item:float(item[1]["time_span"][-1])); return {"kind":"extend","archive_path":path,"metadata":metadata}
    return None

def register_campaign_run(manifest_path: str | Path, *, status: str="submitted") -> dict[str,Any]:
    """Register a JobManager-submitted campaign as lightweight archive history."""
    path=Path(manifest_path).resolve(); payload=json.loads(path.read_text(encoding="utf-8")); campaign_id=str(payload["campaign_id"])
    metadata=dict(payload.get("metadata",{})); now=__import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()
    with _connect() as con:
        con.execute("""INSERT INTO campaign_runs(campaign_id,manifest_path,label,goal,campaign_type,system,preset,created_utc,submitted_utc,last_updated_utc,status,metadata_json)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(campaign_id) DO UPDATE SET
                       manifest_path=excluded.manifest_path,label=excluded.label,goal=excluded.goal,campaign_type=excluded.campaign_type,
                       system=excluded.system,preset=excluded.preset,last_updated_utc=excluded.last_updated_utc,status=excluded.status,metadata_json=excluded.metadata_json""",
                    (campaign_id,str(path),payload.get("label"),payload.get("goal"),metadata.get("campaign_type"),metadata.get("system_name"),metadata.get("preset"),payload.get("created_utc"),now,now,status,json.dumps(payload,sort_keys=True)))
    return campaign_record(campaign_id) or {}

def record_campaign_job(campaign_id: str, job_id: str, *, status: str, job_label: str|None=None, manifest_path: str|Path|None=None) -> None:
    now=__import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()
    with _connect() as con:
        con.execute("""INSERT INTO campaign_jobs(campaign_id,job_id,job_label,status,manifest_path,updated_utc) VALUES(?,?,?,?,?,?)
                       ON CONFLICT(campaign_id,job_id) DO UPDATE SET job_label=excluded.job_label,status=excluded.status,manifest_path=excluded.manifest_path,updated_utc=excluded.updated_utc""",
                    (str(campaign_id),str(job_id),job_label,status,None if manifest_path is None else str(Path(manifest_path).resolve()),now))
        con.execute("UPDATE campaign_runs SET last_updated_utc=?, status=? WHERE campaign_id=?",(now,status,str(campaign_id)))

def record_campaign_solutions(campaign_id: str, job_id: str, paths: list[str|Path]) -> None:
    now=__import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()
    with _connect() as con:
        for path in paths:
            candidate=Path(path)
            if candidate.is_file(): con.execute("INSERT OR IGNORE INTO campaign_solutions(campaign_id,solution_path,job_id,recorded_utc) VALUES(?,?,?,?)",(str(campaign_id),str(candidate.resolve()),str(job_id),now))
        con.execute("UPDATE campaign_runs SET last_updated_utc=?, status='running' WHERE campaign_id=?",(now,str(campaign_id)))

def campaign_record(campaign_id: str) -> dict[str,Any] | None:
    with _connect() as con:
        row=con.execute("SELECT campaign_id,manifest_path,label,goal,campaign_type,system,preset,created_utc,submitted_utc,last_updated_utc,status,metadata_json FROM campaign_runs WHERE campaign_id=?",(str(campaign_id),)).fetchone()
        if row is None: return None
        jobs=con.execute("SELECT job_id,job_label,status,manifest_path,updated_utc FROM campaign_jobs WHERE campaign_id=? ORDER BY job_id",(str(campaign_id),)).fetchall()
        paths=con.execute("SELECT solution_path FROM campaign_solutions WHERE campaign_id=? ORDER BY solution_path",(str(campaign_id),)).fetchall()
    keys=("campaign_id","manifest_path","label","goal","campaign_type","system","preset","created_utc","submitted_utc","last_updated_utc","status","metadata_json")
    result=dict(zip(keys,row)); result["metadata"]=json.loads(result.pop("metadata_json")); result["jobs"]=[dict(zip(("job_id","label","status","manifest_path","updated_utc"),item)) for item in jobs]; result["solution_paths"]=[item[0] for item in paths]; return result

def campaign_solution_paths(campaign_id: str) -> list[Path]:
    with _connect() as con:
        rows=con.execute("SELECT solution_path FROM campaign_solutions WHERE campaign_id=? ORDER BY solution_path",(str(campaign_id),)).fetchall()
    return [Path(row[0]) for row in rows if Path(row[0]).is_file()]

def list_campaign_runs() -> list[dict[str,Any]]:
    with _connect() as con: ids=[row[0] for row in con.execute("SELECT campaign_id FROM campaign_runs ORDER BY submitted_utc DESC").fetchall()]
    return [record for item in ids if (record:=campaign_record(item)) is not None]
