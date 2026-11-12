"""Verified, JobManager-owned input staging; never touches primary archives."""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
from typing import Any
import json, shutil

def _digest(path: Path) -> str:
    digest=sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""): digest.update(block)
    return digest.hexdigest()

def stage_inputs(job: Any, root: str|Path) -> dict[str,Any]:
    """Copy declared files into job-private scratch and verify their hashes."""
    plan=dict(getattr(job,"staging",{})); inputs=[Path(item).expanduser().resolve() for item in plan.get("inputs",[])]
    target=Path(root).resolve()/str(job.job_id); target.mkdir(parents=True,exist_ok=True); copied=[]
    for source in inputs:
        if not source.is_file(): raise FileNotFoundError(f"Staging input is missing: {source}")
        destination=target/source.name
        shutil.copy2(source,destination); source_hash=_digest(source); destination_hash=_digest(destination)
        if source_hash!=destination_hash: raise OSError(f"Checksum mismatch while staging {source}")
        copied.append({"source":str(source),"staged":str(destination),"sha256":source_hash,"bytes":destination.stat().st_size})
    manifest={"job_id":job.job_id,"staging_directory":str(target),"inputs":copied,"output_location":plan.get("output_location","primary"),"stage_out":bool(plan.get("stage_out",True)),"cleanup_scratch":bool(plan.get("cleanup_scratch",True))}
    (target/"staging-manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True),encoding="utf-8")
    return manifest

def cleanup_staging(manifest: dict[str,Any]) -> None:
    if not manifest.get("cleanup_scratch",True): return
    directory=Path(manifest["staging_directory"])
    if directory.is_dir(): shutil.rmtree(directory)
