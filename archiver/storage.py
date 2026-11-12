"""Durable registry and verified copy migration for archive locations."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import os
import shutil
from uuid import uuid4

from archives import ARCHIVE_ROOT

def _stamp() -> str: return datetime.now(timezone.utc).isoformat()
def _atomic_json(path: Path, payload: dict[str,Any]) -> None:
    """Publish registry/relocation changes atomically to concurrent managers."""
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+f".{uuid4().hex}.partial")
    with temporary.open("w",encoding="utf-8",newline="\n") as handle:
        handle.write(json.dumps(payload,indent=2,sort_keys=True)); handle.flush(); os.fsync(handle.fileno())
    os.replace(temporary,path)
def _digest(path: Path) -> str:
    digest=sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""): digest.update(chunk)
    return digest.hexdigest()

@dataclass(frozen=True)
class StorageLocation:
    name: str
    root: str
    writable: bool=True
    purpose: str="solutions"

class StorageRegistry:
    """Registry only; existing canonical archive paths remain unchanged."""
    def __init__(self, root: str|Path|None=None):
        self.root=Path(root or ARCHIVE_ROOT).resolve(); self.path=self.root/"management"/"storage_locations.json"; self.path.parent.mkdir(parents=True,exist_ok=True)
        if not self.path.exists(): self._write({"format":"dynsyslab.storage-registry/v1","locations":{"primary":asdict(StorageLocation("primary",str(self.root)))}})
    def _read(self): return json.loads(self.path.read_text(encoding="utf-8"))
    def _write(self,payload): _atomic_json(self.path,payload)
    def register(self, location: StorageLocation) -> None:
        target=Path(location.root).expanduser().resolve(); target.mkdir(parents=True,exist_ok=True)
        data=self._read(); data["locations"][location.name]=asdict(StorageLocation(location.name,str(target),location.writable,location.purpose)); self._write(data)
    def locations(self) -> dict[str,dict[str,Any]]: return self._read()["locations"]
    def copy_verified(self, source: str|Path, destination_root: str|Path, *, remove_source: bool=False) -> dict[str,Any]:
        source=Path(source).resolve(); destination_root=Path(destination_root).resolve(); target=destination_root/source.name; target.parent.mkdir(parents=True,exist_ok=True)
        if source.name.endswith(".partial.npz"): raise ValueError("Active partial archives may not be migrated.")
        shutil.copy2(source,target)
        if source.stat().st_size!=target.stat().st_size or _digest(source)!=_digest(target):
            target.unlink(missing_ok=True); raise IOError("Migration verification failed; source retained.")
        if remove_source: source.unlink()
        return {"migrated_utc":_stamp(),"source":str(source),"destination":str(target),"verified":True,"source_removed":remove_source}

    def create_relocation_plan(self, *, campaign_id: str, destination: str, source_root: str|Path|None=None) -> Path:
        """Create an inactive migration manifest; it moves no files by itself."""
        locations=self.locations()
        if destination not in locations or not locations[destination].get("writable",False): raise ValueError("Destination must be a registered writable storage location.")
        root=Path(source_root or self.root).resolve(); plans=self.root/"management"/"relocations"; plans.mkdir(parents=True,exist_ok=True)
        path=plans/f"campaign-relocation__{campaign_id}__{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
        _atomic_json(path,{"format":"dynsyslab.campaign-relocation/v1","created_utc":_stamp(),"campaign_id":campaign_id,"source_root":str(root),"destination_location":destination,"destination_root":locations[destination]["root"],"state":"planned","copies":[],"pending_job_ids":[],"finalized":False})
        return path

    def append_completed_copy(self, plan_path: str|Path, source: str|Path) -> dict[str,Any]:
        """Copy a newly completed file for an active relocation plan; source stays intact."""
        plan_path=Path(plan_path); plan=json.loads(plan_path.read_text(encoding="utf-8"))
        if plan.get("finalized"): raise RuntimeError("Relocation plan is finalized.")
        source=Path(source).resolve(); source_root=Path(plan["source_root"]).resolve()
        destination_root=Path(plan["destination_root"]).resolve()
        try:
            source.relative_to(destination_root)
        except ValueError:
            already_routed=False
        else:
            # New work can be written directly to the selected location while
            # a campaign is migrating.  Record it, but never copy it onto
            # itself or create an unnecessary second archive.
            already_routed=True
        if already_routed:
            result={"migrated_utc":_stamp(),"source":str(source),"destination":str(source),
                    "verified":True,"source_removed":False,"already_at_destination":True}
            plan["state"]="copying"; plan.setdefault("copies",[]).append(result)
            _atomic_json(plan_path,plan)
            return result
        try:
            # Keep the archive layout intact: a primary
            # ``solutions/rpe/foo.npz`` becomes the same relative path at the
            # destination, rather than being flattened into one directory.
            relative=source.relative_to(source_root)
        except ValueError:
            # A result may already have been written to a routed secondary
            # location.  It is still copied safely, but placed under a clear
            # external-results namespace to avoid colliding with canonical
            # primary paths.
            relative=Path("external-results")/source.name
        result=self.copy_verified(source,destination_root/relative.parent,remove_source=False)
        result["relative_archive_path"]=str(relative)
        plan["state"]="copying"; plan.setdefault("copies",[]).append(result); _atomic_json(plan_path,plan)
        return result
