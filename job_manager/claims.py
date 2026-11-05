"""Atomic scientific-specification ownership for JobManager."""
from __future__ import annotations
from typing import Any, Callable

def claim_specs(state: dict[str,Any], job_id: str, labels_by_key: dict[str,str], stamp: Callable[[],str]) -> tuple[bool,list[dict[str,Any]]]:
    """Mutate one locked state payload; callers remain responsible for I/O."""
    claims=state.setdefault("spec_claims",{}); jobs=state.setdefault("jobs",{})
    for key,claim in list(claims.items()):
        if jobs.get(claim.get("job_id"),{}).get("status") not in {"claiming","running","quarantined"}: claims.pop(key,None)
    conflicts=[{"label":labels_by_key[key],"owner_job_id":claims[key]["job_id"]} for key in labels_by_key if key in claims and claims[key].get("job_id")!=job_id]
    if conflicts: return False,conflicts
    now=stamp()
    for key,label in labels_by_key.items(): claims[key]={"job_id":job_id,"label":label,"claimed_utc":now}
    return True,[]

def release_specs(state: dict[str,Any], job_id: str, stamp: Callable[[],str]) -> list[str]:
    """Release one job's scientific-specification claims in locked state."""
    claims=state.setdefault("spec_claims",{}); removed=[key for key,value in claims.items() if value.get("job_id")==job_id]
    for key in removed: claims.pop(key,None)
    for row in state.setdefault("jobs",{}).values():
        if row.get("status")=="waiting_for_spec_claim": row.update({"status":"pending","claim_released_utc":stamp()}); row.pop("claim_conflicts",None)
    return removed
