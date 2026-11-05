"""Standalone worker used by local-process and externally scheduled backends."""
from __future__ import annotations
import argparse, json
import os
from pathlib import Path
from .manager import JobManager

def main() -> int:
    parser=argparse.ArgumentParser(); parser.add_argument("--manager-root",required=True); parser.add_argument("--job-id",required=True)
    args=parser.parse_args(); manager=JobManager(root=Path(args.manager_root)); pid=os.getpid(); manager.register_worker(args.job_id,pid=pid)
    result=manager.execute(args.job_id)
    manager.record_worker_exit(args.job_id,pid=pid,exit_code=0 if result.get("status") in {"completed","restructured"} else 1)
    print(json.dumps(result,default=str)); return 0 if result.get("status") in {"completed","restructured"} else 1
if __name__=="__main__": raise SystemExit(main())
