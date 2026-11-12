"""Transactional coordination state; manifests and numerical arrays remain files."""
from __future__ import annotations
from pathlib import Path
import json, sqlite3
from typing import Any

class JobStateStore:
    def __init__(self,path:str|Path):
        self.path=Path(path); self.path.parent.mkdir(parents=True,exist_ok=True)
        db=self._connect()
        try:
            db.execute("CREATE TABLE IF NOT EXISTS state (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            db.execute("INSERT OR IGNORE INTO state(key,value) VALUES('payload',?)",(json.dumps({"format":"dynsyslab.job-manager-state/v2","jobs":{},"spec_claims":{}}),))
            db.commit()
        finally: db.close()
    def _connect(self):
        db=sqlite3.connect(self.path,timeout=30,isolation_level="IMMEDIATE"); db.execute("PRAGMA journal_mode=WAL"); return db
    def read(self)->dict[str,Any]:
        db=self._connect()
        try: row=db.execute("SELECT value FROM state WHERE key='payload'").fetchone()
        finally: db.close()
        return json.loads(row[0])
    def write(self,payload:dict[str,Any])->None:
        db=self._connect()
        try:
            db.execute("UPDATE state SET value=? WHERE key='payload'",(json.dumps(payload,sort_keys=True),)); db.commit()
        finally: db.close()
