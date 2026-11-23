"""Thin compiled phase-space query adapter; all integration stays in native code."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import tempfile
from .workflow import installation_root


def query(request, archive_root, *, executable=None, timeout=600):
    """Explicit lookup/compute through the canonical native metadata cache.

    Use options.section_seeds with dynamics='periodic-forcing' for arbitrary
    stroboscopic returns and refined variational derivatives, without shooting.
    The caller owns the declared domain, forcing phase and tolerances.
    """
    cli=Path(executable).resolve() if executable else installation_root()/'bin'/('dynsyslab-cli.exe' if os.name=='nt' else 'dynsyslab-cli')
    payload={'schema':'dynsyslab.phase-space-query/v1',**dict(request)}
    with tempfile.TemporaryDirectory(prefix='dynsyslab-phase-') as folder:
        path=Path(folder)/'query.json';path.write_text(json.dumps(payload,allow_nan=False),encoding='utf-8')
        result=subprocess.run([str(cli),'phase-space',str(path),'--archive-root',str(Path(archive_root).resolve())],capture_output=True,text=True,encoding='utf-8',timeout=timeout)
    if result.returncode:raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    value=json.loads(result.stdout)
    if value.get('schema')!='dynsyslab.phase-space/v1':raise RuntimeError('Unexpected native phase-space response')
    return value
