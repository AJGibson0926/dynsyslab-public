"""Input split identities and create-only publication for standalone studies."""
from contextlib import contextmanager
import os
from pathlib import Path
import tempfile

def file_id(value):
    """Validate a record or geometry label before using it as a single output path component."""
    if not isinstance(value,str) or not value or value in ('.','..') or any(c in value for c in '/\\:') or Path(value).name!=value:
        raise ValueError('Record and geometry IDs must be filenames')
    return value

def split_paths(config,root):
    """Resolve independent trajectory files and reject shared record IDs or resolved source paths across roles."""
    ids=set();paths=set();result={}
    for role in ('training','validation','test'):
        entries=config[role]
        if not isinstance(entries,dict) or not entries:raise ValueError(f'Missing {role} records')
        result[role]={}
        for name,source in entries.items():
            file_id(name);p=(Path(root)/source).resolve()
            if name in ids or p in paths:raise ValueError('Record ID or source path appears twice')
            result[role][name]=p;ids.add(name);paths.add(p)
    return result

@contextmanager
def publish_directory(destination):
    """Stage a result alongside its destination and publish the complete directory under a new name."""
    dest=Path(destination).resolve();dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():raise FileExistsError(dest)
    with tempfile.TemporaryDirectory(dir=dest.parent,prefix='.misc-study-') as tmp:
        product=Path(tmp)/'result';product.mkdir()
        yield product
        if dest.exists():raise FileExistsError(dest)
        os.rename(product,dest)
