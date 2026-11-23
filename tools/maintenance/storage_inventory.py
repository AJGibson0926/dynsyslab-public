"""Read-only storage inventory for maintenance and archive planning.

Uses metadata only, never follows directory reparse points or opens file
contents. Modification age is NOT last-use age. Results belong outside Git.
"""
import argparse
from datetime import datetime,timezone
import heapq,json,os,stat,time
from pathlib import Path


def inventory(root,minimum_directory_bytes=16*1024**2,largest_count=250):
    root=Path(root).resolve();started=time.time();now=time.time()
    directories=[];largest=[];errors=[];links=[];hardlinks=set();extensions={}
    counts=dict(files=0,directories=0,skipped_reparse_directories=0,errors=0,duplicate_hardlink_names=0,offline_files=0)
    def issue(path,error):
        counts['errors']+=1
        if len(errors)<300:errors.append(dict(path=str(path),error=str(error)))
    def visit(folder):
        total=unique=old=files=0;latest=0;counts['directories']+=1
        try:entries=os.scandir(folder)
        except OSError as error:issue(folder,error);return total,unique,old,files,latest
        with entries:
            for entry in entries:
                try:
                    info=entry.stat(follow_symlinks=False)
                    attrs=getattr(info,'st_file_attributes',0)
                    reparse=entry.is_symlink() or bool(attrs&getattr(stat,'FILE_ATTRIBUTE_REPARSE_POINT',1024))
                    if entry.is_dir(follow_symlinks=False):
                        if reparse:
                            counts['skipped_reparse_directories']+=1
                            if len(links)<200:links.append(entry.path)
                            continue
                        t,u,o,f,m=visit(entry.path);total+=t;unique+=u;old+=o;files+=f;latest=max(latest,m)
                    elif stat.S_ISREG(info.st_mode):
                        size=info.st_size;counts['files']+=1;files+=1;total+=size;latest=max(latest,info.st_mtime)
                        if now-info.st_mtime>90*86400:old+=size
                        duplicate=False
                        if info.st_nlink>1:
                            key=(info.st_dev,info.st_ino);duplicate=key in hardlinks;hardlinks.add(key)
                        if duplicate:counts['duplicate_hardlink_names']+=1
                        else:unique+=size
                        if attrs&getattr(stat,'FILE_ATTRIBUTE_OFFLINE',4096):counts['offline_files']+=1
                        ext=Path(entry.name).suffix.lower() or '(none)';extensions[ext]=extensions.get(ext,0)+size
                        item=(size,entry.path,info.st_mtime)
                        if len(largest)<largest_count:heapq.heappush(largest,item)
                        elif item>largest[0]:heapq.heapreplace(largest,item)
                except OSError as error:issue(entry.path,error)
        if total>=minimum_directory_bytes or Path(folder)==root:
            directories.append(dict(path=str(folder),logical_bytes=total,unique_logical_bytes=unique,
                bytes_modified_over_90_days_ago=old,files=files,newest_file_mtime=latest))
        if counts['directories']%3000==0:print(f"Scanned {counts['directories']} folders / {counts['files']} files",flush=True)
        return total,unique,old,files,latest
    total=visit(root)
    return dict(schema='dynsyslab.storage-inventory/v1',root=str(root),created_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.time()-started,counts=counts,logical_bytes=total[0],unique_logical_bytes=total[1],
        caveats=['Logical sizes, not allocated clusters; sparse/compressed/online-only files can differ.',
                 'Hardlinks counted once in unique totals; ownership attributed to first visited name.',
                 'Reparse directories not followed; inaccessible paths reported without ACL changes.',
                 'Modification timestamps do not establish last use or deletion safety.'],
        directories=sorted(directories,key=lambda d:d['logical_bytes'],reverse=True),
        largest_files=[dict(path=p,logical_bytes=s,mtime=m) for s,p,m in sorted(largest,reverse=True)],
        extension_logical_bytes=sorted(extensions.items(),key=lambda p:p[1],reverse=True),errors=errors,skipped_reparse_examples=links)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root');parser.add_argument('output')
    parser.add_argument('--minimum-directory-mib',type=float,default=16);args=parser.parse_args()
    output=Path(args.output).resolve()
    if output.exists():raise FileExistsError('Use a fresh inventory filename')
    result=inventory(args.root,int(args.minimum_directory_mib*1024**2));output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print(json.dumps({k:result[k] for k in ['root','elapsed_seconds','counts','logical_bytes','unique_logical_bytes']},indent=2))

if __name__=='__main__':main()
