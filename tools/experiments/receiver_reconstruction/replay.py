"""Re-query a fitted receiver lookup on a finer clock within a selected time window."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from study import CASES, load, measurements, rms


def replay(exe, study, output, window, stride=1, case='all'):
    if output.exists():
        raise FileExistsError(output)
    stats = json.loads((study/'results.json').read_text())
    training = load(Path(stats['training']))
    with np.load(stats['held_out']) as held:
        metadata = json.loads(held['metadata_json'].item())
    geometry = dict(metadata['receiver_config'], stride=stride, window_tau=list(window))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as directory:
        tmp = Path(directory)
        (tmp/'geometry.json').write_text(json.dumps(geometry))
        subprocess.run([exe, 'signals', metadata['source'], tmp/'geometry.json', tmp/'signals.npz'], check=True)
        query = load(tmp/'signals.npz')
        F, _ = measurements(training, CASES[case])
        Q, _ = measurements(query, CASES[case])
        scale = rms(F[:, 0])
        np.savez(tmp/'train.npz', features=F, state=training['state'], scales=[scale, scale])
        np.savez(tmp/'query.npz', features=Q, time=query['time'])
        subprocess.run([exe, 'estimate', tmp/'train.npz', tmp/'query.npz', str(stats['k']), tmp/'est.npz'], check=True)
        with np.load(tmp/'est.npz') as est:
            prediction = est['prediction']
        provenance = dict(schema='dynsyslab.receiver-replay/v1', study=str(study.resolve()), source=metadata['source'],
                          training=stats['training'], case=case, neighbors=stats['k'], geometry=geometry)
        np.savez_compressed(tmp/'replay.npz', time=query['time'], state=query['state'], prediction=prediction,
                            pressure=Q/query['ranges_m'], monopole=query['monopole_pressure'],
                            metadata_json=json.dumps(provenance))
        (tmp/'replay.npz').rename(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('exe', 'study', 'out'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--window', nargs=2, type=float, required=True, metavar=('START','STOP'))
    parser.add_argument('--stride', type=int, default=1)
    parser.add_argument('--case', choices=CASES, default='all')
    args = parser.parse_args()
    replay(args.exe.resolve(), args.study, args.out, args.window, args.stride, args.case)
