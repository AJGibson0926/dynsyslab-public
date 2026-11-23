"""Safety integration tests; all mutations are inside a fresh temporary fixture."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

SCRIPT = Path(__file__).with_name('Invoke-StorageCleanup.ps1')
SHELL = shutil.which('pwsh') or shutil.which('powershell')


@unittest.skipUnless(os.name == 'nt' and SHELL, 'Windows PowerShell required')
class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='storage-cleanup-test-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name).resolve()
        self.root = self.base / 'reviewed'
        self.root.mkdir()
        self.plan = self.base / 'plan.json'
        self.policy = self.base / 'policy.json'
        self.result = self.base / 'result.json'
        self.old = self.root / 'old.obj'
        self.old.write_text('generated compiler output')
        os.utime(self.old, (time.time()-86400*20,)*2)
        (self.root / 'results.npz').write_bytes(b'scientific result')
        self.policy.write_text(json.dumps({'roots':[dict(path=str(self.root),
            minAgeDays=7, extensions=['.obj'], reason='Disposable test fixture')]}))

    def run_tool(self, mode, success=True):
        args = [SHELL, '-NoProfile', '-File', str(SCRIPT), '-Mode', mode,
                '-PlanPath', str(self.plan)]
        args += ['-PolicyPath',str(self.policy)] if mode=='Plan' else ['-ResultPath',str(self.result)]
        process = subprocess.run(args, capture_output=True, text=True)
        self.assertEqual(process.returncode == 0, success, process.stdout+process.stderr)

    def test_delete_only_planned_old_extension(self):
        self.run_tool('Plan')
        fresh = self.root / 'new.obj'
        fresh.write_text('new file after plan')
        self.run_tool('Apply')
        self.assertFalse(self.old.exists())
        self.assertTrue(fresh.exists())
        self.assertTrue((self.root/'results.npz').exists())

    def test_changed_file_is_preserved(self):
        self.run_tool('Plan')
        self.old.write_text('changed since planning')
        self.run_tool('Apply')
        self.assertTrue(self.old.exists())
        self.assertEqual(json.loads(self.result.read_text())['skippedFiles'],1)

    def test_outside_path_rejects_entire_plan(self):
        self.run_tool('Plan')
        outside = self.base/'outside.obj'
        outside.write_text('outside reviewed root')
        data=json.loads(self.plan.read_text())
        data['files'].append(dict(data['files'][0],path=str(outside)))
        self.plan.write_text(json.dumps(data))
        self.run_tool('Apply',False)
        self.assertTrue(outside.exists())
        self.assertTrue(self.old.exists())

    def test_junction_is_not_traversed(self):
        outside=self.base/'outside';outside.mkdir()
        (outside/'protected.obj').write_text('protected')
        junction=self.root/'link'
        command="New-Item -ItemType Junction -Path '{}' -Target '{}' | Out-Null".format(
            str(junction).replace("'","''"),str(outside).replace("'","''"))
        subprocess.run([SHELL,'-NoProfile','-Command',command],check=True,capture_output=True)
        self.run_tool('Plan')
        self.run_tool('Apply')
        self.assertTrue((outside/'protected.obj').exists())

    def test_replacing_parent_with_junction_is_refused(self):
        sub=self.root/'sub';sub.mkdir();self.old.rename(sub/'old.obj')
        self.run_tool('Plan')
        outside=self.base/'outside';sub.rename(outside)
        command="New-Item -ItemType Junction -Path '{}' -Target '{}' | Out-Null".format(
            str(sub).replace("'","''"),str(outside).replace("'","''"))
        subprocess.run([SHELL,'-NoProfile','-Command',command],check=True,capture_output=True)
        self.run_tool('Apply')
        self.assertTrue((outside/'old.obj').exists())
        self.assertEqual(json.loads(self.result.read_text())['skippedFiles'],1)


if __name__=='__main__':
    unittest.main()
