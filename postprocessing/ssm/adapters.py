"""Optional adapter boundary for official MATLAB SSM suites."""
from pathlib import Path
ROOT=Path(__file__).parent
def find_ssm_suites(): return {"ssmtool":ROOT/"external-ssmtool","ssmlearn":ROOT/"external-ssmlearn"}
class _MATLABSuiteAdapter:
    suite_name=""
    def __init__(self,matlab_command="matlab"): self.matlab_command=matlab_command; self.path=find_ssm_suites()[self.suite_name]
    def available(self): return self.path.exists() and any(self.path.iterdir())
    def command_hint(self):
        if not self.available(): raise FileNotFoundError(f"{self.suite_name} source snapshot is not installed in {self.path.name}.")
        return f"{self.matlab_command} -batch \"addpath(genpath('{self.path.as_posix()}'))\""
class SSMToolAdapter(_MATLABSuiteAdapter): suite_name="ssmtool"
class SSMLearnAdapter(_MATLABSuiteAdapter): suite_name="ssmlearn"
