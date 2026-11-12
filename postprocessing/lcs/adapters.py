"""Optional adapter boundary for external LCS suites."""
from pathlib import Path
ROOT=Path(__file__).parent
def find_lcs_suites(): return {"lcstool":ROOT/"external-lcstool"}
class LCSToolAdapter:
    """Validates an untouched MATLAB LCStool snapshot before external execution."""
    suite_name="lcstool"
    def __init__(self,matlab_command="matlab"): self.matlab_command=matlab_command; self.path=find_lcs_suites()[self.suite_name]
    def available(self): return self.path.exists() and any(self.path.iterdir())
    def command_hint(self):
        if not self.available(): raise FileNotFoundError("LCStool source snapshot is not installed in external-lcstool.")
        return f"{self.matlab_command} -batch \"addpath('{self.path.as_posix()}')\""
