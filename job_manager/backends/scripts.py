"""Deterministic scheduler script rendering from JobManager manifests."""
from __future__ import annotations
from typing import Any

def slurm_script(job: Any, *, python_executable: str="python") -> str:
    resources=dict(job.resources); deadline=dict(job.deadline); lines=["#!/bin/bash",f"#SBATCH --job-name={job.label}",f"#SBATCH --cpus-per-task={int(resources.get('cores',1))}"]
    if resources.get("memory_mb"): lines.append(f"#SBATCH --mem={int(resources['memory_mb'])}M")
    if resources.get("wall_time_s"): lines.append(f"#SBATCH --time={max(1,int(float(resources['wall_time_s'])//60))}")
    if resources.get("gpus"): lines.append(f"#SBATCH --gpus={int(resources['gpus'])}"+(f":{resources['gpu_type']}" if resources.get("gpu_type") else ""))
    if deadline.get("checkpoint_before_deadline_s"): lines.append(f"#SBATCH --signal=B:USR1@{int(deadline['checkpoint_before_deadline_s'])}")
    if deadline.get("requeue_on_timeout"): lines.append("#SBATCH --requeue")
    lines.append(f"{python_executable} -m job_manager.worker --manager-root $DYNSYSLAB_MANAGER_ROOT --job-id {job.job_id}")
    return "\n".join(lines)+"\n"

def pbs_script(job: Any, *, python_executable: str="python") -> str:
    resources=dict(job.resources); lines=["#!/bin/bash",f"#PBS -N {job.label}",f"#PBS -l select=1:ncpus={int(resources.get('cores',1))}"]
    if resources.get("memory_mb"): lines.append(f"#PBS -l mem={int(resources['memory_mb'])}mb")
    if resources.get("wall_time_s"): lines.append(f"#PBS -l walltime={int(float(resources['wall_time_s'])//3600):02d}:{int(float(resources['wall_time_s'])%3600//60):02d}:00")
    lines.append(f"{python_executable} -m job_manager.worker --manager-root $DYNSYSLAB_MANAGER_ROOT --job-id {job.job_id}")
    return "\n".join(lines)+"\n"
