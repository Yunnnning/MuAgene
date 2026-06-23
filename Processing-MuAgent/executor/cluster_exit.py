"""Force a clean process exit at the end of a scheduler-submitted stage.

After a stage writes its HDF5/h5ad outputs (atomically, via io.write_h5ad_safe:
/tmp staging + fsync + rename), h5py/HDF5 background threads can linger as
non-daemon threads and block the Snakemake child process from exiting. SLURM
then keeps the job RUNNING ("Pid still in cpuset cgroup") even though the output
is complete on disk, hanging the whole pipeline.

`gc.collect()` releases the AnnData/h5py objects so their HDF5 file descriptors
close; `os._exit(0)` skips atexit/finalizers to terminate any remaining threads
immediately. It only fires inside a scheduler-submitted child job (SLURM_JOB_ID
or SLURM_JOB_ID set), so in local mode (job-id env vars unset) it degrades to a
plain gc and never kills a foreground run.

Call this only from non-local `<stage>_execute` rules — never from localrules or
`_propose` rules, which run inside the head-job process and must not exit.
"""
from __future__ import annotations

import gc
import os


def finalize_cluster_exit() -> None:
    """gc.collect(), then os._exit(0) when running as a scheduler child job."""
    gc.collect()
    if os.environ.get("SLURM_JOB_ID"):
        os._exit(0)
