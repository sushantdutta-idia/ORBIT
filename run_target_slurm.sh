#!/bin/bash
#SBATCH --partition=Main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=%x-%j.slurm.log
#
# One target of a run_targets.py run, as a SLURM job on a compute node.
#
# run_targets.py holds each target while it works on it, so one job per
# target can share one work folder: the targets run at the same time
# instead of one after another, and a restart of the Jupyter server no
# longer stops them. A job that is cancelled or killed is taken up again by
# submitting it again: every step keeps its own record, and the claims a
# killed job leaves are taken over once SLURM no longer lists it.
#
# The dataset's own settings follow the target: the pipeline's products
# folder and the work folder, and anything else run_targets.py takes. The
# telescope, band, grid and channel range are read from the products (give
# --channel-range, --cell or --band to overrule them). Submit from the
# pipeline folder, one job per target; --output puts the job's log in the
# work folder, which must exist first. The Band 4 run was:
#
#   mkdir -p targets_run3
#   for T in G0216-0604 G0218-0245; do
#     sbatch --job-name=$T --output=targets_run3/%x-%j.slurm.log run_target_slurm.sh $T \
#       --products uGMRT_Products --workdir targets_run3 --scans G0212-0335=15
#   done
#
# A MeerKAT field needs katbeam for its primary beam, in the Python that runs
# the chain; give that Python as PIPELINE_PYTHON:
#
#   mkdir -p MeerKAT_targets
#   PIPELINE_PYTHON=$HOME/venvs/orbit-imaging/bin/python3 sbatch --job-name=J2344 \
#     --output=MeerKAT_targets/%x-%j.slurm.log run_target_slurm.sh J2344 \
#     --products MeerKAT_Products --workdir MeerKAT_targets \
#     --peel-above 30 --keep-radius 20
#
# WSClean is 3.6 (wsclean-latest.sif) cleaning in sub-images of 2560
# pixels: on G0222-0321's data it made the same image as WSClean 3.1 whole
# in 20 minutes instead of 55, on 16 threads (finding 139). Options for
# run_targets.py may follow the target and overrule these. The job's CPUs
# are WSClean's threads, and three quarters of its memory is WSClean's.

set -euo pipefail

target="${1:?give the target first: sbatch --job-name=T run_target_slurm.sh T --products P --workdir W}"
shift

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# The driver and the scripts it runs need a Python with numpy; WSClean and
# CASA come from their containers.
python=""
for candidate in "${PIPELINE_PYTHON:-}" python3 /opt/jupyterhub/.venv/bin/python3; do
    [ -n "$candidate" ] || continue
    if command -v "$candidate" >/dev/null 2>&1 \
            && "$candidate" -c "import numpy" >/dev/null 2>&1; then
        python="$candidate"
        break
    fi
done
if [ -z "$python" ]; then
    echo "No python3 with numpy on $(hostname). Set PIPELINE_PYTHON to one." >&2
    exit 2
fi
if ! command -v singularity >/dev/null 2>&1; then
    echo "singularity is not on PATH on $(hostname)." >&2
    exit 2
fi

memory_gb=$(( ${SLURM_MEM_PER_NODE:-122880} * 3 / 4 / 1024 ))
threads="${SLURM_CPUS_PER_TASK:-16}"
echo "$(date -u '+%Y-%m-%d %H:%M:%S UTC') ${target} on $(hostname), job ${SLURM_JOB_ID:-?}:" \
     "$(command -v "$python"), ${threads} threads, WSClean memory ${memory_gb} GB"

exec "$python" -u run_targets.py --targets "$target" \
    --threads "$threads" --memory-gb "$memory_gb" \
    --wsclean "singularity exec /software/astro/containers/wsclean-latest.sif wsclean" \
    --parallel-deconvolution 2560 "$@"
