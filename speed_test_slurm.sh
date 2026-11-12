#!/bin/bash
#SBATCH --partition=Main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=speed_tests/%x-%j.slurm.log
#
# One image of finished data, made again with one setting changed, to time
# it against the others: G0222-0321's calibrated data from the last stage
# (targets_run3/G0222-0321/target_whole_field/G0222-0321/round3.ms), on the
# grid and with the cleaning the pipeline uses. The data are read, never
# written to. Submit from the pipeline folder, after `mkdir -p speed_tests`:
#
#   sbatch --job-name=NAME speed_test_slurm.sh NAME [options for experiment_wsclean_image.py]
#
# for example --wsclean 'singularity exec /software/astro/containers/wsclean-latest.sif wsclean'
# or --parallel-deconvolution 0.

set -euo pipefail

name="${1:?give the test a name: sbatch --job-name=NAME speed_test_slurm.sh NAME ...}"
shift

cd "${SLURM_SUBMIT_DIR:-$PWD}"
echo "$(date -u '+%Y-%m-%d %H:%M:%S UTC') ${name} on $(hostname), job ${SLURM_JOB_ID:-?}, ${SLURM_CPUS_PER_TASK:-16} threads"

exec python3 -u experiment_wsclean_image.py \
    --ms targets_run3/G0222-0321/target_whole_field/G0222-0321/round3.ms \
    --workdir "speed_tests/${name}" --outdir uGMRT_Products --target G0222-0321 \
    --imsize 5120 --cell 1 --robust 0 --channels-out 12 \
    --threads "${SLURM_CPUS_PER_TASK:-16}" --memory-gb 90 --channel-range 51 768 \
    --scales 0 5 15 30 --no-local-rms --auto-threshold 1 --nmiter 20 --overwrite \
    --wsclean 'singularity exec /software/astro/containers/oxkat-latest.sif wsclean' \
    "$@"
