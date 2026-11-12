# Self-calibration no-solution recovery patch

This patch addresses the message:

`phase selfcal solve produced no usable calibration table`

## Changes

- Checks that `MODEL_DATA` exists, is finite, unflagged, and non-zero before every phase solve.
- Re-predicts an existing MT-MFS image model into the current MS with a `tclean(niter=0, savemodel="modelcolumn")` restart when scratch model data were lost after `split` or a restart.
- Uses a prioritized list of active antennas with `refantmode="flex"`.
- Runs one configurable long-interval rescue solve (`solint="inf"`, `combine="scan"`, `minsnr=3`, `minblperant=3`, `gaintype="T"`) if the requested solve produces no usable table.
- Validates that calibration solutions are finite and unflagged, not merely that the table has rows.
- Separates the image dynamic-range threshold for starting self-calibration (`phase_start_dynamic_range=20`) from `gaincal` per-solution `minsnr`.
- Treats failure after all solve attempts as a normal low-S/N stop, while retaining the last accepted MS.
- Fixes the candidate-MS order so `CORRECTED_DATA` is checked after `applycal`, not before.

## Restart cleanup

Remove only the failed solve rejection/checkpoint and incomplete phase-1 table before rerunning. Do not delete the initial image/model unless it is visibly wrong.
