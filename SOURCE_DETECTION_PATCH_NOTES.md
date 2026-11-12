# PyBDSF source-detection patch

## Added

- `source_detection.py`: a non-CASA PyBDSF stage that runs after final imaging.
- `source_detection` in the standard pipeline stage order.
- Restart validation and a durable completion report/checkpoint.
- Configuration through `SOURCE_DETECTION_CONFIG` in `config.py` and flat
  `config.yaml` keys.

## Input image selection

The stage reads `imaging/reports/imaging_report.json` and selects the
**PB-corrected** Stokes-I FITS image by default (`image_preference: pbcor` in
both `config.py` and `config.yaml`), so catalogue flux densities are on the
corrected intensity scale off axis; PyBDSF's adaptive local-RMS estimation
handles the position-dependent noise inside the PB cutoff. Set
`source_detection_image_preference: restored` for a uniform-noise detection
experiment, or pass an explicit FITS file directly to `source_detection.py`.

Selection is fail-closed in both directions. If the imaging report names a
product for the requested role but no materialized FITS export of it exists,
the stage raises rather than falling back to filename scoring -- which would
otherwise return the uncorrected restored image for a `pbcor` request. The
role the selected image actually fulfils is recorded as `image_selection` in
the report, and a mismatch against the requested preference is a warning.

## Products

- Grouped-source catalogues: CSV and FITS (`srl`)
- DS9 source regions
- Gaussian-component catalogues: CSV and FITS (`gaul`)
- BBS Gaussian sky model grouped by source
- Gaussian model and residual FITS images
- Local-RMS FITS image
- Island-mask FITS image
- JSON report containing parameters, products, warnings, and source counts

## Defaults

- `thresh_isl = 4.0`
- `thresh_pix = 5.0`
- adaptive RMS boxes enabled
- wavelet (`atrous`) fitting enabled for extended emission
- strict dependency handling disabled

When PyBDSF is unavailable and strict mode is false, source detection is
recorded as skipped and the remainder of the pipeline remains valid. To make a
missing PyBDSF installation fail the run, set
`source_detection_strict: true`.

## Running only source detection

```bash
python3 pipeline.py \
  --outdir /path/to/existing/products \
  --stages source_detection \
  --python /path/to/python-with-bdsf \
  --no-resume
```

## Validation performed

- Compiled all pipeline Python scripts with `python -m py_compile`.
- Confirmed the full pipeline dry-run places source detection after imaging.
- Tested automatic restored-image selection from an imaging report.
- Tested safe handling when `bdsf` is unavailable.
- Tested successful catalogue/image creation with a mock `bdsf` module.
- Tested pipeline completion validation for both successful and skipped runs.

A real scientific PyBDSF run was not possible in the build environment because
PyBDSF was not installed and no real uGMRT FITS image was supplied. Thresholds
and RMS-box settings should be validated on representative fields before final
catalogue production.
