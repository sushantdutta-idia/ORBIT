# Science validation status, acceptance gates, and observing-band recipes

Last audit update: 2026-08-12

## Read this before using a product for science

This document is a validation contract, not a claim that the pipeline is
already observatory-certified. The repository implements and unit-tests many
safety checks, but the audit environment did not contain CASA, real telescope
Measurement Sets, or approved reference reductions. Synthetic and mocked tests
can establish control-flow behavior; they cannot establish an absolute flux
scale, a primary-beam model, image fidelity, astrometry, or RFI completeness.

**No band in the matrix below is publication-ready solely because the pipeline
finishes. Every row remains `NPR` until the real-reference-dataset programme in
this document has passed for the exact software/container versions in use.**

The supported science scope is continuum Stokes I. Full-polarization
calibration, spectral-line processing, transient calibration, tied-array data,
and VLBI/international-baseline reduction are outside this validation scope.

Status terms used below are deliberately narrow:

- **`I` — implemented:** the repository contains a routed code path and
  synthetic/unit coverage. This says nothing about real-data fidelity.
- **`C` — conditional:** the implementation may be used only when the named
  dependency, model, provenance, and all acceptance gates below are satisfied.
- **`E` — external-only:** this repository must not perform that operation;
  use the named observatory/telescope workflow and retain its provenance.
- **`NPR` — not publication-ready:** at least one required scientific
  validation or model is missing. Quick-look output must be labelled as such.

## Honest per-band support matrix

Frequency limits are receiver/profile limits, not a promise that every channel
is sensitive or scientifically usable. The observed SPW coverage, guaranteed
RF range, filter roll-off, persistent RFI, and observing log always take
precedence.

| Telescope / tuning | Frequency scope (MHz) | Flagging, cross-calibration, self-calibration | Imaging, PB, catalogue | DD / 3GC | Audit classification and release blocker |
|---|---:|---|---|---|---|
| uGMRT B2 | Receiver and measured-PB profile 125–250; approximately 175-MHz notch is a special region | `I/C`: raw MS or visibility FITS, CASA bounded flagging, fresh post-calibration reflag/re-solve, CASA phase self-calibration | `I/C`: Taylor-0 correction with the measured 2023 NCRA polynomial; image reference frequency must be 125–250 MHz | `I/C/NPR`: bundled killMS/DDFacet recipe is experimental | `NPR`. CAPTURE's published CASA6 guide did not validate B2; require dedicated B2 references, strong ionospheric/RFI review, and a PB-radial regression |
| uGMRT B3 | 250–500; measured PB model covers 250–500 | `I/C`: same two-round CASA path; band-aware edge/RFI policy | `I/C`: measured 2023 NCRA polynomial Taylor-0 correction | `I/C/NPR`: bundled uGMRT-only experimental recipe | `NPR` until two independent B3 references pass flux, astrometry, noise, and PB-radial tests |
| uGMRT B4 | 550–850 science-valid receiver range; a 400-MHz GWB table through 950 MHz is accepted only with an exact >850-MHz channel mask; measured PB model covers 550–850 | `I/C`: same two-round CASA path | `I/C`: measured 2023 NCRA polynomial Taylor-0 correction | `I/C/NPR` | `NPR`; 850–950 MHz is backend output outside the validated receiver/PB range and must never enter calibration or imaging |
| uGMRT B5 | Receiver profile 1000–1460; measured PB model 1050–1450 | `I/C`: same two-round CASA path | `I/C`: measured 2023 NCRA polynomial Taylor-0 correction; reference frequency must lie inside 1050–1450 MHz | `I/C/NPR` | `NPR`; validate band-edge/RFI behavior, PB, and absolute scale on real B5 data; do not extrapolate the PB model to receiver edges |
| MeerKAT UHF | Digitized profile 544–1088; guaranteed science scope 580–1015 | `I/C`: channels outside 580–1015 are deterministically masked, then CASA parallel-hand Stokes-I uses conservative pre-calibration TFCrop, guarded corrected-data TFCrop/RFlag, and a complete fresh solve; an explicit validated broadband flux model is mandatory | `I/C`: WSClean imaging (the default) with the `katbeam` UHF beam applied to its output by `meerkat_beam.py` (band-averaged beam for the MFS image, each sub-band at its own frequency); CASA `tclean` plus `katbeam` Taylor-0 with `--imager casa`; spectral index is not PB-corrected. The WSClean path has synthetic and unit coverage only: no real-data flux, astrometry or PB-radial reference has been run through it | `E`: use a MeerKAT-specific 3GC workflow | `NPR`. The masked digitizer edges are never part of claimed usable bandwidth |
| MeerKAT L | Digitized profile 856–1712; guaranteed science scope 900–1670 | `I/C`: channels outside 900–1670 are deterministically masked; the persistent masks 933–960, 1163–1299, and 1524–1630 MHz are then scoped only to L; guarded post-calibration flags trigger a complete fresh solve | `I/C`: WSClean (default) with the measured `katbeam` L beam applied to its output (band-averaged for the MFS image, per sub-band otherwise); CASA `tclean` plus Taylor-0 `katbeam` with `--imager casa`; spectral index is not PB-corrected. The WSClean path has synthetic and unit coverage only: no real-data flux, astrometry or PB-radial reference has been run through it | `E`: MeerKAT-specific 3GC | `NPR`. Receiver-edge and persistent masked gaps must be excluded from claimed usable bandwidth |
| MeerKAT S0 | 1750–2625 | `I/C`: explicit `MEERKAT_S0` identity with a shared S-tuning CASA calibration/self-calibration policy | `NPR`: final imaging, PB-required source detection, and publication products are blocked because no validated S-band PB provider is integrated | `E`: validated S-band beam/3GC workflow | `NPR`; validate S0 separately, including tuning-specific model, RFI, beam, and pointing behavior |
| MeerKAT S1 | 1968.75–2843.75 | `I/C`: explicit `MEERKAT_S1` calibration/self-calibration identity | `NPR`: blocked as above | `E` | `NPR`; S1 requires its own real-reference pair |
| MeerKAT S2 | 2187.5–3062.5 | `I/C`: explicit `MEERKAT_S2` calibration/self-calibration identity | `NPR`: blocked as above | `E` | `NPR`; S2 requires its own real-reference pair |
| MeerKAT S3 | 2406.25–3281.25 | `I/C`: explicit `MEERKAT_S3` calibration/self-calibration identity | `NPR`: blocked as above | `E` | `NPR`; S3 requires its own real-reference pair |
| MeerKAT S4 | 2625–3500 | `I/C`: explicit `MEERKAT_S4` calibration/self-calibration identity | `NPR`: blocked as above | `E` | `NPR`; S4 requires its own real-reference pair and high-frequency pointing/DD review |
| LOFAR LBA | Broad profile 10–90; actual station mode and retained subbands come from the observation and LINC | `E`: LINC/DP3/AOFlagger/LoSoTo own raw flagging, clock/TEC, Faraday rotation, bandpass, element beam, DI self-calibration, and bad-station handling. CASA flagging/self-calibration is prohibited here | `I/C`: only a finalized calibrated target MS may enter WSClean/EveryBeam; beam state and data column must be explicit | `E`: LINC/facet/LOFAR-specific DD workflow | `NPR` until a LINC reference product and this exact WSClean/EveryBeam handoff pass LBA flux/PB/astrometry tests |
| LOFAR HBA | Broad profile 110–240; actual station mode and retained subbands come from the observation and LINC | `E`: same upstream-only policy | `I/C`: finalized calibrated target MS, WSClean/EveryBeam, explicit beam state | `E`: LINC HBA selfcal/facet/LOFAR-specific DD | `NPR` until HBA reference validation passes |

The S0–S4 ranges overlap by design because the receiver is wider than the
instantaneously processed bandwidth. A narrow frequency selection can
therefore match several tunings. Pass `--receiver-band S0`, `S1`, `S2`, `S3`,
or `S4` whenever the observed range is ambiguous, and always verify/record the
resolved tuning against SARAO observing metadata. An unambiguous complete
tuning may be detected from frequency coverage; generic `S` is rejected.

## Mandatory release evidence

A science release directory must contain, or link immutably to, all of the
following. A `.done` marker is never sufficient evidence by itself.

| Evidence | Minimum artifact |
|---|---|
| Input identity | Archival URI/observation ID, immutable checksum manifest, byte/file count, observation time, telescope, receiver/tuning, correlator mode, fields, scans, SPWs, correlations, and original data-column meaning |
| Reproducible software | Pipeline commit or source-tree checksum, configuration checksum, CASA/Python versions, container digest, and versions of WSClean, EveryBeam, `katbeam`, PyBDSF, killMS/DDFacet, LINC/DP3/LoSoTo/AOFlagger as applicable |
| Pipeline provenance | `pipeline_run_manifest.json`; for FITS, `ingest/ingest_report.json`; `inspection/inspection_report.json`; the working-copy manifest; and the exact command line |
| Flagging | `flagging/flagging_parameters.json`, iteration and extension reports, flag versions, and antenna/field/SPW/channel/baseline plots; for LOFAR, LINC/AOFlagger reports and parsets instead |
| Calibration | `calibration/calibration_report.json`, every caltable QA sidecar and plot, post-calibration flag report, round-1 archive, round-2 report/tables where required, `split_targets.json`, `split_targets.provenance.json`, and flux-model provenance |
| Self-calibration | `selfcal/reports/selfcal_report.json`, loop history, accepted/rejected candidates, solution plots, and before/after image metrics |
| Imaging / catalogue | `imaging/reports/imaging_report.json`, `imaging_metrics.json`, PB/beam provenance, residual/restored/PB-corrected images and plots, source-detection report, catalogues, and the independent validation measurements below |
| External handoff | For LOFAR or external DD/3GC: upstream workflow version, configuration, input/output checksums, solution files, beam state, applied terms, and data-column semantics |

The pipeline manifest is necessary but not a substitute for an archival
checksum manifest. For an MS directory, hash a deterministic sorted manifest
of every regular file. Keep the raw MS read-only and reduce a verified working
copy. A visibility FITS input must remain immutable. `ingest_report.json` now
performs part of the comparison itself: `structure_comparison` must report
`consistent`, or every `not_comparable` check must be justified by the FITS
convention in use (multi-IF frequency coverage is defined by the FQ table, and
FITS-IDI is outside this cross-check). An `inconsistent` result is a release
blocker. The remaining comparisons -- time range and field structure -- stay
manual.

## Acceptance gates

Every gate below is fail-closed. An expert may issue a project waiver only by
recording the failed criterion, diagnostic evidence, revised selection, and a
comparison reduction. A waiver changes the release decision; it must never
silently change the recorded measurement.

### Gate 1 — input provenance and applicability

Accept only when all of these are true:

1. The immutable input checksum matches the release manifest and the working
   copy is tied to that exact identity. A same-path replacement is a different
   input and requires a new output directory.
2. Telescope identity from `OBSERVATION/TELESCOPE_NAME` agrees with antenna
   names and the explicit CLI selection. Mixed-telescope or conflicting
   metadata is rejected.
3. The complete observed frequency range fits one supported receiver profile.
   MeerKAT UHF/L overlap requires an explicit receiver selection where
   ambiguous. MeerKAT S0–S4 tuning is verified against observing metadata and
   selected explicitly for any narrow overlapping range.
4. The selected flux, bandpass, phase calibrator, target field, reference
   antenna, SPWs, scans, and parallel-hand correlations are explicit and
   present. Calibration never proceeds merely because a catalogue name was
   guessed.
5. Every selected calibrator scan brackets or otherwise supports the stated
   interpolation to target scans. Slews, shadowed scans, hardware failures,
   sub-array changes, and known observatory incident intervals are documented.
   `time_coverage.target_phase_coverage` in the inspection report measures the
   bracketing directly: `bracketed` must be true for the calibrator mapped to
   each target, and `max_visit_gap_minutes` must be short enough that the
   transferred phase remains valid across it for the observing frequency and
   conditions. A leading or trailing extrapolation is a waiver-level finding,
   not a detail.
6. `field_intents` and `calibrator_catalogues` in the inspection report show
   what the field classification was actually inferred from. A catalogue
   recorded as `missing` invalidates any catalogue-based classification,
   because every unmatched field then becomes a science target.
6. The requested science frequencies lie inside the normally guaranteed RF
   range, or an explicit project validation covers the digitized edge. Fully
   masked channels/SPWs are removed from the claimed bandwidth.
7. A calibrated LOFAR input has a complete LINC/DP3 lineage, materialized data
   (not a virtual concatenation unsupported by the downstream tools), a known
   data column, and an unambiguous beam state. Raw LOFAR data is rejected.

### Gate 2 — flagging by field and SPW

The values below are local engineering damage limits, not universal
observatory standards and not proof that RFI is gone. Each automatic operation
must have a restorable flag snapshot and an accepted/rolled-back report.

| Operation | Global newly flagged / total limit | Per-field and per-SPW newly flagged / total limit | Required action on violation |
|---|---:|---:|---|
| uGMRT/MeerKAT pre-calibration automatic pass | 0.20 / 0.70 | 0.35 / 0.85 | Restore the saved flag version and stop |
| Final explicit flag extension | 0.10 newly flagged | 0.20 newly flagged | Restore and stop; do not hide propagation damage in a global total |
| uGMRT corrected-data post-calibration pass | 0.10 / 0.70 | 0.15 / 0.85 | Restore and stop; accepted flags require a complete fresh solve |
| MeerKAT corrected-data post-calibration pass | 0.20 / 0.75 | 0.30 / 0.90 | Restore and stop; accepted flags require a complete fresh solve |
| LOFAR | No CASA limit applies | Use LINC/AOFlagger QA; LINC's standard target chunk gate uses at least 0.5 unflagged fraction | Reject failed chunks/stations or rerun the appropriate external workflow |

For science release, additionally require:

- Separate summaries for every antenna, field, SPW, and their relevant
  field-by-SPW selections. A global fraction cannot certify a calibrator.
- Every **retained** calibrator field/SPW has at most 50% flagged data and
  every retained target field/SPW has at most 70% flagged data. A deliberately
  excised SPW is allowed only when removed consistently from calibrator and
  target selections and from the claimed frequency coverage.
- Every retained calibrator scan has at least 50% unflagged data in the SPWs
  needed for its solve. No required calibrator, target, polarization, or SPW
  may have zero usable samples.
- Antennas above the configured 70% bad-antenna threshold are excluded and
  cannot be the reference antenna. Abrupt antenna/baseline/field/SPW outliers
  below that limit still require manual review.
- The known MeerKAT L masks are applied only to L-band data. They are never
  transferred by frequency coincidence to UHF, S, uGMRT, or LOFAR.
- Before/after waterfalls and amplitude-versus-time/frequency plots show no
  obvious residual broadband bursts, narrow persistent lines, unflagged
  correlator failures, or flagging of credible extended emission. Short- and
  long-baseline behavior must be reviewed separately for extended targets.

The resolved parallel-hand selector is recorded in
`flagging/flagging_parameters.json` and must name both hands of a single feed
basis. A run that could not resolve it fails closed; a recorded selector of
`ABS_ALL` indicates a pre-audit product and invalidates the flagging review,
because cross-hand flags propagate into the science hands through the
extension step.

The recorded reference-antenna choice in `flagging/bad_antenna_report.json`
must be defensible on array geometry as well as flag fraction.

CASA `tfcrop`/`rflag` operate on absolute parallel hands here with internal
extension disabled. The independently checkpointed extension is the one and
only intentional propagation step. Changing that behavior is a new flagging
policy and invalidates the reference result.

### Gate 3 — calibration-table integrity and coverage

For each expected delay, pre-bandpass gain, bandpass, short-interval phase,
scan-averaged transfer phase, amplitude, and flux-scaled table, and for every
self-calibration table:

1. The table is readable, non-empty, tied to the current input/configuration,
   and contains the expected parameter and flag columns.
2. There are no unflagged NaN or infinite values, and the usable solution set
   is not entirely zero. Any such condition is a hard failure.
3. Cross-calibration's total flagged-solution fraction is no greater than 0.50
   and no antenna/SPW/field/observation group exceeds 0.80. Self-calibration
   instead counts flagged, non-finite, and zero gain solutions as unusable and
   applies 0.50 global plus 0.80 antenna/SPW/field/observation/time gates. The
   stricter release rule is
   that every retained antenna/SPW/calibrator selection retains at least 50%
   usable solutions; otherwise exclude it explicitly and derive all dependent
   tables again.
4. Time coverage includes every intended calibrator scan. SPW/channel coverage
   maps exactly to the retained science selection, with no silent extrapolation
   across an uncalibrated SPW or a large masked gap. All target times have the
   intended phase/flux interpolation support.
5. Bandpass amplitude/phase and gain amplitude/phase plots are finite and
   physically smooth at the solution resolution. Discontinuities shared by
   many antennas, reference-antenna jumps, wraps, spikes, and edge divergence
   are explained or cause rejection.
6. If corrected-data flagging was accepted for uGMRT or MeerKAT, the archived
   round-1 tables predate those flags and the active delay-through-flux tables
   are a complete fresh round 2. Mixing table rounds is prohibited.
7. `applycal` produces finite `CORRECTED_DATA` and valid weights for all
   retained selections. The calibrated target MS preserves field, SPW, time,
   UVW, and flag provenance.
8. The applied table list is the intended one. The pre-bandpass gain table is
   an input to the bandpass solve only and must not appear in the later gain
   solves or in `applycal`. Science targets must be corrected with the
   scan-averaged transfer phase table of their mapped phase calibrator, not
   with the short-interval table used for the calibrators themselves. Both
   facts are recorded in `calibration_parameters` in the calibration report.
9. `fluxscale` reference and transfer field sets are disjoint, and the
   effective `fitorder` recorded in `calibration/fluxscale_report.json` is no
   larger than the number of spectral windows minus one.

For LOFAR, apply the analogous tests to the LINC H5Parm/LoSoTo solution sets:
station/subband/time coverage, finite values, flags/weights, clock–TEC and
Faraday-rotation diagnostics, element-beam state, and the finalized target MS.
CASA caltable criteria do not replace LINC QA.

### Gate 4 — flux scale and calibrator model

1. Record the primary calibrator, flux standard, model source, model version,
   valid frequency range, reference epoch if variable, and checksum of every
   external model. Save the positive finite per-SPW fluxes returned by
   `setjy`; for an external CASA image/component list predicted with `ft`,
   retain its content hash and the validated non-zero `MODEL_DATA` evidence.
2. uGMRT defaults to Perley–Butler 2017 in this implementation. MeerKAT uses
   Stevens–Reynolds 2016 where applicable; UHF and any calibrator without a
   validated built-in representation require an explicit broadband model.
   LOFAR keeps the scale and sky-model lineage established by LINC. A named
   standard is not enough if the chosen calibrator/model is unsupported at the
   observed frequencies.
3. `fluxscale` transfer values and uncertainties must be finite, positive, and
   supported by the primary model. A secondary/phase calibrator spectrum must
   be smooth after fitting the expected spectral behavior and consistent
   across scans.
4. Image an independent calibrator or compact, isolated, non-variable field
   sources before accepting the target flux scale. Correct catalog values to
   the image reference frequency using measured or defensible spectral
   indices; do not compare unmatched frequencies as if sources were flat.
5. If amplitude self-calibration is enabled, it must use normalized solutions
   only after phase-only convergence. Median central-field compact-source flux
   may change by no more than 3% from the accepted pre-amplitude image, with no
   trend versus source brightness or distance from pointing centre.
6. Every self-calibration solve, amplitude included, must have run against a
   validated `MODEL_DATA` column re-predicted from the preceding accepted
   image. A solve against CASA's implicit unit point source is void.
7. The image dynamic range that gated each self-calibration mode must be
   recorded and distinct from the `gaincal` minsnr used inside it.

Provisional project-level cross-catalogue gates are below. They are deliberately
conservative local release criteria, **not published observatory guarantees**.
Use at least 20 high-S/N unresolved sources inside PB >= 0.5 when available;
with fewer matches, require an independent calibrator observation and manual
sign-off.

| Band | Allowed absolute median flux-scale offset | Allowed robust source-to-source scatter after uncertainty subtraction |
|---|---:|---:|
| uGMRT B2 | 15% | 20% |
| uGMRT B3/B4/B5 | 10% | 15% |
| MeerKAT UHF/L | 10% | 15% |
| MeerKAT S0–S4 | Central-field calibration check: 10%; no PB-corrected release until a validated S model is integrated | 15% central-field only |
| LOFAR HBA | 15% | 20% |
| LOFAR LBA | 20% | 25% |

### Gate 5 — image noise, fidelity, flux, PB, and astrometry

All final checks use the stated target, data column, weighting, UV selection,
restoring beam, reference frequency, and PB cutoff. Acceptance requires:

- Restored, residual, PSF, PB, and PB-corrected Taylor-0 products are readable,
  finite in their valid region, have consistent WCS/shape/units, and record a
  positive restoring beam. The image is not blank and the CLEAN model is not
  silently empty in the presence of a strong residual.
- Robust blank-sky RMS is measured on the uncorrected residual/restored image,
  with emission and artefacts masked. It is no more than 2 times the larger of
  the independently predicted thermal noise and a documented confusion floor.
  If dynamic-range limited, the product fails this gate until DD calibration,
  peeling, or a documented comparison reduction resolves the excess.
- Blank-sky residual median has absolute value <= 0.2 times its robust RMS.
  Unmodelled isolated residuals above 7 sigma, coherent sidelobe rings/stripes,
  broad negative bowls, or CLEAN components tracking the PSF require rejection
  or an explicit re-imaging/DD investigation.
- Central-field compact-source fluxes pass Gate 4. Within PB >= 0.5, the
  corrected/reference flux ratio has no significant radial trend: a fitted
  change exceeding 5% from centre to the PB=0.5 radius fails the PB model or
  pointing check.
- With at least 20 high-S/N compact catalogue matches, the median radial
  astrometric offset is <= the maximum of 0.1 synthesized-beam FWHM, 0.5 pixel,
  and the combined catalogue/image 1-sigma systematic uncertainty. The 95th
  percentile is <= 0.3 beam and neither RA nor Dec shows a significant bulk
  offset. With fewer matches, validate on a dedicated astrometric calibrator.
- `image_selection.resolved_role` in the source-detection report is the role
  the catalogue image actually fulfils, and
  `image_selection.matches_requested_preference` is true. A catalogue built
  from a restored image must not be published with PB-corrected flux claims.
- Source counts, negative-image detections, and injected-source recovery are
  measured at the intended catalogue threshold. Inside the declared PB/noise
  region, reliability must be >= 95% and completeness at the quoted survey
  limit must be >= 90%, or the catalogue must publish the measured lower
  values rather than claim the nominal limit.
- PB correction is never extrapolated outside the declared cutoff. LOFAR beam
  mode and the meaning of differential/full correction are verified against
  the upstream LINC state. `--lofar-beam-mode none` is quick-look only.
- The declared field of view is the one the image actually delivers. Read
  `products.pbcor_status.effective_aperture` in the imaging report: when
  `image_encloses_pbmin_contour` is false the product is truncated by
  `imsize * cell` at `effective_pb_cutoff`, not at `pbmin`, and that higher
  value is what must be published. Undeconvolved sources between the image
  edge and the main lobe are a dynamic-range risk that must be assessed
  before release.
- `imaging_parameters.flat_imaging_overrides` is empty, or every listed
  override is a deliberate, recorded decision. A flat `config.yaml` imaging key
  replaces that value in every band and discards telescope-specific tuning.
- `products.pbcor_status.pointing_centre_offset` shows the image centre
  coincides with the antenna pointing. The beam models are centred on the
  image reference pixel; an offset mis-centres the beam and biases all
  off-axis photometry.
- Sensitivity is quoted from the residual or restored image, per
  `metrics.rms_reference`, never from the PB-corrected image.
- `imaging_parameters.science_reference_frequency` is recorded. Where a
  receiver digitizes beyond its specified range, the image reference frequency
  must be the centre of the science-valid observed band, not of the full
  digitized span: the measured primary beam is evaluated there, so the wrong
  choice imposes exactly the radial flux trend this gate tests for.
- `pblimit` is recorded in `imaging_parameters`. For a measured-beam provider
  it must be well below `pbmin` so CASA's internal beam model plays no part in
  defining the published aperture; equality with `pbmin`, or CASA's 0.2
  default, means two different beam models are deciding two different things
  in the same image.
- MeerKAT WSClean products: the MFS image is divided by the band-averaged
  `katbeam` beam, exact for a flat-spectrum source and, for a source whose
  spectral index is not zero, biased by an amount that grows with spectral
  index, offset from the pointing centre and the band's fractional width. The
  report records the effective frequencies and weights, the model and package
  version, `pbmin`, and the flux error a single-frequency correction would
  have carried at the half-power radius and at `pbmin`. A MeerKAT image from
  this path must not be released until the same flux, astrometry and
  PB-radial gates as for the CASA path have been run on it, and the two
  backends compared on one data set: they use the same beam model and the
  same image-plane conventions, but different gridders and deconvolvers.
- Only Taylor-0 continuum intensity is within the currently claimed bundled PB
  scope. MeerKAT TT1/alpha and other wideband spectral products are not
  publication-ready under the Taylor-0-only `katbeam` path. uGMRT in-band
  spectral products require independently PB-corrected subbands and their own
  flux/beam validation.

## Operational recipes

These commands show the supported routing. Replace every placeholder with the
inspection-confirmed value, run the dry-run first, and use a new output
directory for a different input, target, model, or validation policy.

### uGMRT B2, B3, B4, or B5

For an archival raw MS:

```bash
python3 run_ugmrt.py /archive/OBS.ms \
  --receiver-band B3 \
  --flux-calibrator FLUX_CAL \
  --bandpass-calibrator BANDPASS_CAL \
  --phase-calibrator PHASE_CAL \
  --target-fields TARGET \
  --target TARGET \
  --reference-antenna REFANT \
  --copy-input-ms \
  --outdir /scratch/OBS_B3 \
  --dry-run
```

For UVFITS/FITS-IDI, replace the positional MS with the visibility FITS file,
add `--fits-format uvfits` or `--fits-format fitsidi`, and omit
`--copy-input-ms`; ingest creates a new MS. A FITS image is not a visibility
input. After checking the plan and inspection report, rerun without
`--dry-run`. Apply the band-specific cautions in the matrix, then execute all
five gates. B2 needs its own manual low-frequency/ionospheric review. The
bundled uGMRT DD stage remains experimental until its separate direction,
solution, image, and flux-degeneracy validation passes.

### MeerKAT UHF or L

```bash
python3 run_meerkat.py /archive/OBS.ms \
  --receiver-band L \
  --flux-calibrator FLUX_CAL \
  --bandpass-calibrator BANDPASS_CAL \
  --phase-calibrator PHASE_CAL \
  --target-fields TARGET \
  --target TARGET \
  --reference-antenna REFANT \
  --copy-input-ms \
  --outdir /scratch/OBS_L \
  --dry-run
```

For UHF, add `--flux-model /models/validated_broadband_model` and archive its
checksum and source. Use an explicit model for any calibrator not represented
validly by the selected `setjy` standard. For L, verify the three persistent
frequency masks against the current observation and exclude the masked
bandwidth from sensitivity and spectral claims. The flagging stage writes an
exact per-SPW channel mask for the digitized channels outside the guaranteed
UHF/L science range; inspect its report and never restore those channels into
the science claim. A high-dynamic-range field still needs an external MeerKAT
3GC workflow.

### MeerKAT S0–S4 calibration-only handoff

Record and pass the exact tuning. Stop before final imaging/source detection:

```bash
python3 run_meerkat.py /archive/OBS_S2.ms \
  --receiver-band S2 \
  --flux-calibrator FLUX_CAL \
  --bandpass-calibrator BANDPASS_CAL \
  --phase-calibrator PHASE_CAL \
  --target-fields TARGET \
  --target TARGET \
  --reference-antenna REFANT \
  --copy-input-ms \
  --stages inspection,flagging,calibration,selfcal \
  --outdir /scratch/OBS_S2_calonly \
  --dry-run
```

The final PB-dependent imaging and source-detection stages must fail closed.
Move the calibrated/self-calibrated MS into a separately validated S-band
imaging/beam/3GC workflow. Do not substitute the UHF/L `katbeam` model, a
Gaussian beam, or a tuning-independent analytic beam for publication without
reference validation.

### LOFAR LBA or HBA

First run the official LINC calibrator and target workflows for the actual
station mode. Preserve their CWL/parset, container digest, calibrator/target
pairing, H5Parm solutions, AOFlagger reports, finalized MS checksums, applied
beam terms, and data-column definition. Only then run this repository:

```bash
python3 run_lofar.py /archive/LINC_finalized_target.ms \
  --data-state calibrated \
  --receiver-band HBA \
  --target TARGET \
  --wsclean-data-column CORRECTED_DATA \
  --lofar-beam-mode full \
  --no-copy-input-ms \
  --outdir /scratch/LOFAR_TARGET \
  --dry-run
```

Use `full` only when WSClean/EveryBeam should derive and apply the complete
station beam from trustworthy MS metadata. Use `differential` only when the
upstream visibilities are already corrected for the phase-centre beam exactly
as WSClean documents. Use `none` only for an explicitly uncorrected quick-look;
it cannot feed a science catalogue. LOFAR flagging, CASA self-calibration, and
the bundled uGMRT DD stage are prohibited.

### Gate 6 -- direction-dependent products

A DD image is released only when all of the following hold:

1. `pb_provider` in `ddcal/reports/dd_image_metrics.json` names the same
   measured beam that final imaging corrected with. A DD product corrected
   with CASA's generic `tclean` beam is on a different off-axis flux scale
   from the continuum image and must not be published beside it.
2. `direction_selection.uncovered_bright_sources` is zero, or every listed
   source is individually assessed and shown not to limit the map. A bright
   source with no direction is not DD-corrected and its artefacts remain.
3. The killMS usable finite/non-zero solution fraction is at or above its
   configured minimum, and the median gain amplitude is close to unity. A
   median far from unity indicates a flux-scale degeneracy between the DD
   gains and the sky model.
4. `DDESolutions-GlobalNorm` is recorded, and the DD and non-DD images agree
   on the flux density of compact sources inside PB >= 0.5 to within the
   Gate 4 limits. DD gains can otherwise absorb a global scale factor.
5. Sensitivity is quoted from `apparent_rms_Jy_beam`, not from the
   PB-corrected RMS.
6. The number of imaging facets and the number of DD solution directions are
   recorded separately and each justified; they are different quantities.

## Real-reference-dataset validation still required

Promotion from `NPR` requires a signed validation manifest for **each matrix
row**, including each of S0, S1, S2, S3, and S4 separately. At minimum:

1. Select two independent, openly retrievable observations per row: one
   calibrator-rich/compact-source field with a trusted comparison catalogue,
   and one difficult field containing realistic RFI plus bright off-axis or
   extended emission. Include low/high elevation or ionospheric conditions
   where relevant. Record archival IDs before claiming support.
2. Produce an independent comparison reduction with the primary telescope
   workflow: current NCRA guidance/CAPTURE-informed manual CASA for uGMRT,
   processMeerKAT or an approved SARAO workflow for MeerKAT, and LINC plus a
   validated LOFAR imaging/facet workflow for LOFAR.
3. Freeze the raw checksum, pipeline source/configuration checksum, CASA and
   all external tool/container versions. Run from a clean output directory;
   never seed the test with products from a different fingerprint.
4. Exercise MS input for every row. Exercise visibility-FITS ingest wherever
   an authoritative equivalent exists, verifying selection, visibility
   statistics, weights/flags, calibration tables, and final fluxes against the
   MS-start reduction.
5. Apply every acceptance gate above and publish the measurements, plots,
   failures, and waivers. Compare per-field/SPW flags, caltable coverage,
   recovered calibrator/field-source fluxes, noise, beam, astrometry, and
   catalogue completeness—not merely file existence.
6. Repeat an identical clean run. Required reproducibility limits are <=0.5
   percentage point change in retained flag fraction, <=1% median compact-source
   flux change, <=3% robust RMS change, and <=0.05 beam median positional change.
   Larger differences must be traced to a documented nondeterministic tool or
   fail validation.
7. Test interrupted-stage recovery and resume against the same input identity.
   The resumed result must pass the same science metrics and must not reuse a
   stale MS, caltable, split target, beam, or image.
8. Sign the resulting manifest with dataset IDs, expected metric ranges,
   reviewer, date, and expiry conditions. Any change to receiver limits,
   calibrator model, flagging policy, PB provider, CASA/tool major version, or
   stage algorithm invalidates the affected row until regression is rerun.

## Primary documents and software references

These links are the basis for the scope and recipes; local thresholds above
are explicitly identified as pipeline/project policy.

### uGMRT

- [Current GMRT system parameters and receiver ranges](https://www.gmrt.ncra.tifr.res.in/doc/GMRT_specs.pdf)
- [NCRA 2023 uGMRT primary-beam polynomial coefficients](https://www.gmrt.ncra.tifr.res.in/doc/beam-shape-v1-29nov2023.pdf)
- [CAPTURE pipeline paper](https://arxiv.org/abs/2010.00196)
- [CAPTURE CASA6 guide](https://www.ncra.tifr.res.in/~ruta/IDAP/capture-cas6-doc.pdf)

### MeerKAT

- [SARAO MeerKAT specifications, guaranteed RF and S0–S4 ranges](https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/277315585/MeerKAT%2Bspecifications.%26)
- [SARAO S-band capability and tuning status](https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1588854789)
- [SARAO MeerKAT primary-beam guidance](https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1481572357/The%2BMeerKAT%2Bprimary%2Bbeam)
- [MeerKAT holography/beam-characterization paper](https://arxiv.org/abs/2301.06752)
- [processMeerKAT cross-calibration sequence](https://idia-pipelines.github.io/docs/processMeerKAT/cross-calibration-in-processmeerkat/)
- [processMeerKAT configuration and persistent L-band masks](https://idia-pipelines.github.io/docs/processMeerKAT/config-files/)
- [processMeerKAT self-calibration guidance](https://idia-pipelines.github.io/docs/processMeerKAT/self-calibration-in-processmeerkat/)
- [processMeerKAT science-imaging guidance](https://idia-pipelines.github.io/docs/processMeerKAT/science-imaging-in-processmeerkat/)
- [`katbeam` measured UHF/L model package](https://github.com/ska-sa/katbeam)

### LOFAR

- [LINC overview](https://linc.readthedocs.io/en/stable/)
- [LINC calibrator workflow and AOFlagger strategy](https://linc.readthedocs.io/en/latest/calibrator.html)
- [LINC target workflow](https://linc.readthedocs.io/en/latest/target.html)
- [LINC standard LBA/HBA workflow variants](https://linc.readthedocs.io/en/latest/running.html)
- [WSClean/EveryBeam primary-beam correction and differential-beam semantics](https://wsclean.readthedocs.io/en/latest/primary_beam_correction.html)

### CASA task semantics used by the acceptance gates

- [`flagdata` TFCrop/RFlag/extension behavior](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.flagging.flagdata.html)
- [`setjy` flux-model behavior](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.setjy.html)
- [`gaincal` solution semantics](https://casadocs.readthedocs.io/en/v6.6.4/api/tt/casatasks.calibration.gaincal.html)
- [`fluxscale` transfer semantics](https://casadocs.readthedocs.io/en/v6.6.0/api/tt/casatasks.calibration.fluxscale.html)
- [`applycal` corrected-data behavior](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.calibration.applycal.html)
- [`tclean` products, restart behavior, and wide-field imaging](https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.imaging.tclean.html)
