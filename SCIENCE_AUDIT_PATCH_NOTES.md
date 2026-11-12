# ORBIT (Observatory-agnostic Radio Band Imaging Toolkit)
Tagline: One pipeline. Multiple arrays. Science-ready images.
“We developed ORBIT (Observatory-agnostic Radio Band Imaging Toolkit), a CASA-based pipeline for the calibration, flagging, self-calibration, and continuum imaging of radio interferometric observations from uGMRT, MeerKAT, and LOFAR.”

# Wideband science-correctness audit and patch

Scope: the uGMRT / MeerKAT / LOFAR continuum path, reviewed for scientific
correctness of the calibration, flagging, self-calibration and imaging
recipes rather than for control flow. Control-flow, provenance and
fail-closed machinery were already extensive and are unchanged except where
they encoded a wrong recipe.

Nothing here is a substitute for the real-reference-dataset programme in
`SCIENCE_VALIDATION.md`. Every band remains `NPR`.

## Defects corrected

### 1. Pre-bandpass phase table was applied to every field and to the target

`prebp.cal` is a `solint='int'`, phase-only table solved on the bandpass
calibrator. Its only purpose is to stop the bandpass calibrator decorrelating
while `bandpass` integrates over the full scan. It was being passed as a
`gaintable` to the final phase solve, the amplitude solve, and `applycal`,
always with `gainfield=<bandpass calibrator>`.

For any field other than the bandpass calibrator that is an interpolation of
one field's per-integration phases onto another field's timestamps. For a
science target observed hours from the bandpass scan, CASA's `nearest`
interpolation in `applycal` supplies a single integration's phase noise as a
constant per-antenna offset, while the solve had used `linear`. The mismatch
does not cancel and shows up as antenna-dependent decorrelation.

Now `prebp.cal` is used only as an input to `bandpass`. It is still solved,
QA-gated, archived between rounds and reported; it is simply no longer
applied downstream. `calibration_parameters.prebp_applied_downstream` in the
calibration report records this explicitly.

### 2. Science targets were corrected with a short-interval phase table

`phase.cal` uses the band's short `phase_solint` (8-60 s). That is the right
table for correcting the calibrators and for the amplitude solve. It is the
wrong table to interpolate across a multi-minute calibrator cycle onto a
target: it transfers one short-interval noise realisation instead of a
scan-averaged solution.

Added `run_phase_transfer_cal()`, producing `phase_scan.cal` with
`solint='inf'` (one solution per calibrator scan, configurable via
`calibration_target_phase_transfer_solint`). `applycal` now gives calibrators
`phase.cal` and targets `phase_scan.cal`. The new table participates in the
per-table QA gates, the round-1 archive, the round-2 recompute, checkpoint
invalidation, resume validation and the QA plots.

### 3. `fluxscale` was given overlapping reference and transfer field sets

With a secondary flux calibrator configured, `reference` was every flux
calibrator while `transfer` was `fluxcals[1:] + phasecals`. The secondary
appeared in both, asking CASA to bootstrap a calibrator from itself.

`transfer` is now strictly the phase calibrators that are not reference
fields. The degenerate case where no transfer field remains (phase calibrator
is itself a modelled flux calibrator) publishes `amp.cal` as `flux.cal`, which
already carries the absolute scale, instead of calling `fluxscale` with an
empty transfer list.

### 4. `fluxscale` spectral order was unbounded for single-SPW wideband data

uGMRT and MeerKAT wideband setups are frequently a single wide spectral
window, where the default `fitorder=1` cannot be constrained. The effective
order is now `min(configured, nspw - 1)` and both the configured and effective
values are written to `calibration/fluxscale_report.json`.

### 5. uGMRT Band 4 received no edge-channel flagging at all

`flag_edge_channels()` returned early whenever the band profile declared a
guaranteed science range. For a B4 observation whose GWB setup lies inside
550-850 MHz the guaranteed-range mask selects nothing, so the function
returned having flagged nothing — the sub-band filter roll-off at every SPW
edge stayed in the data through calibration and imaging.

The two masks address different defects and both now run: the guaranteed-range
mask removes spectrum the receiver does not cover, then the percentage per-SPW
mask removes the roll-off of the windows actually observed. For MeerKAT UHF/L
the percentage cut falls inside the already-masked digitizer edges, so the
change is a no-op there. `flagging/edge_flagging_report.json` now carries both
sections plus the legacy keys.

### 6. Self-calibration imaged with w-projection effectively disabled

All four `tclean` calls in `selfcal.py` requested `gridder='wproject'` without
`wprojplanes`. CASA's default is `wprojplanes=1`, which pays for the
w-projection gridder while performing no w-term correction. Final imaging
passes 96-512 planes; self-calibration was silently using one.

The consequence is not cosmetic: the non-coplanar distortion away from the
phase centre biases the CLEAN model that the self-calibration gains are solved
against, and biases the dynamic-range measurement used to accept or reject
each loop. It also made `predict_model_column()` re-predict `MODEL_DATA` under
a different gridding configuration from the image that produced the model.

`SELFCAL_IMAGING_CONFIG` now carries `wprojplanes` synchronised with
`FINAL_IMAGING_CONFIG` per band, a single `gridder_kwargs()` helper supplies it
to every self-calibration `tclean`, and the `laptop` profile scales it to 64
to match its reduced final-imaging setting.

### 7. Reference-antenna substitution could silently change the phase reference

Every `gaincal` used `refantmode='flex'` with a single reference antenna. If
that antenna is unavailable for a solution interval, CASA substitutes another
antenna, changing the phase reference between the delay, bandpass and gain
tables that are afterwards applied together. That is undetectable and
unrecoverable downstream.

The cross-calibration chain now uses `refantmode='strict'`, which flags such
intervals instead. The existing per-table and per-group flagged-fraction gates
then surface a poor reference-antenna choice as an actionable failure. Revert
with `calibration_refantmode: flex` after reviewing the table QA.

Note that CASA's `bandpass` task has no `refantmode` parameter, so the
bandpass solve retains CASA's internal behaviour. This is a residual
limitation, not an oversight.

`selfcal.py` deliberately keeps `flex`: it passes a ranked list of five
antennas, which is the supported use of that mode.

### 8. LOFAR imaged a wide fractional band as a single channel

`wsclean_channels_out` was 0 for both LOFAR bands, justified in a comment by
the claim that WSClean does not primary-beam-correct the joined-channel MFS
product. That claim is incorrect — with `-apply-primary-beam` and
`-join-channels`, WSClean publishes `<prefix>-MFS-image-pb.fits`.

LBA and HBA are both wide in fractional terms, so the PSF, the station beam
and the sky spectrum all vary appreciably across the band. Deconvolving it as
one channel biases the model and the restored flux scale. The default is now
6 joined channels with the spectral fit the code already builds. The existing
product-detection logic already preferred the MFS outputs, so no change was
needed there. This increases runtime; set `wsclean_channels_out` back to 0 in
the band profile for a quick-look image.

### 9. `--repair-split-targets-only` was passed but never implemented

`pipeline.py` passes this flag to `calibration.py` whenever calibration
solutions are valid but a split-target MS is stale — the recovery path the
README documents. `calibration.py` had no such option, so `argparse` exited
with status 2 and the stage failed. The documented recovery was unreachable.

Implemented as `CalibrationManager.repair_split_targets()`: it requires a
completed `applycal` checkpoint and `CORRECTED_DATA`, re-validates every
calibration table through the existing QA gates, rebuilds only the split
products, and preserves the recorded number of solution rounds so the
orchestrator does not reject its own repaired stage on the next resume.

### 10. `pblimit` let CASA's own beam model define the published aperture

`tclean` was never given `pblimit`, so CASA's 0.2 default applied. That default
masks the restored image using **CASA's internal** beam model, while the flux
correction divides by the pipeline's **measured** beam (the 2023 NCRA uGMRT
polynomial, or katbeam) and cuts at `pbmin=0.3`. Two different, unreconciled
beam models were deciding two different things in the same image.

Whether that mattered depends on how much narrower CASA's model is. The margin
is thin: CASA's 0.2 contour only stays outside the measured 0.3 contour while
CASA's beam is within about 15% of the measured width. The measured uGMRT main
lobe is 7-20% wider than a plain 45 m diffraction estimate (B2 1.11, B3 1.09,
B4 1.07, B5 1.20), so a diffraction-based CASA model would blank part of the
intended science region in B5 and sit within a few per cent of doing so in the
other bands. For MeerKAT the fallback model CASA uses is not the katbeam model
either, and was equally unconstrained.

Measured-beam providers now set `pblimit=0.05` (`pblimit_measured_beam`), far
outside the science cutoff, so `impbcor(cutoff=pbmin)` is the single authority
on the aperture. Where CASA's own beam *is* the correction model,
`pblimit_casa_beam` retains CASA's 0.2. Self-calibration gets the same
treatment for a different reason: blanking the outer field removes real
sources from the model the gains are solved against.

Note the asymmetry this repairs -- the LOFAR/WSClean path already passed
`-primary-beam-limit` explicitly from `pbmin`; only the CASA path was silent.

### 11. `pbmin` was described as the field-of-view limit when it often is not

`pbmin` bounds the published product only if the image is large enough to
contain that contour. Computed against the pipeline's own beam models and its
own `imsize`/`cell` defaults:

| band | image half-width | radius of PB=0.3 | reaches pbmin? |
|---|---:|---:|---|
| B2 @ 250 MHz | 68.3' | 75.8' | no -- truncated at PB=0.38 |
| B2 @ 150 MHz | 68.3' | 126.3' | no -- truncated at PB=0.72 |
| B3 @ 250 MHz | 54.6' | 74.3' | no along the axis |
| B4 @ 850 MHz | 34.1' | 21.4' | yes |
| B5 @ 1450 MHz | 20.5' | 14.1' | yes |

A uGMRT B2 image at the science defaults therefore stops just outside the
half-power point, not at PB=0.3 -- a factor of nearly two in radius. MeerKAT
UHF and L are comfortable at their configured 6144-pixel sizes.

`effective_aperture()` now measures this directly from the beam array the
pipeline just built (so it is exact for both the polynomial and katbeam, not a
model-dependent estimate), prints a warning when the `pbmin` contour is not
enclosed, and records `effective_pb_cutoff` and
`image_encloses_pbmin_contour` in the PB provenance sidecar and in
`products.pbcor_status.effective_aperture`. This is deliberately a warning
rather than a hard failure: a smaller field is often a legitimate choice. What
is not legitimate is describing the result as PB>=0.3-limited when it is not,
and the README and `SCIENCE_VALIDATION.md` have been corrected accordingly.

The undeconvolved-outer-field consequence is the more serious one for dynamic
range and is called out in both documents.

### 12. katbeam's primary beam was validated far more weakly than uGMRT's

The uGMRT path requires every PB pixel to be finite and the peak to be at most
1.05. The katbeam path accepted the plane if `np.isfinite(beam_2d).any()` was
true -- one finite pixel out of millions -- and never checked normalization. A
mostly-NaN beam would have propagated NaNs through `impbcor` into the
published intensity image, and a non-normalized beam would have silently
rescaled the flux-density scale. Both checks are now enforced with the
specific failing pixel count or peak value in the message.

## Self-calibration

### 13. Amplitude self-calibration solved against a model that was not there

`commit_phase_candidate()` publishes each accepted loop with
`split(datacolumn="corrected")`, which does not carry `MODEL_DATA`. The phase
loop handles this -- it calls `ensure_model_column()` before every solve.
`run_amp_phase_selfcal()` did not, and it runs on exactly that split output.

With no `MODEL_DATA` column, CASA supplies an implicit 1 Jy point source at
the phase centre. The amplitude solve would therefore drive every antenna gain
towards making the whole field look like a 1 Jy point source -- not a subtle
bias but a wholesale corruption of the flux scale, applied on top of an
otherwise good phase-only result. It is masked by `enable_amp_selfcal: False`
being the default, so it would first bite whoever opted in.

`run_amp_phase_selfcal()` now re-predicts the model from the last accepted
loop image (or the initial image) before solving.

### 14. One number was doing two unrelated jobs in the amplitude gate

`min_snr_amp` (30-50 depending on band) was used both as the `gaincal`
`minsnr` for the amplitude solve *and* as the image dynamic range above which
amplitude self-calibration is attempted at all. These are different
quantities, and the release notes record this same conflation being fixed for
the phase branch (`phase_start_dynamic_range` vs `min_snr_phase`) while the
amplitude branch kept it.

Read as a DR gate, 30-50 is far too permissive: amplitude self-calibration at
DR=30 absorbs sky-model error straight into the gains. `min_snr_amp` keeps its
symmetric meaning as the `gaincal` minsnr, and a new
`amp_start_dynamic_range` (default 100, `selfcal_amp_start_dynamic_range` in
`config.yaml`) gates the attempt.

### 15. Solution QA was quadratic and would stall on a normal table

`caltable_status()` grouped solutions by building, for every distinct value, a
row mask via a full pure-Python scan with two `repr()` calls per row. `TIME` is
one of the group columns, and a `solint='int'` phase table has as many distinct
times as solution intervals, so the cost is quadratic in the row count.

Measured on this machine: a 20,000-row table takes about 395 s per group
column under the old code and 0.48 s under `np.unique` grouping -- roughly an
800x difference, with identical group statistics. A 100,000-row table across
five group columns would have been hours of apparent hang inside what is meant
to be a validation step.

## Direction-dependent calibration

### 16. The DD image was primary-beam corrected with the wrong beam

`select_primary_beam_image()` matched candidates on a `.pb.tt0` / `.pb`
suffix. The measured uGMRT beam is written as `<imagename>.ugmrt.pb.tt0`,
which also ends in `.pb.tt0`, and the selector returned the first match in the
recorded product inventory. That inventory is built from
`product_candidates()`, where `pb` and `pb.tt0` sit at positions 9 and 10 and
`ugmrt.pb.tt0` at position 20 -- so CASA's generic `tclean` beam always won.

Reproduced directly against the real product ordering: with both beams present
the selector returned `TGT_final_continuum.pb.tt0`.

The consequence is that the DD product -- the final publishable image of the
whole DD workflow -- carried a different off-axis flux scale from the main
continuum image published beside it, while the README stated that the measured
polynomial was used.

Selection now reads `products.pbcor_status.ugmrt_pb_image` (or
`katbeam_pb_image`) recorded by the imaging stage, falls back to an exact
measured-beam suffix, and accepts CASA's own `.pb` only when the imaging
report shows CASA's beam was also the correction model. Otherwise it fails
closed. The provider is recorded as `pb_provider` in the DD metrics.

### 17. DD sensitivity was quoted from the primary-beam-corrected image

`dd_image_metrics.json` recorded RMS and dynamic range measured only on the
PB-corrected image, where the noise is inflated by beam amplification toward
the field edge -- the same mistake the README warns against for the main
imaging product. The metrics now also carry `apparent_rms_Jy_beam`,
`apparent_peak_Jy_beam` and `apparent_dynamic_range` from the flat-noise
image, with a note saying which to quote.

### 18. Direction centroids were wrong across RA = 0

`cluster_sources_into_directions()` merged sources with a flux-weighted linear
mean of right ascension in degrees. Two 1 Jy sources at RA 359.9 and 0.1
produced a centroid at RA 180 -- a direction on the opposite side of the sky,
which would be handed to killMS as a calibration direction. Merging now
averages on the unit circle; the same pair now yields 360.0.

Any field spanning RA=0 was affected. It is silent: the direction file simply
contains a nonsensical position.

### 19. Bright sources were dropped once the direction budget filled

Once `max_directions` was reached, every further source above `min_flux_jy`
was skipped with a bare `continue`. Those are precisely the sources DD
calibration exists to handle, and leaving them uncorrected leaves their
artefacts in the map -- with nothing in any report saying so.

They are now collected, written to
`ddcal/catalogs/uncovered_bright_sources.csv`, counted in
`ddcal/directions/direction_selection.json`, surfaced in the DD report, and
warned about at run time. Verified: 6 bright sources with `max_directions=3`
gives 3 directions and 3 recorded uncovered sources.

### 20. `min_flux_jy` was compared against intrinsic flux

Source detection defaults to the PB-corrected image, so catalogue fluxes are
intrinsic. `min_flux_jy` is a proxy for the signal the DD solver will actually
find, which is the *apparent* flux. A source at PB=0.35 with 0.3 Jy intrinsic
presents 0.1 Jy to the solver but passes a 0.1 Jy threshold as though it were
three times brighter, so directions can be selected that cannot be solved.

Not silently changed -- the threshold semantics are a science choice. The
catalogue's flux scale is now detected and recorded in
`direction_selection.flux_scale`, with an explicit note, so the threshold can
be set knowingly.

## Flagging

### 21. Unreadable correlation metadata silently widened the flagging policy

`parallel_hand_correlation()` caught every exception and set the selector to
`""`, which `absolute_parallel_hand_correlation()` turned into `ABS_ALL`. The
same happened when the `POLARIZATION` table offered neither RR/LL nor XX/YY.
`tfcrop`, `rflag` and `clip` would then run across all correlations.

That is not a safe default. Cross hands are noise-dominated for an unpolarized
field, so their tfcrop polynomial fits and rflag MAD statistics differ from the
parallel hands -- and because `flagdata` flags only the correlations it
selects, while the guarded extension runs with `extendpols=True`, any flags
found in the cross hands then propagate back into RR/LL. A metadata read
failure quietly converted a careful parallel-hand policy into an
all-correlation one whose flags reached the science hands.

`calibration.py` already failed closed on exactly these conditions for its
post-calibration pass. The pre-calibration stage now matches it: unreadable
metadata, mixed circular and linear feeds, or fewer than two parallel hands
all raise with an actionable message.

### 22. The reference antenna was chosen on flag fraction alone

`recommend_refant()` returned the antenna with the lowest flagged fraction and
nothing else. For uGMRT that can be an arm antenna tens of kilometres out --
lower per-baseline S/N and stronger differential ionosphere than the central
square antennas the telescope profile lists, and therefore a poor phase
reference for every solve that follows. The inspection stage does apply the
profile preference, but flagging overrides it and calibration takes flagging's
recommendation first. Finding 7 makes this matter more: with
`refantmode='strict'` a poor reference now costs flagged solutions.

Selection now prefers a profile antenna whose flag fraction is within
`refant_preference_tolerance` (0.05) of the best, and falls back to the
least-flagged antenna otherwise. The ranking, the method and the runner-up are
recorded in `flagging/bad_antenna_report.json`. Verified on three cases: a
central antenna 0.02 behind the best is preferred; one 0.5 behind is not; and a
central antenna that is already cleanest is chosen directly.

## Imaging

### 23. The shipped config.yaml silently discarded telescope-specific tuning

The imaging keys in `config.yaml` are flat: each replaces that value in every
band of both `FINAL_IMAGING_CONFIG` and `SELFCAL_IMAGING_CONFIG`. The file
shipped with `robust: 0.5`, `niter: 300000`, `nterms: 2` and `nsigma: 3.0`
uncommented, so loading the repository as delivered replaced 44 per-band
values, including:

| band | key | config.py | after config.yaml |
|---|---|---:|---:|
| MEERKAT_UHF / L / S0-S4 | robust | -0.5 | 0.5 |
| LOFAR_LBA / HBA | robust | -0.25 | 0.5 |
| all uGMRT | niter | 200000 | 300000 |
| all MeerKAT | niter | 150000 | 300000 |

The robust change is the significant one. MeerKAT has a very dense core, which
is exactly why `config.py` chose -0.5; forcing +0.5 weights those short
baselines heavily and produces a much larger synthesized beam than the profile
intends, with correspondingly worse confusion. The same happened to LOFAR and
to every self-calibration image. Nothing recorded that the band default had
been displaced -- the comment beside the key described +0.5 as "the science
default", which was true for uGMRT and wrong for everything else.

The flat keys are now commented out, so the per-band defaults apply. The
mechanism still works when a user deliberately uncomments one, but every
displaced default is then recorded in `FLAT_IMAGING_OVERRIDES`, warned about at
run time, and published as `imaging_parameters.flat_imaging_overrides`.

### 24. The common-beam product was smoothed after primary-beam correction

`make_common_beam_image()` preferred `self.pbcor_image` as the input to
`imsmooth`. Convolving an already-corrected map mixes pixels whose noise has
been amplified by different primary-beam gains, so the convolved product has no
well-defined noise and its weighting is position-dependent in a way the header
does not describe.

The flat-noise restored image is now convolved first and corrected afterwards
with the same measured beam, publishing both `commonbeam.image` and
`commonbeam.image.pbcor`. The feature is disabled by default
(`common_beam: ""`), so this was a latent defect rather than an active one.

### 25. Nothing checked that the beam was centred on the antenna pointing

Both beam builders evaluate their model about the image reference pixel. That
is the pointing direction for an ordinary single-field image, but a re-phased
Measurement Set moves the image centre away from the antenna boresight, and the
beam would then be applied off centre -- biasing every off-axis flux density
with nothing in the products to show it.

The offset between the image reference direction and the nearest `REFERENCE_DIR`
(falling back to `PHASE_DIR`) is now measured, recorded in the beam provenance
and in `pbcor_status.pointing_centre_offset`, and warned about above 10
arcseconds. The geometry helper handles the RA wrap and the cos(dec) term and
is unit-tested; the Measurement Set read is wrapped so an unreadable FIELD
table records `not_verifiable` rather than failing the stage.

### 26. Dynamic range and RMS were reported from the PB-corrected image

`imaging_metrics.json` reported `rms_Jy_beam` and `dynamic_range` for the
PB-corrected image alongside the restored and residual ones, with nothing
saying which is the map sensitivity. The PB-corrected RMS is inflated by beam
amplification toward the field edge and is not a sensitivity -- the same point
the README makes in prose but the metrics did not encode. A `rms_reference`
block now names the residual image as the sensitivity measurement and explains
why. (Finding 17 is the same defect in the DD stage, where it was the only
metric recorded.)

## Inspection

### 27. CALIBRATE_AMPLI was read as a phase intent

`classify_fields()` folded `*AMPLI*` into the phase-calibrator set alongside
`*PHASE*` and `*GAIN*`. In the CASA intent convention `CALIBRATE_AMPLI` marks
the amplitude/flux-scale calibrator, and a VLA-style Measurement Set commonly
carries `CALIBRATE_AMPLI`, `CALIBRATE_FLUX`, `CALIBRATE_BANDPASS` and
`CALIBRATE_PHASE` on the *same* flux-calibrator field.

The flux calibrator therefore entered the phase list. Downstream that forced
an unnecessary `--target-phase-map` ("Multiple phase calibrators are
available"), and before finding 3 it put one field into both the `fluxscale`
reference and transfer sets.

`*AMPLI*` now contributes to the flux set, which is also the safer direction:
a wrong flux classification fails closed at `setjy` (the 1 Jy fallback and
model-data checks catch it), whereas a wrong phase classification does not. A
field carrying a FLUX or BANDPASS intent is additionally removed from the
inferred phase list whenever a distinct phase calibrator remains -- but only
then, so a snapshot observation using one field for every role still works.

### 28. Nothing recorded what the classification was inferred from

Field roles are inferred from observing intents plus two bundled catalogues,
and none of those inputs appeared in the inspection report. A surprising
calibrator assignment could not be diagnosed after the fact.

Worse, `load_catalogue()` returned an empty set and printed a console warning
when a list was missing. With `flux-cals.list` absent, `is_flux_calibrator()`
returns false for everything and the fallback classifier makes every field a
science target -- recoverable only by noticing one line of stdout.

The report now carries `field_intents` (raw intents per field) and
`calibrator_catalogues` (path, entry count, or `missing`), and a missing
catalogue is a recorded warning.

### 29. The calibrator/target time structure was never measured

`SCIENCE_VALIDATION.md` Gate 1 already required that "every selected
calibrator scan brackets or otherwise supports the stated interpolation to
target scans". Nothing in the pipeline measured it, so that gate could only be
satisfied by hand.

Finding 2 makes it load-bearing: the calibration stage now transfers the phase
calibrator's scan-averaged solutions onto the target by interpolation, which
is only interpolation while calibrator scans exist on both sides of the
target. Outside that range CASA holds the nearest solution constant, so a
target whose calibrator stopped early is corrected with a stale phase for the
remainder -- silently.

`time_coverage` now records the scans and time span of every field, and
`time_coverage.target_phase_coverage` gives, per target and calibrator,
whether it is bracketed, the largest gap between calibrator visits, and the
size of any leading or trailing extrapolation. Unbracketed targets and gaps
beyond 30 minutes become inspection warnings. Verified on mocked metadata: a
normal CAL/TGT cycle reports `bracketed=True` with a 20-minute worst gap and
no warning; a calibrator that stops early reports `bracketed=False` with a
72-minute trailing extrapolation and warns.

## Orchestration

### 30. An explicit stage list ran in the order it was typed

`selected_stages()` validated each name in `--stages` but preserved the typed
order, and `selected_stage_dependencies()` only blocks a stage whose upstream
stage failed *in the same run*. `--stages imaging,calibration` therefore ran
imaging first, against whatever products a previous run had left behind, and
only failed later on an unrelated input-resolution error.

The list is now sorted into pipeline order and de-duplicated, with a printed
note when the order changed.

## Source detection

### 31. An uncorrected image could be substituted for a PB-corrected one

`select_input_image()` guards the `pbcor` preference only against a report
that records *no* PB-corrected product. When the report records one but its
FITS export is missing -- export failed, or the file was removed -- the code
fell through to filename scoring, where the uncorrected restored image scores
+40 for a `pbcor` request and wins.

Reproduced: with `pbcor_image` recorded and its export absent, the selector
returned `TGT_final_continuum.image.tt0.fits`. The resulting catalogue would
carry primary-beam attenuated off-axis flux densities while the report still
recorded `image_preference: pbcor` -- a silent flux-scale error in the
published product. The README states this specific substitution cannot
happen.

Selection now fails closed when the role is decided but unmaterialized, and
the heuristic fallback -- reachable only when the report names no product for
the role at all -- can no longer return a candidate that contradicts it.

Worth noting what was already sound here: the pipeline's own
`validate_stage_completion` binds the recorded `input_image` to the current
imaging role product and checks content hashes, so a *stale* catalogue
surviving an imaging re-run is correctly rejected. The hole was in the
first-pass selection, not in resume.

### 32. The report did not state which image role was actually used

`image_preference` records the intent, not the outcome. Nothing recorded the
role the selected image fulfilled, so a catalogue could claim a corrected flux
scale it did not have -- and the DD stage's `catalogue_flux_scale()` was
reduced to guessing from the filename and the preference.

`image_selection` now records the resolved role, whether the image is
primary-beam corrected, and whether that matches the request; a mismatch is a
warning. `ddcal.py` reads the resolved role in preference to its filename
guess.

### 33. `SOURCE_DETECTION_PATCH_NOTES.md` contradicted the shipped default

The notes stated the stage "selects the restored Stokes-I FITS image by
default", while `config.py`, `config.yaml` and the README all specify
`pbcor`. Anyone reading the notes would have assumed flat-noise input and
attenuated off-axis fluxes. Corrected.

## FITS ingest

### 34. The imported spectral and polarization axes were never validated

`validate_measurement_set()` checked required subtables, row count, required
main columns and the presence of a visibility column -- but nothing about the
axes that every later stage depends on. A Measurement Set with zero-width
spectral windows, non-finite or non-positive channel frequencies, or an empty
correlation setup passed ingest and failed several stages later as an obscure
band-detection or correlation-selection error.

Channel counts, channel frequencies and the correlation setup are now
validated at ingest, and the resulting `spectral`, `polarization` and
`observation` summaries are recorded in the report. The telescope CASA wrote
into the MS is also cross-checked against an explicit `--telescope`, closing
the loop between what the FITS header claimed and what the importer produced.

### 35. The documented FITS-to-MS structural comparison had no implementation

`SCIENCE_VALIDATION.md` required comparing "the imported MS row count, time
range, field/SPW/correlation structure, and frequency coverage with the FITS
metadata before accepting the ingest". Nothing did any of it, so a partial or
restructured import was indistinguishable from a complete one.

`compare_fits_and_ms()` now derives the expected MS row count from the FITS
`GCOUNT` and IF axis (CASA writes one row per group per IF), the expected
correlation count from the STOKES axis, and the expected channel count from
the FREQ axis for a single-IF file. Each check reports `consistent`,
`inconsistent` or `not_comparable`, and the result is recorded as
`structure_comparison` with a warning on inconsistency.

`not_comparable` is used deliberately wherever the convention cannot settle
the question rather than guessing: multi-IF frequency coverage lives in the FQ
table, and FITS-IDI is outside this route. Verified on five constructed cases
-- a matching single-IF import, a matching 4-IF import, a truncated import
(3000 rows where 4000 were expected), an import that dropped correlations, and
a FITS-IDI summary.

Left as a warning rather than a hard failure: the row-count relation holds for
ordinary random-groups UVFITS, but exotic multi-source files were not
available to test against, and a false abort on valid data would be worse than
a loud, recorded warning that the validation contract already requires a human
to read.

## Telescope profiles and CASA compatibility

### 36. A real CASA failure could be reported as a signature mismatch

`compat_task()` caught every exception from a CASA task and inspected the
message for keyword names. When no strict pattern matched, a fallback branch
returned *every* supplied keyword whose name appeared anywhere in the text,
provided the text contained a token such as "invalid parameter".

A genuine failure like `Invalid parameter: field 'XYZ' does not exist. Check
vis, field, spw.` therefore looked like a signature mismatch on `vis`, `field`
and `spw`, and surfaced as "CASA task does not support required science
keyword(s): field, spw, vis. Refusing to run a weakened recipe; use a
supported CASA 6 release" -- pointing at the CASA installation instead of at
the bad field selection.

There was a second hazard behind it. The wrapper retries the task in a loop
after dropping attributed keywords. Because every caller passes
`optional_keywords=()`, `reject_required` currently raises before any retry
can happen, so nothing re-runs today. But a single caller declaring an
optional keyword would arm a path that re-invokes a CASA task after a partial
failure -- for `flagdata` or `applycal` that means applying flags or writing
a calibration table twice.

Attribution is now restricted to `TypeError`, which is what Python raises for
an unexpected keyword argument, plus the strict message patterns. The loose
fallback is gone, and the signature-mismatch error chains the original
exception. Verified: a real CASA-style error propagates unchanged, a genuine
`got an unexpected keyword argument` is still attributed, and a true signature
mismatch still fails closed with the actionable message.

### 37. The uGMRT antenna pattern omitted the south arm

`antenna_pattern` was `^(?:C\d{2}|[EW]\d{2})$`. The array has a central
square plus *three* arms -- east (E02-E06), west (W01-W06) and south
(S01-S06) -- so six of the thirty antennas never matched. Telescope inference
requires at least half the names to match, and 24 of 30 still cleared that,
so a full array was identified correctly by luck. A sub-array dominated by the
south arm was not.

Fixed to `[EWS]`, and checked for collisions: LOFAR's `SE607` and the MeerKAT
`m###` names are unaffected, and a south-arm-only selection now resolves to
uGMRT.

### 38. The primary beam was evaluated at the wrong frequency for B4

With `reffreq` unset, CASA places the image reference frequency at the middle
of the *selected* frequency span. Flagging removes the channels outside a
receiver's specified range, but that does not change the selection, so the
reference frequency still reflects the full digitizer output. The measured
primary beam is then evaluated there -- at a frequency the surviving data do
not represent.

The affected profiles are exactly those where the digitized and science
ranges differ:

| band | digitized | science-valid | tclean default | correct centre |
|---|---|---|---:|---:|
| uGMRT B4 | 550-950 | 550-850 | 750.0 MHz | 700.0 MHz |
| MeerKAT UHF | 544-1088 | 580-1015 | 816.0 MHz | 797.5 MHz |
| MeerKAT L | 856-1712 | 900-1670 | 1284.0 MHz | 1285.0 MHz |

For a B4 observation covering the full 400 MHz GWB output, evaluating the
measured polynomial at 750 MHz instead of 700 MHz gives:

| sky radius | PB at 750 MHz | PB at 700 MHz | flux error |
|---:|---:|---:|---:|
| 10' | 0.829 | 0.850 | 2.5% |
| 15' | 0.650 | 0.689 | 5.9% |
| 20' | 0.454 | 0.506 | 11.4% |
| 24' | 0.310 | 0.365 | 17.8% |

That is a radial flux trend of exactly the kind Gate 5 tests for, and it
exceeds the 5% limit from about 15 arcmin outward -- most of the field.

The reference frequency is now set explicitly to the centre of the
science-valid *observed* band whenever the observation extends outside it, and
both the applied and default values are recorded. Where the observation lies
entirely inside the science range, or the band draws no such distinction, the
CASA default is already correct and nothing changes.

### 39. A sub-band could be planned wholly inside the removed range

`make_subband_plan()` divided the full channel list into equal chunks without
reference to the science range. For a B4 observation spanning 550-950 MHz with
six sub-bands, the top two chunks fall almost entirely above 850 MHz: images
with no unflagged data, and on the uGMRT path a hard failure when the
sub-band's reference frequency falls outside the measured beam model's
validity range.

Planning is now restricted to the science-valid range, with the retained and
total channel counts printed.

### 40. The self-calibration recipe fingerprint described calibration

`stage_recipe_fingerprint()` has explicit branches for `flagging` and
`imaging`; every other stage falls through to a calibration-shaped payload
carrying `CAL_CONFIG`, `CALIBRATION_TASK_CONFIG` and the telescope
calibration overrides. Nothing in it described self-calibration,
source detection or DD calibration.

That fingerprint has a specific job. Completed products are bound to
`stage_fingerprint`, which hashes every implementation file, so any code
change invalidates them. The *recipe* fingerprint is the narrower test that
decides whether a completed stage may be adopted under a new implementation
fingerprint -- the "this was only a bug fix" path. With self-calibration
settings absent from it, a change to solution intervals, imaging weights,
`pblimit`, or the amplitude acceptance gate looked exactly like a code-only
fix, and a completed self-calibration stage would be adopted as valid.

This audit is precisely that case. Findings 6, 13, 14 and the `pblimit` work
all change what self-calibration does; without this fix, resuming an existing
output directory after the upgrade would adopt the old self-calibrated
products and image them as though nothing had changed.

The fall-through payload now carries a `stage_specific` block: `SELFCAL_CONFIG`,
`SELFCAL_IMAGING_CONFIG` and `SELFCAL_TASK_CONFIG` for self-calibration,
`SOURCE_DETECTION_CONFIG` for source detection, and `DD_CONFIG` plus
`DD_TASK_CONFIG` for DD calibration. Verified that the resulting hash moves
when a single self-calibration setting changes.

Two things were already correct and are worth recording: the `flagging` and
`imaging` branches do capture their own configuration tables, so those stages
were never exposed to this; and the strict `stage_fingerprint` still binds
completed products to the source files, so nothing was ever reused *without*
passing through the adoption path.

**Correction (finding 49).** The block described above was never reached for
these three stages. Earlier in the same function, every stage other than
flagging, calibration and imaging returns the implementation fingerprint as
its recipe fingerprint, so the `stage_specific` payload only ever ran for
calibration, where it is empty. The verification above hashed a hand-built
payload rather than calling the function, which is why it did not show this.
The protection this finding claimed was in fact provided another way: an
implementation fingerprint also changes when a setting changes, which is
stricter, not weaker. The real cost was the opposite problem, described in
finding 49.

## Regressions this audit introduced, and how they were caught

### 41. The audit had been verified against 7 of 216 tests

Every stage above was checked with `test_ddcal.py` and
`test_pipeline_ddcal.py` — 7 tests. The repository also carries a `tests/`
directory of 209 tests that had not been run at any point during the audit.
It does not respond to `unittest discover` from the project root, because the
stage scripts deliberately `sys.path.pop(0)` to avoid shadowing the standard
library `inspect` module, and that is enough to hide it from a casual check.

Run against the pre-audit sources it is 209/209 green. Run against the audited
sources it was 209/202 — **seven failures, all mine**. Three were real code
regressions and four were fixtures that had not been told about a new product.

**Report-shape regression (finding 5).** Making `flag_edge_channels()` apply
both the guaranteed-range mask and the per-SPW roll-off cut also changed what
it returns: the combined report nested both sub-reports and dropped `scheme`
and `science_range_MHz` from the top level, where `edge_flagging_report.json`
readers and the validation contract expect them. The combined report now
starts from the guaranteed-range report and layers the roll-off mask on top,
so the published shape is unchanged and both masks are recorded.

**Report robustness regression (finding 2).** `write_report()` gained
`self.band_config["phase_solint"]` and `["amp_solint"]` as hard indexes, while
the two entries immediately above them correctly use `.get()`. A calibration
report is a diagnostic artifact and must not abort a completed run over a
descriptive field; the solve paths index these keys directly and fail far
earlier if a band profile is genuinely incomplete. Both now use `.get()`.

**Fixtures predating `phase_scan.cal` (finding 2).** The scan-averaged phase
table is a new, unconditional calibration product: solved in the main
sequence, re-solved in round 2, archived by `prepare_second_round()`, and
required by resume validation. Four fixtures still described the old
six-table set. They now carry it, and the round-2 sequencing test asserts the
transfer solve is part of the fresh solve rather than reused from round 1 —
which is the behaviour that matters, since mixing rounds is prohibited.

The honest lesson is about the verification, not the code. Every finding above
was individually reproduced and individually tested, and that produced real
confidence in each change while leaving the interactions between them
unmeasured. A regression suite that already existed would have caught all
seven immediately.

One pre-existing quirk was confirmed and left alone: naming both top-level
test modules on a single `unittest` command line fails to import the second,
because the first stage script to load removes the project directory from
`sys.path`. The baseline sources behave identically. Each module passes on its
own and under `discover`.

## Post-audit: flagging killed by SIGKILL

### 42. flagdata loaded a whole scan at a time, with no memory bound

A real run reported `Stage 'flagging' was terminated by SIGKILL (exit code
-9)`. Every automatic flagdata pass -- TFCrop, RFlag, extension, and the
corrected-data pass in calibration -- used `ntime="scan"`, and the remaining
calls used CASA's whole-scan default. flagdata holds one chunk per field and
SPW in memory, so the peak scales with scan length, channel count, correlation
count and array size together. None of that was bounded, and the pattern
predates the audit: it is in the original `flagging.py`.

With an estimated 12 bytes per visibility element (complex visibility, flag
cube, handler flag copy, headroom) against the 232 GB `slurm_mem` in
`config.yaml`:

| setup | whole-scan chunk |
|---|---:|
| uGMRT GWB 4k ch, 2.7 s, 40 min | 76 GB |
| uGMRT GWB 16k ch, 2.7 s, 40 min | 303 GB |
| MeerKAT 4k ch, 2 s, 20 min | 229 GB |
| MeerKAT 32k ch, 8 s, 20 min | 457 GB |

For 16k-channel uGMRT or 32k-channel MeerKAT data, one scan is larger than the
entire allocation. That is enough to explain the kill, but which dataset and
which step were involved is unconfirmed; `sacct` MaxRSS and the last
`[flagging] ... started` line decide it.

`flagging_selectors.py` now holds a CASA-free planner. It measures each
scan/SPW chunk from `msmetadata`, detects the tightest limit (Slurm allocation,
cgroup with ancestors, physical memory, optional configured value), keeps
whole scans wherever they fit within `flagdata_memory_fraction` of it, and
otherwise bounds `ntime` to the longest interval the most demanding scan
allows. When even four integrations cannot fit it raises before any flag is
written and reports the memory needed. Every chunked call in flagging and the
corrected-data pass uses it; quack is deliberately left on whole scans. The
settings live in `PIPELINE_RUN_CONFIG`, not `FLAGGING_CONFIG`, so changing a
memory setting never invalidates a flagging recipe.

`ntime` is not purely an execution choice: TFCrop fits time trends inside a
chunk, and extension's `growtime` is a fraction of the chunk. A bounded value
is therefore never silent -- it is written to
`flagdata_chunk_plan.json` / `postcal_flagdata_chunk_plan.json`, to every
iteration record, to the extension report, and to the post-calibration report.

Also bounded: `active_antenna_names()` read both antenna-ID columns whole and
converted them to Python lists (about 36 bytes per row per column); it now
reads 5M-row blocks and keeps only unique IDs. The opt-in ML flagger still
reads large contiguous DATA blocks and would exhaust memory on real data if
enabled; it is off by default and was left unchanged.

**Caught before shipping.** My first edit to the corrected-data report replaced
the first literal `"ntime": "scan"` in `calibration.py` -- which was inside the
new planner's own fallback branch. Whenever scan metadata could not be read,
the planner would have called itself until Python's recursion limit, while the
report line it was meant to change still said `scan`. The new tests raised a
`RecursionError` immediately; both lines are corrected and an AST check
confirms neither planner calls itself.

Verified: 14 new tests covering chunk planning (fits, bounded against every
constrained scan, insufficient, unknown limit, invalid fraction), limit
detection (smallest source wins, tighter cgroup ancestor, cgroup v1 sentinel,
configured limit only tightens), geometry reading, and both stages failing
before any flag is written. One existing post-calibration test asserted
`ntime == "scan"` and passed only through the no-CASA fallback; its plan is now
pinned. 230 tests pass.

### 43. Calibration read each calibrator's entire MODEL_DATA into memory

The next run passed flagging and was then killed by SIGKILL in calibration.
CASA's calibration tasks read one timestamp at a time, and the corrected-data
flagging pass is bounded by finding 42, so the remaining candidate was the
pipeline's own `validate_selected_model_data()`. After `setjy` or `ft` it read
the whole `MODEL_DATA` column for the calibrator field with a single `getcol`
and then checked every value in a pure-Python loop. It runs for every flux
calibrator and any distinct bandpass calibrator, and again at the start of a
resume once `setjy` has completed. The code predates the audit.

| calibrator data | values | whole read (8-16 B/value) |
|---|---:|---:|
| uGMRT GWB 4k ch, 2.7 s, 30 min | 5.1e9 | 38-76 GiB |
| uGMRT GWB 16k ch, 2.7 s, 30 min | 20.3e9 | 151-302 GiB |
| MeerKAT 32k ch, 8 s, 10 min | 20.4e9 | 152-305 GiB |

For the same high-resolution data that previously killed flagging, one field
alone meets or exceeds the 232 GB allocation; and even where it fits, a
Python loop over twenty billion values would run for days.

Validation now reads `MODEL_DATA` in row blocks of at most 8.4M values
(~128 MiB) and checks at most 67M values per field, counting with NumPy when
available. A field that fits the budget is still checked completely; a larger
field is checked in evenly spaced blocks that always include the first and
last rows, so a model failure confined to either end is still caught. Block
size comes from a 64-row probe using the widest row, so interleaved SPWs of
different widths stay bounded. The pass/fail rule is unchanged -- every value
read must be finite and at least one must be non-zero -- and the QA record now
states `coverage` (`complete` or `sampled`), `rows_checked` and the value
budget.

Nothing in `config.py` changed, so the calibration recipe fingerprint is
unchanged and a resumed run continues from its internal checkpoints.

Verified: the table fake now honours the real `getcol(name, startrow, nrow)`
API; four new tests use CASA-shaped `(ncorr, nchan, nrow)` arrays to check
complete coverage for a small field, sampled coverage with every read inside
the block bound and both end blocks present, a NaN in the final row, and an
all-zero model. 234 tests pass. The orchestrator's own reads were also checked
and are bounded.

### 44. A calibrator scan lost to flagging could not be traced to its cause

A uGMRT Band 4 calibration run logged `Insufficient unflagged antennas` for
the bandpass calibrator in whole channel blocks, and a phase-calibration solve
that lost first most antennas and then all data across one calibrator scan.
The solver is right to refuse; the question is which step removed the data,
and the answer was only partly recoverable.

The flagging stage already saved full CASA summaries -- with per-field,
per-scan and per-SPW breakdowns -- on arrival, after the fixed masks, after
each auto-flag pass and after extension, and inspection recorded which field
each scan belongs to and when it ran. The corrected-data pass in calibration
saved only whole-selection fractions, although its safety limits are field-
and SPW-level (15% for uGMRT), so one calibrator scan among ten can be removed
entirely while the pass is accepted. Nothing tied these products together.

`diagnose_flag_losses.py` is a new read-only tool that does. For every scan it
prints the flagged fraction at each saved step and names the step that removed
most of it; `--time` locates the scan containing a CASA log timestamp. For
every fully flagged channel run it prints the frequency range and classifies
it as a deliberate mask (edge or out-of-band selection, allowing one
neighbouring channel) or data-driven flags. Inside CASA it also counts, per
scan, whether the current and round-1 phase tables contain any usable
solution, which is how a loss introduced by the corrected-data pass shows up.
The corrected-data report now records per-field, per-scan and per-SPW flag
fractions before and after the pass, so future runs can attribute that step
directly.

One scientific note for interpreting an extension loss: CASA's `growfreq`
counts every existing flag, deliberate masks included. With about 30% of a
550-950 MHz window masked on purpose, the 80% rule fires once roughly 71% of
the usable channels in an integration are flagged.

Verified: helper tests (CASA-time round trip, selector expansion, channel
runs, classification with a neighbouring channel, loss attribution) and an
end-to-end fixture of a Band 4-style window in which a flux-calibrator scan is
already flagged on arrival, a phase-calibrator scan is removed by extension, and
the channel runs are two deliberate masks and one RFI block. 241 tests pass.

### 45. The memory fix made a correct calibration fail its completion check

A uGMRT Band 4 calibration ran to completion and was then rejected by the
orchestrator: "corrected-data flagging report does not preserve the
scan-local, non-extending detector configuration". Finding 42 lets that pass
bound flagdata's time chunk when a whole scan does not fit in memory, and
records the chunk actually used. `validate_stage_completion` still required
the report to say `ntime == "scan"`. The report changed; the check that reads it
did not. That is my regression from finding 42, not a data problem.

What the check protects is unchanged: the detectors must never combine data
across scans and must not extend flags. A bounded chunk with
`combinescans=False` never crosses a scan boundary, so it satisfies that. The
check now accepts `scan`, or a positive interval such as `919s` only when the
report's own `flagdata_chunk_plan` records a `bounded` decision with the same
value. It still rejects an interval with no plan, a plan that disagrees, a zero
interval, and `combinescans=True`.

Resubmitting reruns calibration from its internal checkpoints. Every solve,
both rounds, the corrected-data flagging and the split are already
checkpointed, so nothing is re-solved or re-flagged; the bounded MODEL_DATA
check and table QA run again before the stage is accepted.

Two tests added: a recorded bounded chunk is accepted, and each of the rejected
forms fails with the scan-local reason. The existing fixtures only ever wrote
`ntime: "scan"`, which is why the suite did not catch this.

## Post-audit: self-calibration and resuming across code updates

### 46. A partly written MODEL_DATA column passed the self-calibration check

`model_column_status()` decides whether the model column can be reused or
must be predicted again before a gain solve. It read up to 20000 rows in
2000-row chunks from the start of the table and stopped at the first chunk
holding a non-zero unflagged value. `tclean` writes MODEL_DATA in row order,
so a job killed during that write leaves the early rows filled and the rest at
zero. The check passed that column, and the next `gaincal` solved every later
time against a zero model. Confirmed against the previous code: a 100-row
column written to row 60, and a 1000-row column written to row 900, both
returned "usable".

The row limit also bounded memory only by row count. A uGMRT row with 2048
channels and two correlations holds 4096 values, so 20000 rows is about 1.2 GB
of complex data in one read.

The check now sizes reads by value count, as calibration does (finding 43):
at most 2^26 values in total and 2^23 per read. Blocks are spread evenly over
the whole table and always include its first and last rows. Any unflagged
non-finite value fails. Every sampled block that holds unflagged data must
contain a non-zero model value, so an interrupted write is caught at the end
of the table. Flagged rows are ignored. The function returns the same
`(ok, reason)` pair, and the failure message names the row range.

### 47. A killed self-calibration stopped and asked for manual cleanup

Every self-calibration step works on scratch products and marks a checkpoint
only after they are complete. A kill between the two left those products
behind, and the rerun then refused to continue:

- `phase<N>.candidate.ms` or `phase<N>.ms` without the `apply_phase<N>`
  checkpoint raised "Found an uncheckpointed phase-selfcal candidate";
- `final_selfcal.ms` without the `amp_selfcal` checkpoint raised
  "An uncheckpointed final selfcal MS already exists";
- a partial `initial.*` or `final.*` image set had no guard at all. `tclean`
  reuses existing images of the same name by default, so a partial set from
  the killed run would have seeded the new one.

The refusals were safe but needed a person to decide what to delete, and the
image case was not safe. No accepted product is ever changed by these steps: a
phase loop solves against the last accepted MS and applies the solution to a
copy, and the amplitude step's `applycal` rewrites CORRECTED_DATA from DATA
rather than adding to it. Redoing an unfinished step from the last accepted
product therefore cannot apply a solution twice.

Uncheckpointed products are now moved to `<name>.incomplete.<N>`, never
deleted and never overwritten, and the step is redone. This covers the loop's
candidate and split MS plus its `phase_<N>.*` images, `final_selfcal.ms`, and
the initial and final image sets. It follows the pattern calibration tables
and split targets already used. Each move is listed under
`archived_incomplete_products` in the self-calibration report. Checkpointed
products are left alone, as are other loops' images: `phase_1.*` never
matches `phase_10.*`. The guards inside the apply and commit steps remain as
a second line of defence. Archived Measurement Sets are full copies, so
remove them once reviewed if disk space is tight.

### 48. A code update made a correct calibration look like a recipe change

A uGMRT Band 4 run was stopped with "Stage 'calibration' was previously
attempted with a different implementation and its science recipe cannot be
proven unchanged (the stored science-recipe fingerprint changed)". Nothing
about the science had changed. Two faults combined:

1. **Adoption did not stick.** When a completed stage is adopted under a new
   implementation, `pipeline_state.json` is updated but the stage's own report
   keeps the old implementation fingerprint. On the next run the state matched
   the code, the report did not, and the stage was judged stale and re-entered.
2. **Re-entering flagging rewrote a file calibration hashes byte for byte.**
   A completed flagging stage does no flagging when re-entered. It rewrites
   `postflag_summary.json` with the new implementation fingerprint and appends
   a `provenance_migration_history` entry. The calibration recipe fingerprint
   hashes that file's bytes, so a code-only update changed the calibration
   recipe.

Fixes:

- A report's implementation fingerprint is accepted as current when the
  recorded adoption history leads from it to the current fingerprint, and the
  state's product fingerprint is the current one. An adopted stage is now
  reused on the next run instead of being re-entered. A fingerprint the history
  does not reach is still stale.
- When a stored calibration recipe does not match, the orchestrator undoes
  the post-flag report's recorded refreshes, newest first. Each refresh is
  one `json.dump(indent=4)` rewrite that changes only those two keys, so the
  earlier bytes are reproduced exactly. If one of these rebuilt versions, with
  the current configuration, arguments, MS identity and inspection report,
  reproduces the stored SHA-256, the recipe is proven unchanged. Any other
  change, such as a different flag count or calibration setting, still fails.
  The stored value is then used for the rest of the run, so the flux-model
  and solve checkpoints bound to it stay valid and nothing is re-solved.

The calibration, flagging and imaging recipe values are unchanged by this, so
existing state remains valid. Tests reproduce the failure with the real
flagging refresh code and show that a flag-count change or a calibration
setting change made at the same time still fails closed.

### 49. A killed self-calibration could never resume after a code update

Because of the early return described in the correction to finding 40,
self-calibration's recipe fingerprint was its implementation fingerprint.
After a kill, deploying any fix (including finding 47) made the stored and
current values differ. Self-calibration also had no branch in the resume
check, so it stopped with "no matching recipe record" and the accepted phase
loops could not be used.

Self-calibration now has its own recipe fingerprint. It binds the calibrated
target MS instance, telescope and band profile, target, reference antenna,
execution profile, `SELFCAL_CONFIG[band]`, `SELFCAL_IMAGING_CONFIG[band]`,
`SELFCAL_TASK_CONFIG`, and the inspection, calibration and split-target
reports with their implementation-fingerprint keys removed. A failed attempt
resumes from its checkpoints when that recipe matches. An attempt recorded
before this change stored the implementation fingerprint and still fails
closed. Tests check that the recipe differs from the implementation
fingerprint and ignores upstream implementation provenance. They also check
that it changes with an upstream calibration setting, a self-calibration
setting or the target, and that a killed attempt resumes after a code-only
update.

Source detection and DD calibration still use the implementation fingerprint
as their recipe. That is strict, and it means a failed attempt at either
restarts from scratch after a code update rather than resuming.

### 50. A SIGKILL could not be attributed to a cause

The uGMRT Band 4 self-calibration was killed in its first `tclean`, just after
CASA warned "Memory available 3121152 kB is very close to amount of required
memory 9701328 kB". The job had a 232 GB allocation. Three SIGKILLs in this
run (flagging, calibration, self-calibration) all arrived with the same
generic message, and nothing on disk could say whether the job had exceeded
its limit or been killed from outside it.

**What the CASA and casacore sources say.** Read directly from the casacore
and CASA sources:

- casacore's `HostInfo` takes a job memory limit only from a cgroup v1
  `memory.limit_in_bytes` at the process's own cgroup level. Under cgroup v2,
  or when Slurm puts the limit on an ancestor cgroup, it falls back to the
  node's `MemAvailable`. The "Memory available" figure was therefore most
  likely the node's available memory, not the job's.
- W-projection builds its convolution functions on a grid of the padded image
  size rounded up to an even composite: 10000² for an 8192² image. Each OpenMP
  thread holds a complex plane plus an FFTPACK work array, about 1.5 GiB. Planes
  are processed in chunks of at most 100, and `omp_set_num_threads` is
  lowered to the chunk size. The chunk size shrinks when CASA believes little
  memory is free: 5 planes at its 4 GB floor.
- For 128 planes, the chunks are 28 and 100 and the peak is about 62 GiB
  with 64 allocated threads, or about 10 GiB when CASA sees 3 GB free. The
  image store adds about 9 GiB. Neither is close to 232 GB.

So this image, with the settings in use, should not have reached the job's
limit. CASA reported about 3 GB free before `tclean` had allocated anything
large, and the self-calibration process does no large reads before it (checked).
The evidence points to memory being unavailable to the job rather than used
by it. The most likely explanation is other use of a shared node, but the
kill could not be attributed from what the pipeline recorded. That is the
defect fixed here; the ilifu node's cgroup version and the job's `MaxRSS`
still need to be checked.

**Fix.** A new CASA-free module, `memory_diagnostics.py`, reads
`/proc/meminfo` and every cgroup level from the process up to the root, for
both v1 and v2: limit, usage, peak, out-of-memory kill count and limit hits.
The job limit is the tightest limit on that path.

- Before each stage the runner prints the job limit, job usage and the node's
  available memory. It warns when the node has less than half of the job's
  remaining allowance, because a job limit does not keep that memory free on
  the node.
- While the stage runs, the progress loop records the lowest node
  availability and the highest job usage.
- On SIGKILL the runner compares the out-of-memory kill counters of the job's
  own cgroup levels. Levels above the job are excluded, because under v2 they
  also count other jobs. It states one of four causes: the job exceeded its
  limit; the kill came from outside the job (no kill counted against the
  limit); undetermined (the kernel reports no kill counts); or no job limit
  is visible at all. The message includes the node figures and the killed
  process's peak resident memory from `getrusage`, when this stage set it.
  The evidence is written to `diagnostics/<stage>_sigkill_<UTC>.json` and the
  cause is recorded in `pipeline_state.json`.

**Outcome.** The first run with this evidence recorded one out-of-memory kill
in the job's own cgroup, whose limit was **3.0 GiB**, with a job peak of
3.0 GiB. The reading above was wrong: CASA's "Memory available 3121152 kB" was
the job's limit, not the node's free memory, and the job never had 232 GB. A
3 GiB limit is consistent with all three kills. It also explains why
corrected-data flagging had to use bounded chunks (finding 45): with 232 GB,
whole scans would have fitted. Findings 42 and 43 made flagging and
calibration fit into that limit; self-calibration imaging cannot, as finding
51 shows.

Reading these files can never fail a stage. Nothing about the science changes.
Tests cover v1 and v2 layouts, the unlimited sentinel, another job's kills on
the same node, a single kill counted at two v2 levels, an old kernel without
kill counts, the node warning, peak RSS attribution, and the stage runner end
to end with a killed process.

### 51. A stage that could not fit was started anyway

Self-calibration spent its start-up and then died, because nothing compared
its minimum memory with what the process was allowed.

For a serial `tclean` with `gridder='wproject'`, the CASA source fixes a
floor that does not depend on how CASA reads free memory. The padded image
(8192 × 1.2) is rounded up to a 10000-pixel square grid. Planning never
assumes less than 4000 MB, and at most 100 w-planes are computed per chunk.
Chunks run last first, and each lowers the OpenMP thread count to its own
plane count. `wproject_imaging_memory()` evaluates the convolution-function
peak for every chunk size CASA could choose, then adds 4-byte residual, model
and PSF images and one padded complex grid. For the uGMRT B4 self-calibration
image (8192 px, nterms=2, 128 w-planes) the result is:

| Threads | Minimum | Largest |
|--------:|--------:|--------:|
| 1  | 6.4 GiB | about 24 GiB |
| 64 | 8.1 GiB | about 111 GiB |

The largest figure applies when CASA plans from a large memory figure with
many threads. Neither figure counts CASA's other working memory, so the
minimum is a floor, not an estimate of use.

Before every stage the orchestrator now:

- warns when the process's cgroup limit is under half of the Slurm job's
  request (`SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` x
  `SLURM_CPUS_ON_NODE`). That points to a job step or session with a smaller
  share, such as `srun` or an interactive shell. It also warns when a
  cgroup limit applies with no Slurm job at all, as on a login node;
- for self-calibration, prints the minimum and largest memory. It stops
  without starting CASA when the limit (or, with no cgroup limit, the node's
  total memory) is below the minimum, and warns when it is below the largest.
  A dry run reports the same thing without stopping.

Only self-calibration is checked: its imaging is always serial w-projection
from settings the orchestrator reads directly. Final imaging varies by
gridder and MPI layout and is not covered. Tests cover the composite grid
sizes, the B4 floor, a run without w-planes, the three allocation cases, and
a self-calibration refused before CASA starts, which leaves the resume state
unchanged.

### 52. Re-entering calibration made a failed self-calibration unrecoverable

After the memory limit was corrected, the calibration stage ran again (9.8
min) and self-calibration then stopped with "previously attempted with a
different implementation and its science recipe cannot be proven unchanged".
There were two faults:

1. **The self-calibration recipe hashed the whole calibration report.** Much
   of that report is specific to one run: the implementation and recipe
   fingerprints, the working MS's file-stat hash (sizes and modification
   times), checkpoint records, the QA plot path, and the flux model's
   migration history. Re-entering calibration changes them even when every
   solution is identical, so any calibration re-entry changed the
   self-calibration recipe.
2. **A changed recipe after a failed self-calibration was a hard stop.** The
   message asked for a new output directory, which means repeating flagging
   and calibration. Here the failed attempt had no completed checkpoint at
   all.

Fixes:

- The self-calibration recipe now reads only the calibration report's
  science content: MS, instrument, calibrators, targets, reference and bad
  antennas, delay quarantine, Hanning smoothing, calibration parameters,
  table paths and per-table QA, split-target identities, and the flux model
  identity without its pipeline fingerprints. A rebuilt target MS or a
  changed table or setting still changes the recipe; the run-specific fields
  above do not. The payload's schema version is now 2.
- When a failed or interrupted self-calibration meets a changed recipe, the
  orchestrator moves `selfcal/` to `selfcal.superseded.<N>`, records the move
  and both fingerprints under `superseded_attempts`, and runs
  self-calibration from the start. This is safe for self-calibration and not
  for calibration. Calibration writes flags and corrected data into its
  working MS, so an attempt under another recipe has already changed its own
  input. Self-calibration keeps every product under `selfcal/`. The only
  change it makes to the calibrated target is to MODEL_DATA, which the fresh
  initial image rewrites before any gain solve. Nothing from the old attempt
  is reused and nothing is deleted. A dry run reports the move without making
  it. Calibration still stops on a recipe change.

Checked while doing this: self-calibration's MODEL_DATA writes do not make
calibration consider a split target stale. That identity binds the MS
directory inode and its fixed subtables, not MAIN.

Tests: run-specific calibration fields leave the recipe unchanged, while a
changed table path or split-target identity changes it. An attempt from an
older recipe, and one with a changed self-calibration setting, are set aside
and restarted, with earlier archives kept. A dry run moves nothing.

### 53. A correct primary-beam image was rejected for a changed timestamp

uGMRT imaging completed and was then rejected: "uGMRT primary beam no longer
matches PB provenance". Imaging binds its restored, primary-beam and
PB-corrected images with `path_stat_fingerprint`, a hash of each file's name,
size and nanosecond modification time. The beam's fingerprint was taken a
few seconds after `image.fromarray` wrote it to `/idia`. Nothing in imaging
writes to the beam afterwards. It is only read, by `impbcor`, `exportfits` and
the common-beam correction; the geometry descriptor only opens and closes it.
The restored image, written long before its fingerprint was taken, passed the
same check a moment earlier. The likely explanation is that the network
filesystem reported different modification times for the freshly written
files later than it did straight after the write. That is inferred, not
observed; see "still open" below.

The same pattern applied to every beam and PB-corrected image, for uGMRT,
MeerKAT katbeam and uGMRT sub-bands, and to the published FITS inventory,
whose check compared recorded modification times too.

Fixes:

- `imaging_safety.casa_table_content_identity()` hashes the contents of every
  file in a CASA table directory except lock files, per file and overall, and
  ignores timestamps. Imaging now records `source_`, `pb_` and
  `pbcor_content_identity` beside the existing stat fingerprints in every
  beam and PB-correction provenance file.
- `artifact_matches_provenance()` compares by content when a content identity
  is recorded, and by the stat fingerprint only for older provenance. A
  content mismatch names the files added, removed or rewritten. The
  orchestrator's uGMRT and MeerKAT imaging checks and imaging's own reuse
  checks all use it.
- A primary beam that no longer matches its record is rebuilt instead of
  stopping the stage with "Use --force or a fresh output directory". The beam
  is computed from the restored image and a fixed model, so rebuilding it is
  always correct. The old one is moved aside, and every PB-corrected product
  after it is redone, because correction clears its dependent checkpoints.
- Published FITS identities are compared by path, size and SHA-256. The
  recorded modification time is kept for reference only.

Resuming the failed run needs no manual step. The imaging checkpoint keeps the
restored image; the beam's old stat-only record no longer matches, so the
beam is rebuilt and correction and FITS export rerun, now with content
identities.

**Superseded in part.** The timestamp explanation above was wrong; see
finding 55. Measurement Sets are still bound by stat fingerprint: the
calibration working MS, self-calibration input and final MS, and imaging
input. Hashing them is not practical. If modification times on `/idia` also
change for them with no write, the affected stage is judged stale on a later
run. That would explain calibration's unexplained 9.8-minute rerun (finding
52), but it has not been confirmed. The next step is to compare each report's
recorded MS fingerprint with the live one on the existing products.

Tests: new timestamps leave a content identity unchanged; rewritten and added
files are named; lock files are ignored; stat-only provenance is still
compared by stat; FITS identities ignore modification time; a beam with new
timestamps is reused while one with a stale stat-only record is rebuilt; and
the MeerKAT validator accepts moved timestamps but rejects a rewritten beam.

### 54. A stale copy of a stage script failed partway through the stage

After finding 53 the run failed with "Existing uGMRT PB pixels changed after
provenance was recorded", an error that code no longer contains. Every line
number in the traceback matched the previous `imaging.py`, so the cluster was
still running the old file. Files are copied to the cluster one at a time, and
a stale copy shows up only as a failure from old code, after the stage has
already started.

`imaging.py` and `imaging_safety.py` now declare `PIPELINE_INTERFACE_VERSION
= 2`. `pipeline.py` lists the minimum it needs (`SOURCE_INTERFACE_VERSIONS`)
and checks it first in `preflight`, before any stage or dry run. A missing,
unmarked or older file is named, and nothing runs. Pointed at the previous
copies, it reports both files. The version must be raised whenever
`pipeline.py` starts relying on something new in a listed file. Tests cover
the shipped tree, and older, unmarked and missing copies. They also check
that preflight stops before any stage.

### 55. CASA rewrites an image's metadata file, so content identity was too strict

With finding 53 deployed, imaging was rejected again, now with the evidence
the check was built to give: "uGMRT primary beam no longer matches PB
provenance: its contents changed: rewritten table.dat". Only `table.dat`
changed. The pixel files did not. So finding 53's reading, that the network
filesystem reported different modification times, was wrong: the bytes really
change.

`table.dat` holds a casacore table's description and keywords. casacore
rewrites it whenever a keyword changes -- units, image info, a mask entry, the
coordinate record -- and never for a read. Nothing in the pipeline writes to
the beam after it is built; every later use opens it, reads it and closes it.

To find the writer, CASA 6.7.5-18, the same build the cluster runs, was
installed as a standalone Python package here and the beam flow was replayed
against imaging.py's own methods: build the beam from a restored image with
`image.fromarray`, take the geometry descriptor, publish it with `os.replace`,
record the identity, then run the real `run_ugmrt_primary_beam_correction`,
`impbcor` again with the same beam, `exportfits` of the beam and the corrected
image, `image_metrics` (`imstat` and `imhead`), the QA plot read, the
deconvolution diagnostics, and a check from a separate process after the
first exited. The beam stayed byte-identical at every step. The `logtable`
keyword is stored as a relative reference (`././logtable`), so the publish
rename is not the cause either. **The writer was not identified.** Candidates
not reproducible here are MPI (the stage runs under `mpicasa`) and the
cluster's filesystem.

The fix does not depend on identifying it. Provenance now binds the files that
hold the data, not the metadata: every file except `table.dat`, `table.info`,
`table.lock` and the `logtable` subtable. A rewritten `table.dat` is accepted;
rewritten pixels, or an added or removed data file (a mask, for example), are
still rejected and named. A table with no data files at all is compared whole,
so the degenerate case fails closed.

Coordinate metadata is not left unchecked. Provenance records the source and
beam geometry, the orchestrator requires the two to agree, and imaging re-reads
the beam's geometry live and compares it with the restored image before
reusing a beam.

On the next run the beam is reused rather than rebuilt: its pixels are
unchanged, so PB correction reruns from the existing beam.

### 56. Only one science target was ever imaged

Calibration splits every science target, but everything after it ran for one:
`selfcal`, `imaging`, `source_detection` and `ddcal` each wrote to a single
product directory, and with several targets present the orchestrator ranked
them and picked one. An observation of several fields produced one image.

Per-target stages now run once per target:

- `science_target_scopes()` returns the targets and the directory each one
  writes to. The default is every calibrated target; `--target NAME`,
  `--targets A,B`, `--target-selection auto` (rank and image one) and
  `explicit` remain available. An unknown target name fails closed.
- Each target's products go to `Products/targets/<name>/`, with the name made
  filesystem safe. A target already imaged directly under `Products/` keeps
  that location, so a run made before this change is validated and reused
  instead of repeated.
- The orchestrator enters a target scope for the per-target stages:
  report paths, product paths, resume-state keys (`selfcal@<target>`) and
  fingerprint caches all follow the scope, while inspection, flagging and
  calibration stay shared at the run's root.
- `selfcal.py` and `imaging.py` take `--upstream-outdir` for the shared
  inspection and calibration reports; everything else they read and write is
  under their own `--outdir`.

Failure isolation follows from the per-target state: one target failing leaves
the others' completed products validated and reusable.

### 57. Interference survived into the final image

The final PB-corrected image of the Band 4 run showed residual RFI. Two
passes that the recipe never performed explain it:

- **Flags were never grown.** The corrected-data pass in calibration runs
  `tfcrop` and `rflag` with `extendflags=False` and nothing afterwards, so
  the wings of each burst, which sit just under the detection threshold, were
  kept. A separate scan-local `mode='extend'` pass now grows accepted flags
  over neighbouring channels and times (80% by default, `postcal_extend`),
  before the round-2 solutions are solved from the flagged data. The
  detectors keep their scan-local, non-extending configuration; the
  orchestrator now also requires the extension to be scan-local.
- **Nothing looked at the residuals.** `selfcal.py` and `imaging.py` ran no
  flagging at all, so interference comparable to or below the source
  brightness was never separable from the field itself. Self-calibration now
  flags on the residual (CORRECTED or DATA, minus MODEL) once, after the last
  accepted solution and before the final image, using the model from the
  latest self-calibration image. Because the field's own model is subtracted,
  the detector sees interference rather than real emission. The flags are
  then grown within each scan.

The residual pass is bounded like calibration's: the flag state is saved
first, the new fraction is measured, and the whole pass is restored when it
removes more than `residual_flag_max_new_fraction` (10%), which is what a
wrong model or a mis-scaled solution would produce. The outcome, its
parameters and both fractions are written to
`selfcal/reports/residual_flagging.json` and into the self-calibration report,
and the step is checkpointed so a resume does not repeat it.

Both are recipe changes, so they apply to targets that are self-calibrated
and imaged from now on. A completed stage is adopted as it stands; to apply
them to one already imaged, set its `selfcal/` directory aside and rerun.

### 58. Three validators still looked in the run root, and one target could stop the rest

The first multi-target run reached source detection for `G0212-0335`, which
ran correctly in that target's own directory and recorded a skip because
PyBDSF is not installed in the CASA environment. The orchestrator then failed
the stage: "missing completion marker:
.../uGMRT_Products/source_detection/checkpoints/source_detection.done", the
run-root path rather than the target's. Finding 56 scoped the paths written on
one line but missed three built across several lines: the source-detection and
DD-cal completion markers, and the expected imaging FITS product. All three
now use the target's product root, and a test asserts that a validator looks
inside the target directory.

The recorded skip itself was already handled: a source-detection report with
status `skipped` or `disabled` and a reason validates, so a missing optional
dependency does not fail the run once the marker is found where it was
written.

The same run showed a second problem. One target's stage failing raised
immediately, so the remaining targets were never attempted, which defeats the
purpose of imaging them all. Targets are independent: each keeps its own
products and resume state. `run_per_target_stages()` now records a failure,
continues with the next target, and raises at the end naming every target that
did not complete, so the run still fails and resubmitting continues from what
each target has. A single-target run raises directly, as before.

### 59. An optional pass ended the stage, and it exposed that nothing was being deconvolved

Eight of nine targets failed self-calibration at finding 57's residual pass:
`tclean did not produce a usable MODEL_DATA column ... MODEL_DATA is zero
over the unflagged data`. Above the traceback CASA had already said it:
"Restoring with an empty model image. Only residuals will be processed."

Two separate faults.

**The pass could end the stage.** Residual flagging is a quality step; the
images are produced with or without it. It now runs inside a guard that
records any failure in `selfcal/reports/residual_flagging.json` and warns,
instead of raising. It also refuses to run when the image it would subtract
has no clean components: flagging a "residual" that is really the data
itself would let the detector see real emission. `latest_image_base()` only
checked that a model image existed; an empty one exists too.

**Nothing was deconvolved.** The failure is a symptom. For these fields
CLEAN found nothing to clean, so:

- the sky model is empty and self-calibration has nothing to solve against.
  For `G0231-0533`, `current_ms` was still the calibrated input at the end of
  the phase loops, meaning no solution was ever accepted;
- the published image is the dirty map. Sources are not deconvolved and
  sidelobes are not removed, which is what the reported "lots of RFI" in the
  final PB-corrected image looks like.

Self-calibration now says so at the end of the stage, with the residual RMS,
peak and peak/RMS, and records `sky_model_empty` in its report. The stage no
longer completes quietly while having done nothing.

Why CLEAN finds nothing is not settled here and is a property of the data
reaching imaging, not of the deconvolution settings. The candidates, in the
order the evidence should be read: the flagged fraction per field (44.8% on
the one field reported so far); the calibration-table QA sidecars; whether
the phase calibrator brackets each target in time, given the
"Insufficient unflagged antennas" and "Found no unflagged data" messages seen
on phase-calibrator scans; and interference inflating the noise estimate so
that the auto-multithresh mask never opens, which would be self-perpetuating,
since the residual pass of finding 57 needs a model before it can help.

### 60. Self-calibration decided to do nothing and did not say why

The run reported, in full: "Skipping phase selfcal", "Skipping amp+phase
selfcal", and the residual pass skipped for want of a model. Each message was
correct and none of them said what had been measured or why, so a stage that
did nothing looked like a stage that had run.

Three changes, none of which alter what is solved or applied.

**The decisions state their evidence.** Phase self-calibration is skipped when
the initial image's dynamic range is below `phase_start_dynamic_range` (20).
The message now gives the measured dynamic range, peak and RMS against that
threshold. The amplitude message distinguishes "disabled by configuration",
which is the default, from a dynamic range below `amp_start_dynamic_range`.
Every deconvolution now prints whether a model was produced, with the
residual RMS, peak and peak/RMS.

**A field that cannot be masked is now recovered.** With poor uv coverage the
auto-multithresh sidelobe criterion (`sidelobethreshold` x max sidelobe x
peak) can exceed the image peak, so no mask ever opens and CLEAN produces no
model however bright the field. The bounded relaxed-mask retry existed for
this but only fired above 8 sigma, which excludes exactly the marginal fields
that need it. It now fires at 5 sigma, the weakest feature CLEAN can
legitimately deconvolve, and reports whether the retry produced a model.

**A coherence check separates a faint field from incoherent data.**
Calibrated visibilities of a real field add up coherently, so the
vector-averaged amplitude is close to the scalar-averaged one on the samples
a short observation provides. When calibration has not transferred, the
phases are effectively random: the vector average collapses while the scalar
average does not. A faint field lowers both. Self-calibration now measures
both from bounded, spread-out samples of the calibrated target data before
imaging, prints them with the ratio, and warns below
`coherence_warn_below` (0.05) that the evidence points at the calibration
for this target rather than at the field. It is diagnostic only: it reads
nothing but the data, changes nothing, and cannot fail the stage. The result
is recorded as `calibrated_coherence` in the self-calibration report.

What this does not do is make self-calibration solve against a field that is
not there. If the coherence is low the fault is upstream, in the calibration
transferred to these targets, and that is where the next run's evidence will
point.

### 61. The published image was never deconvolved

The PB-corrected image of a target shows several compact sources, each with
the radial spike pattern of the point-spread function, on a background of
arcs. That is a dirty map. The sources are real and detected, so the data are
calibrated well enough to see them: this is not interference and not a
flagging or calibration failure. Nothing was deconvolved, which is what the
"Restoring with an empty model image" warnings of finding 59 were reporting.

The mechanism is auto-multithresh on sparse uv coverage. Its mask threshold is
`max(sidelobethreshold x max_sidelobe x peak_residual, noisethreshold x rms)`.
A snapshot point-spread function like the one in that image has sidelobes at
tens of percent of the peak, so with the configured `sidelobethreshold` of
2.5 the first term lands at or above the image peak. No mask opens, CLEAN
stops immediately, the model stays empty, and the dirty map is published.
Self-calibration is then skipped for want of a model, and the sidelobes that
self-calibration and deconvolution would have removed stay in the image. The
relaxed-mask retry of finding 60 helps only when the relaxed threshold drops
below the peak, which is not guaranteed.

Both self-calibration and final imaging now fall back, after the relaxed
retry has also produced nothing, to cleaning the compact emission inside the
primary beam directly: `usemask='pb'` with `pbmask` 0.2, to
`automask_fallback_nsigma` (6 sigma) and at most `automask_fallback_niter`
(5000) iterations. The depth and the iteration cap are what keep it safe: it
cleans what is plainly above the noise and cannot follow sidelobes down into
it. It runs only when the residual holds a feature at or above
`automask_fallback_min_peak_snr` (5 sigma), only when no mask opened, and can
be turned off with `automask_fallback_enabled`.

The consequences run in the right direction: with the sidelobes deconvolved
the residual RMS falls, so the dynamic range measured for the
self-calibration decision rises, and a field that sat below
`phase_start_dynamic_range` on its dirty map may now pass it. Where the
sources are genuinely too weak, the stage still declines and says so.

### 62. An interrupted stage looked like a crash

A run stopped with a bare `KeyboardInterrupt` traceback from inside
`tclean`'s `makePSF`, printed by casashell from wherever CASA happened to be.
Nothing had gone wrong in the pipeline: an interrupt reaches a batch stage
when the job is cancelled, when it reaches its time limit, or when an
interactive session ends. The traceback said none of that, and the
orchestrator recorded it as "failed with exit code 1", which is what a real
failure looks like.

The long stages (flagging, calibration, self-calibration, imaging) now catch
the interrupt, say what it means and that their checkpoints are kept, and
exit 130, the conventional code for SIGINT. The orchestrator records the
stage as `interrupted` rather than failed and raises `StageInterrupted`.

That exception is deliberately not a `RuntimeError`: the per-target loop of
finding 58 continues past a target that fails, which is right for a failure
of one target and wrong for a cancelled job. An interruption stops the whole
run, and the remaining targets are left untouched for the next submission.

### 63. Every target crashed with SIGSEGV after the new clean passes

All nine targets died with SIGSEGV in self-calibration on the run after
finding 61 added the primary-beam fallback. A segmentation fault is a crash
inside CASA's C++ code, not a Python error, and the pipeline could not say
which of the two plausible causes it was, because it gathered memory evidence
only for SIGKILL.

**Two candidates, both addressed.**

The new code path was exercised for the first time on this run: the relaxed
retry (which finding 60 lowered to 5 sigma, so it now fires) and the
primary-beam fallback both restart tclean over an existing image set with
`calcpsf=False` while changing how the mask is made. A `.mask` written by the
previous pass then meets a pass that builds its mask differently. Both passes
now delete that mask first, so each builds its own, and the fallback runs only
when a primary-beam image exists for `usemask='pb'` to read, which
`calcpsf=False` does not create.

The job's memory limit was 28 GiB, and the convolution-function work of
finding 51 can need more than that at this image size; an allocation failure
inside CASA can surface as a crash rather than as a clean out-of-memory kill.
Memory evidence is therefore now recorded for a stage that dies on any
signal, not only SIGKILL: the job limit, peak and usage, the highest sampled
usage, the node's available memory, the cgroup's out-of-memory counters and
the largest process's resident memory, all written to
`diagnostics/<stage>_sigkill_<UTC>.json` and summarised in the error.

Which of the two it was is answered by that line on the next run: an
out-of-memory count above zero, or a peak at the limit, points at memory; a
crash with the job well inside its limit points at the clean passes, which
can then be turned off with `automask_fallback_enabled` and
`automask_retry_enabled` in `config.yaml`.

### 64. The coherence check was the wrong statistic, and said so about good data

Finding 60's check reported, for a uGMRT Band 4 target: scalar-averaged
amplitude 0.3435 Jy, vector-averaged amplitude 0.0045 Jy, ratio 0.013, and
warned that the calibration had not transferred. The data are fine; the
statistic was wrong.

The scalar average of a target field is set by the per-visibility thermal
noise, here 0.34 Jy, not by source flux. The vector average is the flux at
the phase centre, of order mJy. Their ratio is therefore ~1% for any target
field, well calibrated or not, and only approaches 1 for a bright point
source at the phase centre, which is a calibrator. The 0.05 threshold would
have warned about almost every target ever imaged.

The question the numbers do answer is whether the vector average stands above
its own uncertainty. With N samples of per-component noise sigma, that
uncertainty is sigma/sqrt(N). For this target: sigma = 0.3435/sqrt(pi/2) =
0.274 Jy, N = 80698538, so the uncertainty is 31 uJy, and the measured
4.45 mJy is 146 sigma from zero. Visibilities with random phases cannot do
that. The calibration transferred correctly, and the same arithmetic gives
the thermal noise the data can reach, about 31 uJy/beam.

The check now reports the per-visibility noise, the coherent flux with its
uncertainty and significance, and the thermal noise implied, and judges
against `coherence_min_significance` (5 sigma). It is explicitly one-sided:
a high significance proves coherence, while a low one is reported as
inconclusive, because a field with little flux at its phase centre looks the
same. The ratio is still recorded, labelled as what it is.

### 65. CASA planned its memory from the node, not from the job

The SIGSEGV of finding 63 came with the evidence it was built to produce, the
same on all eight targets: job limit 28.0 GiB, cgroup peak 28.0 GiB, node
available 176.8 GiB, no out-of-memory kill counted, largest process about
15 GiB resident. The job reached its ceiling without the kernel killing it,
which is what an allocation failing inside CASA looks like.

casacore reads a job's cgroup limit only from a v1 hierarchy (finding 50) and
otherwise reports the node's free memory. `WPConvFunc::findConvFunction` sizes
its w-projection work from that figure: it computes convolution functions in
chunks of up to 100 planes, each plane a quarter of a padded grid. For an
8192-pixel image the padded grid is 10000 square, so with 176 GiB apparently
free it plans 100 planes in one cube of 18.6 GiB, on top of the roughly 9 GiB
image store and about 6 GiB of per-thread screens. Against a 28 GiB limit
that cannot fit. Told the job's real limit it plans 18 planes and 3.4 GiB.

The orchestrator now writes `casa/<stage>_resources.rc` for each CASA stage
with `system.resources.memory` set to `casa_memory_fraction` (half) of the
limit this process may actually use, and `system.resources.cores` set to the
stage's thread count, and points `CASARCFILES` at it.
`HostInfo::memoryTotal(true)` reads exactly those keys. casacore takes the
first match across the files named there, so the run's file goes first and
the user's own `.casarc`, `.casa/rc`, `.aipsrc` and the CASAPATH files follow
in casacore's own order, none of them lost. An existing `CASARCFILES` is kept
behind the run's file. Setting the fraction to 0 leaves CASA's view alone.

This bounds what CASA plans; it does not make 28 GiB the right allocation for
imaging at this size. Both are worth fixing.

### 66. The crash was diagnosed three times from outside CASA

The SIGSEGV of findings 63 and 65 recurred after the memory budget was put in
place, with an identical signature. Each round of diagnosis worked from the
orchestrator's own view: exit status, cgroup counters, resident sizes. That
view cannot distinguish an allocation failure in the w-projection cube from
any other fault inside CASA, and it cannot show whether
`system.resources.memory` was honoured, because only CASA reports the figure
it planned from.

CASA writes both to its own log, and that log was being left on the cluster
for someone to go and find. The orchestrator now collects it: when a stage
dies from any signal, it takes the newest `casa*.log` touched since the stage
started, from the working directory or either output directory, keeps its last
200 lines in `diagnostics/<stage>_sigkill_<UTC>.json` beside the memory
evidence, and prints the last 25 to the run's output with the failure. The
line `Memory used in gridding function = X MB from maximum Y MB` states what
CASA budgeted, so whether finding 65's file took effect is answered by the
next failure rather than by another round trip. The capture cannot mask the
failure it documents: a missing or unreadable log yields no record and the
stage still fails as it did.

### 67. Self-calibration could not be made to fit a smaller job

Finding 65 bounds what CASA plans for; it does not reduce what the imaging
actually needs. At uGMRT B4 defaults, self-calibration images 8192 pixels at
0.5" with 128 w-projection planes, and the padded convolution grid, which
carries the largest allocation, grows as the square of imsize. Nothing in
`config.yaml` could change any of it: the only self-cal keys exposed were the
two dynamic-range thresholds and `pblimit`, so fitting the stage into a 28 GiB
job meant editing `config.py` on the cluster.

`selfcal_imsize`, `selfcal_cell`, `selfcal_wprojplanes`, `selfcal_niter` and
`selfcal_robust` now override the self-calibration imaging defaults for every
band. They are validated at load (a size must be a positive integer), and they
deliberately do not touch `FINAL_IMAGING_CONFIG`: self-calibration images only
have to carry a sky model good enough to solve against, so a smaller field
there costs nothing in the published image. Halving imsize to 4096 quarters
the convolution grid.

### 68. The pre-flight check reported a range, not what would happen

Before finding 65 the pre-flight print could only give a range, because CASA
decided its own w-projection chunk size from a figure the orchestrator did not
control: at uGMRT B4 defaults, at least 8.1 GiB and up to about 39 GiB. Now
that the run sets `system.resources.memory`, the chunk size is determined, so
the check computes it. For the observed job it reports 18 w-planes at a time
and about 10.3 GiB, against a 28 GiB limit.

That figure is what makes finding 66 conclusive: if the casarc took effect,
w-projection is using a third of the job's memory and the SIGSEGV is not the
convolution allocation. The warning that follows a plan too large for the job
now names both levers with their arithmetic, `selfcal_imsize` and `--threads`,
rather than only threads.

The model shows the peak is not monotonic in the budget, and it is worth
knowing why before reading these numbers: chunks are built last first, and
each sets the OpenMP thread count to at most its own plane count for the rest
of the process. 128 planes in chunks of 18 leaves a final chunk of 2, which
clamps the build to 2 threads and 10.3 GiB; chunks of 13 leave a final chunk
of 11 and cost 22.8 GiB. A larger allowance can therefore plan a larger peak,
which is one more reason to state the figure rather than a bound.

### 69. auto-multithresh cannot open a mask on the map it is shown

A published PB-corrected image came back with every source wearing its full
point-spread function: radial spokes, concentric rings and negative bowls,
at the same brightness as the sources. That is a dirty map. Findings 59, 61
and 63 each treated a symptom of it; this is the mechanism.

auto-multithresh opens a mask at

    max(sidelobethreshold * max_sidelobe * peak, noisethreshold * rms)

and measures both terms on the residual, which on the first pass is the dirty
map. In a field with several comparable sources that map's rms is sidelobe
confusion, not noise. For the observed image, sources near 1.5 mJy/beam over a
confusion rms near 0.3 mJy/beam give `4.25 * rms` = 1.3 mJy, essentially the
peak itself. No mask opens, so `niter` is never spent, the model stays empty,
and the restored image is the dirty map. Self-calibration then has nothing to
solve against and skips, which is the report the observer actually saw.

Every recovery this pipeline had ran *after* that measurement and inherited
its premise. The fix runs before it. A bounded first pass cleans what is
plainly there, `usemask='pb'` inside pb>0.2 at `nsigma=5` with a capped
iteration count, writing no model column. It cannot run away, because nsigma
is re-measured from the residual at each major cycle: the first cycle reaches
only the brightest peaks, and as their sidelobes come out the threshold
follows the noise down. The auto-multithresh pass then continues from that
residual (`restart=True`, `calcpsf=False`, `calcres=False`), its components
adding to the ones already found, and does the deep structure-following clean
on an rms that now means what it says. The relaxed retry and the
primary-beam fallback are kept behind it, and a run can switch the bootstrap
off with `clean_bootstrap_enabled: false`.

### 70. The pipeline published a dirty map and called it a final image

The defect above is a deconvolution failure. That it reached the observer
three times is a separate defect, and the worse of the two: with no CLEAN
model the stage printed a warning and published anyway.

An undeconvolved image is not a faint version of the right answer. It is the
sky convolved with the array's point-spread function, so each source carries
sidelobes at the level of its neighbours and no flux measured from it means
anything. Final imaging now refuses to publish one when the residual holds a
peak above `undeconvolved_max_peak_snr` (5 sigma): the emission is real and
every clean path failed to reach it. A field whose residual holds nothing
above that threshold is consistent with empty, so the image still stands, but
it is labelled as the dirty map it is. An image kept from an earlier run is
held to the same standard on resume rather than passing because it exists.

The report now also carries `model_flux_Jy`, the flux in the CLEAN model,
which answers "did the deconvolution do anything" in one number instead of
requiring the image to be looked at.

`selfcal.py` now declares a `PIPELINE_INTERFACE_VERSION` as `imaging.py` does,
and both are at 3. Files reach the cluster one at a time, and two earlier
rounds were spent on a stale copy (finding 54); a half-applied release of this
change would produce exactly the dirty map it fixes, so preflight names the
file instead.

### 71. The preview could only ever show noise

The next image after finding 69 still looked like static, and part of the
reason was the plot rather than the data. `plot_casa_image` stretched every
image linearly between its 1st and 99.5th percentiles. Well over 99% of a
radio image's pixels are noise, so that range is about -2.3 to +2.6 sigma:
the noise fills the full colour range, every source saturates to white, and a
perfectly thermal-noise-limited image looks the same as a broken one.

The preview now measures the noise robustly (median absolute deviation) in the
middle quarter of each axis, where a primary-beam correction has barely raised
it, and displays from -3 sigma to the peak, linear within 3 sigma and
logarithmic above. The measured noise, peak and their ratio are printed in the
title, and the colour bar is labelled in Jy/beam. `qa_products.py` now
declares an interface version, so a stale copy cannot silently keep the old
plots.

The colour bars of the two images the observer sent still carried numbers. A
1-99.5 percentile stretch spans about 4.9 sigma, so the noise of the
PB-corrected image went from about 0.5 to about 0.28 mJy/beam between them:
the bootstrap deconvolution of finding 69 did help. What remained is set out
in the next two findings.

### 72. The thermal noise the pipeline reported was the sample's, not the data's

Finding 64 reported "the thermal noise this data can reach" as
sigma/sqrt(N), with N the number of samples in the coherence check. That check
reads a bounded sample of about 20000 rows, so the figure, 31 uJy/beam for the
observed target, was the noise of the sample. The full dataset holds
proportionally more samples and reaches lower. The significance test is
unaffected, since it is a statement about the sample, but the floor was
overstated, and it was used to judge how far the image was from its limit.

The check now scales the unflagged sample count by the fraction of rows read
and reports the naturally weighted thermal noise of all the data, labelled as
such. It also records what a gain solve will see: antennas with unflagged
data, spectral windows, scans, correlations, the integration time, and the
on-source time from the number of distinct timestamps, which it counts by
reading the TIME column a block at a time.

### 73. Self-calibration was gated on the number that gain errors hold down

Phase self-calibration started only when the initial image reached a peak/RMS
dynamic range of 20. For a field whose noise is set by antenna gain errors,
the spokes around every bright source that the observer's image shows, that
is circular: the errors that self-calibration corrects are what keep the
dynamic range below 20, so it never runs.

What decides whether phase solutions can be trusted is their own S/N. One
antenna's phase is solved from every baseline to it, over one solution
interval, so its S/N is

    A_model * sqrt(n) / sigma

where sigma is the per-component visibility noise from finding 72, n the
number of unflagged samples on that antenna in one interval (one spectral
window, both hands for gaintype T, one for G, the interval capped at a scan
since solutions do not span scans), and A_model the model's rms visibility
amplitude. By Parseval that last is the square root of the sum of the squared
model pixels, which `imstat` returns as `sumsq`; components spread over
several pixels make it an underestimate, so the estimate errs low. For a field
like the observer's, a few mJy of model against 30 antennas and a few hours on
source gives an S/N near 100 at a 4-minute interval.

Phase self-calibration now starts when either the dynamic range reaches
`phase_start_dynamic_range` (20) or that S/N, for the longest configured
interval, reaches `phase_start_min_solution_snr` (10, twice gaincal's minsnr
for these bands). Neither condition weakens a safeguard: gaincal still
discards solutions below minsnr, the unusable-solution gates still reject a
poor table, and each loop is still accepted only if the image improves. The
skip message and the self-calibration report give both numbers.

### 74. The final image now says whether it is noise-limited

Final imaging compares its residual RMS with the thermal floor measured by
self-calibration. Within about a factor of two, Briggs weighting and
w-projection account for the difference and the image is noise-limited. Above
`noise_over_thermal_warn` (3) it is limited by errors in the data, which more
cleaning cannot remove, and the stage says so, naming whether phase
self-calibration ran and, if not, the solution S/N that kept it from running.
Both figures are in the imaging report as `thermal_noise_Jy_beam` and
`noise_over_thermal`.

### 75. Every peak the pipeline measured was a clipped statistic

Every peak, dynamic range and residual significance in the self-calibration
and imaging stages came from `imstat(algorithm='chauvenet', maxiter=5)`,
taking `max` alongside `rms`. Clipping algorithms exclude outliers before
computing *every* statistic they return, the extremes and their positions
included. Checked directly against `casatools.image.statistics` in CASA: a
256x256 image of 1e-4 Jy noise with one 6 mJy pixel gives, classic, max 6.0
mJy at (200, 50); chauvenet, max 0.41 mJy at (149, 101). A 60-sigma source
is reported as a 4-sigma maximum somewhere else.

So every "peak" was about 4 sigma whatever the image held, and everything
built on it was inert:

- the self-calibration start gate (`phase_start_dynamic_range` 20) compared
  a dynamic range of about 4 with 20 and could never pass; the observer's
  field, at a true dynamic range of 126, was never self-calibrated by it;
- the per-loop acceptance test compared two clipped maxima, so loops were
  accepted or rejected on the noise in a 4-sigma statistic;
- the refusal to publish a dirty map (finding 70) and the relaxed-retry and
  primary-beam-fallback cleans (findings 59, 61) all trigger on a residual
  peak of 5 sigma or more, which a clipped maximum never reaches. That is why
  a dirty map was still published after finding 70 was in place.

Findings 60, 69, 70 and 73 each explained one of these symptoms without
seeing the common cause; the explanations were incomplete, and 73 in
particular attributed to gain errors a dynamic range that the statistic
itself was holding near 4.

`casa_compat.peak_and_noise` now takes the extremes and their positions from
classic statistics, which see every pixel, and the noise from Chauvenet
clipping, over the same pixels. Self-calibration's initial and per-loop
measurements, both stages' residual diagnostics, the imaging metrics and the
DD-calibration metrics script all use it.

### 76. A significant peak must be one that noise does not reach

Once peaks are true maxima, a fixed 5-sigma trigger fires on noise. An image
of N independent beams has its largest |value| above t with probability
about 2N Q(t); a 2048-pixel image at about 10 pixels per beam holds roughly
400,000 beams, and noise alone exceeds 5 sigma in most such images. The
triggers had been set against a statistic that never exceeded about 4.

`casa_compat.noise_peak_snr` gives the level noise alone exceeds with 1%
probability for the image's size: 5.1 sigma at 512 pixels, 5.6 at 2048, 6.1
at 8192. The relaxed retry, the primary-beam fallback, the refusal to
publish an undeconvolved image and the continuation of CLEAN (finding 77)
each use the larger of their configured level and that floor. A Monte Carlo
test confirms the floor: over 200 noise images it is exceeded about 1% of
the time, where 5 sigma is exceeded in more than 10%.

### 77. CLEAN's own account of why it stopped was discarded

The observer's final residual still held a 6.07 mJy/beam feature, 45 times
its noise, with no record of why CLEAN had stopped: every `tclean` return
value was thrown away. That value carries the stop code, the iterations and
major cycles run, and the PSF's largest sidelobe. Its codes were checked
against CASA's `imager_base.getStopDescription`: 1 iteration limit, 2
threshold, 3 force stop, 4 no change over two major cycles, 5 and 6 a
diverging residual, 7 empty mask, 8 nsigma, 9 nmajor.

Each pass is now summarised, printed and kept under `clean_passes` in both
reports. When the auto-multithresh pass stops at its iteration or
major-cycle limit with a residual peak above the significance floor, it is
continued from where it stopped (`restart=True`, `calcpsf=False`,
`calcres=False`, keeping its mask and model) up to `clean_max_continuations`
(3) further runs of `niter`, and a warning names what to raise if that is
still not enough. In self-calibration this matters twice: the flux an
unfinished model lacks is absorbed into the gain solutions.

### 78. The bootstrap stopped at the primary beam's 0.2 level

Finding 69's bootstrap cleaned only where CASA's primary beam exceeds 0.2,
while final imaging deconvolves to `pblimit` 0.05. A bright source between
the two, in the outer primary beam, is exactly the source whose sidelobes
most inflate the noise that auto-multithresh then measures, and the
observer's residual shows its two strongest features outside the
PB-corrected circle. The bootstrap now reaches as far as the deconvolution
(`clean_bootstrap_pbmask` defaults to the run's `pblimit`); an explicit
setting still wins.

### 79. The largest residual is located and its cause named

A strong residual has three causes that need opposite fixes: CLEAN was
never allowed there (outside its mask), CLEAN stopped first, or CLEAN reached
it and could not model it. Final imaging now reports the largest residual's
pixel, distance from the pointing centre, CASA primary-beam level and mask
membership, and names which of the three it is from the mask and the recorded
stop code. When CASA did not report a stop code it says the cause is unknown
rather than presuming convergence.

In the third case, for a source off-axis, it also computes how much the
primary beam alone bends the apparent spectrum. A beam whose width scales as
1/frequency multiplies a spectrum by pb**((nu/nu0)**2), which in log-log is
a curvature of 2 ln(pb); a linear (nterms=2) fit misses it by up to
(2/3)|2 ln pb| h**2 of the apparent flux over a band of half-width h in
ln(nu). At pb 0.2 across 550-750 MHz that is about 5%, and at 550-850 MHz
about 10%: a residual of several hundred microjansky around a 10 mJy
apparent source, whatever the calibration. nterms=3 models it.

### 80. Pixel sampling is checked against the restoring beam

Too few pixels across the beam leave every point source between pixel
centres, which CLEAN's delta-function components cannot represent; the
error sits around each bright source at a fixed fraction of its flux and no
cleaning removes it. Final imaging reports pixels across the restoring
beam's minor axis and warns below `min_pixels_per_beam` (3), giving the cell
that would sample it at four pixels per beam.

`casa_compat.py` and `ddcal.py` now declare interface versions, and
`imaging.py` and `selfcal.py` are at 5, so preflight names any file of this
release that did not reach the cluster.

### 81. Every channel of a calibrated uGMRT target had the same weight

GMRT data imported with importgmrt or importuvfits carry one weight per row
and correlation, shared by every channel; nothing in the pipeline created
WEIGHT_SPECTRUM, ran statwt or initweights. A channel's calibrated noise in
Jy follows the SEFD there and whatever interference remains in it, and both
vary across Band 4, so imaging gave noisier channels the same say as better
ones. For per-channel noise sigma_c, equal weighting leaves the image noise
sqrt(mean(sigma_c^2) * mean(1/sigma_c^2)) times the inverse-variance
weighted value.

This finding was first written with a wrong premise, and the correction
matters for what to expect from it. It said the calibrated noise scales as
1/|B_i B_j|, so that the observer's bandpass, falling about twentyfold from
600 to 850 MHz, implied a noise ratio of several hundred across the band and
an equal-weighting penalty near 6. That holds only if the raw correlator
noise is the same in every channel. For unnormalised correlator output, the
electronics' gain multiplies signal and noise alike and divides out in
calibration; the calibrated noise is set by the SEFD, which the bandpass
amplitude need not track. The size of the benefit is therefore not predicted
here but measured: the run reports the spread of the statwt weights and the
penalty they remove.

Each target is now reweighted right after it is split: statwt on DATA (where
the split puts the calibrated visibilities) in bins of 16 channels by 120 s,
per baseline and correlation, writing 1/variance in Jy^-2 to WEIGHT_SPECTRUM.
Per CASA's own documentation of statwt, a channel bin other than "spw" is
what creates and fills that column. Bins with fewer than 16 unflagged values
are flagged rather than weighted from a handful of numbers, which would let
one bin dominate the image. A marker is written only when the run completes;
a target MS without one is preserved as incomplete and split again, never
reused half-reweighted. Targets split by an earlier release are left alone,
with a warning, so their downstream products stay consistent.

The run reports the spread of the resulting per-channel weights and the
noise penalty equal weighting would have carried. It is on by default for
uGMRT (`ugmrt_statwt_enabled`), off for MeerKAT and LOFAR, whose data carry
channel weights of their own.

statwt counts all scatter as noise, including the sky's. That is negligible
when the field's visibility amplitude is small against the per-visibility
noise (for the observer's field, about 7 mJy against 0.27 Jy). Once
self-calibration has a model it compares the two and warns above
`statwt_max_sky_to_noise` (0.1), when long baselines would have been
under-weighted.

### 82. The thermal floor now comes from the weights

With inverse-variance weights the naturally weighted image noise is
1/sqrt(sum of weights), directly. The coherence check now sums the weights of
the samples it reads, scales the sum to the whole dataset as it does the
sample count, and uses 1/sqrt of it as the thermal floor whenever calibration
reports statistical weights for the target. The equal-weight estimate of
finding 72 is kept alongside, and without statistical weights it remains the
floor, since a row weight of 1 has no unit.

### 83. The calibration plots drew flagged solutions as real ones

`plot_caltable_time_series` and `plot_caltable_bandpass` read each table's
FLAG column and never used it. CASA stores a flagged gain as 1+0j, so the
amplitude plots carried a false line at exactly 1.0; after fluxscale scaled
those placeholders by 1/sqrt(S) it moved to 0.48 in every phase-calibrator
scan, and the phase plots carried a false line at 0 degrees on top of the
reference antenna's. Flagged bandpass channels appeared as flat plateaus at a
fill value. Only the first polarization was drawn, every calibrator field was
pooled, and the bandpass had no frequency axis or spread.

Flagged solutions are now omitted and flagged channels are gaps. Both
polarizations are drawn, each calibrator field has its own colour, the
bandpass shows the median with the 16th-84th percentile spread across
antennas (circular spread for phase) against frequency, and the phase axis
says it is relative to the reference antenna.

### 84. What the gain tables say is now stated in numbers

A scatter of every antenna's phase cannot show whether phase transfer is
good; the change in each antenna's phase between consecutive visits to the
phase calibrator can, because the target's phases are interpolated across
exactly those gaps. The calibration stage now reports, per calibrator field,
the median gap and the median and 90th-percentile absolute phase change
(wrapped to +/-180 degrees, reference antenna and flagged solutions
excluded), with the five antennas that move most. It also lists antenna and
polarization chains whose median flux-scaled gain lies outside half to twice
the array median, the weak chains visible as low outliers in every scan of
the observer's plots. Both are in the calibration report under
`calibration_stability`.

`calibration.py` and `config.py` now declare interface versions: without the
new `config.py` the uGMRT default that turns statistical weights on would be
absent, and the weights would silently not be written.

### 85. A reused working copy had its flags stacked run after run

The observer's flagging plots read 50.5% flagged. The progression shows the
split: 42.7% was already flagged when the flagging stage began ("pre" is
measured before any of this pipeline's own flagging, the deterministic cuts
included), and the stage itself added 7.8 points, or 13.6% of the data it
received unflagged. That is an ordinary interference loss for Band 4.

The 42.7% is the anomaly, and the plots show where it came from. The
observation's 400-MHz GWB setup digitizes 550-950 MHz, and this pipeline
cuts everything above Band 4's specified 850 MHz: about a quarter of the
channels. Applied to unflagged data in this run, that cut alone would have
raised the total by about 14 points during the first iteration. The whole
first iteration, band cut included, raised it by 6.7. So the band cut had
already been applied: the Measurement Set arrived carrying an earlier run's
flags.

That happens whenever a run in a new --outdir reuses a working copy: a fixed
`working_ms` path (the README's own example sets one), a FITS run's reused
`--imported-ms`, or `copy_input_ms: false`. The managed copy inside the
output directory is fresh for each new directory; these are not. Each such
run counts the earlier flags as the input's own and flags interference
again on top, so the data erode run after run.

The flagging stage now reads `<ms>.flagversions/FLAG_VERSION_LIST` before it
does anything. Every run of this pipeline saves flag versions whose comments
identify it, and at that point this run has saved none, so a version of ours
can only come from an earlier run. It stops and says how to start from a
fresh copy, without ever suggesting that an MS that may be the observer's
original data be deleted. `flagging_allow_previously_flagged_input: true`
permits it deliberately. Resuming in the same output directory is unaffected:
the check runs only before this run's pre-flag summary exists.

### 86. The flagging report separates what arrived from what it flagged

A cumulative flagged fraction cannot answer "was this over-flagged",
because it mixes the input's flags with the run's. The flagging stage now
reports the fraction that arrived flagged and the share of the data it
received unflagged that it went on to flag, (post - pre) / (1 - pre), in the
log, in `flagging/flagging_attribution.json`, and on the progression plot.

### 87. A dead antenna stopped calibration once its flags were no longer inherited

With a fresh working copy, calibration failed at the pre-bandpass phase
solve: `prebp.cal: ANTENNA1=3: group has no usable solution`. Antenna 3 is
C03, which the earlier (stacked) runs showed 100% flagged. A dead or
disconnected receiver produces noise, and interference flagging leaves noise
alone, so on a clean copy C03 reached calibration only partly flagged. The
delay solve, which averages a whole calibrator scan, found a few solutions for
it; the pre-bandpass solve, at the integration time, found none. In the
earlier runs C03 had been fully flagged by a previous run's flags, so the
flagging stage listed it as bad and calibration excluded it: finding 85's
flag stacking had been hiding this.

The calibration-table QA exempts only antennas absent from the data or
excluded from the solves, and the only automatic exclusion was the delay
stage's quarantine. That quarantine is now general: when a solve fails
because, and only because, some antennas have no usable solution at all
(ANTENNA1 the sole failed grouping, as before), the pre-bandpass stage
quarantines them under the same bounds (one per event by default, the
array's total exclusions within 10% of its antennas, at least 8 left, never
the reference antenna), then re-solves the delay table without them, since
"excluded from every solve" includes the one already made, and re-solves the
pre-bandpass table. The quarantine record now says which table and which
failure caused each exclusion. When the bounds prevent it, the error names
the antennas and the settings that would allow it, instead of repeating the
table's failure.

### 88. Raising the quarantine limit would have orphaned the calibration

The next run found two dead antennas in the pre-bandpass solve, C03 and S03,
and stopped at the per-solve limit of one, naming both and the setting to
raise. Raising it in `config.yaml` would then have been refused on resume:
`config.yaml` is part of every stage's implementation fingerprint, and the
whole calibration task configuration is part of its science recipe, so the
failed calibration would have looked like a different recipe and required a
new `--outdir`, and with it a fresh copy and a full re-flag.

The four quarantine settings decide only what happens when a solve cannot use
an antenna. They cannot affect a finished checkpoint: the solves done before a
quarantine either did not involve it, or are redone without the quarantined
antenna when it happens, and what the quarantine decided is recorded on its
own. So they are now failure-handling policy, not recipe. When a recorded
calibration recipe, of a failed attempt or of finished products, does not
match, the pipeline reproduces it with the quarantine settings the attempt
ran under (recorded from now on) or, for attempts recorded earlier, the
values the pipeline has always shipped with, optionally together with the
existing post-flag provenance refresh. If that reproduces the recorded
SHA-256, the recorded value is pinned, as for a provenance refresh, so the
flux-model checkpoint bound to it stays valid. Any other change, a science
setting included, still fails closed. The state records, for every recipe
fingerprint, the policy it was computed under.

### 89. The WSClean stage could not start at all

The orchestrator builds the WSClean stage's command with `--upstream-outdir`,
which finding 56's multi-target work added so that a per-target stage can read
the run's shared reports. `wsclean_imaging.py` had no such option, so its
argument parser refused the command (`unrecognized arguments:
--upstream-outdir`) before doing anything. That affected every WSClean run,
LOFAR included, from the multi-target change onwards. The existing tests
checked the command's contents and never fed it to the script's own parser.

`wsclean_imaging.py` now accepts `--upstream-outdir` and reads the inspection
and calibration reports and the split targets from it, and its own
self-calibration report from `--outdir`, as `imaging.py` does. The new tests
build the real command with the real orchestrator, for LOFAR, MeerKAT and a
per-target scope, and parse it with the stage's real parser, which is the
check that was missing.

### 90. MeerKAT now images with WSClean; uGMRT is pinned to tclean

`telescope_profiles.py` makes WSClean the default imager for MeerKAT, as it
already was for LOFAR. uGMRT stays on CASA `tclean`, and this is enforced in
four places so no configuration can send uGMRT to WSClean: the profile default;
the capability check, which refuses `--imager wsclean` for uGMRT at preflight;
the WSClean stage itself, which refuses a uGMRT MS; and `run_ugmrt.py`, which
pins `--imager casa` so that an `imager: wsclean` meant for MeerKAT in a
`config.yaml` shared between launchers cannot reach it. MeerKAT can still use
`--imager casa`; LOFAR still requires WSClean.

The capability check had forbidden WSClean for MeerKAT, with the reason that
the generic WSClean path had no telescope-specific primary-beam correction.
That reason was right, and the default was not changed by removing the check:
the next finding supplies the correction.

### 91. A MeerKAT primary-beam correction for WSClean's output

WSClean's own correction (EveryBeam) is validated here only for LOFAR, and the
stage refused every non-LOFAR beam except "none". `meerkat_beam.py` applies the
same measured `katbeam` model the CASA backend uses, to WSClean's FITS output.
The orientation is the CASA path's exactly (offsets from the reference pixel
with the image's own increments, `CDELT1`'s sign included), so a pixel gets the
same beam value from either backend, and it keeps the CASA path's safeguards:
a beam that is not finite or whose peak is not unity is refused, pixels below
`pbmin` are blanked, and an unversioned `katbeam` is refused.

It differs in how the band is treated, deliberately. A MeerKAT band spans a
factor of two and a primary beam's width goes as 1/frequency, so one
frequency's beam is wrong across most of an MFS image. The sub-band images are
each divided by the beam at their own effective frequency, using only the part
of each sub-band inside the receiver's science range, since the flagging stage
removed the rest, and a sub-band wholly outside it is skipped. The MFS image
is divided by the band-averaged beam: the MFS image of a flat-spectrum source
is the weighted mean of its sub-band images, so it is attenuated by the
weighted mean of their beams, and dividing by that recovers the flux. This is
tested on a uniform flat-spectrum sky seen through a frequency-dependent beam
(the corrected image equals the true flux to 2e-4), and a single-frequency
correction is shown to leave exactly the pixelwise error the report predicts.
Spectral index is not corrected and the report records that. The report also
states how large the single-frequency error would have been at the half-power
radius and at `pbmin`, so the choice is visible rather than assumed.
`meerkat_mfs_beam_mode: reference_frequency` restores the CASA backend's
single-frequency division.

WSClean is run with joined-channel deconvolution for MeerKAT (UHF 6, L 8 output
channels, `-fit-spectral-pol` from `nterms`): without `-channels-out` WSClean
makes one flat-spectrum image and ignores `nterms`. The same channels give the
sub-band images.

### 92. The WSClean stage accepted only LOFAR-shaped input

The stage required an inspection report bound to the exact imaging MS, a fixed
data column checked against that report, and a field id taken from it. None of
that holds for MeerKAT, whose imaging MS is a self-calibration product that
inspection never saw. For a LOFAR MS the checks and their messages are
unchanged. For MeerKAT the MS must be one the pipeline recorded, the
self-calibrated MS in `selfcal_report.json` or an entry of `split_targets.json`,
and then its target is known, the single field is not selected by a raw-MS
field id (which would select the wrong field, or none), and the data column is
chosen as the CASA backend chooses it, so both read the same data: DATA for
`final_selfcal.ms`, else CORRECTED_DATA when the MS declares one. Declared
columns are read from the plain strings in `table.dat`, which needs no CASA.
An MS that is neither inspected nor a pipeline product is refused.

The product classifier looked for "model", "beam" and the like anywhere in a
filename, so a target called `model_galaxy` could hide its own image. It now
reads only what follows WSClean's prefix.

### 93. The orchestrator validates a MeerKAT WSClean report from its provenance

`expected_backend` was "wsclean if LOFAR else casa". It is now what this run
uses, and an old CASA-imaged MeerKAT directory is told to pass `--imager casa`.
The CASA contract for MeerKAT reports is unchanged. A WSClean report is held to
the FITS counterpart of it: the model and package version, the reference
frequency, that only the flux scale was corrected, and that the restored image,
the beam and the corrected image are still the files the correction was made
from (size and SHA-256 recorded in a provenance file written with them). The
tests produce a report with the real stage and validate it with the real
orchestrator, and each tampering (a changed pixel in the corrected image or the
beam, a provenance file naming other images, a changed model, version,
frequency, mode or method) is rejected.

### 94. A missing tool now stops the run at the start

Imaging is the last stage, so a missing `wsclean`, or a Python without
`katbeam`, would have been found after calibration and self-calibration had
used the job. With the telescope given explicitly the preflight now checks for
WSClean and, for MeerKAT, imports `katbeam` and `astropy` in the `--python`
that will run the stage, saying how to proceed (for MeerKAT, `--imager casa`;
for LOFAR, which has no CASA path, only to install it). The stage's own
preflight repeats the check before WSClean starts. Required files are
stage-dependent: `meerkat_beam.py` is demanded only of a run that can be
MeerKAT and use WSClean, so a uGMRT run never asks for it.
`telescope_profiles.py` now declares an interface version, since a stale copy
would quietly keep MeerKAT on CASA.

When a job has a memory limit, WSClean is passed `-abs-mem` at
`wsclean_memory_fraction` (0.8) of it. WSClean sizes itself from the node, not
the Slurm job, the same mistake CASA made in finding 65. The cap is not added
when the node is the limit, and it is excluded from the comparison that decides
whether a finished image is reused, so a different allocation does not discard
hours of imaging.

### 95. A new setting must not change a recipe it cannot affect

Adding keys to the imaging configuration would have changed the science-recipe
fingerprint of every existing imaging run, uGMRT's included, and made a resume
in an old output directory refuse. The WSClean memory cap is now never part of
the recipe, the MeerKAT beam mode only for a MeerKAT run on WSClean, and
`wsclean_channels_out` only for a WSClean run. Checked directly against the
code from before this change: the imaging recipe for uGMRT B4, LOFAR HBA,
MeerKAT L with CASA and MeerKAT UHF with CASA is byte-identical, so those
directories resume as before. MeerKAT on the new default differs, as it should,
being a different backend.

**What has not been verified.** No `wsclean` or `katbeam` was available where
this was written, so the WSClean command line was checked against the LOFAR
path's existing flags (`-join-channels`, `-channels-out`, `-fit-spectral-pol`,
`-multiscale`) and the new `-abs-mem`, which is documented but not exercised;
the beam correction was tested on real FITS files with a stand-in beam model of
the right frequency dependence; and neither has been run against a real
MeerKAT data set or compared with the CASA backend's image of the same data.
Self-calibration still uses `tclean` for MeerKAT, because its gain solves read
the `MODEL_DATA` that `tclean` writes.

## Post-audit: why a deconvolved image was still seven times thermal

Findings 96 to 101 come from one run read end to end: uGMRT Band 4, target
`G0212-0335`, the reports under `targets/G0212-0335/` and the CASA logs of its
calibration, self-calibration and imaging. The published image had a residual
rms of 87 uJy/beam against a thermal floor of 11 uJy/beam (`noise_over_thermal`
7.7), with the point-spread function visible around every source. The numbers
quoted are that run's.

### 96. One flagged antenna rejected every self-calibration solution

Phase self-calibration was attempted and its solve succeeded: `gaincal` at
`solint=4min` found solutions in all seven intervals that held data. The table
was then rejected, the "rescue" solve was rejected the same way, and the stage
carried on with `phase_tables: []`. The published image was made from data that
had never been self-calibrated. The reason recorded was `ANTENNA1=2 has
unusable solution fraction 0.857`.

Antenna 2 is C02, which arrived 97.8% flagged. Eight of the thirty antennas
were in that state: three quarantined during calibration (C03, S03, W04, 100%
flagged) and five that had lost their calibration for the target's main scan
(C02, C04, C12, S02, W06, 97-98%). `gaincal` can only return flagged solutions
for an antenna with no data, and `caltable_status` failed the whole table when
any one antenna exceeded 80% unusable. With a flagged antenna present, which
is the normal condition of an array, no target could ever pass. Calibration's
own table check (schema 2) already left out quarantined and placeholder
antennas; the self-calibration copy of it did not.

Self-calibration now counts unflagged data per antenna before the first solve
(`measure_antenna_data`, one `flagdata` summary). An antenna with less than
`min_antenna_data_fraction` (0.20) of the median antenna's unflagged samples is
recorded as unsolvable and its rows are left out of every coverage figure.
Among the antennas that do have data, up to
`max_caltable_failed_antenna_fraction` (0.25) may fail without rejecting the
table, and likewise `max_caltable_failed_time_fraction` of the solution
intervals; they are listed in the table's QA sidecar as
`tolerated_failed_groups`. A spectral window, field or observation without
solutions still rejects the table, as do non-finite or zero solutions and a
global unusable fraction above 0.5. The census is in the report as
`antenna_data`, and the reference-antenna list is ranked by unflagged data
instead of row count, which had offered C02, C03 and C04.

### 97. A failed self-calibration solution flagged the data it could not improve

Solutions were applied with `calflagstrict`. Each phase loop solves an
increment on data already corrected by the loops before it, so where an
increment cannot be solved the right outcome is the correction those data
already have. `calflagstrict` flags them instead, and at short solution
intervals on a faint field that can remove a large part of an observation
while the dynamic range, the acceptance test, goes up. Self-calibration now
applies with `calonly` (`selfcal_applymode`). Calibrator-based calibration
keeps `calflagstrict`, where a missing solution really does mean uncalibrated
data.

The rescue solve combined scans (`solint='inf'`, `combine='scan'`): one
solution per antenna across visits to the target that were hours apart, over
which antenna phases are unrelated. It now solves one solution per scan.

### 98. auto-multithresh could not open a mask on this point-spread function

After the bootstrap of finding 69, the masked pass stopped at once: CASA's
`SidelobeLevel = 0.440183`, "No regions found by automasking", stop code 7,
zero iterations of 200000. The mask criterion is `sidelobethreshold x
max_sidelobe x peak_residual`; with 2.5 and a sidelobe of 0.44 that is 1.1
times the peak, which no pixel can exceed. Two usable scans and 33 minutes on
source give that PSF. So every image of this run, self-calibration and final,
was cleaned only by the 5 sigma unmasked bootstrap, on an rms that sidelobes
had inflated, and each source kept up to 5 sigma of its PSF.

Three changes:

- The sidelobe threshold is lowered per image, from the largest sidelobe that
  the first pass on that image reports, so that the criterion sits at no more
  than `automask_max_sidelobe_fraction` (0.6) of the residual peak
  (`casa_compat.adaptive_sidelobethreshold`). A PSF with sidelobes near 0.1 is
  unaffected. What was used is recorded as `automask_adaptations`.
- The bootstrap is shallower, `clean_bootstrap_nsigma` 8 instead of 5. Its job
  is to make the rms mean something. At 5 sigma it left nothing above
  auto-multithresh's own noise threshold (4.25 sigma, less the area a region
  must have), so even a mask that could open had nothing to open on.
- If auto-multithresh still opens no mask, what remains above
  `clean_floor_nsigma` (5) is cleaned without one (`floor_clean`), so an image
  is never left shallower than before.

### 99. CASA's uGMRT beam is narrower than the real one, and nothing outside it was cleaned

`pblimit` was 0.05. CASA masks, and does not deconvolve, anything below that
level of its own beam model, which for uGMRT is the Airy disk of a uniformly
illuminated 45 m dish. The bootstrap's mask covered 66.45% of the 8192-pixel
image, a radius of 31.4 arcmin, which is where that model reaches 0.05 at 700
MHz. The measured beam is 0.16 there, and at 550 MHz it is 0.34. The dirty map
peaked at 30.9 mJy/beam over the full image against 20 mJy/beam inside the
published circle: the brightest source in the data lay outside it. Sources
between 31 arcmin and the image edge were left with their full PSF, which is
the streaks that cross the published region without converging on anything in
it.

`pblimit` is now -0.01 for self-calibration and for final imaging with a
measured beam, which switches CASA's blanking off. The unmasked passes use
`usemask='user'` with no mask in that case, so they reach the whole image. The
published aperture is unchanged: it is set by `pbmin` on this pipeline's own
beam.

### 100. The Band 4 image was four times too fine and still too small

8192 pixels of 0.5 arcsec put 9 pixels across a 4.8 x 4.5 arcsec beam and
reached 34 arcmin from the pointing centre, where the measured beam is 0.11
at 700 MHz and 0.28 at 550 MHz. Final imaging took 20 hours for 17 major
cycles and 15261 multi-scale iterations, in an 8-thread session.

Band 4 now images 5120 pixels of 1.0 arcsec for both self-calibration and
final imaging: 4.5 pixels across that beam, 43 arcmin from the centre, and 39%
of the pixels. The multi-scale list drops its 135-pixel scale, which at 1.0
arcsec would be 2.25 arcmin; cost grows as the square of the number of scales.
Other bands are unchanged here, but their defaults have the same shape of
problem and should be reviewed against the beam before use: Band 2 is 1.0
arcsec cells for a beam near 15 arcsec, over a field less than half the
primary beam's width.

### 101. `strict` reference-antenna mode removed a whole science scan

The target's first scan (scan 3, 28% of its data and a third of its hour-angle
coverage) arrived 100% flagged. It had been 35.6% flagged after the flagging
stage. The calibration log has the cause:

    At 2020/11/22/14:08:02.4 (Spw=0, Fld=1), refant (id=0) was flagged;
    flagging all antennas strictly.

The reference antenna was flagged during the first phase-calibrator scan.
`calibration_refantmode: strict` then flagged every antenna's amplitude
solution for that scan, and linear interpolation of a time-dependent table
flags data bracketed by a flagged solution, so the science scan next to it
was flagged in full by `applycal`. Scans 4 to 7 show the same decay, 75% to
61%, against about 42% later in the observation.

For parallel-hand, total-intensity calibration a change of reference antenna
adds the same phase to every antenna at that time and cancels on every
baseline, so `strict` buys nothing here. The default is now `flex`, CASA's
own; `strict` remains available for work that needs one physical reference
throughout. This changes the calibration recipe, so it takes effect only in a
new output directory.

The same interpolation rule is why the five antennas of finding 96 were lost
for the main scan: each had a failed solution at one of the two calibrator
scans bracketing it. That is not changed here.

### 102. Flags written by provisional solutions were never taken back

The re-run with `flex` still showed the first phase-calibrator scan 83%
flagged and the science scan after it 86%, so finding 101 was at most a small
part of that loss. The run's own flag summaries show where it came from:

| after | flagged |
|---|---|
| the flagging stage | 32.7% |
| the round-1 applycal | 49.2% |
| corrected-data tfcrop/rflag/extend | 51.1% |

The detectors, which are what the two-round design exists for, added 1.9%.
The first applycal added 16.5%. About 9 points of that are the two
quarantined antennas (C03 and S03, 100% flagged), which are excluded from
every solve: two of thirty antennas are on 13% of the baselines. That part is
correct and round 2 flags it again. The other 7 to 8 points, more than a
tenth of the data that was usable going in, are visibilities whose
provisional solution failed.

With corrected-data flagging enabled, calibration is solved twice. The first
set of solutions is provisional: it is solved on data that still holds the
interference the detectors are about to remove, so it fails often, and its
only purpose is to give those detectors a CORRECTED_DATA column. But applycal
runs with `calflagstrict`, which writes a flag for every visibility whose
solution failed, and `run_postcal_flagging` saved its flag version *after*
that applycal. Nothing between there and the second solve took those flags
back, so every failed provisional solution became a permanent loss. Linear
interpolation in time makes each one expensive: a time-dependent solution
flagged at one calibrator scan flags that antenna in the science scans on
both sides. The second solve then failed more often still, for want of the
baselines the first had removed.

`calibration.py` now saves the flags before the round-1 applycal
(`save_flags_before_round1_applycal`) and, once the detectors have run and
before the round-2 tables are solved, writes

    A | (C & ~B)

into the working MS, where A is the state before the first applycal, B the
state after it, and C the state after the detectors
(`release_provisional_solution_flags`). That is what the raw-data flagging
found plus what the corrected-data detectors found, and nothing else. CASA's
flagmanager can only AND, OR or replace, so the combination is done with the
table tool on the saved flag-version tables, a block of rows at a time, and
is idempotent. The round-2 applycal still flags whatever its own solutions
cannot calibrate.

It is built so that it cannot leave the run worse off than before: the
current flags are saved first and restored in full on any failure, after
which round 2 proceeds on the old flags with a warning. A run whose first
applycal predates this change has no saved state A and is left alone. The
outcome is written to `calibration/round1_flag_release_report.json`
(`fraction_released` is the figure to read).

Two things this does not do. Visibilities released here were hidden from the
corrected-data detectors, so they carry only the raw-data flagging; the
residual pass in self-calibration still sees them. And the round-2 applycal
keeps the interpolation rule, so an antenna whose *final* solution fails at
one calibrator scan is still lost in the neighbouring science scans. How
much of those 7 to 8 points comes back is therefore not known until it is
run; an antenna that is really unusable for part of the observation will be flagged
again, correctly.

No configuration key was added for this, deliberately: a new key changes the
calibration recipe fingerprint, and every existing output directory would
then refuse to resume (finding 95). `release_provisional_flags: false` in
the calibration task configuration switches it off if that is ever wanted.
Because the recipe is unchanged, a calibration already completed is adopted
as it stands; the change takes effect when calibration is next run from the
start.

**What has not been verified.** CASA does not start on the machine these
changes were written on, so none of them has been run on data. The logic is
covered by unit tests against stubbed CASA tasks (618 passing), and the
`tclean` arguments involved (`pblimit < 0`, `usemask='user'` with an empty
mask, a per-call `sidelobethreshold`) are documented CASA usage. The first
real run should be read the way this one was: `phase_history` should show
accepted loops with rising dynamic range, `clean_passes` should show the
auto-multithresh pass doing iterations, and `noise_over_thermal` should fall.

### 103. Every image gridded all 4096 correlator channels

The first run with findings 96 to 101 in place showed them working on
`G0212-0335`: three phase loops accepted, the dirty-image peak up from 30.4 to
53.8 mJy/beam as the phases came together, the restored dynamic range from 110
to 256, auto-multithresh opening a mask of 0.34% of the image, and the residual
from 103 to 75 uJy/beam against a thermal 20.6 uJy/beam. `calonly` held the
flagged fraction at 73.8% through every loop.

It was also slow. The initial image took 5.0 hours, phase loop 1 5.1, loop 2
4.2; each major cycle took 6 to 7 minutes and an image needed about 24 of
them, 11 in the 8-sigma bootstrap and 13 in the masked pass. At that rate
self-calibration alone would take more than a day for this 33-minute target
and nine targets about two weeks. The cause is in the data, not the
deconvolution: calibration split each target with all 4096 channels of 97.7
kHz (`applycal` reported 662,814,720 visibilities for 80,910 rows, 2
correlations x 4096 channels), and nothing downstream averaged them. Every
major cycle of every image read and gridded all of them.

Calibration now averages a target to the band's
`target_max_channel_width_khz` when it is split. Band 4 sets 400 kHz, which
averages four channels into 390.6 kHz: a quarter of the data in each major
cycle, and a quarter of the disk for every self-calibration candidate MS. The
cost is bandwidth smearing. For a Gaussian beam the loss of peak brightness
depends only on `channel width x radius / (beam x reference frequency)`,
because smear and beam both scale as 1/frequency; with the 4.9 x 4.5 arcsec
beam at 750 MHz it is 1.8% at the edge of the 5120 x 1 arcsec field (43
arcmin), 0.9% at 30 arcmin and 3.6% in the image corners. Integrated flux
densities are unchanged. The data are Hanning-smoothed before calibration, so
neighbouring native channels were already correlated.

The order of operations keeps statwt as it was: the target is split at full
resolution, statwt weights the native channels in its usual 16-channel bins,
and a second `split` with `width` averages them with those weights, summing
them into the averaged channel's WEIGHT_SPECTRUM. An averaged channel is
flagged only if all four of its channels were. The statwt marker is written
for the averaged MS last, so a run stopped anywhere before it splits the
target again rather than reusing it, and the full-resolution copy is removed
once the averaged one exists. The factor is the largest that keeps channels
within the limit and divides every window's channel count; a band without
the setting, or channel widths that cannot be read, keeps native channels
(with a warning in the second case). What was done is written to
`calibration/target_channel_average.json` and the calibration report, only
when a band was averaged, so nothing changes for other bands.
`calibration_target_max_channel_width_khz` in config.yaml sets one value for
every band, 0 to switch it off.

Only Band 4 is changed. Its fields and beams were reviewed in finding 100;
the other bands should be checked the same way before they are given a
width.

This adds a key to the Band 4 calibration configuration, so it changes the
Band 4 calibration recipe. It is meant for a new output directory, the one
that the calibration fix of finding 102 also needs. Uploaded into the
directory of a run already past calibration, it would make that run's
completed calibration refuse to resume.

**What has not been verified.** It has not run on data. The tests (636
passing) cover the averaging factor, the order of the two splits and statwt,
interruption and reuse, and the configuration. The smearing figures are the
standard Gaussian-beam formula, not a measurement. The speed-up is expected to
be about four times per major cycle, because gridding and reading scale with
the number of channels, but CASA's fixed costs per cycle (the PSF patch, the
minor cycle, the w-kernels) are unchanged, so the real gain should be measured
in the new run's CASA log.

### 104. The masked self-calibration clean had no bound but its depth

Every auto-multithresh pass of that run oscillated around its threshold
instead of settling onto it. The peak residual inside the mask at the start of
each major cycle was:

| image | major cycles | peak residual in mask (uJy/beam) |
|---|---|---|
| initial | 21 | 1289, 911, 1098, 700, 775, 795, 508, 768, 736, 456, 751, ... 279, 272 |
| phase 1 | 19 | 896, 636, 717, 489, 610, 406, 497, 353, 445, 300, ... 230, 222 |
| phase 2 | 13 | 831, 592, 485, 367, 639, 383, 353, 597, 358, 240, 219, 223, 211, 207 |
| phase 3 | 16, not finished | 816, 558, 643, 453, 523, 370, 426, 337, 655, 410, 265, 320, 520, 313, 531, 321 |

A minor cycle that cleans to within about 1.3 times the nsigma threshold
(about 210 uJy/beam) is followed by a major cycle that finds a peak of 520 to
655 uJy/beam inside the mask, which five to nine iterations then remove. CASA
reported "Possible divergence. Peak residual increased by 10% from minimum" in
the minor cycles of three of the four passes (1, 5 and 4 times). At 6 to 16
minutes a major cycle, half of each pass's 2.4 to 3.6 hours went to undoing
the cycle before. Where the returning peak sits cannot be read from the log;
the minor cycle uses an 80-pixel PSF patch against a PSF whose largest
sidelobe is 0.43, and every component it adds below about 5 sigma carries a
noise-dominated spectral term (the Taylor-1 sensitivity is eight times the
Taylor-0), so both are candidates.

Such a pass ends only when a cycle happens to land below the threshold, or at
`niter` (100000). A pass that reaches `niter` was then continued up to
`clean_max_continuations` (3) more times by `continue_to_convergence`, whose
test is the whole image's residual peak: in this field a 10-sigma negative
artefact beside the brightest source, outside the mask, which cleaning does
not remove. That would have been tens of hours on one image.

Self-calibration's masked pass now carries `nmajor` =
`clean_max_major_cycles` (20; 0 for no limit; `selfcal_clean_max_major_cycles`
in config.yaml), and a pass stopped by that limit (stop code 9) is not
continued. Of the passes above that converged (13, 19 and 21 cycles), twenty
would have changed only the initial image's, stopping it one cycle early at
279 instead of 272 uJy/beam. The fourth was still cycling at 22 (212, 563,
332, 208, 480, 297, 212 uJy/beam against a threshold of 201) after 4.5 hours,
with about a hundred iterations per three cycles: at that rate `niter` was
hundreds of hours away. The bootstrap is not bounded. A pass stopped by
`niter` is continued as before, now with the same bound. Final imaging
(`imaging.py`) is unchanged until its own masked pass has been seen on this
data.

The default lives in `selfcal.py` (`DEFAULT_CLEAN_MAX_MAJOR_CYCLES`), not in
`SELFCAL_TASK_CONFIG`, deliberately: the self-calibration recipe, which decides
whether an interrupted self-calibration may keep its accepted loops, is
unchanged, so `selfcal.py` alone can be put into a run already under way.
Setting the key in config.yaml does change the recipe. Before this was
offered for the running directory, the cluster's `selfcal.py` was checked:
its SHA-256 equals this file with finding 104 taken out, so the update
carries nothing else.

A smaller limit would save more time. With averaged data (finding 103) each
cycle should cost about a quarter as much, so the value is better set from
the new run's passes than from these.

**Seen on data since.** The next target of that run, `G0213-0251`, was
self-calibrated with this change in place. Its masked passes were called with
`nmajor=20` and its bootstraps with `nmajor=-1`; the initial image and phase
loops 1 and 2 each stopped with CASA's "reached nmajor" after 20 major cycles
and about three hours, with the peak in the mask at 534, 286 and 596
uJy/beam against thresholds of 498, 279 and 241, and none was continued.
Loops 1 and 2 were accepted (restored peak 62 to 112 to 126 mJy/beam,
residual rms 200 to 106 to 92 uJy/beam). So on this field the limit, not the
threshold, ends every masked pass. Tests: 645 passing when this was written.

### 105. The self-calibration report never recorded its final image

`G0212-0335` finished self-calibration with three accepted phase loops, and
its report read `"final_metrics": {}`. The peak, noise and dynamic range of
the final image were measured at the end of `generate_qa_plots`, after the
plotms figures. When those figures moved to the CASA-native path the function
began with `return self.generate_native_qa_plots()`, and nothing after it has
run since: no `final_metrics`, no `selfcal_metrics.json`. The measurement is
now its own step (`measure_final_image`), run before the figures, and the
unreachable plotms code is gone. A failure to read the image is a warning, not
a failed stage. This is reporting only and changes no recipe.

What that run's self-calibration did, from its report and CASA log:

| step | solint | peak (mJy/beam) | rms (uJy/beam) | dynamic range | residual rms (uJy/beam) | masked-pass major cycles |
|---|---|---|---|---|---|---|
| initial | - | 12.75 | 115.6 | 110 | 102.6 | 21 |
| phase 1, accepted | 4min | 21.68 | 93.7 | 231 | 80.9 | 19 |
| phase 2, accepted | 2min | 22.57 | 88.0 | 256 | 75.2 | 13 |
| phase 3, accepted | 1min | 23.06 | 86.4 | 267 | 73.1 | 23 |
| phase 4, rejected | int | 23.59 | 86.7 | 272 | 73.4 | 18 |
| final (phase 3 + residual flagging) | - | not recorded | not recorded | not recorded | 72.7 | 14 |

Loop 4 was rejected for improving the dynamic range by 1.9% against a
threshold of 2%. Residual flagging took the flagged fraction from 73.80% to
74.34%. The thermal floor from the statistical weights is 11.2 uJy/beam for
natural weighting and CASA's figure for the robust 0.5 weights actually imaged
is 20.6 uJy/beam, so the final residual is 6.5 and 3.5 times thermal on those
two measures. The whole stage took 32.7 hours.

### 106. The final image's masked clean stalled in the same way

Finding 104 left `imaging.py` alone until its masked pass had been seen. On
`G0212-0335` it ran multi-scale (scales 0, 5, 15, 45 pixels): the 8-sigma
bootstrap took 14 major cycles and 3.0 hours, and the masked pass took the
peak residual inside its mask through 581, 500, 424, 360, 318, 276, 276, 193
uJy/beam in its first eight cycles. From there, with the nsigma threshold at
136, the peak at the start of each cycle was 248, 239, 234, 234, 229, 151,
239: each minor cycle cleaned to its cycle threshold (0.65 of the peak, about
150) and the next major cycle put the peak back. Seven cycles and two hours
bought nothing, at 14 to 35 minutes a cycle, and nothing but a lucky cycle or
`niter` (200000, of which 22000 were used) would end it. Its
`continue_to_convergence` had the same flaw as self-calibration's: it would
have continued a pass stopped by `nmajor`, on the strength of the whole
image's residual peak.

`imaging.py` now bounds its masked pass with the same
`clean_max_major_cycles` (default 20 in `DEFAULT_CLEAN_MAX_MAJOR_CYCLES`, 0
for no limit) and does not continue a pass stopped by it; a pass stopped by
`niter` is continued as before. A pass stopped by the bound is reported as
not converged, with CASA's reason, and where its largest residual lies inside
the mask `describe_residual_peak` says CLEAN stopped before reaching it. In
config.yaml, `clean_max_major_cycles` sets both stages and
`selfcal_clean_max_major_cycles` / `imaging_clean_max_major_cycles` one. As
in finding 104 the default is in the code, so no recipe changes.

Two things from that pass that are not changed here. The multi-scale model
holds 0.772 Jy against 0.57 Jy for self-calibration's point-source model of
the same data, and the residual's robust rms is 46 uJy/beam against 69, so a
third of the field's flux is extended and self-calibration solves against a
model without it. And the 45-pixel scale has a net flux of -0.038 Jy, which is
a negative bowl being modelled. Both are worth a decision once this image and
the averaged-channel run exist to compare against.

**What has not been verified.** As finding 104: a stop by `nmajor` has not
been seen on this data. Tests: 655 passing.

### 107. Residual flagging subtracted a rejected loop's model

`latest_image_base` is documented as "the image whose model describes the
current Measurement Set", and returned the image of the highest-numbered loop
in the history with a model on disk. A rejected loop has one: on
`G0212-0335` loops 1 to 3 were accepted and loop 4 rejected, the data kept
were `phase3.ms`, and the residual flagging before the final image recorded
`"image": ".../images/phase_4"`. It subtracted the model of a candidate whose
solutions had been discarded from data that never had them.

Here the two models differ by under 1% (loop 4 missed acceptance by a tenth
of a percent of dynamic range) and the pass flagged 0.54% of the data. But a
loop is also rejected when its solutions make the image worse, and then its
model is wrong for the data kept: real emission stays in the "residual" and
the detector may flag it, up to the 10% the pass allows before it reverts.
Only accepted loops are considered now; with none, the initial image is used
as before.

### The first image from findings 96 to 101

`G0212-0335`, final imaging (multi-scale 0, 5, 15, 45 pixels, robust 0.5,
reference 700 MHz), from `imaging_report.json`:

| | this run | the run that prompted finding 96 |
|---|---|---|
| restored rms | 46.9 uJy/beam | - |
| residual rms | 47.1 uJy/beam | 87 uJy/beam |
| restored peak | 33.9 mJy/beam | - |
| dynamic range | 724 | - |
| `noise_over_thermal` | 4.19 | 7.7 |
| clean passes | bootstrap 14 major cycles, auto-multithresh 21, both stopped at nsigma | bootstrap only; auto-multithresh found no mask |
| model flux | 0.779 Jy | - |
| largest residual | 347 uJy/beam (7.4 sigma) at (2271, 3923) | - |

The thermal figure in `noise_over_thermal` is the natural-weighting floor of
the statistical weights, 11.2 uJy/beam; against CASA's sensitivity for the
robust 0.5 weights imaged (20.6 uJy/beam) the residual is 2.3 times thermal.
Imaging took 12.0 hours and the target 44.7 hours in all.

What is still wrong with it, from the preview images: each bright source
carries radial spikes and negative patches (the restored minimum is -4.8
mJy/beam, a seventh of the peak), and long linear stripes cross the field.
Those are gain errors that phase-only, direction-independent
self-calibration leaves behind. Three things bear on them and none is changed
here: amplitude self-calibration is off by default (`enable_amp_selfcal`),
and when on it is applied with no before-and-after test; the
self-calibration model is point sources only and holds 0.57 Jy of the 0.78
Jy the multi-scale image finds; and 74% of the target's data are flagged,
which finding 102 addresses in part. Source detection was skipped on this
target because PyBDSF is not installed in the pipeline's Python
(`No module named 'bdsf'`).

## Why the first corrected image was still limited

`G0212-0335` came out at 47 uJy/beam with spikes and negative patches round
every bright source. Findings 108 to 111 come from asking of that run's own
reports where the data and the calibration fell short. The answer is mostly
that the target was imaged from less than half of the data it should have
had.

About 40 to 45% of every scan is flagged for structural reasons: three dead
antennas (C03, S03, W04; 19% of baselines) and the part of the band above 850
MHz that the bandpass rejects. `G0213-0251` sits at that level, 49.3%.
`G0212-0335` was 73.8% flagged. It was observed in three scans (3, 15 and
29) at the start, middle and end of an eight-hour run, which is the best
hour-angle coverage a three-scan target can have, and it was imaged from two
of them with 22 antennas.

### 108. One failed calibrator solution cost an antenna the scans on both sides

The working MS's flag summary after the first calibration round reads, scan
by scan: 0.83 (scan 2, the first visit to the phase calibrator), 0.86 (scan
3, the target), 0.70, 0.67, 0.67, 0.61, 0.53, 0.53, 0.50, 0.45, and then 0.40
to 0.45 for the rest of the night. It steps down at each calibrator visit.
That is the signature of how `applycal` treats a time-dependent table: a
visibility is interpolated between the two solutions that bracket it and is
flagged when *either* is flagged. An antenna whose solution failed at one
visit loses the science scans on both sides of it.

The final gain tables show how little it takes. `phase_scan.cal` holds
twelve solutions per antenna; C02 had three flagged and C04 one. Both kept 6
to 7% of their data on this target, as did C12, S02 and W06
(`antenna_data.json`: about 1.0e6 unflagged samples each against 15.4e6 for
a working antenna). Five antennas were lost to a handful of failed
solutions, and with them a third of the baselines the target could have had.

A failed solution says that the calibrator's data for that antenna at that
visit could not be solved, and those data stay flagged. It says nothing about
the antenna half an hour later on the target, when the same antenna has good
solutions on either side. Interpolating across such a gap is what calibration
transfer ordinarily means.

For uGMRT, the two time-dependent tables a target receives (`phase_scan.cal`
and `flux.cal`) are now applied from copies, `phase_scan.target.cal` and
`flux.target.cal`, in which a failed solution is replaced when the same
antenna, in the same field and spectral window, has good solutions close by:

- between two good solutions no more than
  `target_transfer_bridge_max_gap_minutes` apart (120, which spans one missed
  visit of a calibrator observed about hourly), the linear interpolation of
  the two, amplitude and phase separately, phase along the shorter arc;
- at the start or end of the series, the nearest good solution if it is
  within half of that;
- never from another filled solution, and not at all for an antenna with
  fewer than `target_transfer_bridge_min_valid_fraction` (0.5) of its
  solutions good, which is a failing antenna and not a failed solution.

The solved tables are not modified and remain what the calibrators receive.
`applycal` still runs with `calflagstrict`, so whatever is left unbridged is
flagged as before. Every fill is listed in
`calibration/target_transfer_bridge.json` and summarised in the calibration
report. If a copy cannot be made or read, the targets get the solved table
and a warning. `ugmrt_target_transfer_bridge: false` in config.yaml switches
it off. Other telescopes are unchanged.

The phases this gives a bridged antenna are approximate, as any transferred
phase is; phase self-calibration then solves for that antenna against a model
built from the whole array. Together with finding 102, which stops the first
round's failures being inherited by the second, this is what should return
the missing scan and antennas. How much comes back is not known until it has
run: the selfcal report's `calibrated_coherence` (`antennas_with_data`,
`scans`) and the flagged fraction in the first self-calibration `applycal`
are the numbers to read.

### 109. Amplitude self-calibration was applied without looking at the result

Phase self-calibration stopped improving this field at a dynamic range of
267, with the artefacts that amplitude errors leave. Amplitude
self-calibration was off (`enable_amp_selfcal`), and as written could not
safely be turned on: the solve was applied to the accepted Measurement Set
and committed, with no image and no test. A phase loop has to improve the
dynamic range by 2% to be kept; the one step that can also move the flux
scale had no condition at all.

It is now a loop like the others. The `ap` solution (per scan for Band 4,
`solnorm=True` with `normtype='median'`) goes onto a copy of the
phase-calibrated MS, the copy is imaged exactly as a phase loop's is, and
`judge_amp_selfcal` keeps it only if

- the dynamic range rose by at least `amp_improvement_stop` (the phase
  loops' 2% unless set), and
- the peak moved by no more than `amp_max_peak_change` (10%). The solutions
  are normalised, so a larger shift means the solve has moved the flux scale
  to suit its model.

A rejected candidate is deleted and the phase-only MS stays in use; the
decision and the numbers are in the report as `amp_history`, and a resumed
run does not solve again. An accepted one becomes `final_selfcal.ms` as
before, and residual flagging then subtracts its model. The phase-calibrated
MS is no longer written to at any point.

Two smaller faults in the same place. Whether to attempt it was decided from
the *initial* image's dynamic range, before any phase loop: a field starting
at 60 and reaching 250 after phase self-calibration was never tried. The test
is now made on the image the solve starts from. And on resume after an
accepted solve the report lost `amp_table`.

It remains off by default. `selfcal_enable_amp: true` in config.yaml turns
it on, which changes the self-calibration recipe.

### 110. Antennas that could not be self-calibrated were imaged anyway

The five antennas of finding 108 had too little data for a self-calibration
solution, and with `calonly` (finding 97) their visibilities passed through
every loop with the phases the calibrator transfer gave them. Those are the
phases the first loop improved on by 70% for the rest of the array. About 3%
of the visibilities in the final image were therefore calibrated differently,
and worse, than everything around them, which is a direct source of the
structure round bright sources.

Finding 97's reasoning holds for data whose *increment* cannot be solved:
they keep the calibration of the loops before. It does not hold for an
antenna that never had a solution. When a loop's candidate MS is made, the
remaining data of antennas the census marks unsolvable (under 20% of the
median antenna's unflagged data) are now flagged in it. If the loop is
rejected the candidate is discarded and the flags with it.
`selfcal_flag_unsolvable_antennas: false` switches it off. With finding 108
these antennas should mostly have their data back and be solved normally;
this is for the ones that do not.

### 111. Self-calibration can be given a multi-scale model

Not a fault, an option. The self-calibration images are cleaned with point
components, and on this field that model held 0.57 Jy of the 0.78 Jy the
multi-scale final image found. `selfcal_scales: [0, 5, 15]` in config.yaml
passes those scales to every self-calibration clean. It is off by default
because a multi-scale minor cycle took about four times as long per
iteration in the final imaging of this target; it is the thing to turn on if
extended sources still carry artefacts once the data and the amplitudes are
right.

**What has not been verified.** None of findings 108 to 111 has run on data.
The table edits in 108 use the same table-tool calls as finding 102, on
copies; the gaincal and tclean arguments in 109 and 111 are ones the
pipeline already uses elsewhere (`normtype` is in CASA 6.7's `gaincal`
signature as logged). 696 tests pass, covering the bridging arithmetic and
its limits, which tables each field receives, acceptance and rejection of an
amplitude candidate on each condition, resume after either, and the
interruption paths.

### Findings 102, 103 and 108 on data

The run started on 2026-10-05 with findings 102 to 111 in place (same data,
new output directory, `selfcal_enable_amp: true`) is the first time any of
them ran. Through calibration, from its reports and CASA log:

| | this run | the run before |
|---|---|---|
| flagged after the round-1 applycal (whole MS) | 44.4% | 49.2% |
| flagged after corrected-data flagging | 46.6% | 51.1% |
| released before round 2 (finding 102) | 11.4% of all visibilities, to 35.2% flagged | not done |
| failed solutions bridged, round 2 (finding 108) | 7 of 61 in `phase_scan.cal`, 75 of 224 in `flux.cal` | not done |
| `G0212-0335` after the round-2 applycal | 48.3% flagged | 73.8% |
| `G0213-0251` | 41.6% | 49.3% |
| the other seven targets | 43.8 to 45.0% | not reached |
| antennas excluded | C03, S03 | C03, S03, W04 |
| calibration stage | 7.1 hours | 7.9 hours |

The release handled all 1,054,875 rows in blocks, none singly, and round 2
solved on the released data: W04, quarantined before for want of a delay
solution, was kept. The number of failed solutions barely changed between
rounds (60 and 213 in round 1), so the failures are in the calibrator data
themselves and not a product of the provisional flags; it is the bridging
that stops them reaching the targets. With the bridged tables the two gain
tables added nothing to seven of the nine targets' flags, 2.2 points to
`G0231-0533` and 4.1 to `G0212-0335`; on the phase calibrator, which is not
bridged, the same two tables added 13.9 points.

Targets were split, weighted at 97.7 kHz and averaged to 1024 channels of
390.6 kHz (`target_channel_average.json`), about 14 minutes each. The first
self-calibration image of `G0212-0335` then ran its major cycles in about
two minutes against six to seven before, and CASA's sensitivity for the
weights imaged is 15.0 uJy/beam against 20.6.

### What twice the data did, and did not do

`G0212-0335` in that run entered self-calibration 47.4% flagged, with 28
antennas and all three scans (`calibrated_coherence`), against 73.8%, 22 and
two before. Self-calibration took 12.5 hours for six images against 32.7:

| step | solint | dynamic range | residual rms (uJy/beam) | masked pass |
|---|---|---|---|---|
| initial | - | 78 | 96.7 | 8 major cycles, 23456 iterations |
| phase 1, accepted | 4min | 198 (+156%) | 92.5 | 8, 26237 |
| phase 2, accepted | 2min | 224 (+12.9%) | 89.5 | 9, 29494 |
| phase 3, accepted | 1min | 232 (+3.5%) | 88.4 | 11, 29219 |
| phase 4, accepted | int | 240 (+3.3%) | 88.0 | 10, 29532 |
| final | - | 241 (peak 25.2 mJy, rms 104.6 uJy) | 87.5 | 10, 28458 |

Every masked pass now ends at its nsigma threshold in 8 to 11 major cycles,
where the earlier run needed 13 to 23 and the 20-cycle limit: with the third
scan the PSF's sidelobes fell and the bootstrap halves the peak each cycle
instead of removing a third. The first-loop phase corrections span the full
plus or minus 180 degrees (`phase1_phase_vs_time.png`): the calibrator is
20 degrees from the targets and its phases transfer poorly.

The noise did not fall. The final residual is 87.5 uJy/beam against 72.7
before, and the final imaging's bootstrap ended at 80 against 72, while
CASA's thermal figure went from 20.6 to 15.0. So the image is not limited by
thermal noise or by the amount of data but by errors that scale with the
sources, and the recovered scan and antennas, whose amplitudes are bridged
or transferred across longer gaps, carry more of them. For N antennas with
amplitude errors e that are constant over each of M scans the dynamic range
is of order sqrt(M) N / e; with 28 antennas, 3 scans and e of 5 to 7% that
is 700 to 1000, where the images sit (724 before). Per-scan amplitude
self-calibration is the step that addresses exactly that, and it did not
run:

### 112. The amplitude solve could not succeed on a field like this

With `selfcal_enable_amp: true` the solve ran and its table was rejected:
91% of the solutions were flagged. Three things in the solve, none of them
the data's fault:

- `minsnr` was the band's `min_snr_amp`, 40. The model visibility of this
  field is about 8 mJy rms; one baseline, scan and polarization has about
  0.9 mJy of noise, so an antenna's per-scan amplitude solution reaches a
  signal-to-noise near 40. gaincal kept 7 of 81.
- One solution per polarization (`gaintype='G'`), which halves the data in
  each.
- A single reference antenna. Its own solution fell below the threshold in
  two of the three scans, and with no alternate CASA logged "refant (id=0)
  was flagged; flagging all antennas strictly" and discarded both scans
  whole. The phase loops have used a ranked list since finding 96.

The solve now uses that ranked list with `refantmode` as for phase,
`gaintype='T'`, and `minsnr` 10 (`DEFAULT_AMP_SELFCAL_MINSNR`,
`DEFAULT_AMP_SELFCAL_GAINTYPE`; `selfcal_amp_minsnr` and
`selfcal_amp_gaintype` in config.yaml). A solution good to 10% is admitted;
whether the set of them improves the image is finding 109's test, which is
the right place for that judgement. As with the major-cycle limit the
defaults are in `selfcal.py`, so the recipe of a self-calibration already
run is unchanged and its accepted phase loops can be kept: only the
amplitude step and what follows it need redoing. `min_snr_amp` in the band
tables and `amp_gaintype` in the task configuration are no longer read by
the solve. A solve that produces no usable table is now recorded in
`amp_history` (the report said `null`), and an accepted solve on a resumed
stage clears the `final_image` checkpoint, so the final image is remade from
the data now in use.

### 113. The gain-outlier summary was lost on every run

`calibration_stability` in the calibration report carried
`"antenna_gains_error": "The truth value of an array with more than one
element is ambiguous"`. `caltable_stability_summary` computed the outliers and
then evaluated `data.get("antenna_names") or []` on the numpy array of names
the table tool returns, which raises; the result was dropped before it was
stored. The same expression for field names lost nothing only because it
came after the assignment. Names are now converted without a truth test.
This is the summary that names antennas whose gains sit far from the array's,
which is what one wants when asking whether recovered antennas are sound.

**What has not been verified.** The new solve settings have not run. 704
tests pass.

### Finding 112 on data

The run was stopped on 2026-10-06, the files uploaded, and G0212-0335's
self-calibration re-entered. The orchestrator adopted ingest, flagging and
calibration, the four phase loops were restored from their checkpoints, and
only the amplitude solve ran (27 s). gaincal used `T`, `minsnr=10` and the
ranked reference antennas, as intended. Of the three per-scan intervals the
first gave 24 of 27 antennas a solution and the second all 27. The third gave
none, and the table was rejected for it. Findings 114 and 115 are what that
showed.

### 114. An 86-second scan counted as a third of the observation

The table test tolerates a quarter of a table's `TIME` groups failing, since
data without a solution keep the calibration they have (`calonly`). It
counted groups. G0212-0335's three scans hold 53, 125 and 8 integrations of
10.7 s: the last, at 22:08 UT, is what flagging left of the final scan, 86 s
in which only the central-square antennas and two of the south arm have
data. No antenna reaches a signal-to-noise of 10 on it in amplitude, CASA
then reports "refant (id=0) was flagged; flagging all antennas strictly" for
that interval, and one failed group of three is 33%. The solutions for the
other 96% of the observation were discarded with it.

`caltable_status` now weighs a failed `TIME` group by the observing time it
covers. `solution_time_shares` reads the `TIME` column of the Measurement Set
in use, gives each integration to the solution nearest to it, and takes a
solution's share as the fraction of rows it receives; the 25% tolerance
(`max_caltable_failed_time_fraction`, unchanged) applies to the summed share
of the failed groups. Here that is 8/186, 4.3%. The share of each group and
`failed_observing_time_share` are written to the table's `.qa.json`. The
test is stricter where it should be: one failed interval in five passed by
count even when it was the scan holding most of the data. When the
Measurement Set cannot be read, or was not observed at the solutions' times,
groups are counted as before. Antenna groups are still counted.

### 115. Per-scan amplitude solutions were interpolated between scans

The amplitude table was applied with `interp='linear'`, as the phase tables
are. With one solution per scan that draws a line from one scan's gain to the
next scan's, hours later. Worse, applycal does not correct a visibility that
lies between a good solution and a failed one (the behaviour finding 108
dealt with for the calibrators). With the last interval failed, every
visibility in the second half of the 22-minute scan lies between that scan's
solution and the failed one: a third of the data would have kept its old
amplitude, with a step in every antenna's gain at the scan's midpoint. The
three antennas without a solution in the first scan would have lost the
first half of the second scan as well, and so all of it.

The table is now applied with `interp='nearest'` (`AMP_SELFCAL_INTERP`), so
each visibility takes the solution of its own interval, which is what the
solve assumed. Data whose own solution failed stay as phase self-calibration
left them. The phase tables are unchanged: their solutions are seconds to
minutes apart within a scan, where linear interpolation is the right choice.

**Findings 114 and 115 on data.** Uploaded and restarted at 11:28 UT on
2026-10-06. The same solve (24, 27 and 0 antennas of 27 in the three
intervals) now passes the table test with `failed_observing_time_share`
0.043, and applycal logged `interp=nearest`. It also logged that 3.9% of
all visibilities, 7.7% of the unflagged ones, have no amplitude solution:
the last scan's remnant and the first-scan data of C13, S04 and W05. With
`calonly` they are imaged as phase self-calibration left them. What 115
describes for linear interpolation was not observed, since it was never
applied that way.

### 116. A resumed self-calibration reported none of what it restored

The report written by the resumed stage above had `clean_passes: {}`,
`residual_flagging: null` and `unsolvable_antenna_flagging: null`, although
five images had been cleaned and both flagging steps had run in the first
attempt. Those records were held in memory by the process that did the work;
a resumed stage restores accepted loops from checkpoints and repeats none of
it. Each clean pass is now written to `reports/clean_passes.json` as it
finishes and read back when a stage starts, an image made again replaces its
earlier entry, the residual-flagging record is read from the file it was
already saved to, and the unsolvable-antenna record is saved to
`reports/unsolvable_antenna_flagging.json`. Failing to write one of these
does not stop the stage. Reports of stages resumed before this change stay
incomplete; the numbers are in the CASA logs of the first attempt.

**What has not been verified.** Finding 116 has not run on data. 721 tests
pass.

### Amplitude self-calibration on data, and where that leaves the image

The amplitude candidate for G0212-0335 was accepted at 14:38 UT on
2026-10-06: dynamic range 240 to 287 (+20%), image rms 104.5 to 87.5
µJy/beam, peak +0.4%. That returns the self-calibration image to where the
first run's was (dynamic range 267) before it went on to a 47 µJy/beam final
image, so the final image of this run was expected near 45-55 µJy/beam
against a thermal figure of about 15. Sushant's requirement is an
artefact-free image at about 20 µJy/beam, and he stopped the run on that
expectation. Direction-independent self-calibration has now been taken as
far as it goes on this field without reaching it.

What the reports say about the remainder: the brightest source in the image
(25 mJy/beam apparent, at 02:12:05.8 -03:07:24) lies 29.5 arcmin from the
pointing centre. The uGMRT Band 4 primary beam there falls from about 0.42
at 550 MHz to 0.12 at 850 MHz, so the source is about 86 mJy seen through a
beam that changes by a factor of 3.5 across the band and steeply with
pointing; the deepest residual in the self-calibration image (-0.87 mJy)
is beside it. Errors towards such a source are not ones a single gain per
antenna can remove.

### 117. A measurement before the next run: `diagnose_noise_floor.py`

Every test of a change so far has been a pipeline run of many hours on a
target. `diagnose_noise_floor.py` answers the prior question in minutes, on
the self-calibrated Measurement Set that already exists and without changing
it. It makes dirty Stokes V images, of all the data and of each scan, on the
grid and weighting of the science image: with RR and LL the sky nearly
cancels, so their rms is the floor a perfect Stokes I image of these
visibilities could reach (an upper limit on it, since R and L gain
differences leak some Stokes I). It also measures the existing Stokes I
image's robust rms in a 10 by 10 grid of boxes, which shows whether the
excess is concentrated around particular sources or uniform. The report is
`noise_floor_test/<target>/noise_floor.json`.

If Stokes V is at or below the goal, the data are sound and the work is
towards the off-axis sources (peeling, or a spectral model that follows the
beam) and deconvolution. If it is well above, the limit is in the flagging
and weights, and no imaging change reaches the goal.

**First run, G0212-0335, 18:15-18:31 UT on 2026-10-06 (robust 0.5, 16
minutes with the stopped run's imaging process still sharing the node).**
In µJy/beam, robust rms:

| | all data | scan 3 (9.5 min) | scan 15 (22 min) | scan 29 (86 s) |
|---|---|---|---|---|
| dirty Stokes V | 23.6 | 74.8 | 25.0 | 92.1 |

The self-calibration Stokes I image has 80 overall, 59.5 in the quietest
tenth of its boxes and 141 in the box holding the brightest source; the
centre of the field is at 105-134 and the corners at 49-61. Its residual
has 64, and 47 in its quietest tenth.

What this shows. The data's own noise is at most 23.6, a third of what the
Stokes I image has anywhere, so the excess is sky, and it follows the sky:
highest where the sources are. What it does not show is the floor itself.
The dirty V image has a peak of 2.4 mJy/beam, a tenth of the brightest
source's Stokes I: at 29.5 arcmin the R and L beams differ by that much,
and the sidelobes of that leaked source are in the 23.6. The per-scan
values are affected more, their beams being worse, so scan 3's 74.8 against
the 39-44 its length predicts is a question, not a finding. The tool's own
sentence for this case said that no Stokes I image could reach the goal,
which a dirty image cannot support.

**Second version.** The tool now cleans the V image shallowly (5000
iterations, 5 sigma, four major cycles, continuing from the dirty products)
and reports the residual's rms as the floor, calling a dirty value an upper
limit where the clean is not made. It measures the point-spread function
(width, largest sidelobe, far-sidelobe rms) and takes several `--robust`
values in one run, which is the measurement behind a choice of weighting:
the sidelobes are what a lower `robust` buys and the V floor is what it
costs. The report is saved after each step.

**Second run, 19:15-19:39 UT (three weightings, all data).** The shallow
clean and the point-spread-function measurement ran on CASA without change.

| robust | beam | largest sidelobe | far-sidelobe rms | Stokes V dirty | Stokes V cleaned |
|---|---|---|---|---|---|
| 0.5 | 4.9" | 0.33 | 0.0019 | 23.6 | 22.5 |
| 0 | 3.7" | 0.28 | 0.0011 | 22.7 | 22.2 |
| -0.3 | 3.7" | 0.26 | 0.0009 | 25.9 | 25.7 |

(µJy/beam; the width is that of a circle with the beam's half-power area
and is quantised by the 1" pixel, so the last two are not distinguished.)

- Cleaning the leaked source out of V changes the floor by 5%: its far
  sidelobes, 2.4 mJy times 0.0019, add only 4.6 µJy in quadrature. The
  worry recorded above was larger than the effect. The floor of these data
  is 22 µJy/beam, half as much again as the 15 the weights predict; the
  goal of about 20 is the floor itself, not something below it.
- Robust 0 has the same noise per beam as 0.5, a beam a quarter narrower
  and far sidelobes 40% lower. Robust -0.3 costs 15% in noise for little
  more. That the noise does not rise from 0.5 to 0 says the short
  baselines 0.5 favours are worth less than their weights claim.
- The Stokes I image is at 80 (60 in its quietest tenth) against that 22.
  Far sidelobes of uncleaned faint sources account for perhaps 15-25 of
  it on any plausible source count, so most of the excess is not a
  deconvolution depth that a lower threshold would cure. A deep unmasked
  clean to three times the floor was considered as the next test and
  dropped for that reason: its threshold would be at the one-sigma level
  of the present image, where CLEAN removes anything, and it would take
  hours to say nothing about the cause.

**Third version: `--model-residual`.** For each weighting and each scan the
tool now images what the pipeline's own model leaves in the data: `tclean`
with `niter=0`, the self-calibration model images as `startmodel`, the
w-projection gridder the model was made with, and `savemodel='none'`. One
major cycle each, minutes not hours. Set beside the Stokes V of the same
scan and weighting it gives, per scan, how many times its noise the
residual is. A scan that carries calibration errors stands far above its
floor while the others sit near theirs; if the first scan, which is short,
was recovered by bridged calibrator solutions and has three antennas
without self-calibration solutions, turns out to hold the excess, leaving
it out costs little (scan 15 alone has a floor of 25.0 against 22.5 for
everything). The `startmodel` call has not run on CASA. 741 tests pass.

**Third run, 20:09-22:02 UT (G0212-0335, `--robust 0 0.5 --model-residual`).**
It took 1 h 53 min against the hour expected. The `startmodel` call ran
unchanged, and the robust 0.5 residual of all the data came out at 64.3
µJy/beam, the figure the pipeline's own residual image has, so the method
measures what tclean left. Robust rms in µJy/beam:

| all data | Stokes V (cleaned) | Stokes I model residual | ratio | quietest tenth | worst box |
|---|---|---|---|---|---|
| robust 0.5 | 22.5 | 64.3 | 2.9 | 46.6 | 122 |
| robust 0 | 22.2 | 45.0 | 2.0 | 31.6 | 94 |

| robust 0, per scan | Stokes V (dirty) | Stokes I model residual | ratio | quietest tenth | worst box | largest residual |
|---|---|---|---|---|---|---|
| scan 3 (9.5 min) | 61.0 | 199.9 | 3.3 | 131.6 | 604 | 31.5 mJy |
| scan 15 (22 min) | 23.9 | 38.0 | 1.6 | 28.0 | 79 | 1.1 mJy |
| scan 29 (86 s) | 110.9 | 214.2 | 1.9 | 145.8 | 389 | 19.1 mJy |

- The weighting alone takes the residual of the same model and data from
  64.3 to 45.0.
- The scans are not alike. Scan 15 by itself sits at 1.6 times its noise
  and at 1.2 times in its quietest tenth: nearly what its data allow, with
  the existing calibration. Scan 3 sits at 3.3 times, and in the box that
  holds the off-axis source its residual is 604 against scan 15's 79, far
  more than the ratio of their noise (2.6). Scan 3 is also 1.7 times
  noisier than its length predicts from scan 15 (61.0 against 36.7).
- The three add up: weighting the scans by their noise (86%, 13% and 4%)
  gives 33, 26 and 9 in quadrature, 43, against the 45.0 measured. A third
  of the variance of the combined residual comes from a scan with an eighth
  of the weight.
- Leaving scans 3 and 29 out would give 38 for a floor of 23.9 in place of
  45 for 22.2, with the model as it is; a model made from scan 15 alone no
  longer has to average over what scan 3 sees and may do better.
- Not understood: scan 3's largest residual, 31.5 mJy/beam, exceeds the
  brightest source in the image (25 mJy/beam) and cannot be reconciled
  with the 1.1 mJy/beam of the combined residual at a 13% weight. Its
  position was not recorded; the tool should report it.

### 118. A WSClean image of the same data, as a test: `experiment_wsclean_image.py`

Sushant asked whether WSClean would do better than tclean here. For the
imaging, probably: it deconvolves sub-bands jointly with a polynomial
through them, which suits a source at 29.5 arcmin whose apparent spectrum is
the primary beam's (a factor of 3.5 across the band) better than two Taylor
terms; its mask and depth follow a local noise estimate; and WGridder is
faster than w-projection. For the excess noise, only in part: the gains are
the same whichever imager reads them. The pipeline's WSClean backend serves
MeerKAT and LOFAR and refuses uGMRT, so before any of that is rebuilt the
question is put to the data.

The script images one target's self-calibrated Measurement Set with WSClean
3.6 (on ilifu with WGridder, EveryBeam and IDG) on the pipeline's grid:
`-gridder wgridder`, `-join-channels -channels-out 6 -fit-spectral-pol 3`,
`-auto-mask 4 -auto-threshold 1.5 -local-rms`, multi-scale 0, 5, 15 pixels,
robust 0 by default, no primary-beam correction. It first makes a dirty
Stokes V image with the same weighting, which is the noise of the data as
WSClean weights them, and reports the Stokes I image against it and against
the tclean image, in the same grid of boxes. Three things it takes care of:

- `-no-update-model-required` and `-data-column DATA`: the Measurement Set
  is a pipeline product, its `DATA` holds the self-calibrated visibilities
  and its `MODEL_DATA` the pipeline's model; WSClean writes to neither.
- `-channel-range`, measured from a flag summary: Band 4 is recorded to
  950 MHz and usable to 850, and equal sub-bands over all 1024 channels
  would leave one and a half of six empty.
- `-abs-mem 20 -j 8`: WSClean sizes itself from the node's memory, not the
  session's 29 GB.

Not yet run. The WSClean command line is tested as text and the run against
stand-ins; whether WSClean 3.6 accepts every option as written is for the
first run to show. 760 tests pass.

**First attempt, 2026-10-07 09:08 UT.** It did not start. On ilifu WSClean
is not a program on the session's PATH: it is reached with `singularity exec
/software/astro/containers/oxkat-latest.sif wsclean`, and that container
holds WSClean 3.1.0, not the 3.6 whose version text Sushant had pasted. The
script looked for an executable called `wsclean`, found none, and said so;
he then passed the script to WSClean as if it were a Measurement Set. Two
changes: `--wsclean` takes a whole command, split into words, so a container
wrapper is given in quotes; and before launching, the script asks that
WSClean for its version and help text and builds the command from what it
offers (`-use-wgridder` where there is no `-gridder`; `-local-rms` and
`-parallel-deconvolution` only where listed). A command that does not answer
as WSClean stops the run before any work. 769 tests pass.

**Second attempt, 09:50 UT.** `ModuleNotFoundError: No module named
'casatools'`, from a Python 3.13. The terminal that can run the container
is not the one whose `python3` is CASA's: one has CASA and no WSClean, the
other the reverse, and the script asked for both. It now needs neither CASA
nor anything outside the standard library to launch WSClean:

- the grid and the channels that hold data are taken from `--imsize`,
  `--cell` and `--channel-range` where CASA cannot be imported, with a
  message naming the option if one is missing. For this dataset the range
  is 51 to 768: the flagging report masks native channels 0-203 (5% edge)
  and 3072-4095 (above 850 MHz) of 4096 ascending from 550 MHz, and the
  targets are averaged by four;
- the WSClean images are read by a FITS reader of some twenty lines that
  uses numpy alone; without numpy the images are still made, and
  `--measure-only` measures them later from a Python that has it (CASA's),
  which also adds the comparison with the tclean image;
- both tools drop their own folder from the import path before importing
  anything, as the pipeline's modules do: the folder holds the pipeline's
  `inspect.py`, which numpy would otherwise import in place of the standard
  library's.

Run end to end here without CASA against a stand-in that answers as WSClean
3.1.0 does. 778 tests pass.

**Third attempt, 10:00 UT: the first WSClean numbers, and a crash.** The
script found WSClean 3.1 in the container, chose `-use-wgridder` from its
help, and the dirty Stokes V image of scan 15 at robust 0 came out at 23.4
µJy/beam; CASA had measured 23.9 for the same scan and weighting, so the two
imagers agree on what the data hold. One inversion of the six sub-bands took
77 s. The Stokes I run then wrote its six beams (5.2 by 3.6 arcsec in the
top sub-band) and six dirty images and aborted as deconvolution began:
`double free or corruption (!prev)`, exit status -6. Nothing in the command
was rejected; this is a fault inside WSClean 3.1.0. The options that first
act at that moment are multi-scale, the local-rms mask and the splitting of
the image for parallel deconvolution, and the last is the one with a history
of memory faults in that release. It is now left out before WSClean 3.2
unless a size is given (`--parallel-deconvolution`), keeping multi-scale and
the local-rms mask, which are the point of the test. Whether that is the
cause is for the next run to show. The log is now flushed line by line so
progress can be followed from outside. 779 tests pass.

**Fourth attempt, 10:10 UT.** It was not the cause: without
`-parallel-deconvolution` WSClean aborted at the same point, after the sixth
dirty image and before the first line of deconvolution, this time with
`munmap_chunk(): invalid pointer`. Two different complaints from the
allocator at one place say the heap is being overwritten there, by code
that runs as deconvolution is set up: multi-scale, the local-rms image or
the spectral fit. oxkat, whose container this is, runs WSClean with
`-local-rms`, joined channels and a spectral fit and without multi-scale,
which points at multi-scale but proves nothing.

Guessing one option per round trip costs Sushant an upload and a run each
time, so the script now finds out for itself. `feature_ladder` lists the
deconvolution set-ups from fullest to plainest, each dropping one of
multi-scale, the local-rms mask and the spectral fit; when WSClean dies on a
signal (a negative exit status, or 134/139 through a wrapper) the crash's
log is kept as `...I.crash<n>.log`, what it left of the image is removed,
and the next rung runs. An exit that is not a signal, such as a rejected
option, is not retried. The report lists every attempt and says which
set-up made the image. Run here against a stand-in that aborts when asked
for multi-scale: the second rung completed. 787 tests pass.

**Fifth attempt, 10:20-10:54 UT: the first WSClean image.** Started without
multi-scale (`--scales` with no values), it did not abort: the fault in this
build is in multi-scale. Scan 15, robust 0, local-rms mask, three-term
spectral fit; five major iterations, 3366 components, 34 minutes (inversion
22, prediction 6, deconvolution 2). Robust rms in µJy/beam:

| | whole image | quietest tenth | worst box |
|---|---|---|---|
| WSClean dirty Stokes V | 23.4 | | |
| WSClean Stokes I, this run | 64.3 | 42.6 | 171 (off-axis source) |
| tclean's model, same scan and weighting (residual) | 38.0 | 28.0 | 79 |

It is worse than what tclean's model leaves, and the log says why. With
`-local-rms` the first major iteration measures the noise on the dirty
image, where sidelobes and the sources themselves raise it (lowest local
value 39 µJy, global estimate 109). Scaled by that, everything was under the
4-sigma mask threshold after 542 components, in one major iteration, and the
mask was frozen there. The remaining four iterations cleaned only those
pixels. The residual still peaks at 13.7 mJy/beam and falls to -6.4: the
off-axis source has a mask a few pixels wide and cannot be deconvolved
inside it. WSClean's own estimate of the residual noise went 109, 79, 71,
66, 64 and stopped.

Two incidental observations. The source peaks at 41.9 mJy/beam here against
25.1 in the tclean image, which is not a disagreement: WSClean's joined
image is the band average and tclean's first Taylor term is the flux at its
750 MHz reference, and for a spectrum this steep (the beam's) the average is
about 1.6 times that. And the last major iteration began on a peak of
-0.48 mJy/beam at pixel 2560,2560, the phase centre exactly, which is where
signals that do not fringe end up; it should be looked at in the tclean
residual too.

The next run turns the local-rms mask off (`--no-local-rms`), so the mask
threshold follows a global noise estimate that falls as cleaning proceeds,
with `--auto-threshold 1.0`. No change to the script.

**Sixth run, 11:00-11:32 UT: 30.8 µJy/beam.** Scan 15, robust 0, global
auto-mask at 4 sigma, cleaning to 1 sigma inside it, three-term spectral
fit, no multi-scale, six threads. Six major iterations, 141 205 components,
30 minutes (inversion 11, prediction 7, deconvolution 12). WSClean's noise
estimate fell 109, 85, 58, 42, 34, 31 as the mask grew over four major
iterations before it was frozen. Robust rms in µJy/beam:

| | whole image | quietest tenth | worst box |
|---|---|---|---|
| WSClean dirty Stokes V | 23.4 | | |
| WSClean Stokes I | 30.8 | 24.1 | 47 (off-axis source) |
| its residual | 30.6 | 24.1 | 42 |
| tclean's model, same scan and weighting (residual) | 38.0 | 28.0 | 79 |
| the pipeline's image (all scans, robust 0.5) | 80 | 60 | 141 |

The quietest tenth is at the noise of the data. On the same visibilities and
the same calibration WSClean leaves a fifth less than tclean's model does
overall and half as much around the off-axis source. Residual extremes are
+0.24 and -0.21 mJy/beam.

**What the image shows.** Read into the browser and displayed (whole field
at +-0.3 mJy/beam; 160-pixel cut-outs at +-5% of peak around the four
brightest sources). The background is noise. Every bright source carries
the same pattern: an hourglass of positive wings with negatives between
them, reaching about -10% of the source's peak 10 to 35 arcsec away:

| source (pixel) | peak (mJy/beam, band average) | deepest negative beside it |
|---|---|---|
| 3162,4226 (29.5 arcmin off axis) | 37.8 | -3.74 |
| 2562,2566 (phase centre) | 23.1 | -2.15 |
| 1641,2102 | 20.6 | -2.03 |
| 3226,1859 | 15.6 | -1.87 |
| 383,2079 | 11.5 | -1.19 |
| 4148,2906 | 11.4 | -1.16 |

These are in the model, not the residual: inside its mask CLEAN fitted the
pattern with negative components. The same fraction and shape at the phase
centre as at 29.5 arcmin rules out the primary beam, pointing, and
smearing, all of which vanish on axis. It is a direction-independent error
in the visibilities: the response of the data to a point source is not the
point-spread function their weights imply. The sub-band images at 593 and
640 MHz show it separately (8-13% each) with source positions steady to
0.4 pixel, so it is not a disagreement between sub-bands either.

This withdraws the reading recorded under "Amplitude self-calibration on
data": the off-axis source is not what limits this image, and peeling it is
not the next step. Candidates are gains solved against tclean's incomplete
model (0.70 Jy, 28 000 components, against 141 000 here), which biases
amplitudes by baseline length, and errors within the scan that one
amplitude per scan cannot follow.

Solving again against the WSClean model is not safe as it stands: that
model contains the pattern, as negative components, and gains solved
against it would be solved to keep it. Before any solve, a cheaper
question: does the amplitude scale move within the scan? If every source
were, say, 15% brighter at the end of the 22 minutes than at the start (a
ramp is what linear interpolation between two calibrator scans gives when
their gains differ), the image of the whole would carry the difference
between the early and late point-spread functions around every source, on
axis or off, in each sub-band, which is what is seen. `--time-split N`
makes N dirty Stokes I images of consecutive stretches (WSClean
`-intervals-out`, no deconvolution), finds the brightest sources in their
mean and reports each one's flux per stretch against its own mean, with the
median over sources. A common change of more than 5% is called a moving
scale. Run here against a stand-in with a 24% ramp, which it reported as
24%. 793 tests pass.

**Time split on data, 12:21-12:23 UT.** Four stretches of scan 15, eight
sources of 12 to 46 mJy/beam, interval noise 225 µJy/beam: relative fluxes
1.008, 0.975, 0.998, 1.006, with a source-to-source scatter of 1.3 to 3.5%.
The amplitude scale does not move through the scan. That candidate is out.

What remains is an error that is the same for all directions, steady in
time and present in each sub-band: gains that differ between antennas in a
fixed way. If the central-square antennas sat on a lower amplitude scale
than the arm antennas, a point source would be seen with more flux on long
baselines than on short ones, and its image would be the sharp long-baseline
beam with too little of the broad short-baseline one: a core with negatives
around it, shaped by the array's three arms. `--uv-split` tests that the
same way as the time split: one natural-weighted dirty image per range of
baseline length (WSClean `-minuv-l`, `-maxuv-l`), the brightest sources
located in the longest range and measured in each, their fluxes compared.
The shortest range is reported and not judged, since at its resolution
neighbouring sources blend. Run here against a stand-in with a 20% rise.
799 tests pass.

**uv split on data, 12:29-12:31 UT: the test answered a different question.**
Peak flux in mJy/beam, natural weighting, no deconvolution:

| source (pixel) | 0.6-3 kλ | 3-8 | 8-20 | 20-35 | 35-70 |
|---|---|---|---|---|---|
| 3162,4226 (off axis) | 112.5 | 86.6 | 59.5 | 51.1 | 28.4 |
| 2562,2567 (phase centre) | 104.8 | 80.1 | 44.8 | 24.7 | 13.4 |
| 1641,2102 | 42.1 | 38.1 | 31.6 | 28.1 | 20.7 |
| 383,2080 | 19.3 | 27.2 | 20.7 | 13.1 | 9.5 |
| 4148,2906 | 22.5 | 25.9 | 18.6 | 14.1 | 9.0 |

Every source falls with baseline length, by factors of 2 to 6 between 3-8
and 35-70 kλ and by a different factor each. The amplitude scale cannot do
that, since it would move them together; source structure can. The
brightest sources of this field are resolved: the target at the phase
centre goes from 105 mJy at arcminute resolution to 13 mJy in a 4 arcsec
beam. The premise of the test, a field of point sources, was wrong, and so
was the sentence the tool printed ("antennas on different scales"): a fall
that differs between sources is now read as resolved sources, and only a
rise as a fault of the scale. The same explains the sub-band peaks recorded
above (0.67 to 0.84 between 593 and 640 MHz): a resolved source loses peak
flux as the beam shrinks, on top of its spectrum.

**What that does to the reading of the image.** The 30.8 µJy/beam image was
cleaned with point components only, because multi-scale is what aborted in
WSClean 3.1.0, and cleaned to one sigma inside its mask. Point-component
CLEAN driven deep into extended emission is the textbook source of negative
holes and of ripples that follow the beam's sidelobes, the same around
every extended source and in each sub-band, on axis or off. That fits what
was seen at least as well as a calibration error does, and unlike one it
is known to be present. So "direction-independent error in the
visibilities" above is withdrawn as a conclusion: the pattern is first of
all a property of the deconvolution used, and whether anything remains once
the sources are cleaned with scales is not yet known.

Multi-scale has so far been tried in this build only together with the
local-rms mask (both aborts). The run that worked had neither. The next run
is multi-scale without the local-rms mask; if that aborts too, this WSClean
cannot do it and a newer one, or tclean, has to. 800 tests pass.

**Seventh run, 12:35-13:25 UT: multi-scale works, and is not the answer.**
Scales 0, 5, 15, 30 pixels with the global mask: no abort, so the fault in
WSClean 3.1.0 is multi-scale together with the local-rms mask. 29.0
µJy/beam (30.8 with point components), quietest tenth 23.5, worst box 42,
residual 28.7 with extremes of +0.18 and -0.18 mJy/beam; 49 minutes, 34 of
them deconvolution. The negative beside the brightest source is -3.62
mJy/beam where it was -3.74. Point-component cleaning was not its cause,
and that reading, entered above an hour earlier, is withdrawn in turn.

**The pattern, measured.** Cut-outs of 73 pixels around six sources away
from the centre, each divided by its peak and averaged, then split into the
part that is unchanged by a half-turn about the peak and the part that
changes sign:

- the part that changes sign reaches +13% of the peak at (+8, -9) pixels
  and -13% at (-8, +9); rms over the annulus 2.3%, the same as the
  symmetric part's 2.2%;
- each of the sixteen brightest sources (3.1 to 37.4 mJy/beam, 0.1 to 41
  arcmin from the phase centre, all position angles) has that same
  antisymmetric pattern, with coefficients 0.7 to 1.2 of the average (two
  faint outliers at 1.6 and 2.7); the source at the phase centre has 1.02;
- it is in the dirty images, which no deconvolution has touched, at 1.2
  times the strength it has in the restored one; in the dirty images of the
  first and last quarter of the scan it is 1.21 and 1.15;
- in the six sub-band images it has one sign throughout: 0.97, 1.07, 1.08,
  1.11, 0.85, 0.54 from 593 to 827 MHz.

A point-spread function is unchanged by a half-turn whatever the weights,
so an antisymmetric response to every source is imaginary visibility, that
is phase. Identical on axis and off: not the beam, pointing or smearing.
One sign across the band: not delays, which change sign about the band
centre. Constant through the scan. What is left is a fixed phase per
antenna, some tenths of a radian on a few antennas or about 0.1 rad rms
over all of them.

Four loops of phase self-calibration at decreasing intervals did not remove
it, and the likely reason is the oldest one. A response common to every
source is exactly what antenna phases produce, so the data cannot tell "the
antennas are offset" from "every source has this shape"; only the model
can, by insisting that sources are compact and positive. The pipeline's
first model is cleaned to 3 sigma inside a mask, far below lobes that are
13% of 25 mJy, so the lobes went into the model, as positive and negative
components, and each solve reproduced them. That is inference from the
measurements above; the solve that would confirm it has not been made.

The test: a WSClean model of the cores only (no negative components, an
absolute threshold above the lobes of the brightest source), one phase per
antenna for the scan against it, applied, and the field imaged again as in
the seventh run.

### 119. The test of fixed phase offsets, in three commands

- `experiment_wsclean_image.py --core-model 6`: WSClean with
  `-threshold 0.006 -no-negative`, joined sub-bands and the spectral fit, no
  mask and no scales, and without `-no-update-model-required`, so that it
  writes the model to `MODEL_DATA`. Six mJy/beam keeps the cores of the
  eight sources above it and is above the 4.9 mJy/beam lobe of the
  brightest. It is the only mode that writes to a Measurement Set and
  refuses one inside the pipeline's products.
- `experiment_phase_offsets.py` (CASA): `gaincal` with `calmode='p'`,
  `solint='inf'`, `gaintype='G'` and the pipeline's reference antennas
  against that model; prints the phase of each antenna and polarization and
  their rms; `applycal` (`calonly`, `nearest`) and `split` to
  `<name>_ph.ms`. An rms of 3 degrees or more is called fixed offsets.
- the image of `<name>_ph.ms` as in the seventh run.

The image report now carries `source_asymmetry`, the measurement made by
hand above: cut-outs of the brightest sources divided by their peaks and
averaged, the largest value and the rms of the part that changes sign under
a half-turn, and the median over sources of the deepest negative beside
them. Before and after are then two numbers.

None of the three has run on data. 820 tests pass.

**The three steps on data, 13:44-15:03 UT.** All ran unchanged.

*Cores-only model.* 170 components in two major iterations, 196 mJy in 71
pixels, none negative.

*Phases against it* (degrees, R/L; C00 is the reference):

| | |
|---|---|
| central square | C01 3/3, C02 -1/-2, C04 -2/-2, C05 -4/-4, C06 1/0, C08 0/0, C09 -3/-2, C10 0/0, C11 -4/-4, C12 5/6, C13 0/0, C14 -2/-2 |
| east arm | E02 -7/-6, E03 -7/-7, E04 -4/-4, E05 2/2, E06 -4/-3 |
| west arm | W01 -9/-9, W02 -14/-14, W03 -20/-20, W04 -12/-13, W05 9/8 |
| south arm | S01 -5/-5, **S02 -82/-82**, **S04 26/26**, S06 -10/-10 |

Rms 17 degrees, twelve antennas beyond 5, the two polarizations within a
degree of each other throughout. The offsets are real, they sit on the
arms, and one antenna was 82 degrees out through four loops of phase
self-calibration. The inference of the previous section is confirmed as
far as the offsets go.

*The corrected data, imaged as in the seventh run:*

| | before | after |
|---|---|---|
| peak of the brightest source (mJy/beam) | 37.4 | 46.2 |
| peaks of the next seven | | +25% to +37% |
| image robust rms (µJy/beam) | 29.0 | 28.2 |
| worst box | 42 | 38 |
| the antisymmetric pattern of the seventh run, per source | 1.0 | 0.05 to 0.45 |
| negative beside the brightest source | -9.7% | -4.1% |
| negative beside the others | about -10% | -7% to -10%, at (+3, +8) pixels |
| antisymmetric rms in the dirty image, eight sources | 3.9% | 2.4% |

A quarter of every source's peak had been scattered out of it. The wide
lobes at 13 arcsec are a third of what they were. What is left is tighter
and was there before: in the dirty images of both the uncorrected and the
corrected data the antisymmetric part peaks at 14-15% some 7 pixels from
the source, nearly north-south, the same for most sources. North-south is
the south arm, whose antennas had the largest offsets, and S02 and S04 are
the two that gave no solution at 4, 2 and 1 minute in the pipeline's loops
and did at one integration, which is how an antenna behaves whose phase
does not hold still. One correction for the scan takes out a wandering
antenna's mean and nothing else. This is a reading of the numbers; the
solve that shows each antenna's phase through the scan is the next step.

**For that step.** `experiment_phase_offsets.py` now reports, for any
`--solint`, each antenna's mean phase and the rms of its phase about that
mean (on the circle, so wraps do not count), names those offset by more
than 5 degrees and those moving by more than 10, joins solutions linearly
when there are several, and solves amplitudes with `--calmode ap`
(normalised to a median of one). The cores-only model step measures the
antisymmetric part in its own dirty image, which reads the state of the
calibration in ten minutes without a fifty-minute image. 828 tests pass.

**Round 2 on data, 15:10-18:46 UT.** Phases at one minute against a new
cores-only model of the once-corrected data (210 mJy in 59 pixels), 23
solutions per antenna:

- no antenna wanders. The rms of every antenna's phase about its own mean
  is 2 degrees or less. The reading offered above, that S02 and S04 move
  within the scan, is wrong, and one correction per scan is the right kind;
- the means are not zero. Central square within 2 degrees, east arm +1,
  west arm 2, 5, 5, -2, 0, south arm S01 -10, S02 -16, S04 +12, S06 -5. The
  first solve had left the south arm 10 to 16 degrees short.

The antisymmetric part in the dirty image, eight brightest sources:

| data | peak | rms |
|---|---|---|
| as the pipeline left them (measured by hand, slightly different sources) | 15% | 3.9% |
| after round 1, one phase per antenna | 11.2% | 2.2% |
| after round 2 | 6.5% | 1.7% |

It converges, a round at a time, as self-calibration against a model made
from still-distorted data does: the first model's cores were themselves
displaced by an 82-degree error, so the first solve could not be exact.
The peak stays at (-1, -7) pixels as it shrinks, so it is one error being
worn down, not a new one appearing.

**Round 3 and its image, 18:54-19:36 UT.** One more phase per antenna
against the model left by round 2: rms 2.4 degrees, largest 6.6; S02 -7,
S01 -6, W02 and W03 +5 (as in round 2), S04 +3, the rest within 2. Then the
field as in the seventh run:

| | uncorrected | round 1 | round 3 |
|---|---|---|---|
| image robust rms (µJy/beam) | 29.0 | 28.2 | 27.5 |
| quietest tenth | 23.5 | 23.3 | 23.0 |
| worst box | 42 | 38 | 37 |
| peak, brightest source (mJy/beam) | 37.4 | 46.2 | 45.0 |
| peak, target at the phase centre | 23.6 | 30.8 | 33.3 |
| antisymmetric part, restored image, eight sources | about 13% | 9.0% | 3.8% |
| deepest negative beside them, median | about -10% | -8.9% | -4.7% |
| most negative pixel (mJy/beam) | -3.62 | -3.04 | -1.84 |

Sources are 20 to 43% brighter at their peaks than the pipeline's
calibration left them, and the artefacts beside them have halved twice.
They are not gone: -4.7% of 45 mJy/beam is 2 mJy/beam, seventy times the
noise.

**What is left differs across the band.** In the sub-band images of round
3 the antisymmetric part is 0.9% rms (peak +4.0%) at 593 MHz and 1.8% rms
(peak -11.8%) at 827 MHz, and the two patterns, compared at radii scaled by
the frequency ratio as the beam scales, correlate at -0.43. Before any
correction the six sub-bands agreed in sign (0.5 to 1.1 of the average
pattern): the fixed offsets dominated. With those removed, what remains
reverses between the ends of the band, which is what a delay does once a
single phase has been taken out of it: zero at the centre, opposite at the
ends. A 0.4 ns delay runs the phase through 40 degrees across 280 MHz.

`experiment_phase_offsets.py --delays` solves a delay per antenna first
(`gaincal`, `gaintype='K'`, one per scan), prints each in nanoseconds with
the phase it runs through across the band, solves the gains with the delays
applied, and applies both. 834 tests pass.

**Round 4 on data, 19:52-20:40 UT: delays.** Before it, the dirty-image
reading on the round-3 data: 5.2% peak, 1.7% rms, against 6.5% and 1.7%
after round 2. A third phase-only round had bought nothing, as it should
not if what was left was not a single phase.

Delays against a cores-only model of 216 mJy (ns, R/L):

| | |
|---|---|
| S04 | 2.26 / 2.21 |
| S02 | 0.71 / 0.61 |
| W03 | 0.44 / 0.37 |
| C09 | 0.42 / 0.43 |
| W02, C06, C12, W01, C01, C11 | 0.20 to 0.37 |
| W04 | -0.13 / 0.75 |
| the other sixteen | under 0.2 in magnitude |

2.26 ns turns S04's phase through 228 degrees across 280 MHz: over the band
its contribution to every baseline it forms averaged to little. The two
polarizations, solved separately, agree to about 0.1 ns except on W04,
which is either a fault in one hand of that antenna or a poor solve. The
phases solved on top are large (S04 83 degrees) only because CASA refers a
delay to the spectral window's reference frequency and the phase takes up
the difference; it is the pair that is applied.

The field, as in the seventh run:

| | uncorrected | round 1 | round 3 | round 4 |
|---|---|---|---|---|
| image robust rms (µJy/beam) | 29.0 | 28.2 | 27.5 | 27.1 |
| worst box | 42 | 38 | 37 | 35 |
| peak, brightest source (mJy/beam) | 37.4 | 46.2 | 45.0 | 49.1 |
| peak, target at the phase centre | 23.6 | 30.8 | 33.3 | 37.3 |
| antisymmetric part, eight sources: peak | about 13% | 9.0% | 3.8% | 2.1% |
| antisymmetric part: rms | about 2.3% | 1.0% | 0.65% | 0.27% |
| deepest negative beside them, median | about -10% | -8.9% | -4.7% | -3.9% |

The target is 58% brighter at its peak than the pipeline's calibration made
it. The antisymmetric part, which is phase, is an eighth of what it was.
The negatives have not followed it down: about -4% of the peak beside each
of the ten brightest sources, now 6 to 9 pixels away and with a negative or
nothing at the opposite point. In the averaged cut-out the symmetric part
(0.70% rms) is now twice the antisymmetric (0.33%). A symmetric error is
amplitude. The pipeline's amplitude self-calibration ran while S04 and S02
were decorrelated by these delays and offsets, and an amplitude solve
answers lost coherence by raising the antenna's gain; with the coherence
restored those gains are too high. That is the next thing to measure.

`experiment_phase_offsets.py --uvrange` restricts the gain solve to given
baselines, for an amplitude solve against a cores-only model, which lacks
the extended emission that short baselines see. 836 tests pass.

**Round 5 on data, 20:51-21:32 UT: amplitudes, and a mistake in applying
them.** Before it, the dirty-image reading on the delay-corrected data:
2.2% peak, 0.60% rms (5.2% and 1.7% before the delays).

Amplitude and phase per antenna and polarization against a model of
everything above 3 mJy/beam (336 mJy in 164 pixels), baselines beyond 3 kλ.
Phases 5 degrees or less. Amplitudes, R/L, where 1 is the median:

| | |
|---|---|
| S02 | 1.21 / 1.35 |
| S04 | 1.10 / 1.29 |
| W04 | 0.75 / 0.68 |
| S06 | 0.74 / 0.84 |
| central square | 0.90 to 1.17, most above 1 |
| east arm | 0.86 to 1.08, most below 1 |
| west arm (W01-W03, W05) | 1.12, 1.01-1.06, 0.90-0.96, 0.83-0.90 |

S02 and S04 are high, as expected of antennas whose gains were raised while
they were decorrelated. Eighteen antennas are more than 10% from the median,
and the two hands differ by up to 0.2 (C12, E03, S04), which a
polarization-averaged solve, the pipeline's, cannot see. There is also a
trend that should not be trusted: the central square mostly above one and
the outer arms below. A model of cores lacks the extended emission that
shorter baselines see, and that is the sign such a lack gives the solve,
even beyond 3 kλ.

The image:

| | round 4 | round 5 |
|---|---|---|
| dirty Stokes V (µJy/beam) | 23.6 | 26.5 |
| image robust rms | 27.1 | 30.6 |
| image over Stokes V | 1.15 | 1.15 |
| deepest negative beside the sources, median | -3.9% | -2.6% |
| most negative pixel (mJy/beam) | -1.91 | -1.18 |
| antisymmetric part: peak, rms | 2.1%, 0.27% | 3.0%, 0.31% |

The negatives fell by a third, so the amplitudes were part of what was
left. The noise of the data rose by 12%, and that is the tool's fault: it
applied the amplitudes with `calwt=False`. Dividing W04 by 0.7 raises its
noise by 1.4 and its weight has to drop by half to match; left unchanged,
the antennas that were turned up count as much as quiet ones. The ratio of
image to data noise is unchanged at 1.15, so the loss is all in the floor.

Corrected: the weights follow whenever amplitudes are applied
(`calwt=True` for that table only). And for the bias, the image run can now
make the model itself: `--no-negative --write-model` cleans the whole field
as usual but with positive components and leaves the result in
`MODEL_DATA`, so the amplitude solve has the extended emission in its model
and can use every baseline. 841 tests pass.

**Round 5 redone with weights, and imaged with positive components only,
2026-10-08 05:55-07:20 UT.** The dirty Stokes V is 23.4 µJy/beam: with the
weights following the amplitudes the floor is where it was (23.6), where
the mistaken application had it at 26.5.

| | round 4, usual clean | round 5, positive components |
|---|---|---|
| image robust rms (µJy/beam) | 27.1 | 30.8 |
| quietest tenth | 23.1 | 25.2 |
| worst box | 35 | 43 |
| residual: robust rms, most negative | 26.9, -0.16 mJy | 30.5, -0.91 mJy |
| most negative pixel of the image (mJy/beam) | -1.91 | -0.51 |
| deepest negative beside the sources, median | -3.9% | -1.3% |
| antisymmetric part: peak, rms | 2.1%, 0.27% | 2.7%, 0.31% |

The two columns differ in two things at once and must be read with that in
mind. The amplitudes took the negatives from about 4% to 1.3% of the peak
and the deepest pixel from -1.9 to -0.5 mJy/beam. Forbidding negative
components took 27.1 to 30.8: what the usual clean had been fitting with
negative components beside every source, and so removing together with its
sidelobes across the field, is now left in the residual, and the quietest
tenth of the image rises from 23.1 to 25.2 with it. The 27.1 was in part
the look of an image whose errors had been modelled away. 30.8 against a
floor of 23.4 is the more honest statement of where the calibration is:
errors worth 20 µJy/beam in quadrature remain.

The run took 75 minutes, against 45: writing the model column at every
major iteration more than tripled inversion and prediction. The model is
now in `MODEL_DATA` of the round-5 copy for the amplitude solve on all
baselines.

**Round 6, 07:32-08:52 UT: amplitudes against the whole field.** Amplitude
and phase per antenna and polarization against the positive model of the
whole field, every baseline: all amplitudes between 0.93 and 1.04, none
more than 10% from the median where eighteen had been; phases within 2
degrees. The central square now reads slightly low (0.93 to 1.01) and the
arms slightly high (1.00 to 1.04): the reverse, at a third of the size, of
what the cores-only model had given. The first amplitude solve had
over-corrected the central antennas by about 5%, the bias suspected of it,
and this one takes that back. Phase, delay and amplitude per antenna have
each stopped moving.

Positive-component image, as in the round before:

| | round 5 | round 6 |
|---|---|---|
| image robust rms (µJy/beam) | 30.8 | 30.1 |
| quietest tenth | 25.2 | 24.9 |
| worst box | 43 | 42 |
| most negative pixel of the image (mJy/beam) | -0.51 | -0.48 |
| deepest negative beside the sources, median | -1.3% | -1.2% |
| antisymmetric part: peak, rms | 2.7%, 0.31% | 3.2%, 0.33% |

Looked at with yesterday's display (cut-outs at +-5% of peak, the field at
+-0.3 mJy/beam): the hourglass around every bright source is gone. The four
brightest are compact, the target at the phase centre with what looks like
real structure to its south. What remains is faint radial structure around
the brightest, at the level of 1% of their peaks, and a noise in the
source-free corners of the image of 24 to 25 µJy/beam against 23.5 for the
data.

**Where the scan stands against where it began (G0212-0335, scan 15):**

| | pipeline (tclean, robust 0.5, all scans) | WSClean, pipeline's calibration | after rounds 1-6 |
|---|---|---|---|
| image noise (µJy/beam) | 80 (60 quietest) | 29.0 (23.5) | 30.1 with positive components (24.9); 27.1 with the usual clean before the amplitudes |
| peak of the target (mJy/beam) | | 23.6 | 39.2 |
| negative beside bright sources | | about -10% | about -1% |

What was wrong, in the order found: the imaging (tclean's mask and depth,
robust 0.5, two Taylor terms); a first scan worse than its share of the
data; fixed phase offsets of up to 82 degrees that four self-calibration
loops had frozen into their own models; delays of up to 2.3 ns that a
single phase per antenna cannot follow; amplitudes solved while those were
present, one per antenna for both polarizations, and applied by me at first
without their weights. None of these is specific to this target.

### 120. The six rounds as one command: `selfcal_recipe.py`

Rounds 1 to 6 took two days because each was two commands in two terminals
(WSClean in one container, CASA's Python in another session) with a person
carrying the result between them. Nine targets cannot be done that way.
ilifu has CASA containers beside the WSClean one
(`/software/astro/containers/casa-modular-*.sif`), so one script that needs
neither can run both.

`selfcal_recipe.py` is that script. It is plain Python, started from the
terminal that has `singularity`. It copies the Measurement Set it is given
(never writes to it, and refuses a work folder inside the pipeline's
products), then runs the sequence that worked by hand, each model with
`experiment_wsclean_image.py` and each solve with
`experiment_phase_offsets.py` under the CASA container:

| stage | model | solve | repeated |
|---|---|---|---|
| `phase` | cores above 6 mJy/beam, positive | phase per antenna per scan | up to 4 times, until the largest phase found is 10 degrees or less |
| `delay` | the same | delay per antenna, phase on top | once |
| `amplitude_cores` | cores above 3 mJy/beam | amplitude and phase per polarization, baselines beyond 3 kλ | once |
| `amplitude` | whole field, positive components, written by the image run | amplitude and phase, every baseline, weights following | up to 2 times, until amplitudes are within 10% of one and phases within 5 degrees |
| `final` | | image with positive components; nothing written to the data | |

The limits are what the hand-run showed: phases found of 82, 16 and 7
degrees, then amplitudes against the whole field all within 7% of one. A
round is ended on a small correction, and never on a number its solve did
not report. On scan 15 the plan therefore runs exactly the hand-run's
steps, with one difference: round 2 by hand solved a phase every minute (to
see whether any antenna wandered; none did), and here it is one phase for
the scan like the others.

What it does besides chaining:

- **Asks first.** Before anything is copied it runs `wsclean --version`
  and has the CASA command import `casatasks` and open the Measurement Set
  read-only, which also shows the folder is visible from inside the
  container. Without `--casa` it tries the four CASA containers in turn
  and uses the first that answers. A tool that does not answer stops the
  run in its first minute, with what it said.
- **Continues.** `recipe_report.json` in the work folder records every
  finished step, what is running and since when, and what each round found
  (lopsided part of the dirty image, largest phase, delay and amplitude
  with their antennas, each image against the noise of its data). After a
  failure or an interruption the same command starts at the step that did
  not finish. A work folder that holds a run on other data is refused;
  `--restart` discards an earlier run, and only a folder that this script
  made.
- **Cleans up after itself.** Each round's copy of the data, with its model
  and corrected columns, is three times the size of the data; it is
  removed when the next copy exists (`--keep-intermediate` keeps them).
  The gain tables and each solve's report stay.

Two things in it are this field's and not yet general. The thresholds of
the cores-only models (6 and 3 mJy/beam, `--core-threshold` and
`--amplitude-threshold`) are a quarter and an eighth of G0212-0335's
brightest peak before correction; another field needs its own, and they
should come from that field's dirty image. And the image grid and channel
range are given on the command line because the terminal that has
`singularity` has no CASA to read them from the pipeline's products.

**One change to the solve for more than one scan.**
`experiment_phase_offsets.py` joined solutions linearly whenever an antenna
had more than one. That was written for several solutions within a scan.
With one solution per scan and several scans it would have drawn a line
from one scan's offset to the next one's across the hour between them, and
corrected each scan with a mixture of the two (finding 115, in another
place). One solution per scan (`--solint inf`) is now applied to its own
scan only. For scan 15 alone nothing changes.

887 tests pass. The stand-ins in `tests/test_selfcal_recipe.py` leave what
each script leaves, so the order of the steps, the Measurement Set each
acts on, the repetition limits, continuing after a failure, and the guards
are tested; the commands the driver builds were also given to the two real
scripts' own argument parsers and to `experiment_wsclean_image.py
--dry-run`, which printed the WSClean commands of the hand-run. What is
not tested here is the CASA container itself: whether its `python3` has
`casatasks` and sees `/idia` is what the first minute of the run finds out.

**Expected on scan 15**, since the steps are the hand-run's: phases found
near 82, 16 and 7 degrees (S02), a largest delay near 2.3 ns (S04), a
final image near 30 µJy/beam with positive components (25 in its quietest
tenth) against 23.5 for the data, and negatives beside the bright sources
near 1% of their peaks. A result far from these is a fault in the
automation, not in the data. About three and a half hours.

**On data, 2026-10-08 09:20-14:30 UT (`recipe_scan15`).** The check before
starting found WSClean 3.1 and CASA 6.7.6.14 in
`casa-modular-py3.12-latest.sif`, the first container tried. The run then
went through without a failure and without a hand on it:

| round | automated | by hand |
|---|---|---|
| `phase1` | 82.3 deg (S02); S04 +26, W03 -20, W02 -14, W04 -12, S06 -10, W01 -9 | the same to the decimal |
| `phase2` | 16.0 deg (S02) | 16 |
| `phase3` | 6.6 deg (S02), so no fourth round | 7 |
| lopsided part of the dirty image before each round | 13.4, 11.2, 6.5, 5.2, 2.1% | 11.2 and 6.5 at rounds 2 and 3 |
| `delay` | 2.26 ns (S04) | 2.26 |
| `amplitude_cores` | 0.68 to 1.35 | 0.68 to 1.35 |
| `amplitude1` model (positive, whole field) | 30.7 µJy/beam, 25.2 quietest tenth, data 23.4, negatives -1.3% | 30.8, 25.2, 23.4, -1.3% |
| `amplitude1` solve | 0.93 to 1.04, phases within 2.1 deg; not repeated | 0.93 to 1.04, within 2 |
| `final` | 30.0 µJy/beam, 24.9 quietest tenth, data 23.4, negatives -1.2%, antisymmetric 2.8% / 0.30% | 30.1, 24.9, 23.5, -1.2%, 3.2% / 0.33% |

WSClean took the same path through the first model as by hand (81, then
170 iterations to 6 mJy). The automation is the hand-run. (The 83 degrees
reported for S04 in the delay round is the constant that goes with its
delay, not a new offset; and "peak" in the report's reading is the
brightest pixel of the image, an off-axis source at 50.7 mJy/beam, where
the 39 mJy/beam quoted before is the target.)

It took 5 h 10 min, not three and a half hours. The cores-only models took
20 to 27 minutes each against about ten by hand (first model: inversion
11:28 against 3:50, deconvolution 3:42 against 0:41, so the processor and
not only the disk); the whole-field model 70 minutes, as by hand; the
final image 69. And the JupyterLab server of the session answered only at
09:20, 09:35, 10:24-10:25, 13:16-13:29 and from 14:44: it was silent
whenever WSClean was at work and came back within minutes of each WSClean
step ending. Yesterday's runs by hand, with the same six threads, cost it
an occasional time-out. What differs is not known (the node's load, or
how many cores the session has today); `nice` on the command is the first
thing to try, since a server starved of processor time by its own
session's job is what this looks like.

**What is left in the image, measured on its residual.** Around the seven
brightest sources (11 to 51 mJy/beam), in rings about each, the residual's
rms is:

| ring (pixels = arcsec) | 8-15 | 15-30 | 30-60 |
|---|---|---|---|
| residual rms (µJy/beam), brightest to faintest | 131, 131, 113, 89, 56, 99, 55 | 95, 83, 89, 74, 55, 81, 55 | 72, 67, 63, 64, 49, 56, 42 |
| point-spread function, rms sidelobe | 9.6% | 6.3% | 3.2% |

against 24.9 µJy/beam in the quiet parts of the image and 23.4 for the
data. Three things follow.

- The structure does not scale with the source: a source four times
  fainter has more than half of it. Fitted as a part proportional to the
  peak plus a constant, it is 0.24% of the peak in the inner ring and
  0.12% in the outer, on about 40 µJy/beam that every bright source has
  around it whatever its flux. Divided by the sidelobe level, the
  proportional part is what gain errors of 2.5 to 4% would leave.
- It is only partly the same from source to source. The residuals around
  the five brightest, each divided by its peak, correlate at 0.1 to 0.4
  where one error common to the field would give about 0.7. Part is
  therefore common, which a direction-independent solve can still remove
  (gains varying within the 22 minutes are the candidate: the four
  stretches of the scan differed by up to 2.5% on 2026-10-07, which was
  "steady" when the errors were 10% and is not at 0.2%), and part is each
  source's own. The brightest source is 29.5 arcmin from the pointing
  centre, the fifth 27 and the sixth 37, on the flank of the primary beam
  and beyond, where one gain per antenna for the whole field cannot be
  right.
- The "deepest negative beside the sources" (-1.2%) is no longer a
  measure of calibration. Every one of those minima lies within 7 arcsec
  of its source's peak, under the source, and the residual there
  correlates only weakly with the point-spread function (0 to -0.2). All
  seven sources are resolved (the model holds 1.5 to 3.3 times the peak
  within an arcminute), and a multi-scale clean that may not use negative
  components cannot take back what a large component put in too much.
  That constraint was right for the models that calibration solves
  against. It is not right for the image that is kept.

So the scan stands at 24.9 µJy/beam away from bright sources, 6% above its
data, with ripples of 50 to 130 µJy/beam within an arcminute of seven
sources. And the data's own noise is 23.4 µJy/beam for this scan and about
22 for all three: an image of this observation at 20 µJy/beam does not
exist to be made.

**Do the gains move within the scan? 2026-10-08, about 15:30 UT.**
Amplitude and phase per antenna and polarization every minute, against the
positive whole-field model in the round-5 copy
(`experiment_phase_offsets.py --calmode ap --solint 1min --suffix _t1`,
CASA from the container; 22 solutions per antenna):

| | amplitude, rms about its mean over the 22 minutes | phase, rms |
|---|---|---|
| central square (12 antennas without C12) | 0.7 to 2.6% | 0.3 to 1.7 deg |
| arms (E02-E06, S01, S04, S06, W01, W02, W04, W05) | 1.9 to 3.3% | 0.3 to 1.1 deg |
| W03; S02 right-hand | 0.8 to 1.3%; 1.6% | 0.3 to 0.8 deg |
| **C12**, both polarizations | **12.4% and 9.2%** | 0.6 to 0.7 deg |
| **S02**, left-hand only | **11.9%** | 0.9 deg |

The phases are steady, as the one-minute phase solve of round 2 had found.
The amplitudes are not. Noise in a solution moves its amplitude, as a
fraction, as much as its phase in radians; the arm antennas' phases move
by 0.3 to 1.1 degrees, which is 0.6 to 1.9%, and that includes whatever
their phases really do, so most of their 1.9 to 3.3% in amplitude is the
antennas and not the noise. Both polarizations of an antenna move by the
same amount. Two are of another size: C12 in both polarizations and S02 in
one move by about a tenth.

That is the size the residual asked for. A scale that is off by 2 to 3%
for a few minutes at a time leaves, around every source, that fraction of
the point-spread function of those minutes, whose sidelobes are several
times those of the whole scan: of order 0.1 to 0.2% of the peak, where the
part of the residual that scales with the source was measured at 0.12 to
0.24%. One amplitude per antenna for the scan, which is what the recipe
solves, cannot follow it. Whether these variations are the antennas' own
or are the off-axis sources seen through a primary beam that turns on the
sky, taken up by a solve that has one gain for the whole field, the solve
cannot say; the image of the corrected copy
(`wsclean_test/G0212-0335_scan15_ph_r2_r3_r4_r5_t1.ms`) will, by whether
the source at the pointing centre improves with the others or at their
expense.

### 121. Two runs in one work folder destroyed each other

The image of that corrected copy was started on 2026-10-08 at 20:03 UT and
was dead nine minutes later, in its first step, with

    Error opening temporary data file
    'wsclean_test_scan15_t1/G0212-0335/tmp-part0005-V-b0.tmp'

The log of that step holds the output of two WSClean processes laid over
each other, one sentence cut off by the other's ("Cleaning up temporary f
== Constructing image =="), and two exceptions. The command had been
started twice. Both copies were given the same work folder, so they wrote
the same image names and the same temporary files, and the one that
finished its sub-band first removed the files the other was reading. Each
then failed on a file the other had taken. Nothing was running from 20:12
until the next morning.

Nothing in either script noticed. For the imaging script that costs a run;
for `selfcal_recipe.py` it would cost more, since each round removes the
previous round's copy of the data, and a second copy of the driver would
remove it from under the first.

A work folder is now held by the run that is using it.
`experiment_wsclean_image.claimed(workdir)` writes `running.lock` (process,
machine, time) with an exclusive create before anything is made in the
folder, and removes it when the run ends, however it ends. A second run
given the same folder stops at once and says who holds it. A claim left by
a run that was killed is taken over when its process no longer exists on
this machine; a claim from another machine cannot be checked, so the
message names the file to remove. Every branch of the imaging script that
runs WSClean takes the claim, and the driver takes one on its own folder
before it touches anything in it, including under `--restart`, which now
empties the folder only once it holds it. 895 tests pass.

On data, 2026-10-09 10:54 UT: the command started once ran through, with
one process holding the folder, and under `nice -n 19` the JupyterLab
server answered in under two seconds throughout, where the day before it
had been silent for hours. A session starved of processor time by its own
job is what it was.

### 122. Two images of the minute-by-minute copy, and what they overturn

**Positive components, like for like (11:29 UT).** Amplitude and phase
every minute instead of once for the scan:

| | one amplitude per scan | one per minute |
|---|---|---|
| image robust rms (µJy/beam) | 30.0 | 29.7 |
| quietest tenth | 24.9 | 24.6 |
| noisiest box | 42.2 | 40.6 |
| residual rms in rings of 8-15 / 15-30 / 30-60 arcsec about the 51 mJy source | 131 / 95 / 72 | 121 / 93 / 70 |
| the same about the target (39 mJy) | 131 / 83 / 67 | 131 / 85 / 67 |
| about the 33, 24, 19, 18 and 11 mJy sources | | each within 2 µJy/beam of before |

I had expected the rings to fall by a third to a half. They did not move.
The amplitudes do vary within the scan as measured, and correcting that
is worth one per cent of the noise; it is not what is left around the
sources, and a stage for it does not belong in the recipe.

**The usual clean of the same data (12:08 UT).** 27.9 µJy/beam, 23.9 in
the quietest tenth against 23.4 for the data, median box 28.2, noisiest
35.4; residual within +-0.16 mJy/beam. It shows two things.

*The positive models hold flux that is not there.* The usual model sums to
1.02 Jy (1.34 positive, 0.32 negative); the positive one to 1.60 Jy.
Within an arcminute of the seven brightest sources the positive model
holds 128, 131, 59, 80, 28, 39 and 19 mJy where the usual one nets 106,
107, 42, 65, 25, 30 and 15: 15 to 40% more. That excess is the negative
bowl measured in the positive image's residual (about -36 µJy/beam over
the arcminute around the brightest, which is -20 mJy). Positive components
cannot take back what a large one overshot, so the model ends above the
data, and it does so in extended emission, which the short baselines see.
The amplitude solve against that model then found the central square at
0.93 to 0.99 and the arms at 1.00 to 1.04, and applying it scaled the
shortest baselines up by some 8% and the longest down by 4%. That is the
"trend reversed" of round 6: the cores-only model had too little extended
flux and pushed the centre one way, the positive whole-field model has too
much and pushed it the other. Neither is the antennas. The whole-field
model for amplitudes has to be a usual clean, whose flux is the data's;
positive components belong to the cores-only models of the phase and
delay rounds, where they are what breaks the degeneracy and where a
threshold of several mJy keeps the overshoot out. `selfcal_recipe.py`
still has the positive whole-field model and is to be changed, with the
final image, which should likewise be a usual clean.

*There are holes beside the bright sources that are not noise.* In the
restored image the deepest negative within an arcminute of each of the
seven is -1.1 to -2.9% of its peak (-0.27 to -1.10 mJy/beam), and 4 to 11%
of the pixels in that arcminute are below -3 sigma where noise would put
0.1%. Displayed at +-0.3 mJy/beam, the target and the off-axis source
each sit in spokes of alternating sign, in the same orientation in both.
Spokes through a source are the fan beams of the array's arms; a spoke
left behind means the baselines along one arm are not subtracted as the
model says, and the same orientation in two sources half a degree apart
says the error is in the antennas rather than in the direction.

**What has not been solved for.** Each antenna now has, per scan, a delay
and a phase, which together are a straight line in phase across the band.
The ionosphere's part of a phase goes as one over the frequency. Over 570
to 850 MHz that is not a straight line: if the 228 degrees that S04's
2.26 ns delay took out across the band were ionospheric, the line through
it leaves a bow of 22 degrees at the band centre, some 8 degrees rms, on
every baseline to that antenna; S02, W03 and C09 in proportion. No step so
far could see it, because every solve averaged over the band.

`experiment_phase_offsets.py --bandpass N` is the test: a gain per antenna
for every N channels against the model (`bandpass`, `solint='inf,Nch'`,
one set per scan), reported per antenna as the phase's rms across the
band, the same with the best straight line removed (what a delay could
not have taken), and the amplitude's rms. Phases only are applied unless
`--bandpass-amplitude` is given: amplitudes across the band depend on how
well the model has each source's spectrum through the primary beam, so
they are reported and left. `--prior` applies an earlier solve's table
first and together with the result. 909 tests pass.

### 123. The gains against frequency: the model is a staircase, and W04's upper band is wrong

`--bandpass 30` on the round-5 copy with the minute-by-minute table applied
first (2026-10-09 13:16 UT): 25 bins of 30 channels (11.7 MHz) hold data,
for 27 antennas.

**The bow is not there.** With the best straight line removed, S04's phase
moves by 1.8 and 3.3 degrees rms across the band, not eight; most antennas
by 1 to 2. Its 2.26 ns was a delay. The hypothesis of finding 122 is
wrong.

**What is there, per antenna.**

- Slopes that a delay should have taken and did not, on the central
  square: C09 -22 and -24 degrees across the band, C06 -10 and -20, C12
  -13 and -16, C10 -11 in one polarization, S06 -12 in one. The delays
  were solved against the cores-only model, which lacks what the short
  baselines of these antennas see.
- **W04.** In both polarizations its amplitude is right below about
  735 MHz and then climbs: 1.13, 1.25, 1.42, 1.60, 1.93, 2.21, 2.20 from
  740 to 815 MHz, with the phase swinging by 16 degrees either way in one
  polarization and no solution at all above 820 MHz. Its visibilities in
  the upper third of the band are up to 2.2 times what the sky gives on
  every other antenna.
- S04, S02 and E06, mostly in one polarization, fall to 0.5 to 0.75 in the
  last three bins (above 820 MHz).

**What is there in every antenna alike.** The median amplitude over the
antennas, bin by bin, is a sawtooth: 1.049 falling to 0.982, 1.019 to
0.973, 1.007 to 0.965, and so on, the high points in bins 1, 6, 10, 14, 18
and 22. WSClean's six sub-bands begin at channels 51, 170, 290, 409, 529
and 648, which are bins 1, 6, 10, 14, 18 and 22. A thing all antennas
share is not an antenna. The model is one image for each sub-band, a
staircase in frequency with steps 47 MHz wide, and the sky falls smoothly
across each step: about 2% either way in a gain, so 4% either way in a
visibility, which is an apparent spectral index near -1 for the field as a
whole and much steeper for a source on the flank of the primary beam.

That is the residual that would not move. It is a fixed fraction of each
source (the part of the rings that scaled with the peak, 0.12 to 0.24%);
it differs from source to source, because each has its own apparent
spectrum (correlations of 0.1 to 0.4); it is worst for the sources
farthest out (holes of -2.2 and -2.9% at 29.5 and 37 arcmin, at the same
offset from both); it has the form of the point-spread function stretched
along the radius, which is spokes of alternating sign; and no calibration
touches it, because it is in the model that every calibration and every
image so far has used (rings unchanged to 2 µJy/beam through three
calibrations). Six sub-bands across 570 to 850 MHz are too few.

**Changes.** `--bandpass-amplitude` now corrects each antenna's amplitude
as its ratio to the median of all antennas in that bin
(`common_amplitude`), so that the model's steps are never written into the
data; the shared part is reported with its size. With amplitudes applied,
channels of an antenna that the whole scan could not solve are flagged
(`applymode='calflag'`, `nearest,linearflag`) rather than given the
nearest gain: W04 above 820 MHz is not to be corrected by what it needed
at 815. The phase-only correction is as it was. 912 tests pass.

**Next, in one command:** the solve with amplitudes
(`--suffix _bpa`), a usual clean of the result with six sub-bands, which
isolates W04 and the others against today's 27.9 / 23.9, and the same
with twelve (`--channels-out 12`), which halves the steps.

**The same accident, in the one script without the guard (16:12 UT).** The
command was started at 16:12:41 and again at 16:12:55, before the new
script had been uploaded. The two `bandpass` solves wrote the same table:
one failed with "Directory::remove: ...bandpass.cal/OBSERVATION is not
empty", the other with "error in open or create of file", both before
anything was applied, so the data are untouched. A third start a minute
later, with the new script, stopped on the broken table they left. The
claim of finding 121 covered the imaging script and the driver; the solve
script now has it too (`claimed`, a file `<stem><suffix>.lock` beside the
products, taken in `main` before anything is removed or written), and the
commands given from now on append to their log (`>>`) so that a second
start cannot empty the first one's. 916 tests pass.

**On data, 16:20-17:59 UT.** One run, started once. The solve with
amplitudes, then two usual cleans of the corrected copy
(`..._r5_bpa.ms`):

| | before (finding 122) | bandpass corrections, 6 sub-bands | the same, 12 sub-bands |
|---|---|---|---|
| noise of the data, dirty Stokes V (µJy/beam) | 23.4 | 22.5 | 22.5 |
| image robust rms | 27.9 | 27.3 | 26.7 |
| quietest tenth | 23.9 | 23.7 | 23.2 |
| median box, noisiest box | 28.2, 35.4 | 27.6, 33.7 | 26.9, 32.3 |
| most negative pixel (mJy/beam) | -1.10 | -0.84 | -0.83 |

Within an arcminute of the seven brightest sources, deepest pixel as a
percentage of the peak, and the number of pixels below -75 µJy/beam:

| source | before | 6 sub-bands | 12 sub-bands |
|---|---|---|---|
| 51 mJy, 29.5 arcmin off-axis | -2.15, 1268 | -1.65, 916 | -1.63, 1001 |
| 39 mJy, the target | -1.09, 965 | -0.80, 669 | -0.89, 676 |
| 33 mJy | -1.57, 1205 | -1.49, 638 | -1.56, 644 |
| 24 mJy | -1.32, 1168 | -1.25, 1090 | -1.15, 992 |
| 19 mJy | -1.41, 528 | -1.17, 338 | -1.19, 313 |
| 18 mJy, 37 arcmin off-axis | -2.94, 991 | -2.69, 894 | -2.72, 845 |
| 11 mJy | -2.73, 473 | -2.01, 268 | -1.69, 206 |

The per-antenna corrections did what one expects of them: the noise of
the data itself fell by 4%, the first step to lower the floor rather than
the image above it (W04's upper band had been adding its doubled noise to
everything), and a third of the negative area around the sources went.

Twelve sub-bands did not do what I said they would. The image is 2%
quieter everywhere, quiet parts included, which is what a smoother model
is worth; but the holes beside the sources are where they were, to a
tenth of a per cent of the peak. The staircase is real and is not what
makes them. That is the second prediction about these holes to fail, and
the reasoning of finding 123 that tied them to the model's steps is
withdrawn: a model error that changes sign across each sub-band spreads
thinly and far from the source, it does not dig next to it.

**The holes, sub-band by sub-band (18:20 UT).** From the twelve sub-band
images of that run, around the same seven sources.

- The sub-bands are far from equal. Their noise is 88, 76, 66, 63, 65, 66,
  73, 82, 95, 108, 128 and 226 µJy/beam from 581 to 835 MHz: the top one
  has a twelfth of the weight of a good one, and the upper third of the
  band a third of the weight of the lower two.
- In the top sub-band every source's peak is displaced by 0.2 to 1.1
  arcsec, each in its own direction. But leaving it out of the average
  changes no hole by more than 0.05% of a peak. It is too light to matter.
- Each hole is at one position in every sub-band where it is seen ((-4,
  -19) pixels from the 51 mJy source from 651 to 766 MHz, (2, 8) from the
  33 and 18 mJy sources). A sidelobe would move by three pixels over that
  range. These are negative components of the model, which is one set of
  positions for all sub-bands.
- As a fraction of the source the holes grow with frequency, and most for
  the sources farthest from the pointing centre:

| source, distance off-axis | deepest hole, 581-697 MHz | 720-812 MHz |
|---|---|---|
| 39 mJy, the target, 0 | -0.7% | -1.5% |
| 24 mJy, 16 arcmin | -1.0% | -2.4% |
| 33 mJy, 17 arcmin | -1.1% | -2.8% |
| 19 mJy, 27 arcmin | -1.0% | -3.2% |
| 51 mJy, 29.5 arcmin | -1.4% | -3.0% |
| 18 mJy, 37 arcmin | -2.2% | -5.1% |

  In the single sub-bands the 18 mJy source goes from -1.5% at 581 MHz to
  -8% at 789, while its apparent flux falls from 35 to 10 mJy/beam as the
  primary beam narrows past it.

An error that is larger the farther a source is from the pointing centre
and the higher the frequency, that no antenna-based solve for the field
as a whole has touched, is what antennas that do not all point alike, or
whose beams are not all alike, give: on the flank of the primary beam each
antenna then sees the source with its own gain, the more so the steeper
the flank. That is the reading that fits; it has not been tested. The
test is to solve gains towards one such source alone, with the rest of
the field subtracted, and see whether they differ from antenna to
antenna by the several per cent this needs. The target's own holes (at
45 arcsec, -0.9%) are not explained by this and may be the faint
extended emission of the source itself, left under the 4-sigma mask.

### 124. One source alone: `experiment_extract_source.py`, and the recipe brought up to date

Sushant's decision, 2026-10-09: the direction-dependent work is to be done
on this field before any production run.

**The test.** Gains towards one source can only be solved where that
source is alone. `experiment_extract_source.py` makes that Measurement
Set: from the per-sub-band model images of a WSClean run it writes a set
with a disc about the source set to zero (the rest of the field, with its
negative components, which is what the data hold); WSClean predicts that
into the model column of a copy of the data; and CASA, running this same
script with `--casa-part`, subtracts it (`uvsub`), keeps the result, and
moves the phase centre to the source (`phaseshift`). The source's position
comes from the images' own sine projection. The data given are copied and
never written to; the folder is claimed as in finding 121; the copy with
its model column and the models of the rest are removed when they have
served.

`selfcal_recipe.py` on that Measurement Set, with a small image, is then
self-calibration in one direction, with everything the recipe has learnt:
a cores-only positive model first, then delays, amplitudes, and the gains
across the band, which is where an error that grows with frequency should
show. `--image-first` images the data before any round, so that the holes
before and after are measured on the same footing.

**The recipe brought up to date**, since this test runs it and the
production run will:

- the whole-field models and the final image are usual cleans; positive
  components are kept for the cores-only models (finding 122);
- a `bandpass` round after the amplitudes, `--bandpass 30
  --bandpass-amplitude` on a whole-field model (finding 123);
- `--channels-out`, default 12, passed to every model and image
  (finding 123: worth 2% of the noise);
- the report reads a bandpass solve by what runs across the band.

945 tests pass. Not tested here: WSClean's `-predict` from the
per-sub-band images, and CASA's `uvsub` and `phaseshift`, which the unit
tests stand in for.

**What would count as an answer.** For the 51 mJy source 29.5 arcmin
out: before, its deepest hole is -1.6% of the peak in the whole band and
-3.0% in the upper half. If antennas see it differently, the solves
towards it will find gains that differ from antenna to antenna by several
per cent, more at the top of the band, and the holes will shrink to a
fraction. If the gains come out at one within their noise (1 to 2%) and
the holes stay, the error is not antenna-based in this direction either,
and what is left is the baselines themselves or the deconvolution of this
source.

### 125. The first source alone: the error is in its direction, and one round of the recipe was biasing amplitudes

**On data, 2026-10-09 19:16-19:39 UT** (`extract_test/source1`,
`recipe_source1`). One run. The three steps that had not run on data
worked: WSClean predicted the rest of the field (860 mJy) from the twelve
model images, CASA subtracted it and moved the phase centre, and the
source's model alone holds 97 mJy. The recipe on the result took 21
minutes, the image being 1024 pixels.

| round, towards this source alone | found |
|---|---|
| before any round | peak 50.0 mJy/beam; eight holes within 80 arcsec between -0.6 and -1.5% of the peak |
| `phase1` | up to 8 degrees (S02), W02 -6, W03 -4 |
| `delay` | up to 0.28 ns, differing between the two polarizations of an antenna (C02 +0.09 and -0.14, S02 +0.16 and -0.10, S06 0 and +0.21); 14 antennas more than 5 degrees from the reference |
| `amplitude_cores` (cores, beyond 3 kλ) | 0.86 to 1.08: W03 0.86/0.88, S06 0.87/0.93, W04 0.89, W05 0.92/0.94, S02 0.94/0.92, the central square 1.00 to 1.08 |
| `amplitude1` (usual clean of the whole image) | 0.99 to 1.02 |
| `bandpass` | per antenna, phases up to 13 degrees rms across the band, amplitudes 5 to 13% rms about what the antennas share |
| final image | peak 55.4 mJy/beam; one hole at -1.15%, four between -0.4 and -0.5% |

**The holes went, all but one.** The hole 19 arcsec south that had been the
deepest in every image of the field (-1.5%) is at -0.4%; the one 7 arcsec
away went from -1.25 to -0.45%. One to the north-east stays (-1.3 to
-1.15%). The report's single number for "negative beside the sources"
(-1.5 to -1.15%) hides this, being only the deepest.

**The gains towards this source are the antennas' own.** Across the band,
each antenna's amplitude over the median of all antennas:

- below 700 MHz it moves by 1.7% rms, and the two polarizations of an
  antenna agree (correlation 0.55 over 324 pairs): what an antenna
  pointing a little off does to a source on the flank of its beam, the
  same in both polarizations;
- above 700 MHz the scatter grows to 6 to 7% and the polarizations no
  longer agree overall (-0.06), much of which is noise, the source having
  fallen from 90 to 27 mJy/beam and the band's noise having trebled. But
  single antennas run smoothly and far, and oppositely in the two
  polarizations: S02 +19% in one and -21% in the other around 760 MHz,
  W03 -23% in one. Noise does not run smoothly over eight bins. Opposite
  signs in the two polarizations on the flank of the beam is what beams
  that are displaced from each other give.

So the reading of finding 123 holds where it could be tested: towards a
source 29.5 arcmin out the antennas do not have one gain, and giving each
its own takes most of the holes away.

**But the peak's rise of 11% is not to be believed, and the recipe is
changed for it.** Phases of a few degrees cannot raise a peak by a tenth;
the amplitudes of the cores-only round did, by scaling the arms up by
their 0.86 to 0.94 and the central square down. That round's model lacks
the 30 mJy this source has beyond 8 arcsec, which the baselines of the
central square see: its antennas come out above one, all thirteen of
them, and the normalisation puts the arms below. The same thing was seen
in the whole field in round 5. There the next round, against a model of
positive components, pushed the other way and seemed to take it back.
Against a usual clean it does not: a usual clean reproduces whatever the
data are, a tilt between short and long baselines looks to it like the
shape of the source, and the solve returns one (0.99 to 1.02 here). The
bias stays in the data. So the scatter among the arm antennas above is
partly real and partly this, and the 55 mJy/beam is too high.

`selfcal_recipe.py` no longer has the cores-only amplitude round.
Amplitudes are solved only against the whole field, on every baseline, up
to three rounds until they are within 10% of one. What differs from
antenna to antenna is found that way; how short and long baselines stand
to each other is left as the calibrators set it, which no model made from
the target can improve. `--time-solint` adds a last round at a shorter
interval within each scan, for the hole that is left: a pointing that
drifts through 22 minutes is the next thing an off-axis source would
show. 947 tests pass.

**The same source with the recipe as changed, 19:47-20:10 UT**
(`recipe_source1b`, with `--time-solint 2min`). Holes within 80 arcsec as a
percentage of the peak, and the number of pixels there below
-100 µJy/beam:

| image | peak (mJy/beam) | pixels below | 19 arcsec S | 7 arcsec | 17 arcsec NE | the next ones |
|---|---|---|---|---|---|---|
| before | 50.0 | 1010 | -1.52 | -1.26 | -1.31 | -0.90, -0.89, -0.72 |
| after phases and delays | 51.7 | 968 | -0.77 | -0.56 | -1.58 | -0.75, -0.72, -0.68 |
| after amplitudes (0.97 to 1.05) | 51.8 | 836 | -0.71 | | -1.56 | -0.70, -0.65, -0.63 |
| after the gains across the band | 52.2 | 684 | -0.50 | -0.40 | -1.41 | -0.62, -0.54, -0.51 |
| after the round within the scan | 52.5 | 604 | -0.52 | -0.49 | -1.32 | -0.61, -0.52, -0.48 |

Without the biased round the peak rises by 5%, not 11, and that is to be
believed: it comes in steps with the phases, the delays and the bandpass.
The amplitudes towards the source, solved against its whole image, are
within 5% of one. Two of the three deep holes go to a third, and the
negative area by 40%. Within the scan nothing moves: amplitudes by 2% rms
at most, phases by under a degree, and the round for it changes little.

The hole 17 arcsec to the north-east does not respond to anything: not to
the field's calibration, not to gains towards the source in phase, delay,
amplitude, frequency or time. In the sub-band images it is the hole of the
lowest three sub-bands, where the source is brightest and highest on the
beam and where pointing matters least. It is not an antenna's error. It
has not been explained.

**What taking the source out leaves behind.** The residual images of the
two cleans of the isolated source say what a subtraction would leave in
the field: with the field's calibration (the "before" image) and with the
source's own (the final one).

| within 100 arcsec of the source | field's calibration | the source's own |
|---|---|---|
| residual rms inside 30 arcsec, over the noise far away | 1.23 | 1.10 |
| the same between 30 and 100 arcsec | 1.26 | 1.19 |
| extremes (µJy/beam) | -149, +155 | -133, +123 |
| pixels beyond 3 sigma, where noise alone gives 85 | 468 | 315 |

Either way the holes go with the source: a usual clean's model holds
them, and subtracting the model subtracts them. The source's own
calibration makes what is left smaller, from a quarter above the noise to
a tenth near the source. So the field can be rid of the holes beside its
bright off-axis sources by taking those sources out, each with the gains
the antennas had towards it; what it cannot have is those sources in the
image and no holes beside them, one hole in each having resisted
everything.

### 126. Taking a calibrated source out of the field: `--peel`

`experiment_extract_source.py --peel <recipe folder>` is the way back.
From the recipe's report it takes the source's calibrated data and the
model of its final image. WSClean predicts that model into a copy of
those data; CASA subtracts it and moves the phase centre back to the
field's; WSClean predicts the model of the rest of the field (the discs
of this source and of any already taken out left blank) into the result,
and CASA adds it (`uvsub(reverse=True)`). What comes out is the field's
data without the source.

The calibration towards the source is not undone on what is left. What is
left of the isolated data once the source is subtracted is noise and the
little the field's model did not describe; gains within a few per cent of
one applied to that change nothing that can be seen. The rest of the
field comes back as its model, which is in the field's own calibration.

`--also-without X Y` names sources already taken out, for extracting the
next one from the peeled data and for peeling it: the rest-of-the-field
model then leaves their discs blank too. 957 tests pass. Not yet run on
data: the predict of a small-field model into the shifted data, the shift
back, and `uvsub` in reverse.

**On data, 2026-10-10 09:45-10:46 UT.** One run. The peel took under three
minutes and each step worked the first time; the field without the source
(`extract_test/source1/source1_peeled.ms`) was then imaged as before
(`wsclean_test_scan15_peel1`, twelve sub-bands, usual clean).

| | with the source | without it |
|---|---|---|
| image robust rms (µJy/beam) | 26.7 | 26.3 |
| quietest tenth; noise of the data | 23.2; 22.5 | 23.2; 22.4 |
| brightest pixel (mJy/beam) | 50.9, the source | 39.2, the target |
| most negative pixel (mJy/beam) | -0.83, the source's hole | -0.53, another source's |
| the box that held the source | 32.3, joint noisiest of the image | no longer among the five noisiest (below 30.7) |

Where the source stood, within 100 arcsec: the brightest pixel went from
50.9 mJy/beam to 0.18, the deepest from -0.83 to -0.20, and the pixels
below -75 µJy/beam within an arcminute from 1001 to 186. The patch is
not at the noise: its rms is 36 µJy/beam where the patch opposite it
across the pointing centre has 24 to 28, a third to a half above. Elsewhere
nothing moved: the target is at 39.16 mJy/beam (39.05) with the same
-0.89% hole in the same place, the other five sources keep their peaks to
0.03 mJy/beam and their holes, and the opposite patch is the same to
1 µJy/beam. So the model of the rest, subtracted and added back, returns
the rest as it was.

### 127. The remaining sources in one run: `peel_sources.py`

Three commands per source by hand, for five more sources, is the
situation of finding 120 again. `peel_sources.py` runs them: for each
source in the order given it extracts, calibrates towards it
(`selfcal_recipe.py`, 1024 pixels, an image first), and peels, the next
source being taken from what the last one left with `--also-without`
naming all that are gone; then it images the field. Sources taken out in
an earlier run are named with `--already-out`. It holds its folder,
records each source's peak and deepest pixel before and after from the
recipe's report, and continues with the first source not yet out when the
same command is given again.

Each source has its own two numbers. The threshold of its cores-only model
is about an eighth of its peak, as 6 mJy/beam was for the first. The
channels per gain of its bandpass round grow as it gets fainter (60 for
the 33 and 24 mJy sources, 120 for the 19 and 18, 240 for the 11): the
bandpass round flags what it cannot solve, and a faint source must not
cost the field data because one of its bins was too noisy to solve.

969 tests pass, and the sixteen commands the driver builds for the five
sources were given to the three scripts' own argument parsers.

**On data, 2026-10-10 13:12-14:39 UT (`peel_field`): four sources out, a
failure on the fifth, and a fault of mine that the logs showed.** One run.
The loop works, at 17 to 20 minutes a source.

| source (off-axis) | peak before, after (mJy/beam) | deepest pixel before, after | towards it |
|---|---|---|---|
| 33 mJy (17 arcmin) | 33.2, 33.6 | -0.47, -0.43 | amplitudes 0.91 to 1.09; across the band 8 deg and 7% |
| 24 mJy (16 arcmin) | 24.2, 24.1 | -0.33, -0.25 | 0.94 to 1.09; 10 deg and 8% |
| 19 mJy (27 arcmin) | 19.0, 19.6 | -0.25, -0.27 | delay 0.98 ns on a central antenna; 0.94 to 1.18; 34 deg and 22% |
| 18 mJy (37 arcmin) | 18.7, 18.9 | -0.42, -0.27 | 0.80 to 1.12; 20 deg and 29% |
| 11 mJy (16.5 arcmin) | 11.6, 11.8 | | delay 17.9 ns on C01; amplitudes across the band by 53%; then the peel failed |

The peel of the last stopped in WSClean: "the input image contains
non-finite values". Its final model had no data in the top two sub-bands.
Counting the visibilities WSClean gridded in each sub-band, at the start
of each source's recipe, says why:

| sub-band | at the start | after the 33 mJy | 24 mJy | 19 mJy | 18 mJy source |
|---|---|---|---|---|---|
| 9 (789 MHz) | 2.57 M | 2.57 | 2.51 | 2.42 | 2.42 |
| 10 (812 MHz) | 2.41 M | 2.41 | 2.30 | 0.94 | 0.65 |
| 11 (835 MHz) | 1.39 M | 1.12 | 1.04 | 0.28 | 0.10 |

Each bandpass round flagged the channels it could not solve, as finding
123 had it do for the whole field, where an antenna without a solution
has nothing usable there. Towards one faint source a bin without a
solution says only that the source was too faint in that bin; but the
data are the field's, and the flags went back into the field with what
was left. Five sources between them flagged the top of the band away.
Wider bins for fainter sources, which was my precaution, were not enough.
And the numbers of the last three sources are not gains: a delay of 18 ns
on an antenna of the central square and amplitudes running by a half
are a solve with more unknowns than the source can carry.

Three changes. `experiment_phase_offsets.py --keep-unsolved` leaves what a
bandpass solve could not solve unflagged and uncorrected; the recipe
passes it on, and `peel_sources.py` always gives it. `selfcal_recipe.py
--phase-amplitude-only` keeps only the cores-only phase rounds and the
whole-field amplitude rounds, a phase and an amplitude per antenna and
scan, which an 11 mJy source determines to under a per cent. And
`peel_sources.py` takes `BINS` of 0 to mean that. The field data of this
run (`peel_field`) carry the flags and are not to be used; the run is to
be made again from the data without the first source, in another folder,
with the 33 and 24 mJy sources in full and the three fainter ones with
phases and amplitudes only. 974 tests pass.

**The run made again, 2026-10-10 14:49-16:54 UT (`peel_field2`).** One
run, no failure. Each sub-band holds at the end what it held at the start
(2.41 and 1.39 million visibilities in the top two, after every source).
The gains towards the three faint sources are gains now: amplitudes of
0.94 to 1.18, 0.79 to 1.12 and 0.80 to 1.07 in the first round, within 4%
of one in the second.

The field without its six bright off-axis sources
(`peel_field2/source6/source6_peeled.ms`, image in `peel_field2/image`):

| | all sources in (finding 123) | six sources out |
|---|---|---|
| image robust rms (µJy/beam) | 26.7 | 25.0 |
| quietest tenth; noise of the data | 23.2; 22.5 | 22.5; 21.7 |
| median box; noisiest box | 26.9; 32.3 | 25.1; 29.0 |
| brightest pixel (mJy/beam) | 50.9, an off-axis source | 38.5, the target |
| most negative pixel (mJy/beam) | -0.83 | -0.35, the hole beside the target |

Where each source stood, within 100 arcsec, and the patch opposite it
across the pointing centre for comparison:

| source that stood there | brightest, deepest pixel (µJy/beam) | rms inside 30 arcsec | pixels below -75 µJy/beam | the opposite patch: rms, pixels |
|---|---|---|---|---|
| 51 mJy | 172, -149 | 34 | 150 | 24, 30 |
| 33 mJy | 178, -159 | 32 | 123 | 25, 8 |
| 24 mJy | 120, -125 | 33 | 149 | 26, 48 |
| 19 mJy | 125, -115 | 30 | 87 | 25, 16 |
| 18 mJy | 138, -117 | 33 | 85 | 23, 4 |
| 11 mJy | 110, -101 | 30 | 53 | 27, 85 |

No source and no hole is left at any of the six: holes of 0.2 to 0.8
mJy/beam have become patches a fifth to two fifths above the noise, with
nothing beyond 0.18 mJy/beam. The noise of the data fell from 22.4 to
21.7 µJy/beam and the beam from 4.76 by 3.67 to 4.64 by 3.55 arcsec; both
are the weights, which each source's amplitude corrections rescaled in
what was left, and not a gain in the data. The target's peak is 38.5
mJy/beam where it was 39.1, which is that smaller beam on a resolved
source; its flux within an arcminute is 105 mJy.

What is left of the old kind is beside the target: -0.92% of its peak
45 arcsec to the south-west, where it has been in every image, and five
more at half a per cent. The target was calibrated with the field, and
the field's one gain per antenna was a compromise between it and six
brighter directions that are now gone.

`selfcal_recipe.py --final-imsize` lets the final image be larger than
the models of the rounds, for the next step: the recipe on these data,
which the target now dominates, with models of the inner 51 arcmin and
the image of the whole field. 976 tests pass.

### 128. One delay of 17.85 ns: a solve is no longer applied whatever it says

**On data, 2026-10-10 17:24 UT (`recipe_peeled`), stopped as useless.**
The recipe on the peeled field found phases of up to 12 degrees towards
the target and, in the delay round, delays within half a nanosecond for
every antenna and polarization but one: C10 in one polarization, 17.85 ns
(0.32 ns in its other). The same number, 17.9 ns, had come out for C01
towards the 11 mJy source. A delay is found on the baseline to the
reference antenna alone; C00 to C10 is a few hundred metres, the model of
the cores there is one source of 108 mJy, and the search settled on a
ripple of the data. Applied, it turns that polarization of C10 five
times round across the band. The amplitude rounds then gave C10 0.27 and
1.31, then 1.54 and 0.48, while every other antenna stayed within 5% of
one; the image went from 25.0 to 43.5 µJy/beam in the small field's
terms, and its holes from -0.9% to -3.3%.

Nothing in the scripts looked at a solution before applying it. Three
limits are now in `experiment_phase_offsets.py`:

- a delay beyond `--max-delay-ns` (5; the largest real one in this
  observation was 2.3) is set to zero in the table, so the antenna keeps
  the delay it has and still gets its phase (`wild_delays`);
- in an amplitude solve, a solution outside `--amplitude-limits` (0.4 to
  2.5) is flagged in the table and so not applied;
- in a bandpass solve, a bin whose amplitude over what the antennas share
  is outside the same limits is flagged in the table.

What is rejected is printed, kept in the solve's report
(`rejected_delays`, `rejected_amplitudes`, `rejected_bins`), and carried
into the recipe's reading as "not believed and not applied". The limits
are wide: W04's real 2.2 across the band (finding 123) passes, and so do
the 0.68 to 1.35 of the first amplitude round on the whole field.
982 tests pass.

**On data, 2026-10-10 20:33-22:29 SAST (`recipe_peeled2`).** The fixed
scripts reached the cluster four minutes into a run started without them;
it was stopped and started again in the same folder, and went on from the
phase round it had finished. In the delay round: "Delays beyond 5 ns are
not believed and are not applied: C10 (+17.85 ns)". Every other delay is
within 0.53 ns. The amplitudes that followed are 0.98 to 1.05; the
bandpass round rejected nothing (phases to 10.7 degrees across the band,
mostly a slope; amplitudes to 7.8%), and what the antennas share across
the band is now 0.6% where it was 2.4% with six sub-bands.

The field without its six bright off-axis sources, calibrated again on
what is left (`recipe_peeled2/G0212-0335/round5.ms`, image in
`images/final`):

| | before this round (finding 127) | after |
|---|---|---|
| image robust rms; quietest tenth; data (µJy/beam) | 25.0; 22.5; 21.7 | 25.2; 22.6; 21.5 |
| the target's peak (mJy/beam); flux within an arcminute (mJy) | 38.5; 105 | 39.0; 102 |
| the hole 45 arcsec south-west of the target | -0.92% (-0.35 mJy/beam) | not among its seven deepest (shallower than -0.17) |
| deepest pixel within 80 arcsec of the target | -0.35 mJy/beam | -0.26, 7 arcsec from the core |
| pixels below -75 µJy/beam within an arcminute of the target | 582 | 488 |
| the six patches where sources stood | rms 30 to 34 | 29 to 35, as before |

The hole that had stood beside the target in every image is gone, as the
one of the upper band beside the first off-axis source went when that
source had its own gains. The image as a whole is as it was.

**What the image still holds.** Its most negative pixel is -0.375
mJy/beam, 7 arcsec from the 11.6 mJy source that lies 103 arcsec south of
the target; then -0.35 beside an 8 mJy source 38 arcmin out, -0.32 beside
a 5.8 mJy source 39 arcmin out, the -0.26 beside the target, -0.23 beside
a 9.8 mJy source, and one of -0.23 with nothing bright near it. 2021
pixels of 26 million are below -5 sigma where noise would put 8, and
11 182 below -4 sigma where noise would put 831. The remaining sources of
5 to 12 mJy/beam each keep a hole of 3 to 4% of their peak, a larger
fraction than the bright ones had; in mJy/beam the deepest is now 0.375
where it was 1.10 three days ago and 0.83 before the sources were taken
out.

**Where scan 15 of G0212-0335 stands.**

| | pipeline (all scans, tclean) | WSClean, its calibration | recipe, all sources | six sources out, calibrated on the rest |
|---|---|---|---|---|
| image robust rms (µJy/beam) | 80 | 29.0 | 26.7 | 25.2 |
| quietest tenth | 60 | 23.5 | 23.2 | 22.6 |
| noise of the data | about 23.5 | 23.4 | 22.5 | 21.5 |
| deepest negative (mJy/beam) | | about -2.4 (10% of 23.6) | -0.83 | -0.375 |

### 129. All the targets in one command: `run_targets.py`

Sushant, 2026-10-10: "Run the entire pipeline (without subbands) using
WSClean, peeling, DD calibration for all the targets in one go."

`run_targets.py` takes each target the pipeline calibrated (the nine in
`uGMRT_Products/calibration/split_targets.json`) through the sequence of
findings 118 to 128, one after another: the recipe on the whole field,
peeling, the recipe on what is left, and the primary-beam correction. It
runs the other scripts and needs neither CASA nor WSClean itself. Each
step keeps its own record, so the same command given again goes on where
it stopped. Until one target has finished, a failure stops the run, since
what stopped the first would as likely stop the next hours in; after
that, a target that fails is left and the others go on.

"Without subbands" is read as: no sub-band images as products. WSClean
still deconvolves in twelve frequency sub-bands, without which a band
this wide cannot be imaged (finding 123); the product of a target is one
image of the whole band, corrected for the primary beam.

What was done by hand on scan 15 of G0212-0335 and had to become a rule:

- **Where the data start.** The calibrated data of the pipeline, all
  scans, which have had no self-calibration; the hand-run started from the
  pipeline's self-calibrated copy. One phase per scan cannot follow a
  phase that drifts within a scan, so after the rounds per scan the recipe
  now solves phases every minute against the cores-only model
  (`--time-phase-solint`, stage `phase_time`, up to three rounds).
- **How deep the cores-only models go.** Not 6 mJy/beam, which was
  G0212-0335's, but a fraction of each field's own brightest pixel
  (`--core-fraction 0.13`, which is what 6 was of 46.8). WSClean ends its
  first major iteration when the peak has fallen to (1 - mgain) of where
  it began, so one major iteration with mgain 0.87 stops there without
  the peak having to be known (`experiment_wsclean_image.py
  --core-model-fraction`).
- **When a round has settled, with several scans.** The solve's report
  now gives each antenna's largest single solution beside its mean, and
  the recipe judges by that: +80 degrees in one scan and -80 in the next
  have a mean of nothing.
- **Which sources are peeled** (`peel_sources.py --auto-sources`): those of
  at least 10 mJy/beam as the telescope saw them, more than 4 arcmin from
  the pointing centre, at most eight, brightest first, none within two
  discs of a brighter one. For G0212-0335 that is the six that were
  peeled by hand, and it leaves the 11.6 mJy source 103 arcsec from the
  target where it is. Each gets a cores threshold of an eighth of its
  peak, and a bandpass round in 30 channels from 40 mJy/beam, in 60 from
  20, and none (nor delays) below.
- **The primary beam** (`correct_primary_beam.py`): the image is divided
  by the mean of its twelve sub-bands' beams, weighted as WSClean weighted
  the sub-bands (`WSCIMGWG`), each beam the NCRA polynomial the pipeline
  already holds for Band 4 at that sub-band's frequency. It is blank
  beyond the radius at which the narrowest beam has fallen to 0.3 and
  wherever the mean beam is below 0.3. The image itself is not changed.
- **Disk.** A step at 5120 pixels and twelve sub-bands writes 6 GB of
  images, and each peeled source leaves copies of the data. `--tidy` in
  the recipe removes the sub-band images of finished rounds, and of the
  final image all but the sub-bands' images and models that peeling and
  the beam correction read; in `peel_sources.py` it removes a source's own
  data, the copies its calibration worked on, and the field as the source
  before it left it. About 15 GB a target is left.

1028 tests pass, and the commands the driver builds for a target were
given to the three scripts' own argument parsers. **Not run on data:**
the start from data without self-calibration and the phases every minute;
any solve with more than one scan; the cores-only model by fraction; the
choice of sources; the tidying; the beam correction; and the driver
itself. G0212-0335 is first in the pipeline's record, so the first target
through is the one whose scan 15 is known: 25 µJy/beam with six sources
out, a target of 39 mJy/beam. With all three scans the data's noise
should be a little lower than scan 15's 21.5 µJy/beam, unless the first
scan, whose own noise was 1.66 times its share, costs more than it gives.

At the pace of scan 15 a target takes about fourteen hours: eight or nine
for the field (six or seven cores-only models, three or four images of
the whole field), two or three to peel six sources, three or four for
what is left. Nine targets are five to six days.

### 130. Two scans that cost more than they add: `select_scans.py`

**On data, 2026-10-10 23:22 to 2026-10-11 07:47 SAST: the first run of
`run_targets.py`, G0212-0335 with all three of its scans.** It was
stopped in the field's third amplitude round, eight and a half hours in,
because of where it was going. In µJy/beam:

| | scan 15 alone (findings 120 to 125) | all three scans |
|---|---|---|
| image of the whole field, phases only | 29.0 (23.5 in the quietest tenth) | 41.0 (30.4) |
| after one amplitude round per scan | 27.9 (23.9) | 39.8 (29.6) |
| Stokes V | 23.4 | 22.9, then 21.9 |
| negative beside the sources | -2% of their peak | -7%, then -6% |
| an image of the whole field takes | about 1 h | 3 h 50 min, then 3 h 03 min |

The two other scans lower the data's noise by 4% and raise the image's by
40%. Scan 15 by itself, taken through to the end, is at 25.2.

What the solves show of the two:

- **Scan 29** (22:08 UT, 1.3 min, what flagging left of the last scan).
  In every solve 25 to 32 of its 56 solutions fail at a signal-to-noise
  of 3, and the rest are near it: they fit noise. Its delays are -1187 to
  +1190 ns for every antenna, all refused by the 5 ns limit of finding
  128; its amplitudes for seven antennas are outside 0.4 to 2.5 and
  refused; its two one-minute phase solutions are what kept the "largest
  phase" of the phases-every-minute rounds at 173, 94 and 84 degrees, so
  those rounds never settled and all three were run.
- **Scan 3** (14:12 UT, 9.3 min, cut short: the next scan starts eight
  minutes after it ends, against two and a half between whole scans).
  Its solutions do not fail, 2 of 52. But its delays are not scan 15's:
  C09 -5.0 ns, which is 502 degrees across the band, S04 +3.1, S06 +1.2,
  S02 +1.0, W04 -1.0 in one hand only; scan 15 has 2.3 ns on S04 and
  nothing else above 0.5. Its amplitudes per antenna and hand differ from
  scan 15's by up to 60%, and a second round still finds 0.84 to 1.16
  where scan 15 alone gave 0.93 to 1.04. Finding 117 had already measured
  it: 1.66 times as noisy as its length allows in Stokes V, and a
  residual 3.3 times its own noise against 1.6 for scan 15.

This was known and was run anyway. Finding 117 says that leaving scans 3
and 29 out costs little, and finding 129 names scan 3 as the risk of the
first target. The recipe was worked out on scan 15 and nothing in it
repairs a scan that is noisy rather than miscalibrated. The cost was a
night of the node.

**What the other targets have** (`inspection_report.json`,
`time_coverage`): G0213-0251 one scan of 22 min; the seven others two
scans of 22 min each, four and a half hours apart. G0212-0335 is the only
target with pieces of scans. So the other targets need both their scans,
half their data each, and what has to be decided by rule is whether a
scan is a whole one.

**The measure.** `select_scans.py` reads a target's calibrated data, one
scan at a time, and sums over the cross-correlations and the channels
that hold data the difference of the two parallel hands, (RR - LL) / 2,
which holds no sky to speak of, with the weights the data carry:

- the scatter of one visibility, `rms_jy`: how noisy the scan is for the
  amount of data it has;
- the noise the scan alone gives an image with natural weights,
  `noise_jy` = sqrt(sum w² |V|² / 2) / sum w, and from it the scan's
  share of the target's weight;
- the plain scatter over a robust one (1.4826 times the median distance
  from zero): above about 1.2 there is interference left in.

It also gives the noise of the scans together as an imager adds them, by
the weights they carry and not by their noise, for every scan and for the
kept ones.

**The rules.** A scan is left out if it holds less than 10% of the weight
(`--least-share`): too little for solutions of its own, and the image
loses at most 5% without it. Or if it is more than 1.4 times as noisy,
for the same amount of data, as the scan with the largest share
(`--most-excess`): it carries more than receiver noise. The scan with the
largest share is always kept. `--scans` names the scans to use and
overrules both; what the rules say of each scan is still recorded. If
scans are left out, CASA splits the kept ones into
`<workdir>/<target>.ms`; if none are, the data given are what goes on.
The data given are opened for reading only.

The limits rest on one target. From the images of finding 117 (61.0,
23.9 and 110.9 µJy/beam for 9.3, 22.2 and 1.3 min) scan 3 should come
out 1.65 times as noisy as scan 15 with 13% of the weight, and scan 29
1.1 times with 4%: both left out, each by one rule. **That the scatter of
the visibilities reproduces those figures from the images has not been
seen**, so the first run names G0212-0335's scan outright
(`--scans G0212-0335=15`) and what the rules say of it is the test of the
measure.

**In `run_targets.py`.** The scans are a stage of their own and it is
done for every target before the first image is started: a few minutes
each, and the table of all the targets' scans is printed and in
`targets_report.json` (`targets.<T>.scans`) half an hour into a run of
days. The field's recipe is given the kept scans. A work folder whose
field was made from other scans is not gone on with and not removed: the
target stops with "the scans to use have changed". `--scans TARGET=3,15`,
`--least-share` and `--most-excess` are passed on. The summary line of a
target says which scans it is without.

1080 tests pass: the sums against arrays worked through the long way on a
stand-in for `casatools.table` (flags, row flags, autocorrelations, four
correlations, weights per row, reading in pieces), the rules on
G0212-0335's own figures, and the driver's flow. **Not run on data:**
`select_scans.py` itself, and still any solve on two whole scans.

**On data, 2026-10-11 08:09 to 08:15 SAST (`targets_run2/`, one run,
`--scans G0212-0335=15`).** The nine targets' scans were measured in six
minutes, 26 to 47 s a target, and the split of scan 15 is in that. Per
scan: minutes with data, antennas, flagged part of the channels in use,
the weighted scatter of one visibility of (RR - LL) / 2, the scatter of
its bulk (unweighted, robust), the plain scatter over the bulk's, the
noise the scan alone gives an image with natural weights, its share, and
its noise against the main scan's for the same amount of data.

| target | scan | min | ant. | flagged | scatter (mJy) | bulk (mJy) | plain/bulk | alone (µJy/beam) | share | excess |
|---|---|---|---|---|---|---|---|---|---|---|
| G0212-0335 | 3 | 9.5 | 26 | 31% | 131.8 | 452.1 | 3.35 | 91.2 | 3% | 1.79 |
| | 15 | 22.4 | 27 | 22% | 73.5 | 91.6 | 1.34 | 16.7 | 90% | 1.00 |
| | 29 | 1.4 | 28 | 19% | 60.7 | 93.9 | 1.71 | 61.4 | 7% | 0.83 |
| G0213-0251 | 17 | 22.2 | 28 | 15% | 131.3 | 123.9 | 1.31 | 28.2 | 100% | 1.00 |
| G0216-0604 | 5 | 22.2 | 27 | 23% | 119.8 | 147.4 | 1.42 | 27.9 | 37% | 1.19 |
| | 18 | 22.2 | 28 | 16% | 100.3 | 112.8 | 1.18 | 21.3 | 63% | 1.00 |
| G0218-0245 | 6 | 22.2 | 27 | 22% | 170.8 | 186.8 | 1.20 | 38.8 | 39% | 1.15 |
| | 20 | 22.4 | 28 | 15% | 148.1 | 155.5 | 1.11 | 31.0 | 61% | 1.00 |
| G0220-0318 | 8 | 22.2 | 27 | 23% | 77.8 | 98.6 | 1.51 | 17.7 | 41% | 1.12 |
| | 21 | 22.2 | 28 | 16% | 69.7 | 83.9 | 1.30 | 14.8 | 59% | 1.00 |
| G0222-0321 | 9 | 22.2 | 27 | 24% | 65.7 | 92.9 | 1.57 | 15.8 | 45% | 1.00 |
| | 23 | 22.2 | 28 | 16% | 65.5 | 84.3 | 1.34 | 14.3 | 55% | 1.00 |
| G0224-0620 | 11 | 22.6 | 27 | 23% | 77.5 | 97.7 | 1.35 | 18.0 | 53% | 1.00 |
| | 24 | 22.2 | 28 | 16% | 84.5 | 106.8 | 1.40 | 19.1 | 47% | 1.09 |
| G0229-0556 | 12 | 22.6 | 27 | 23% | 94.5 | 105.5 | 1.36 | 21.5 | 52% | 1.00 |
| | 26 | 22.2 | 28 | 16% | 99.4 | 112.8 | 1.40 | 22.2 | 48% | 1.05 |
| G0231-0533 | 14 | 22.2 | 27 | 23% | 58.4 | 79.3 | 1.39 | 13.4 | 47% | 1.01 |
| | 27 | 22.4 | 28 | 16% | 57.8 | 78.6 | 1.37 | 12.7 | 53% | 1.00 |

- **The seven targets with two scans keep both.** Their scans are alike:
  at most 1.19 times as noisy as each other, shares of 37 to 63%.
- **G0212-0335: the rules agree with what was named**, scan 15 only, but
  not for the reason predicted. Scan 3 was expected at 13% of the weight
  and 1.65 times as noisy. It holds 3%, and is left out by that rule
  before the other is asked; it is 1.79 times as noisy weighted, and the
  bulk of its visibilities scatters five times as much as scan 15's (452
  against 92 mJy) under weights twelve times lower (22 against 268). Why
  the images of finding 117 gave it 13% (61.0 against 23.9 µJy/beam) is
  not known: those were of the pipeline's self-calibrated data, at robust
  0, and their rms a robust one, any of which may be the difference.
  Scan 29 holds 7%, expected 4%. Without the two the data allow 16.7
  µJy/beam with natural weights against 15.9 with them, 5% more.
- **Wrong in the script as uploaded:** every scan whose plain scatter is
  over 1.2 times its bulk's is said to have "interference left in". All
  eighteen scans are over 1.1 and whole ones reach 1.57: channels and
  antennas differ in noise, and that alone spreads the distribution. Only
  scan 3's 3.35 says anything. The words are now given only to a scan
  1.5 times above the main one. Nothing is decided by them; the copy on
  the cluster still prints them.
- **Not receiver noise.** The noise of a whole scan alone runs from 12.7
  (G0231-0533) to 38.8 µJy/beam (G0218-0245), a factor of three for the
  same 22 minutes, antennas and flags, and the bulk scatter from 79 to
  187 mJy. It goes with the field and not with the time: both scans of a
  field, four and a half hours apart, are alike, and scans next to each
  other on different fields are not (17, 18, 20, 21 at 18:27 to 19:46 UT:
  124, 113, 156, 84 mJy). Receivers do not know which field they look
  at. What does is the sky: the two hands' beams differ away from the
  centre (in the dirty Stokes V of G0212-0335 the source at 29.5 arcmin
  is a tenth of its Stokes I, finding 117), so (RR - LL) / 2 holds a part
  of every bright off-axis source. **Reading, not tested:** G0218-0245,
  G0213-0251 and G0216-0604 have much brighter sources on the flank of
  the beam or beyond it than G0212-0335 has; their first images will
  say. Two consequences. The comparison of a target's scans stands, since
  both see the same sky, though a sky that dominates the scatter would
  hide a noisy scan (G0218-0245's 1.15 is of scatter that is mostly not
  noise). And the figure "the data allow" is above the receivers' noise
  in such fields: the receivers' own is at most what G0231-0533 has, 58
  mJy per visibility, which for two whole scans is 9 µJy/beam with
  natural weights and about 12 at robust 0 if the ratio of G0212-0335's
  scan 15 holds (21.5 in its Stokes V image against 16.7 here).

1081 tests pass.

### 131. A hole beside every source: the amplitude round that was taken out

**On data, 2026-10-11 08:15 to 10:55 SAST (`targets_run2/`, G0212-0335,
scan 15, the field's recipe from the pipeline's calibrated data).**

| round | found |
|---|---|
| `phase1` to `phase3` | largest phase 161 (C01), 23, 6 degrees (S02) |
| `phase_time1`, `2` (every minute) | 12, 4 degrees; within the scan phases moved by 5.6 degrees rms, then 0.4 |
| `delay` | 2.34 ns (S04); none refused |
| `amplitude1`: image of the whole field | 23.8 µJy/beam (20.1 in the quietest tenth) against 20.8 for the data; peak 45.8 mJy/beam; 63 min |
| `amplitude1`: solve | 0.89 to 1.11, almost all of it between the two hands of an antenna (E03 0.89/1.11, C12 0.95/1.11, S04 0.94/1.10); the mean of the two hands within 3% of one for every antenna but W04 (0.94) |
| `amplitude2`: image | 23.7 (20.1) against 20.0; peak 45.6 |
| `amplitude2`: solve | 0.95 to 1.03 |

The noise is the best this field has had: 23.8 where the same stage gave
29.0 before, on data whose own noise is 20.8 and not 23.4, with 27
antennas solved and not 25. Those earlier data were the pipeline's
self-calibrated ones (`experiment_wsclean_image.py` takes them when no
`--ms` is given, and the scan-15 file all of findings 118 to 128 worked
on was split from them); these are its calibrated ones, before its own
self-calibration, which flagged what it could not solve.

**But every source has a hole beside it.** In the `amplitude2` image the
deepest negative is -4.01 mJy/beam, six pixels south of the brightest
source (45.6 mJy/beam): -8.8% of it. The same at the others, peak and
hole in mJy/beam and the hole's place:

| source (pixel) | peak | hole | of the peak | where |
|---|---|---|---|---|
| 3162, 4227 | 45.6 | -4.01 | -8.8% | 6 px south |
| 2562, 2567 (the target) | 34.7 | -2.60 | -7.5% | 6 px north |
| 1641, 2103 | 29.2 | -1.77 | -6.0% | 6 px north |
| 3227, 1866 | 21.5 | -1.67 | -7.8% | 6 px north |
| 4148, 2906 | 16.3 | -1.17 | -7.2% | 6 px north |
| 2554, 2464 | 9.8 | -0.84 | -8.6% | 6 px south |
| 1820, 3216 | 10.1 | -0.82 | -8.1% | 5 px south |
| 1254, 3342 | 4.4 | -0.42 | -9.5% | 6 px south |

The same part of the peak for a source of 4 mJy/beam and one of 46, on
the axis and half a degree off it, so one error for the whole field and
not a direction's. Always north or south, at a beam and a half. Negative
on both sides where a source is compact (-6.0 and -3.7% at 1641, 2103;
-7.8 and -4.0% at 3227, 1866), and on the side away from its extension
where it is not, the other being filled by the source: an even pattern,
which is amplitude. The residual holds nothing of it (within +-0.15
mJy/beam): the holes are negative components of the model, which sums to
1.32 Jy positive and -0.52 negative. The amplitude round did not move
them (-7.8 to -7.6%, -4.010 to -4.009 mJy/beam).

This is the hourglass of finding 119, at twice the depth. There, on the
pipeline's self-calibrated data, the negatives were at -3.9% after the
phases and delays, and the solve against a cores-only model took a third
of that away at once (amplitudes 0.68 to 1.35). Finding 125 took that
round out of the recipe, for the trend between short and long baselines
it leaves, and wrote down why the round that stayed cannot do its work:
"a usual clean reproduces whatever the data are ... and the solve returns
one". That was true and was the reason to keep a round that does not. On
data that had been through the cores round already it did not show. The
recipe as it stood had never been run from data without amplitude
self-calibration (finding 129 says so), and the first run from such data
shows it.

**The round is back, as an option.** `selfcal_recipe.py
--amplitude-cores` puts, after the delays and before the whole field's
amplitude rounds, `amplitude_cores`: amplitude and phase per antenna and
hand against a cores-only model of positive components, on baselines
beyond 3 kλ (`--amplitude-cores-uvrange`), up to three times until the
amplitudes found are within 10% of one. Its cores go twice as deep as the
phases' (half the threshold, or half the fraction of the field's peak):
by hand they were 3 and 6 mJy/beam. `run_targets.py` asks for it in the
field's recipe, whose data are the pipeline's calibrated ones, and not in
the recipe on what is left after peeling, nor do the solves towards
single sources have it: those data have been through it.

What is not settled by this. The trend of finding 119 (central square
above one, outer arms below, the sign a model without extended emission
gives) comes back with the round, and nothing after it takes it out.
Whether a third of the holes goes, as before, or more, since more of the
amplitude error is still in these data, is not known; the first whole
field's image after the round says, 1 h 40 min into a run. If the holes
stay near -7% the cause is not the antennas' amplitudes.

1085 tests pass.

**On data, 2026-10-11 11:23 to 13:12 SAST (`targets_run3/`, one run, the
same data and command with the round in).** Phases and delays came out as
before to the digit. Then:

| `amplitude_cores1` (R/L) | |
|---|---|
| S02 | 1.59 / 1.76 |
| S04 | 1.25 / 1.47 |
| W04 | 0.71 / 0.64 |
| S06 | 0.73 / 0.83 |
| central square | 0.91 to 1.20 |
| east arm | 0.86 to 1.08 |
| W01, W02, W03, W05 | 1.11, 1.02/1.07, 0.91/0.96, 0.89/0.96 |

The same antennas in the same sense as on 8 October (S02 1.21/1.35, S04
1.10/1.29, W04 0.75/0.68, S06 0.74/0.84), S02 and S04 further out. They
were put down then to the pipeline's amplitude self-calibration having
raised antennas that delays and phase offsets had decorrelated. These
data have not had that self-calibration and are worse, so it is the
calibrators' amplitudes that did it, and the pipeline's self-calibration
had taken part of it back. Two antennas of the south arm at one and a
half times their amplitude is what puts a hole north and south of every
source.

`amplitude_cores2` found 0.91 to 1.03 and the round ended. The two hands
agree to the third decimal in it: what is left is the model's, not the
antennas'. It is the trend: the central square at 1.00, the outermost
antennas below (S06 0.91, W05 0.92, W04 0.95, E06 0.97).

The whole field's image after it, against the one without the round:

| | without (`targets_run2`) | with |
|---|---|---|
| image (µJy/beam), quietest tenth | 23.8, 20.1 | 24.1, 20.4 |
| Stokes V | 20.8 | 20.5 |
| deepest negative (mJy/beam) | -4.01 | -1.49 |
| pixels below -0.3 mJy/beam | 2249 | 261 |
| beside the sources, median | -7.8% | -3.0% |
| the target, on the axis (35 mJy/beam) | -7.5% | -1.3% |
| the brightest, 29.5 arcmin out (45) | -8.8% | -3.3% |
| the other eight measured | -5.3 to -9.5% | -1.8 to -4.1% |
| restoring beam (arcsec, PA) | 4.44 x 3.49, 39 | 4.56 x 3.94, 22 |
| solve against it | 0.89 to 1.11 | 0.98 to 1.03 |

The hole beside the target is gone to within what its own structure
leaves. What stays goes with distance from the axis, 2 to 4% of the
sources away from it, which is where the scan stood before its gains
across the band and the peeling (-1.1 to -2.9%, finding 122). The beam is
a little wider: the weights follow the amplitudes, and the south arm's
two antennas count for less than they did.

**The first target through the whole chain, 2026-10-11 11:23 to 18:28
SAST (`targets_run3/`, G0212-0335, scan 15; 7 h 05 min).** Scans 6 min,
field 3 h 53 min, peeling 1 h 16 min, what is left 1 h 55 min, the beam
under a minute.

| | field, all sources | five sources out, calibrated on the rest |
|---|---|---|
| image (µJy/beam); quietest tenth; data | 23.5; 20.2; 19.5 | 21.9; 19.5; 18.9 |
| deepest negative (mJy/beam) | -1.38 | -0.55 |
| restoring beam (arcsec, PA) | 4.60 x 3.84, 27 | 4.36 x 3.43, 40 |
| the target: peak (mJy/beam); flux within an arcminute (mJy) | 34.7; 79.1 | 32.5; 79.6 |
| the source 103 arcsec south of it: peak; flux within 25 arcsec | 10.5; 37.7 | 9.2; 38.0 |

- **The sources chosen** by rule were five: 44.6, 29.2, 21.5, 16.2 and
  15.9 mJy/beam at 29.5, 17.1, 16.0, 27.1 and 37.2 arcmin. By hand there
  were six; the sixth is at 9.85 mJy/beam now, under the limit. Where the
  five stood the image has an rms of 27 to 31 µJy/beam and nothing beyond
  +-0.15 mJy/beam.
- **The peaks fell and the fluxes did not.** Every source is 5 to 13%
  lower at its peak in the last image than in the field's, and has the
  same flux to within a few per cent (eight sources measured). The beam
  is 15% smaller in area, the weights having followed the gains across
  the band once more, and these sources are resolved.
- **The flux scale is the calibrators'.** The amplitude solves are
  normalised to their median. NVSS has the target at 43.5 +- 1.4 mJy and
  TGSS at 314 +- 32: a power law between them gives 79.0 mJy at 710 MHz,
  and the image has 79.6 within an arcminute, on the axis. For the source
  to the south they give 41 and the image has 38 within 25 arcsec.
  **The fluxes of findings 119 to 128 were too high by about a fifth**
  (the target at 39.0 mJy/beam and 102 mJy): those data carried the
  pipeline's amplitude self-calibration and the amplitude rounds of the
  hand chain, none of them normalised. On this scale that image's 25.2
  µJy/beam was about 21 and its deepest negative -0.31, so the two
  results are level in noise and this one has the deeper holes.
- **What the image still holds.** -0.55 mJy/beam seven arcsec from the
  target's core (-1.7% of its peak; -1.1% in the field's image), -0.55
  beside the source to its south, -0.51 beside a 9 mJy source 16 arcmin
  out, and 3 to 6% of the peak beside each remaining source of 4 to 10
  mJy/beam, still north or south of it. 1609 pixels are below -5 sigma
  where noise would put 8. The recipe on what is left made the holes near
  the axis deeper in mJy/beam (-0.38 to -0.55 at the target, -0.29 to
  -0.55 at its neighbour) while lowering the noise: not understood.
- **The beam correction** goes out to 21.7 arcmin, where the narrowest
  sub-band's beam is 0.3 and the mean beam 0.46: a fifth of the image.
  Checked on the image: 1.000 at the target, 0.996 at 1.6 arcmin, 0.876 at
  9.2, 0.647 at 16.5 (the polynomial at 685 MHz gives 0.648), blank
  beyond.

### 132. What is left after peeling was calibrated against cores it no longer has

The last image of G0212-0335 has deeper holes near the axis than the
field's image had before peeling (finding 131). Where they grew, from the
images of each step, deepest negative within 14 arcsec in mJy/beam:

| | field, final | after peeling, phases and delays | after the amplitudes | final |
|---|---|---|---|---|
| the target (33 mJy/beam) | -0.38 | -0.57 | -0.55 | -0.55 |
| 2554, 2464 (9) | -0.29 | -0.78 | -0.63 | -0.55 |
| 1820, 3216 (9) | -0.36 | -0.67 | -0.57 | -0.51 |
| 1712, 1956 (8) | -0.34 | -0.40 | -0.32 | -0.27 |
| 1254, 3342 (4) | -0.16 | -0.31 | -0.25 | -0.23 |
| 3000, 2895 (5) | -0.20 | -0.39 | -0.29 | -0.25 |

They doubled between the field's image and the first image of the recipe
on what is left, and the rounds against the whole field then took part of
it back. Between those two images lie the peeling and three solves
against cores-only models, two of phase and one of delay. Two things are
known of them.

*The solves.* With the five sources out, the cores-only model holds 76
mJy where the field's held 180 to 370. `phase1` found, against it, -8 to
-14 degrees on ten antennas of the central square and the two inner ones
of the south arm (C05, C08, C09, C10, C11, C14, S01, S02 and others), the
same in both hands to a degree, and next to nothing on the arms; `delay`
then found 0.2 to 0.5 ns on the same central antennas (C01, C09, C10,
C12, C14). The same in both hands is not noise and not the antennas: it
is the model. The central square makes the short baselines, which see
everything the cores-only model leaves out, and in this field that is now
most of the flux. These data had phases good to 2 degrees when the field's
recipe ended. A cores-only model is for data whose phases are tens of
degrees out, where it breaks the degeneracy between a lopsided source and
an offset one; here it moved antennas that were right.

*The peeling.* It changes the weights. The restoring beam is 4.60 x 3.84
arcsec in the field's last image and 4.39 x 3.50 in the first one after,
which is what it was before the amplitude round against cores had put
S02 and S04 down and W04 and S06 up; through the five sources the minor
axis goes 3.72, 3.69, 3.61, 3.52, 3.49, 3.45, a step at each source's own
calibration. A source is subtracted from data that have its own gains
applied, and the rest of the field's model is added back to that: what is
not in a model, the noise and the faint sources, keeps the gains of five
directions that are not its own, and the weights follow them. That is
finding 126's design and it was not seen as a fault then.

Which of the two deepened the holes the images made so far cannot say:
there is none of the peeled data before the first solve. **A test that
does:** `selfcal_recipe.py --whole-field-only` leaves out every round
against a cores-only model, so the first thing it makes is an image of
the data as they came, and what follows is amplitude and phase against
the whole field, and the bandpass. On
`targets_run3/G0212-0335/peel/source5/source5_peeled.ms`, in a folder of
its own beside the run: if its first image already has the holes at
-0.57 and -0.78, it is the peeling; if it has them near -0.38 and -0.29,
it was the solves, and its last image is what the target should have
been. 1087 tests pass.

**On data, 2026-10-11 19:32 SAST (`target_test/`, beside the main run):
it was the solves.** The image of the peeled data before any solve, with
the beam the peeling left (4.39 x 3.50 arcsec):

| deepest negative within 14 arcsec (mJy/beam) | field, final | peeled, before any solve | after the cores-only phases and delays |
|---|---|---|---|
| the target | -0.38 | -0.37 | -0.57 |
| 2554, 2464 | -0.29 | -0.36 | -0.78 |
| 1820, 3216 | -0.36 | -0.38 | -0.67 |
| 1712, 1956 | -0.34 | -0.28 | -0.40 |
| 1254, 3342 | -0.16 | -0.14 | -0.31 |
| 3000, 2895 | -0.20 | -0.25 | -0.39 |

Peeling left the sources near the axis as they were, change of weights
and all, and the deepest negative of the inner field is -0.38 where the
three cores-only solves had it at -0.78. Against the whole field the
first solve finds 2 degrees and amplitudes of 0.98 to 1.04.

`run_targets.py` now gives the recipe on what is left
`--whole-field-only`, and has it work in `<target>/target_whole_field`:
a `target` folder that an earlier run finished is neither used nor
removed, and the stage and the beam correction are done again from the
peeled data, which are kept. 1088 tests pass.

The fault in the peeling's design stands as described and is not what
hurt here; it is left as it is.

**The test's last image, 21:20 SAST (1 h 48 min beside the main run).**
G0212-0335 with five sources out, solved against the whole field only,
beside the same with the cores-only rounds first:

| | with the cores-only rounds | whole field only |
|---|---|---|
| image (µJy/beam); quietest tenth; data | 21.9; 19.5; 18.9 | 21.8; 19.5; 18.9 |
| deepest negative (mJy/beam) | -0.55 | -0.375 |
| pixels below -5 sigma; -4 sigma | 1609; 5351 | 1285; 4946 |
| across the band, phases up to | 12.5 degrees rms | 5.4 |
| beside the target (32.7 mJy/beam, 79.4 mJy) | -0.55 (-1.7%) | -0.375 (-1.1%) |
| 2554, 2464 | -0.55 | -0.24 |
| 1820, 3216 | -0.51 | -0.29 |
| 1712, 1956 | -0.27 | -0.28 |
| 1254, 3342 | -0.23 | -0.17 |
| 3000, 2895 | -0.25 | -0.17 |

The bandpass round of the earlier run had spent itself undoing the
delays. What is left beside the sources of 4 to 10 mJy/beam is 2.4 to
4.4% of their peak, still north or south, and 1.1% beside the target: the
field's own, which peeling and this stage neither add to nor take away.

### 133. The same rounds in the solve towards each peeled source

`peel_sources.py` calibrates towards each source with the recipe as it
stood: a cores-only phase round, a cores-only delay round, then amplitude
and bandpass against the source's whole image. The data it is given are
the field's, calibrated already, which is the case of finding 132. What
the image of each source read before any round and after the two
cores-only ones, in µJy/beam (quietest tenth):

| | source | before | after phase and delay | final |
|---|---|---|---|---|
| G0212-0335 | 44.6 mJy/beam, 29.5 arcmin | 29.6 (25.2) | 31.3 (25.6) | 29.5 (25.1) |
| | 29.2, 17.1 | 29.0 (25.8) | 28.9 (25.7) | 28.1 (25.7) |
| | 21.5, 16.0 | 29.0 (26.2) | 29.2 (26.3) | 27.9 (25.8) |
| G0213-0251 | 131 and 67 mJy/beam, a double of 300 mJy, 9 arcmin | 39.6 (33.8) | 53.8 (39.1) | 47.9 (35.6) |

No source was better for the two rounds, and the bright double was a
third worse and did not come back: its peak fell from 129.8 to 126.9
mJy/beam, the deepest negative beside it stayed at -3.2. A cores-only
model of a double holds its two peaks and not what lies between and
around them, and the phases solved against it are wrong by what that is
worth. On G0212-0335 the sources are compact and it cost little.

`peel_sources.py` now gives each source's recipe `--whole-field-only`:
amplitude and phase against the source's whole image, then the bandpass,
or the first alone for a source too faint for more. 1091 tests pass, one
of them that the recipe takes each source's command and makes those
rounds of it. **Not run on data.** A source already peeled is not peeled
again by this; G0213-0251's four were taken out the old way.

### 134. Four targets through, two of them with two scans; an amplitude that was real and was refused

**Read on 2026-10-13 at 23:15 SAST, 44 hours after the last reading**
(the session was idle and the browser signed out). `targets_run3/`, one
run since its restart on 11 October at 22:00 SAST, no step failed.
In µJy/beam; "near the centre" is a robust rms of the image within 2
arcmin of the pointing centre, and in the ring from 2 to 5 arcmin.

| target (scans) | image | quietest tenth | data | within 2 arcmin | 2 to 5 arcmin | peeled | deepest (mJy/beam) | hours |
|---|---|---|---|---|---|---|---|---|
| G0212-0335 (15) | 21.8 | 19.5 | 18.9 | 31.7 | 26.2 | 5 | -0.36 | 7 |
| G0213-0251 (17) | 25.7 | 20.8 | 16.9 | 34.1 | 31.2 | 4 | -0.63 | 6 |
| G0216-0604 (5, 18) | 18.8 | 15.2 | 13.5 | 31.4 | 27.4 | 8 | -1.58 | 17 |
| G0218-0245 (6, 20) | 20.2 | 15.8 | 13.0 | 31.6 | 28.5 | 8 | -0.65 | 14 |

G0220-0318 (8, 21) has its field's image, 14.5 (12.7) against 12.0, and
was at the seventh of its sources.

**Two scans work.** Phases per scan and within them settle (175, 37, 7
degrees on G0216-0604; 162, 28, 5 on G0218-0245), delays are found per
scan (S04 3.4 and 3.2 ns, against 2.3 in the single scans later in the
day), and an image of the whole field takes about as long as for one
scan. Their data's noise is 13 against 17 to 19 for one scan.

**The noise near the pointing centres is not the image's.** 26 to 34
within 5 arcmin of the centres of the three fields with a source of 131
to 257 mJy/beam 9 to 11 arcmin out (G0213-0251, G0216-0604, G0218-0245):
twice the quietest tenth, which lies in the corners. The figures a run
prints are of the whole image and of its quietest tenth, and neither is
where the targets are.

*Corrected on 2026-10-14 at 08:05 SAST.* Written first as "every centre,
a bright source in the field or none", on four fields. The fifth says
otherwise: G0220-0318, whose brightest source is 44.5 mJy/beam, has 16.8
within 2 arcmin, 16.4 from 2 to 5 and 16.4 from 5 to 10, against 15.0
for the image and 12.7 for the data. On G0212-0335 (brightest 45) the 2
to 5 arcmin ring is 26.2 against 21.8, and within 2 arcmin the target's
own structure is in the figure. So the excess goes with the very bright
sources that were peeled, whose rings (below) and sidelobes are spread
over the inner field: it is the peeling of those sources that has to do
better, not the field's calibration.

**An amplitude that was real, refused (G0216-0604, and G0220-0318).**
The round against cores found S04 beyond 2.5 in one of G0216-0604's two
scans, and the limit of finding 128 answered "not believed and not
applied". It was the antenna: S04 is the one found at 1.25 to 1.47 on
G0212-0335 and its delay is larger early in the day, so its calibrator
amplitudes are further out there. Left as it was, the field's images
have -8.8 to -7.2% beside the sources where the other fields have -3%,
and the last image -9 and -10% beside two sources of 8 and 9 mJy/beam and
-3.9% beside the 22 mJy/beam source 2.5 arcmin from the centre. On
G0220-0318 the same happened to S02 in one scan; there the rounds against
the whole field then moved it (0.80 to 1.32) and the field ends at -3.6%.

`experiment_phase_offsets.py`: the limits are 0.25 to 4 (`--amplitude-limits`),
and an antenna whose amplitude is still outside them, in a solve with one
solution per scan, is flagged for that scan in the corrected data, not
left in at an amplitude that is not believed. `--keep-unsolved` keeps
that from happening, and `selfcal_recipe.py` now passes it to the gain
solves as well as to the bandpass one, so a solve towards one source or
on what is left of a field still flags nothing. The bins of a bandpass
keep their own limits, 0.4 to 2.5 (`--bin-limits`). 1096 tests pass.
**Not run on data.** G0216-0604 has to be done again from its field.

**A ring left around each very bright source that was peeled.** The
deepest negatives of three images lie 96 to 180 pixels from where a
source of 131, 176 or 257 mJy/beam stood, just outside the disc of 90
pixels that was taken out with it: -0.6 mJy/beam on G0213-0251 and
G0218-0245, -1.6 on G0216-0604. They are components of the field's model
that described that source's errors out there; the model of "the rest"
keeps them, and they are given back after the source is gone. Inside the
disc the image is at 50 to 84 µJy/beam where the faint sources' patches
are at 25 to 35. A disc that grows with the source would take them; not
written.

### 135. All nine targets through once

`targets_run3/`, from 11 October 11:23 SAST to 16 October 04:04 SAST, one
command, restarted twice: on 11 October to take in finding 132, and on 15
October after the Jupyter server it ran under was restarted at 11:42
SAST and took the run with it (the new server was on another machine, so
four locks had to be removed by hand). Fixes were uploaded while it ran:
the per-source solves of finding 133 from G0222-0321's peeling on, the
amplitude limits of finding 134 from G0224-0620's field on. In µJy/beam;
the rings are robust rms within 2 arcmin of the pointing centre, from 2
to 5 and from 5 to 10.

| target | scans | image | quietest tenth | data | centre: 0-2, 2-5, 5-10 | peeled | at the centre (mJy/beam) |
|---|---|---|---|---|---|---|---|
| G0212-0335 | 15 | 21.8 | 19.5 | 18.9 | 31.7*, 26.2, 24.0 | 5 | 32.7 |
| G0213-0251 | 17 | 25.7 | 20.8 | 16.9 | 34.1, 31.2, 32.0 | 4 | 0.75 |
| G0216-0604 | 5, 18 | 18.8 | 15.2 | 13.5 | 31.4, 27.4, 26.0 | 8 | 1.2 |
| G0218-0245 | 6, 20 | 20.2 | 15.8 | 13.0 | 31.6, 28.5, 27.8 | 8 | 4.3 |
| G0220-0318 | 8, 21 | 15.0 | 13.2 | 12.7 | 16.8, 16.4, 16.4 | 7 | |
| G0222-0321 | 9, 23 | 14.1 | 12.7 | 12.3 | 16.0, 15.1, 15.2 | 8 | 0.37 |
| G0224-0620 | 11, 24 | 15.5 | 13.1 | 11.9 | 26.1*, 19.2, 17.1 | 7 | 29.2 |
| G0229-0556 | 12, 26 | 14.4 | 12.5 | 11.6 | 16.8, 15.7, 15.6 | 7 | 2.35 |
| G0231-0533 | 14, 27 | 14.1 | 13.1 | 12.3 | not measured | 4 | |

\* the target itself, resolved, is in the inner circle.

- **Five fields are near the data's noise to the centre**: G0220-0318,
  G0222-0321, G0224-0620 outside its target, G0229-0556 and, by its
  whole-image figures, G0231-0533.
- **Three are not**: G0213-0251, G0216-0604 and G0218-0245, the fields
  with a source of 131 to 257 mJy/beam 9 to 11 arcmin out. Their rings
  of -0.4 to -1.6 mJy/beam just outside the peeling disc, and the excess
  within 10 arcmin of the centre, are those sources' (finding 134).
- **G0216-0604 is to be done again** from its field: S04's amplitude was
  refused in one of its scans (finding 134). G0220-0318 and G0222-0321
  had S02 refused in one scan, but their whole-field rounds then moved it
  and their holes are at -3.6 and -3.0%, as in the other fields.
- Hours per target: 7 and 6 for the two single scans; 12 to 17 for two
  scans, the field stage 7 to 9 of them.

### 136. A disc that grows with the source

What finding 134 found: around the three sources of 131, 176 and 257
mJy/beam that were peeled, rings of -0.4 to -1.6 mJy/beam 96 to 180
pixels out, just beyond the disc of 90 pixels, and an excess of noise
over the inner field of those three targets (26 to 34 µJy/beam within 5
arcmin of the centre, against 15 to 17 in the fields without such a
source). The field's model held that source's errors out there; the model
of "the rest" kept them, and the peeling gave them back.

`peel_sources.py` now gives each source its own disc: 90 pixels
(`--radius`) up to 30 mJy/beam, larger by the square root of how much
brighter, at most three times as large, and never reaching closer to the
pointing centre than the kept radius (4 arcmin): the target stays out of
every disc. A source whose peak lies inside the disc of a brighter one
goes out with it and is not taken on its own. The disc goes to both the
extraction and the peel, and each source already out is passed on with
its own (`--also-without X Y RADIUS`; `experiment_extract_source.py`
still takes X Y alone and gives it the present source's disc).

| field | source | at | disc (pixels) | rings seen at | disc's edge from the centre |
|---|---|---|---|---|---|
| G0213-0251 | 131 mJy/beam | 9.0 arcmin | 188 | 96 to 179 | 5.8 arcmin |
| G0216-0604 | 177 | 11.0 | 218 | 97 to 103 | 7.4 |
| G0218-0245 | 257 | 10.5 | 264 | 96 to 177 | 6.1 |
| G0213-0251 | 87 | 24.0 | 153 | | 21.4 |
| G0224-0620 | 99 | 21.5 | 163 | | 18.7 |
| G0229-0556 | 71 | 14.6 | 139 | | 12.3 |
| G0212-0335 | 45 | 29.5 | 110 | | 27.7 |

What a larger disc costs: whatever else lies in it is taken out too and
is not in the last image. Around the three bright sources that is the
rest of each source's own structure and the faint sources nearest it.

1100 tests pass. **Not run on data.** The square root is a guess fitted
to where three rings ended, not a measured law; if a ring is left just
beyond the new disc the scale is wrong, and the images will show it.
G0213-0251 and G0218-0245 are to be peeled again from their fields, and
G0216-0604 done again from the start.

### 137. Where the hours go, and the targets side by side on SLURM

Sushant, 2026-10-16: "Can't these runs be faster? Because I want this
pipeline to be as fast and reliable as possible."

**Where the hours go**, G0222-0321 (two scans, 12.3 h on the Jupyter
node, 6 threads at the lowest priority):

| | hours |
|---|---|
| nine rounds against cores-only models (model and solve, about 7 min each) | 1.1 |
| five images of the whole field, 5120 pixels (71 to 76 min each) | 6.0 |
| peeling eight sources (about 21 min each) | 2.8 |
| three models of the inner field (24 min each) | 1.2 |
| the last image, 5120 pixels | 1.15 |

Of an image of the whole field, WSClean's own clock: inversion 9 min,
prediction 6 min, deconvolution 52 to 57 min. Three quarters of every
large image, and so about two thirds of a target, is the minor cycles of
the multi-scale clean, run whole on one image.

**Parallel deconvolution was never tested after the crash it was
blamed for.** On 2026-10-07 WSClean 3.1 aborted as deconvolution began,
`-parallel-deconvolution` was taken out first, the abort came back
without it, and the cause was the local-rms mask (finding 118); the
option stayed off for WSClean before 3.2. A test now images
G0222-0321's finished data (`target_whole_field/.../round3.ms`) with
`--parallel-deconvolution 2560` and nothing else changed, beside its own
last image.

**The targets side by side.** `run_targets.py` held its whole work
folder, so one run did the targets one after another on one Jupyter
node; nine took four and a half days, and the restart of that node's
Jupyter server on 15 October stopped everything. Now:

- each target is held by the run working on it, for the length of that
  run, in its own folder (`<workdir>/<target>/running.lock`); a target
  another run holds is left to it. One run for all the targets works as
  before; so do several runs of one target each, in one folder;
- each target's record is its own file (`<target>/target_report.json`),
  written only by the run that holds it; `targets_report.json` is made
  from them at every save, written whole under a name of its own and then
  moved, so several runs can write it;
- a run of an earlier version, which held the whole folder, is waited
  for: a new run refuses while that claim is held by a run still going;
- a claim made inside a SLURM job records the job; it is taken over once
  `squeue` no longer lists that job, wherever it ran. A job that is
  killed no longer leaves claims to be removed by hand (on 15 October
  four had to be);
- threads (default: the CPUs a SLURM job was given), WSClean's memory
  (`--memory-gb`) and `--parallel-deconvolution` reach the field's
  recipe, the peeling and its recipes, and the recipe on what is left;
- `run_target_slurm.sh`: one target per job on the Main partition, 16
  CPUs, 120 GB, three days; it checks for a Python with numpy and for
  singularity before starting, and gives WSClean three quarters of the
  job's memory.

1111 tests pass. **Not run on ilifu's compute nodes.** Whether their
Python has numpy, and whether singularity is on their PATH, the job's
first lines will say.

**The test, 2026-10-16 11:00 to 12:01 SAST** (`speed_test/`, on the
Jupyter node beside the run, 6 threads), against the image the run made
of the same data on 14 October (on another Jupyter node, alone):

| | whole (14 Oct) | in sub-images of 2560 |
|---|---|---|
| WSClean 3.1 | ran | ran: the abort of 7 October does not come back |
| major iterations | 6 | 7 |
| inversion, prediction, deconvolution (min) | 9.5, 6.3, 49.9 | 6.3, 4.3, 48.1 |
| image; quietest tenth (µJy/beam) | 14.14; 12.66 | 14.14; 12.65 |
| brightest; deepest (mJy/beam) | 12.44; -0.358 | 12.44; -0.368 |
| beside the sources, median | -1.95% | -1.98% |

The image is the same. The cleaning is not faster on 6 threads that are
shared with another run: the minor cycles are bound by the threads they
have, whole or split. What splitting is for is many threads, which a
whole multi-scale clean does not use well; on a node of its own with 16
it is worth having, and costs nothing in the image. The gain has to come
from more threads and from the targets running side by side.

### 138. First targets on SLURM, and what a node of its own does and does not buy

**2026-10-16 12:15 SAST**, three jobs from `run_target_slurm.sh` with
`--parallel-deconvolution 2560`: G0213-0251 (compute-104, its last
image), G0216-0604 (compute-240, from its field) and G0218-0245
(compute-240, its peeling). The compute nodes have `/usr/bin/python3`
with numpy and singularity on their PATH; WSClean got `-j 16 -abs-mem 90`.
No claim had to be removed by hand.

| | Jupyter node, 6 threads | compute node, 16 threads |
|---|---|---|
| a round against a cores-only model | about 7 min | 4 to 5 min |
| G0216-0604's eleven such rounds | (63 min on G0222-0321) | 43 min |
| an image of the whole field, one scan (G0213-0251's last) | about 60 min | 68 min |
| its major iterations | 5 to 6 | 20, the most allowed |
| its inversion, prediction, deconvolution (min) | | 8.8, 7.0, 50.7 |

The cleaning does not get faster with the threads: 50 minutes on 6
shared ones whole, 48 on 6 split, 51 on 16 split. Split, the image takes
up to the most major iterations it is allowed instead of six, because
each sub-image stops at its own threshold; the image is not the worse for
it (G0222-0321 in finding 137; G0213-0251 below). So what a target gains
on SLURM is the small rounds and, above all, running beside the others;
not the large images. Those are bound by the multi-scale minor cycles of
WSClean 3.1. A newer WSClean (`wsclean-latest.sif`) and the whole image
on 16 threads are what `speed_test_slurm.sh` compares next, each on a node
of its own, on G0222-0321's finished data.

**G0213-0251 with the disc of 188 pixels** (finding 136):

| | disc of 90 | disc of 188 |
|---|---|---|
| image; quietest tenth; data (µJy/beam) | 25.7; 20.8; 16.9 | 24.1; 19.7; 16.0 |
| 96 to 180 pixels from the double: deepest | -0.63 | -0.16 |
| inside the disc: rms; extremes | | 37 µJy/beam; +-0.2 mJy/beam |
| pixels below -5 sigma | 9641 | 5287 |
| within 2 arcmin of the centre; 2 to 5; 5 to 10 | 34.1; 31.2; 32.0 | 32.7; 29.3; 30.1 |

The ring is gone. The excess near the centre is not: within 5 arcmin the
image is at 29.7 µJy/beam where its corners are at 21.2, and Stokes V
there at 18.0 where its corners are at 15.0. Within 5 arcmin, 2.1% of
the pixels are below -3 times the quietest tenth's 19.7 and 5.0% above:
symmetric about a scatter of 30, as noise would be, with real sources on
top. So it was not the bright double's ring; it is something spread over
the inner field and as much negative as positive. Not understood. One
candidate: sources too faint for the 4-sigma mask are never cleaned, and
where the beam holds most of them, near the centre, their sidelobes add
up.

### 139. WSClean 3.6, in sub-images: the same image in a third of the time

**2026-10-16 13:59 to 14:57 SAST**, `speed_test_slurm.sh`: G0222-0321's
calibrated data from its last stage, imaged on the pipeline's grid and
cleaning, each on a compute node of its own with 16 threads.

| | WSClean 3.1, whole | WSClean 3.6, in sub-images of 2560 |
|---|---|---|
| node | compute-003 | (queued later) |
| major iterations | 6 | 5 |
| inversion, prediction, deconvolution (min) | 5.3, 3.6, 46.0 | 3.2, 1.7, 14.3 |
| about | 55 min | 19 min |
| image; quietest tenth (µJy/beam) | 14.14; 12.66 | 14.11; 12.64 |
| brightest; deepest (mJy/beam) | 12.44; -0.352 | 12.44; -0.344 |
| beside the sources, median | -1.9% | -1.9% |

The image is the same to the third figure; it is made three times as
fast. WSClean 3.6 splits the cleaning well where 3.1 needed up to twenty
major iterations to do it (finding 138), and its gridding and prediction
are faster as well. WSClean 3.6 whole was still cleaning after 55
minutes. The container is `wsclean-latest.sif` (WSClean 3.6, 2025-02-07);
every image's report records the version that made it.

All of the pipeline's calls of WSClean take their gridder option from the
version (`gridder_arguments`), so 3.6 needs no other change.
`run_target_slurm.sh` now runs WSClean 3.6 with
`--parallel-deconvolution 2560` unless told otherwise. A large image of
two scans was about 75 minutes; it should be about 25. A target of two
scans should take five to six hours where it took twelve to fourteen.

### 140. GMRT LTA files as input

The GMRT Online Archive delivers raw data in GMRT's own LTA format as
well as FITS. The pipeline took only a Measurement Set, UVFITS or
FITS-IDI, so an LTA file had to be converted by hand first. It now takes
an LTA file as its input (`ms: /path/obs_gwb.lta` in `config.yaml`, or
the positional argument) and runs the ingest stage on it as follows:

1. `lta_convert.py` (new) runs NCRA's `listscan` on the file, which writes
   `<name>.log` naming `TEST.FITS` as the output, then `gvfits` on that
   log, which writes the FITS. Both run in a new, empty staging folder on
   a link to the LTA file, with no terminal input. So nothing is written
   beside the archive file, a failed conversion leaves no partial FITS, and
   two conversions cannot share a folder. The log is used exactly as
   listscan wrote it: every scan, antenna and channel goes into the FITS,
   because the pipeline flags data and chooses scans itself. These are the
   same steps CAPTURE takes (`fromlta`), except that CAPTURE renames the
   FITS by editing the log; here the FITS is renamed after gvfits has
   finished.
2. The FITS goes to `<outdir>/ingest/lta/<name>.FITS`, with
   `<name>.FITS.conversion.json` recording the LTA file it came from (path,
   size, time, SHA-256) and the two programs (command, SHA-256). The
   listscan log and both programs' output are kept beside it.
3. The FITS is inspected like any input (it must be random-groups UVFITS,
   and its telescope must agree with `--telescope`) and imported with
   `importgmrt`, as uGMRT FITS always was. The structure comparison of
   finding 35 applies to it.
4. Once the Measurement Set is validated and the report and completion
   marker are written, the converted FITS is removed, because it is as
   large as the data and the MS holds everything it held. `--keep-lta-fits`
   (`keep_lta_fits: true`) keeps it. If the import fails, the FITS is kept,
   and the next attempt reuses it instead of converting again.

The ingest report's `source_fits` is the LTA file itself, so a resumed run
is bound to the archive file's SHA-256 exactly as it is bound to a FITS
file's. The report also gives `source_format: lta` and the conversion
record. Reports written before this change have no `source_format`; they
are read as FITS, so finished FITS runs resume unchanged. The ingest
command for FITS input is unchanged.

Checks: an LTA input is recognised by its name (`.lta`; `.ltb`, the GSB
lower sideband; `.lta_1` etc., parts of a split recording) or by
`--fits-format lta`. With any `--telescope` other than `auto` or `ugmrt`
it is refused. Before the run starts, `pipeline.py` checks that listscan
and gvfits can be found and stops with directions if not (a dry run only
warns). When it plans the ingest of an LTA file, it needs `fits_ingest.py`
and `lta_convert.py` of interface 1, so an old `fits_ingest.py` stops the
run at the start and not hours in. The conversion is refused when the
output folder has less free space than the size of the LTA file. It warns
below three times that size, because the FITS and the MS each need about
that much. Where the programs live (`--listscan`, `--gvfits`) and whether
the FITS is kept are not science: the resume comparison of stage commands
leaves them out, so moving the programs does not redo a finished ingest.

The programs are not part of CASA. Linux builds of both are in
https://github.com/ruta-k/CAPTURE-CASA6. Each option takes a path, a name
on PATH, or a whole command, for example one that runs the program in a
container. They run inside the staging folder, so a relative program path
is resolved to an absolute one first; paths inside a wrapper command must
already be absolute. `lta_convert.py IN.lta OUT.FITS` also converts a file
on its own.

Not covered: legacy GSB dual-frequency recordings, which CAPTURE also
leaves to be converted by hand one polarisation at a time. Joining a split
recording (GMRT's `ltacomb`); each part is converted on its own. Neither
program has been run on real data here: there is no LTA file on this
machine and no GMRT tools, so the tests use stand-in programs that behave
as CAPTURE's code expects listscan and gvfits to (a `<name>.log` naming
`TEST.FITS` in the working folder, then the FITS written there). The first
real LTA run should be read at `<outdir>/ingest/lta/` (log and program
output) before calibration results are trusted. 38 new tests; 1149 pass.

Addendum (2026-10-19): `run_target_slurm.sh` had the Band 4 run's
settings built in (`--workdir targets_run3 --channel-range 51 768
--scans G0212-0335=15`, log in `targets_run3/`). With any other dataset
`run_targets.py` stops at once, because `--scans` names a target that is
not among the calibrated ones, and that option cannot be removed by one
given later. The script now passes only the job's resources and WSClean;
the products folder, work folder, channel range, band and any scans
follow the target on the `sbatch` line, and `sbatch --output=` places
the log. The header gives the Band 4 run's command in the new form.

### 141. One command to run on the Jupyter node or as a SLURM job

The first Band 5 start on ilifu (2026-10-19) failed at once with "CASA
executable not found". On ilifu CASA is an environment module
(`module load casa/6.7.0`), and that terminal had none loaded. The next
two starts as SLURM jobs were refused with "Requested node configuration
is not available", which doesn't say what is wrong. The cause was
`config.yaml`'s job settings (16 tasks of 4 CPUs, then 32 tasks of 16
from the default `ntasks`), larger than any node. Running in a Jupyter
terminal also needed `nohup nice -n 19 ... &` to outlive the terminal.

`pipeline.py` now has three ways to start, all equivalent in what they
run:

- `python3 pipeline.py` runs in the foreground, as before. A SLURM job
  uses this form.
- `python3 pipeline.py --background` runs on this machine (the Jupyter
  node), in a new session with no terminal and at niceness 19. Its output
  is appended to `<outdir>/pipeline_run.log`. The command waits 20 seconds
  and, if the run has stopped by then (no CASA, a missing input, another
  run in the folder), shows the end of the log and returns its status;
  otherwise it prints the process number, the log and how to stop it.
  Closing the terminal or the tab doesn't stop the run; a restart of the
  Jupyter server still does, and the same command then resumes it.
- `python3 pipeline.py --submit` writes the SLURM job from `config.yaml`
  (or to `--write-sbatch PATH`) and submits it with `sbatch --parsable`.
  First it asks `sinfo` for the partition's node sizes and refuses a job
  no node can take, naming the CPUs and memory asked for and available.
  It then prints the job number, its log files and `scancel`. It notes
  when `--mpi off` leaves all but one of several tasks idle.

`casa_module` (`--casa-module`) names the module that provides CASA. The
run asks a login shell for the environment after `module load` (with
`env -0`) and makes it its own before any stage starts, so CASA is found
both on the Jupyter node and in a SLURM job. The PATH entries the run
had are kept after the module's own. A module that doesn't load stops
the run at once (a dry run only warns). The stage commands are
unchanged: they still name `casa`, which now resolves through the
module's PATH. So the resume comparison of commands sees no change.

A run may now be started on the Jupyter node and as a SLURM job on a
compute node. flock on a shared filesystem may not be seen from another
node, so the run lock now also records the host and its time is
refreshed every 60 seconds while the run lives. A run (and `--background`
or `--submit` beforehand) refuses an output folder whose lock names a
SLURM job `squeue` still lists, a live process on the same host, or
another host's process whose lock was refreshed within 5 minutes. A run
that ends clears its record. Records from before this change carry no
host and are left to flock, as before.

Neither launcher is used inside a SLURM job, and `--background` refuses
`--dry-run`, `--list-stages` and `--resume-status`. 26 new tests (the
module script runs in a real bash with a stand-in `module`); 1175 pass.
Checked on this Mac: `--background` with no CASA reports the error within
a second; `--submit` says sbatch is not available here.

Addendum (2026-10-20), a modular CASA. ilifu's `casa/6.7.0` module puts
no `casa` program on the PATH. It configures a modular CASA container
(`/software/astro/containers/casa-modular/casa-modular-v6.7.0-py3.10-2025-12-01.sif`)
and defines `casapy` as a shell function (`singularity exec IMAGE
xvfb-run -a ipython "$@"`) and the GUIs as aliases. A program can't call
a shell function, so the first run with `casa_module: casa/6.7.0` loaded
the module and still found no CASA.

A new setting, `casa_python` (`--casa-python`), is the command that
starts a Python with casatasks and casatools. Each CASA stage then runs
as `<command> script args` in place of `casa --nologger --nogui -c script
args`. The stage scripts already import CASA explicitly, draw with Agg
and don't use plotms, and they drop their own folder from `sys.path`
(where `inspect.py` would stand in for the standard module), so they run
unchanged. `fits_ingest.py`, the one without that guard, now has it.
When `casa_module` gives no `casa` program, the run asks the shell what
`casapy` and `casa` are after loading the module. If one of them runs a
container (`singularity` or `apptainer exec IMAGE.sif|.simg|.img`), the
stages run with `singularity exec IMAGE python3`, and the run says so.
`casa_python` given explicitly wins. Before the first CASA stage, the run
starts that Python once with `import casatasks, casatools` and prints
the CASA version, so a wrong image stops it at once. A modular CASA has
no mpicasa: CASA MPI is off with it, and `--mpi on` is refused. The
generated SLURM job carries an explicit `--casa-python`; one found from
the module is found again inside the job. 8 more tests (the container
detection runs in a real bash against a stand-in `module` that defines
`casapy` as ilifu's does); 1183 pass.

Addendum (2026-10-20), a Jupyter session is a SLURM job. On ilifu the
Jupyter server itself runs as a SLURM job, so `SLURM_JOB_ID` is set in
its terminals. `--background` refused to start there ("This is already
a SLURM job"), and `--submit` would have done the same. The lock rule
was wrong for the same reason: it held a lock while the job named in it
was listed, and a Jupyter session's job outlives a run started in it.
Both launchers now start from inside a SLURM job; `--background` notes
that the run lasts as long as that session's job. A lock is now judged
like this: on the same machine, by whether its process lives; from
another machine, it is free once its SLURM job has left the queue, and
otherwise held only while its time is being refreshed (within 5 minutes).
1186 tests pass.

Addendum to finding 140 (2026-10-20), first LTA run on ilifu. The
conversion stopped at once: `listscan` exited with status 127, "error
while loading shared libraries: libgfortran.so.3". The NCRA builds (from
CAPTURE-CASA6) link the gcc 5/6 Fortran runtime, which neither the CASA
6.7 modular container they now run in nor current systems carry. The
converters now look first for libraries in the program's own folder and
in a `lib/` folder beside it: `LD_LIBRARY_PATH` gets those folders ahead
of its own entries, for the two programs only. A failed program whose
output names a missing shared library gets an error that says which
library, where to put it, and where libgfortran.so.3 can be had
(conda-forge `libgfortran` 3.0.0). `lta_convert.py` is interface 2, and
`pipeline.py` requires it for LTA input, so a stale converter stops the
run at the start. 2 more tests; 1188 pass.

### 142. Antennas with no signal at the bandpass calibrator

The first Band 5 run (GWB, 4096 channels, observed 2021-02-22, the same
nine fields as Band 4) went from the LTA file through ingest, inspection
and flagging in about an hour. It stopped in calibration: "bandpass.cal:
ANTENNA1=3: group has no usable solution". C12 (no pre-bandpass phase
solution) and W06 (no delay solution) had already been excluded
automatically. In the bandpass solve, CASA reported "4 of 56 solutions
flagged due to SNR < 3" in every channel: those were C03 and C05, both
polarisations, on 3C48. Every other antenna kept 94% of its solutions,
the rest being the flagged band edges. C03 kept none and C05 1%. Both
passed the delay and pre-bandpass solves, which combine the whole band
and so have about 8 times the signal-to-noise per solution. Only 9% and
12% of their data had been flagged: a receiver with no signal produces
noise, and noise doesn't look like interference. These two antennas are
not working in this Band 5 observation.

The pipeline excluded such antennas after the delay and pre-bandpass
solves but not after the bandpass. Now it does the same there: an
antenna whose bandpass solutions are all, or all but 5%, unusable is
excluded. That covers C05's 1%, which the rule for the other solves,
exactly none, would miss. The usual limits apply (per event, share of
the array, antennas left, never the reference antenna). The delay,
pre-bandpass and bandpass tables are then solved again without it. An
antenna with more than 5% of its bandpass left stays a hard error for a
person to look at.

When a limit stops an exclusion, the error now names the limit and what
would pass. Here: "with the 2 excluded before (C12, W06) that makes 4 of
30 antennas (13%), more than calibration_max_auto_quarantine_delay_fraction
(0.1); 0.14 would allow it". The default 10% stays: losing a seventh of
the array should be a decision, made in `config.yaml`.

`lta_convert.py` also keeps gvfits's own `gvfits.log` (as
`<name>.FITS.gvfits.log`), which its output points to and which was
discarded. 6 new tests; 1194 pass.

### 143. When the two hands disagree: per-hand calibrator gains

The second Band 5 run (2026-10-21) passed the bandpass at 12:41 SAST.
C03 and C05 were excluded by finding 142, after C12 and W06, and the
delay, pre-bandpass and bandpass tables were solved again without them.
The phase solve on the calibrators then failed at 12:43 SAST:
"phase.cal: ANTENNA1=10: flagged fraction 0.914530 exceeds the configured
group maximum". CASA's per-antenna count of good solutions, out of 351
ten-second intervals (151 on 3C48, 200 on 0323+055): C11 30, S02 60, S03
60, the reference C00 100, C02 226, E03 246, the rest 308-351. The
flagging stage had flagged these antennas no more than the others
(11-14% on both calibrators). In the per-hand (G) pre-bandpass solve on
3C48, at the same cadence, C11, S02 and S03 had kept 100% and C00 88%.

The phase solve uses `gaintype='T'`, one phase for both hands, chosen
for uGMRT for its signal-to-noise. That holds only while each antenna's
R-L phase stays where the time-independent bandpass left it. Where it
drifts, the hands cancel, and a solve that each hand would pass alone
fails SNR 3. Every antenna that failed T had passed G on the same
calibrator, which no other cause on the list explains.

`run_phasecal` now solves again with `gaintype='G'` when the T table
fails its check. If G passes, `calibration/calibrator_gaintype.json`
records the switch, and every calibrator gain solve of the run uses G:
the phase, the scan-averaged phase transferred to the targets, and the
amplitude. A resumed run reads the record and solves with G at once. If
G fails too, the run stops with both reasons. The calibration report
records the switch (`calibrator_gaintype_switch`). Datasets on which T
passes, Band 4 among them, are solved as before. G costs a factor of
sqrt(2) in signal-to-noise per solution, which bright calibrators can
spare. 6 new tests; 1200 pass.

Confirmed on data (resumed run, 2026-10-21 13:08-13:14 SAST). T gave the
same counts again. G, on the same 351 intervals, gave good solutions per
hand (RR, LL) of C11 332/334, S02 331/330, S03 328/332, C00 298/303, C02
331/333, E03 320/322, and 298-335 for every antenna. phase.cal passed its
check (20% flagged, the excluded antennas included), the switch was
recorded, and the scan-averaged phase solve followed with G.

### 144. Averaging channels before the per-target chain (Band 5)

Band 5's calibration finished at 18:50 SAST on 2026-10-21, after 5 h 40
min. The split was slow because Band 5 targets keep all 4096 channels of
97.7 kHz: only Band 4 has `target_max_channel_width_khz`. `run_targets.py`
would have imaged them at that resolution, four times the visibilities
of Band 4's 390.6 kHz channels in every WSClean run and every copy of a
target's data. Smearing doesn't call for it: 390.6 kHz at 1250 MHz, with
Band 5's 2" beam, smears a source at the half-power radius (12.6') by
about 11%, less than the 16% Band 4 had at its own half-power radius.

`run_targets.py` and `select_scans.py` take `--average-channels N`. The
scans step's `split` then averages every N channels (`width=N`, weighted
by the data's per-channel weights), and makes that copy even when every
scan is kept. `--channel-range` still counts the calibrated data's
channels. Every later step gets the range in averaged channels, keeping
only channels made wholly of data: [ceil(first/N), floor(end/N)). For
Band 5 the data are channels 122-3973 (the flagging stage's 3% band
edges, 1061.9-1438.1 MHz, no gap fully flagged inside), so `122 3974`
with N=4 gives 31-993. With N=1 (the default) nothing changes, and the
scans rule and split record are as before, so finished Band 4 runs are
not redone. 6 new tests; 1206 pass.

Band 5 settings for `run_targets.py`, scaled from Band 4's: `--band B5`
(the NCRA polynomial; FWHM 25.3' at 1250 MHz, the 0.3 level at 14.2' at
the top of the band), `--cell 0.5` (about 4 pixels across the 2" beam,
as 1" was for Band 4's 4.4" beam), `--imsize 6144` (51', about twice the
half-power radius each way, as 5120 x 1" was for Band 4), and
`--inner-imsize 4096` (34', the corrected area with the margin Band 4's
3072 had). To be checked on the first target's images before the other
eight are started.

### 145. Band 5 delays do not drift: a hypothesis tested and dropped

The Band 5 test target (G0222-0321) had, after its amplitude rounds, an
image at 21.2 µJy/beam against 19.8 for the data, but a negative of
-7.4% of the peak beside a source. Its self-calibration showed C01 and
C08 with per-scan phases 130-175° apart between its two scans (3 h
apart), delays of 3.5-15 ns in one hand, and 1-minute phases moving by
up to 46° within a scan. I took this for antenna delays drifting during
the night: calibration solves K once, on 3C48 in the first half hour,
and applies it to the whole 5.7 h. That would also have explained the
T/G failure of finding 143. I did not change calibration on that
reading. `diagnose_delay_drift.py` (new; tests in
`test_diagnose_delay_drift.py`) solves one delay per scan on the
calibrators with the pipeline's delay and bandpass tables applied, and
prints what is left per antenna, hand and scan.

Result (2026-10-22 13:10-13:15 SAST, 3C48 and seven scans of 0323+055
through the night): every antenna stays within +-0.4 ns of its 3C48
delay, except W03-W05, which rise slowly to +0.6 ns, and S04, which sits
at +1.0 to +1.3 ns on 0323+055 against 0 on 3C48. C01 and C08 are steady
to 0.1 ns. So there is no delay drift worth tracking (0.6 ns is less
than a quarter turn across 400 MHz), and calibration stays as it is.
S04's offset depends on the field, not on the time, which points to an
error in its recorded position. Self-calibration, which solves delays
per target, absorbs it (G0222-0321's first round found S04 at 163°, and
the delay round found 1.6 ns). What C01 and C08 show on the target is
phase jumps in one hand, not drift; the 1-minute per-hand self-cal
rounds follow them. The -7.4% stays open, to be judged on the image.
3 new tests; 1209 pass.

Band 5 test target (2026-10-22, SLURM job 14420459, G0222-0321, the
settings of finding 144): through all stages in 4 h 02 min (scans 1 min,
field 2 h 37 min, peel 11 min, last stage 1 h 12 min, primary beam 1
min). The beam is 2.26" x 2.13" (2.29" x 1.82" after peeling), 4.4
pixels at 0.5", as planned. Final image 22.5 µJy/beam (21.4 in its
quietest tenth) against 21.2 for the data. Rings around the centre
(0-2', 2-5', 5-10', 10-14', 14-20'): 24.3, 23.3, 23.1, 23.0, 22.7
µJy/beam, with 0, 1, 5, 116 and 77 pixels below -5 sigma (the field's
own image had 239 and 285 beyond 10'). The one peeled source (11.0
mJy/beam at 17.5') left 0.13 mJy/beam. Artefacts remain only around the
6.4 and 5.9 mJy sources at 13.6-13.7', outside the half-power radius
(12.6') and just inside the corrected area (14.2'): 34 and 69 pixels
below -5 sigma, deepest -0.23 mJy/beam. A Band 5 `--peel-above 5` would
take those out, and out of the image. Not changed: `run_targets.py`
takes a target's peeling as done once `peel_report.json` lists its
choice, so a changed `--peel-above` or `--keep-radius` silently keeps
the old peeling; the peel and target folders must be moved aside by
hand. To be fixed.

### 146. A field too faint for the rounds within the scan

Band 5, 2026-10-22, the eight remaining targets on SLURM. In every field
but one, the brightest pixel rose from round to round as the phases came
right (G0213-0251 16.6 -> 44.0 mJy/beam, G0218-0245 56 -> 109, G0222-0321
6.4 -> 10.4). G0231-0533 has one 12-minute scan, and its brightest source
is 3-4 mJy/beam (natural-weight floor 24.1 µJy/beam). It rose through its
per-scan phase rounds (3.27 -> 4.06) and then fell through the 1-minute
rounds and the delay round: 4.08, 3.54, 2.79, 2.02 mJy/beam. In the delay
round nearly every antenna's delay was refused as beyond 5 ns. The
cores-only model of the next round (6.5% of the brightest pixel, about
130 µJy/beam against a noise of 33) then sat at about 4 sigma. WSClean,
asked for one major iteration of up to 50 000 components per sub-image,
cleaned noise for more than 25 minutes without a line of output, and
would have left a model full of noise peaks.

Two changes. `experiment_wsclean_image.py`: a cores-only model
(`--core-model-fraction`) also stops at 5 sigma (`-auto-threshold 5`,
`CORE_MODEL_MIN_SIGMA`); WSClean stops at whichever comes first. The
brighter fields' cores thresholds were 7-15 sigma, so nothing changes
for them. `run_targets.py`: `--time-phase-solint none` leaves the rounds
within the scan out of the field stage, for a field too faint for them.
2 new tests; 1210 pass. Not built: a check that undoes a round after
which the field's brightest pixel has fallen; the rounds are kept
however they turn out. To do.

### 147. The cores-only amplitude round and resolved targets

G0212-0335 (Band 5) finished at 20:08 SAST on 2026-10-22, with nothing
peeled, at 25.5 µJy/beam against 22.9 for the data. Its deepest pixel was
-0.86 mJy/beam, 4.7" from the target's 16.3 mJy/beam peak: a bowl of
negatives around the resolved target (330 pixels below -10 sigma; within
30" the positive flux outweighed the negative only 3.2 to 1). Within 60"
of the target the image holds 27 mJy (39 mJy counting only pixels above
3 sigma). NVSS's 43.5 mJy at 1.4 GHz, which takes in all the extended
emission, and TGSS's 314 mJy at 150 MHz (spectral index -0.89) predict
48 mJy at 1250 MHz. So about 40% of the target's flux is missing.

The cause is the cores-only amplitude round of finding 131. Its first
solve put the central square at 1.03-1.38 and the arms at 0.33-0.89, the
far arm antennas lowest (E06 0.36/0.57, W05 0.33/0.38, S06 0.36/0.46).
The target's peak rose 84% at that round (9.0 -> 16.6 mJy/beam). An
antenna fault stands alone; a split by distance from the array's middle
is a model that holds the wrong share of a resolved source's flux on
short baselines against long ones. The amplitudes then make the source
as compact as the model. At Band 4's 4.4" the targets were compact
enough for this to stay small (G0212-0335's flux checked out at 79.6
mJy). At Band 5's 2" it is not. Across the targets, the ratio of the
central square's median amplitude to the arms' in that round, against
the peak's rise: G0212-0335 1.63 (+84%), G0220-0318 1.39 (+55%),
G0229-0556 1.17 (+34%), G0213-0251 1.13 (+15%), G0222-0321 1.09 (+8%),
G0216-0604 1.04, G0218-0245 0.99, G0224-0620 0.95. The last three are
fields where compact off-axis sources dominate.

Refusing the round would also refuse what it is for: on Band 4 it found
S02 at 1.59 and shrank the holes beside every source from -7.8% to -3%.
So `experiment_phase_offsets.py --detrend-radius F` sets each antenna's
median amplitude against log10 of its distance from the median position
of the antennas solved, and draws a Theil-Sen line, which single faulty
antennas cannot pull. When the inner and outer halves of the array differ
by more than F, it divides the line out of every solution (phases
untouched; the median antenna keeps its amplitude). It does this before
the amplitude limits are judged, so the far antennas are not flagged for
the model's bias. The report gives `radial_trend` (inner and outer
medians, ratio, the line per antenna). `selfcal_recipe.py` passes F =
0.10 (`AMPLITUDE_CORES_DETREND`) to the cores-only amplitude rounds only;
the whole-field rounds solve against usual cleans, which hold the
extended emission. The round's reading says when the run was divided out.

Predicted for G0212-0335 redone: the amplitude_cores1 reading reports a
run divided out (inner about 1.15, outer about 0.70); the target's peak
stays near its pre-amplitude 9-10 mJy/beam instead of 16.6; the bowl
beside it falls from -0.86 to about -0.3 mJy/beam or less; and the flux
within 60" rises from 27 mJy towards 45-50. 7 new tests (a GMRT-like
array with a model trend and one faulty antenna; the solve end to end);
1217 pass.

On data (G0212-0335 redone, 2026-10-22 20:30-23:10 SAST). Every cores-only
amplitude round of the five redone targets reported its run divided out
(G0212-0335 inner/outer halves 1.11/0.64, then 1.25/0.87 and 1.25/0.83;
G0220-0318 1.07/0.65, then 1.16/0.88). The run comes back each round,
because the cores model still lacks the extended emission, and each time
it is taken out instead of applied. Same stage (the first whole-field
model), before -> with the fix: beam 3.47" x 2.72" -> 2.58" x 2.04" (the
biased round had weighted the far arm antennas down, so the far arms
count again); noise 26.8 -> 30.7 µJy/beam; target peak 16.7 -> 15.1
mJy/beam; flux within 60" 31.1 -> 61.0 mJy; deepest within 30" -0.61 ->
-0.45 mJy/beam; pixels below -10 sigma within 30" 134 -> 8. The peak
prediction was wrong: the jump at the first cores round (9.0 -> 15.2
mJy/beam) stays, and at the restored resolution it is real; the bowl and
the flux came right. Open: 61 mJy is 27% above the 48 mJy predicted from
NVSS and TGSS; the Band 5 flux scale is to be checked on compact sources
against Band 4 before any flux is used.

### 148. Every amplitude solve, not only the cores-only one

All nine Band 5 targets finished overnight (2026-10-22/23). G0212-0335,
redone with finding 147, ended where it began: beam 3.55" x 2.49", flux
within 60" 26.4 mJy, the bowl beside the target back (-0.84 mJy/beam, 287
pixels below -10 sigma). At its first whole-field model it had been at
2.58" x 2.04", 61 mJy and 8 pixels. The whole-field amplitude rounds,
solved against usual cleans and not detrended, brought the bias back:
central square over arms 1.32, 1.18, 1.09, so 1.6 together, as much as
the cores round had had. The data noise fell meanwhile (28.6, 25.9, 23.0
µJy/beam), because the far arm antennas were weighted down again: a
better noise figure bought with resolution and flux. The bandpass round
ran the other way (0.93: E06 1.19, W05 1.12, W04 1.11). A usual clean of
a field whose amplitudes are wrong reproduces the wrong structure, as the
Band 4 notes had warned; at 2" no model of a field with a bright extended
source is complete enough.

The products of the central-square/arm ratios over the whole-field
amplitude and bandpass rounds of both self-calibration stages, per target:
G0212-0335 1.59, G0231-0533 1.41, G0218-0245 1.31 (1.30 in the first
round after peeling), G0216-0604 0.80, G0229-0556 1.11, G0224-0620 1.05,
G0213-0251 1.01, G0220-0318 0.98 (1.27 in the field stage, 0.76 after
peeling, cancelling), G0222-0321 0.97. Every target has rounds beyond 5%.

An antenna's gain has no reason to depend on its distance from the
array's middle. What does depend on baseline length, decorrelation for
one, is baseline-based, and antenna gains cannot carry it correctly
anyway. So `selfcal_recipe.py` now passes `--detrend-radius` to every
amplitude solve (calmode ap: cores-only, whole-field and time rounds) and
to the bandpass round, at 5% (`AMPLITUDE_DETREND`; `--amplitude-detrend
0` switches it off). `experiment_phase_offsets.py` divides the run out
of the bandpass round's amplitudes too, after each antenna's amplitude
has been divided by what all antennas share in each bin, and reports
`radial_trend` there. The round readings say when it was divided out.
4 new tests (every amplitude round and the bandpass pass it, phase rounds
do not, it can be switched off; the bandpass round end to end); 1221
pass. All nine targets to be redone from their field stage.

### 149. Dividing the run out moved the flux scale

Finding 148's rerun (2026-10-23) gave G0212-0335 everything predicted for
it: beam 2.50" x 1.98", 52.4 mJy within 60", deepest pixel beside the
target -0.40 mJy/beam, 13 pixels below -10 sigma (287 before). But its
Stokes V noise, which measures nothing but the data's scale, rose from 22
to 28 µJy/beam, and in G0216-0604's target stage it fell round by round
(17.0, 14.1, 11.9, 10.8) while the brightest peak fell with it (8.4, 7.1,
6.2, 5.8 mJy/beam).

`detrend_amplitudes` scaled every antenna to the median of the line
`radial_trend` drew, meaning to keep the median antenna's amplitude. The
antennas crowd at the array's middle, the line is drawn against the
logarithm of distance, and its median is not the antennas' median: on
G0212-0335's first cores-only solve the antennas' median was 1.07 and the
line's 0.85, so the median antenna came out at 0.83 and every visibility
was turned up by 1/0.83^2 = 1.44. The bandpass round did the same.

Per target, the product of 1/median^2 over the amplitude rounds of both
stages: G0220-0318 1.52, G0212-0335 1.44, G0231-0533 1.41, G0218-0245
1.17, G0213-0251 1.06, G0222-0321 0.97, G0229-0556 0.95, G0224-0620 0.88,
G0216-0604 0.83. Divided by it, the ratio of the Stokes V noise to what
the weights allow with natural weighting comes together: six targets at
1.41-1.47, against 0.87-2.21 before. G0216-0604 (1.05), G0229-0556
(1.83), G0231-0533 (1.30) and G0222-0321 (1.59) stay apart; G0216-0604
loses more per round than its gains explain and is still open.

Against NVSS, with this factor taken out: G0224-0620's brightest compact
field source (13.8 mJy at 1.4 GHz, 12.9 mJy/beam peak here) comes within
2% of NVSS scaled by (1.25/1.4)^-0.7. The phase calibrator 0323+055
(J0323+0534) was bootstrapped at 2.16 Jy against 2.9-3.0 Jy expected from
NVSS 2.79-2.91 Jy, FIRST 2.80 Jy (1.1" x 0.7", compact) and alpha -0.34
from TGSS: 25% low. A compact source that matches NVSS while the phase
calibrator does not points to the calibrator having faded since the
1990s, not to a wrong scale on the targets; to be checked on every field
after the rerun.

`detrend_amplitudes` now divides by the line only, and the new
`normalise_amplitudes` sets each hand's median amplitude, over the
unflagged solutions of every antenna, time and channel, to one. Every
amplitude solve passes through it, detrended or not, before its limits
are judged (CASA's own solnorm is undone by the line, and left G0216-0604
at a median of 1.03 per round); the bandpass round after its line is
divided out. The reports carry `amplitude_medians`, and a round whose
median was more than 2% from one says so in its reading ("median
amplitude 0.83/0.84 set back to 1"). 4 new tests, among them G0212-0335's
own 26 antennas; 1225 pass. All nine targets to be redone from their field
stage: the intermediate rounds are not kept, so no round can be resumed.

### 150. How low the Band 5 noise can go, and an image nearer natural weights

Asked for about 10 µJy/beam in Band 5 on the grounds that Band 4 reached
about 20, the data were set against the radiometer equation first. Each
Band 5 target has two scans of 12 minutes (G0231-0533 one, G0229-0556's
second 3.4 minutes), 26 antennas and about 330 MHz once the edges and
interference are out: SEFD 341 Jy (Tsys 75 K, 0.22 K/Jy) gives 13.7
µJy/beam with natural weights, as select_scans.py measured from the data
(13.3-14.1 for the two-scan targets, 17.6 for G0229-0556, 24.1 for
G0231-0533). Band 4 had two scans of 22 minutes and 27-28 antennas:
WSClean gridded 5.1-5.3 million visibilities per sub-band there against
2.2-3.2 million here, 0.64 million in the sub-band at 1156-1187 MHz (GPS
L5 and Galileo E5a). The flagging is not the cause: of the targets' 34%
before the post-calibration flagging, 26 points are the four excluded
antennas and 6 the band edges. Once finding 149 is redone the robust-0
images should come to about 21 µJy/beam for the two-scan targets, which
is what half Band 4's data give at its noise. 10 µJy/beam would need
about four times the time on each target at robust 0.

What can be had from these data is part of what robust 0 pays: its
Stokes V noise is 1.55 times the natural figure (21.2 against 13.7 on
G0222-0321), the GMRT's central square holding most of the
visibilities. `run_targets.py --extra-robust R ...` now images each
target's calibrated data again at those robust values, cleaned as the
recipe's final image is (scales 0 5 15 30, auto-threshold 1, 20 major
cycles; a test holds these to the recipe's defaults), corrected for the
primary beam, and tidied to what the correction needs. It goes in
`images/final_robust_<R>/` beside the robust-0 image, which stays the
target's product; the report carries each image's noise and restoring
beam, and the summary line names them. A target already finished makes
only the new images when given the option later. Expected for robust
0.5 before any is made: about 16-18 µJy/beam for the two-scan targets
and a beam near 3" instead of 2.3". 3 new tests; 1228 pass.

On the rerun (2026-10-23) every finished target met it: robust 0 at
20.6-21.8 µJy/beam for the two-scan targets, robust 0.5 at 17.5-18.2 with
beams of 2.6"-3.5" by 2.2"-2.4". G0212-0335: 20.5 and 17.5 µJy/beam, 36.1
and 35.2 mJy within 60", peak 10.1 mJy/beam, deepest pixel -0.28 and
-0.12 mJy/beam (10 and 0 pixels below -10 sigma): finding 149's
predictions (36 mJy, about 10 mJy/beam) held.

### 151. A peeled source's gains were handed to the whole field

`experiment_extract_source.py --peel` subtracted the source's model from
the data its recipe had calibrated towards the source, and gave the field
back from those: every gain solved towards the source, amplitudes and the
weights CASA moves with them, went to every other source and to the
noise. For a bright source with clean solutions this costs little; for a
faint one it is the whole field's scale. G0216-0604 in Band 5: its second
peeled source, 10 mJy/beam at 10', had amplitudes of 0.44 to 1.41 moving
by 56% within a scan (W04); over its three rounds the source's own peak
fell from 10.0 to 8.6 mJy/beam and the hole beside it deepened from -4%
to -16% of the peak, and the field's Stokes V noise fell from 16.5 to 13.4
µJy/beam across the peel while G0218-0245 (110 mJy) and G0224-0620 moved
by under 4%.

The peel now takes the source out of the data the recipe started from,
as those data see it: model times gain, the gain of each visibility being
the recipe's input over its output, which calibration makes exact. The
rest of the field keeps its own calibration and weights. The two copies
are checked to be the same data row for row (rows, antennas, times,
shapes); a visibility the calibration left at zero takes the model as it
is, and the share of those is printed. Where the recipe's input is gone
the old way is used, with a warning, and `peel.json` says whose
calibration the field was given back with (`calibration_kept`). 5 new
tests, among them the arithmetic on gains of 0.4 to 1.4; 1233 pass. The
five targets that peeled (G0213-0251, G0216-0604, G0218-0245,
G0224-0620, G0229-0556) to be redone from their peel; their fields stand.

Redone (2026-10-23/24): every `peel.json` says the field's calibration
was kept, the gains towards the sources at medians of 0.99-1.01, and
G0216-0604's target stage started at a Stokes V noise of 16.9 µJy/beam
instead of 13.4. Its compact sources kept their peaks across the peel
(3.98, 2.20, 5.88 mJy/beam before, 3.95, 2.20, 5.72 after).

### 152. Amplitude rounds in the faintest fields, and a phase-only test

The target stage of G0216-0604 still lost Stokes V noise round by round
(16.9, 15.8, 14.3, 13.5 µJy/beam), and its three compact sources lost 5
to 8% of their peaks across it. In G0231-0533's field (one 12-minute
scan, nothing brighter than 4 mJy/beam) the brightest compact source
went from 2.86 mJy/beam in the data as calibrated to 4.03 after the
phase and delay rounds, then down through the amplitude rounds to 3.45
(-14%; the next two -16% and -12%), while its flux within 15" stayed
within the aperture's noise (5.28 then 5.22 mJy). G0229-0556's compact
sources lost 6-13% of their peaks between its field and target images
and its beam narrowed from 2.02" to 1.73" in its minor axis. After
finding 149 every round's median amplitude is one, but in these fields
the same antennas come out high round after round (G0216-0604: E04
1.24, 1.13, 1.08; E06, S02, E05, C13), and CASA's weights follow the
gains: the weighted scale moved by 0.924, 0.948, 0.972 for the measured
0.935, 0.905, 0.944 of the noise. The rounds then fit the noise and the
model's shortcomings as much as the antennas. The brighter fields hold
their noise within 5%.

Against NVSS (aperture sums in the primary-beam-corrected images, alpha
-0.7) the sources compact at 2" in the bright fields come out at 0.87 to
1.05; most of the faint fields' sources are 25-40" across in NVSS and
resolved here, and were short of NVSS already before any amplitude
round, so this comparison cannot judge the faint fields to 10%.

`selfcal_recipe.py --phase-only` solves no amplitude: the whole-field
rounds solve phases only (`phase_field`, until 5 degrees, at most 3), and
the cores-only amplitude rounds and the bandpass round are left out;
`run_targets.py --phase-only` gives it to both stages and leaves out
`--amplitude-cores`. 6 new tests; 1239 pass. The test: G0231-0533 redone
with it. Expected before the run: its three compact sources at about
4.0, 1.55 and 1.2 mJy/beam (as in the last phase-only round; 3.45, 1.31,
1.05 now) and the robust-0 noise within 5% of the present 33.1 µJy/beam.

The test (2026-10-24, job 14468748): phase-only came out better on every
count, though the peaks fell short of the prediction. Robust 0: 32.5
µJy/beam against 33.1, beam 2.15" x 1.84" against 2.64" x 1.56", the three
compact sources at 3.65, 1.47 and 1.09 mJy/beam against 3.45, 1.31 and
1.05 (predicted about 4.0, 1.55 and 1.2: 5-9% too high; the reference had
been a shallow image made before the delay round). Robust 0.5: 26.4
against 28.2 µJy/beam, 2.79" x 1.94" against 3.29" x 1.87", peaks 3.75,
1.47, 1.07 against 3.61, 1.34, 1.02. The Stokes V noise stayed at 30.3
µJy/beam through every round. G0231-0533 keeps this result.

The other peeled targets' compact sources across their target stages
(field image, then target image): G0213-0251 2.26/2.22 and 4.21/4.20,
G0218-0245 0.76/0.81 and 0.77/0.80, G0224-0620 12.40/12.47, 1.35/1.55
and 0.81/0.84 mJy/beam; they keep their amplitude rounds. G0216-0604
(-5 to -8%) and G0229-0556 (-6 to -13%) lose: `run_targets.py
--phase-only-after-peel` solves their target stage in phase only and
leaves their field stage, whose bright sources carry the amplitudes, as
it is. 1 more test; 1240 pass.

Redone that way (2026-10-24): G0216-0604's three compact sources at
3.96, 2.21 and 5.79 mJy/beam (predicted 3.95, 2.20, 5.72; 3.71, 2.02,
5.46 with amplitudes); G0229-0556's twelve compact sources at 0.97-1.03
of their field-stage peaks (0.87-0.94 with amplitudes). Their Stokes V
noise held through the rounds (16.9 to 16.5 and 30.8 to 30.3 µJy/beam),
and their robust-0 images came to 22.3 µJy/beam (18.7 in the quietest
tenth; 19-21 predicted) and 31.9 (30-31 predicted). The lower figures of
the amplitude runs, 19.5 and 29.7, had come with the shift of weight.

### 153. The Band 5 flux scale against Band 4 and FIRST

Finished set (2026-10-25), robust 0 and 0.5 in µJy/beam: G0212-0335
20.6/17.5, G0213-0251 20.2/17.9, G0216-0604 22.3/22.6 (17.8 in the
quietest tenth), G0218-0245 21.2/18.3, G0220-0318 21.6/17.8, G0222-0321
21.8/18.2, G0224-0620 18.7/16.1, G0229-0556 31.9/26.1, G0231-0533
32.5/26.4. G0216-0604's noise is 28-36 µJy/beam in boxes of 3' around the
peeled 163 mJy (NVSS) source 11' out, against 18.5 in its quietest boxes.

Band 4 (710 MHz) against Band 5 (1250 MHz), 12" apertures in the
primary-beam-corrected images, sources compact at 2" (peak over aperture
above 0.6) and within 12': median spectral index -0.67 for the 22 above
2 mJy at 1250 MHz, -0.90 for all 62 (the fainter steepened by Band 5's
aperture noise and resolution); G0229-0556 and G0231-0533 -1.68 on 4 and
2 sources. Band 4 was made before findings 147-149 and is not a reference
by itself.

FIRST (1.4 GHz, 5", Fint) against 8" apertures in Band 5, where
(1250/1400)^alpha makes 1.04-1.12 the expectation for alpha -0.35 to
-1: the median over 46 sources is 1.02, but the bright compact sources
(FIRST at least 2.5 mJy, compact at 2") come to about 0.78, and those
within 3.3' of the pointing centres, where the primary beam does not
matter, to about 0.80 (G0212-0335 0.71, G0213-0251 0.82, G0218-0245
0.80, G0220-0318 0.91, G0216-0604 0.76). G0231-0533, self-calibrated in
phase only and so on the calibrators' scale alone, gives 0.50-0.60 on its
two bright compact sources. The field medians scatter from 0.55 to 1.35
(G0224-0620), more than variability over the decades since FIRST
explains for 5-8 sources a field. With the phase calibrator's bootstrap
(2.16 Jy against 2.9-3.0 Jy expected, 0.72-0.75), this points to the
targets' scale being about 25% low, not to 0323+055 having faded as
finding 149's notes supposed: a faded calibrator does not lower the
targets.

`diagnose_flux_scale.py` (new) reads the calibrated data of the
calibrators scan by scan and prints the vector average of Stokes I on
baselines up to 5 klambda (on which both are points), for the band and
four sub-bands, with 3C48 set against Perley-Butler 2017, and the median
of each scan's uncalibrated autocorrelations, which follow the system
temperature unless the correlator normalises them. 4 new tests; 1244
pass. Its reading decides where the scale is lost: 3C48 below its model
would put it in the calibration of the flux calibrator, 3C48 right and
0323+055 low in the transfer between them.

Its reading (2026-10-25, Band5.ms CORRECTED_DATA, baselines within 5
klambda): 3C48 at 17.89, 18.46, 19.84 and 20.66 Jy at 1400, 1300, 1200 and
1100 MHz, 1.10, 1.07, 1.08 and 1.05 of Perley-Butler 2017; 0323+055 at
2.26-2.35 Jy on six of its seven scans (2.14 Jy at 1400 MHz, 2.57 at 1100:
alpha -0.76), scan 10 corrupted (1.22 Jy, 0.33-3.05 across the band); the
data hold no autocorrelations. The calibration is consistent with itself:
3C48 comes out as its model, and 0323+055 as fluxscale made it.

The calibrator against surveys: RACS-low (887.5 MHz, April 2019) 4.05 Jy,
NVSS 2.79-2.91, FIRST 2.80 (1.1" x 0.7"), RACS-mid (1367.5 MHz, 28 July
2021, five months after this observation) 2.88 Jy. RACS-low to NVSS gives
alpha -0.77, as these data do within the band, and 25 years of surveys
agree to a few percent: it is not variable. At 1250 MHz it should be 3.08
Jy; these data give 2.30. The targets against RACS-mid (Duchesne et al.
2024, VizieR J/other/PASA/41.3; same epoch, so no variability), 12"
apertures, sources with Fpeak/Ftot above 0.8 at RACS-mid's 11": median
0.767 over 46 sources in 7 fields (G0222-0321 and G0229-0556 have none),
0.791 for the 20 of at least 4 mJy, 0.835 within 5' of the centres, field
medians 0.63 (G0231-0533) to 0.91 (G0213-0251), against 1.03-1.07 for
alpha -0.35 to -0.75. So the Band 5 flux scale is low by a factor of
about 1.3-1.35 on axis (0323+055 3.08/2.30 = 1.34; sources near the
centres 1.07/0.835 = 1.28), and perhaps 10% more at 8-12', which would be
the primary-beam model's. The cause is in the 3C48 scan or its transfer
to 0323+055: the calibration reproduces 3C48 and then finds 0323+055 25%
fainter than it is. Everything after the transfer is linear in it (the
self-calibration keeps the median amplitude), so the images and their
noise are all low by the same factor: robust 0 at about 21 µJy/beam here
is about 28 on the true scale.

### 154. A calibrator given its flux density outright

Chosen fix for finding 153: recalibrate with 0323+055 as an absolute
reference of its own. `config.yaml` key
`calibration_calibrator_flux_models` (parsed by `config.py`'s
`parse_calibrator_flux_models`) names calibrators with a stated Stokes I
flux density, reference frequency and spectral index, as
`FIELD=FLUX_JY@REFFREQ_MHZ[,SPIX[,CURVATURE]]`, several separated by `;`
(a flat string, which the fallback parser keeps when PyYAML is absent),
or as a mapping. It is set in `CALIBRATION_TASK_CONFIG` only when given,
so other runs keep their recorded recipe. In `calibration.py`,
`setjy_field` models such a flux calibrator with setjy's manual standard
(fluxdensity, spix, reffreq, scalebychan) and every other one as before;
the stated models enter the flux-model identity, so a change of them
invalidates setjy and everything after it. Named as
`secondary_flux_calibrator`, the phase calibrator is then a reference
field, fluxscale has nothing to transfer, and amp.cal is published as
flux.cal: the targets take their gains from 0323+055 at its stated scale,
and 3C48 keeps Perley-Butler 2017 for itself and for the delay and
bandpass solves (its in-band shape is within 5% of the model). 7 new
tests; 1251 pass.

For Band 5: `secondary_flux_calibrator: 0323+055` and
`calibration_calibrator_flux_models: "0323+055=3.08@1250,-0.76"`
(RACS-mid 2.88 Jy at 1367.5 MHz, alpha -0.76 from these data and from
RACS-low to NVSS). Expected: fluxscale_report status not_required;
diagnose_flux_scale giving 0323+055 about 3.1-3.3 Jy and 3C48 as before;
the targets' fluxes and noise up by about 1.34 (robust 0 near 28
µJy/beam, G0212-0335 near 48 mJy within 60", as NVSS and TGSS
predicted), and their ratio to RACS-mid near 1.03.

### 155. A recipe changed in config.yaml was adopted as a code upgrade

With finding 154's two lines in `config.yaml`, `pipeline.py --to-stage
calibration` (2026-10-25 12:30 SAST) reported every stage "completed by an
older pipeline implementation; its unchanged science command and all
durable products were validated and adopted", and calibrated nothing.
`completed_stage_implementation_is_compatible` accepted a completed stage
whose recorded science-recipe fingerprint no longer matched as long as
its command line did; settings taken from `config.yaml` are not on the
command line. A stage with a recipe of its own (flagging, calibration,
self-calibration, imaging) whose recorded recipe differs from its
implementation fingerprint and from the current recipe is now refused:
"the recorded science-recipe fingerprint changed", and the pipeline then
asks for a new output directory, as it is meant to for a science change
of a stage that writes into the working Measurement Set. Records without
a recipe of their own, and legacy ones whose recipe was their
implementation fingerprint, are still compared by command. 4 new tests;
1255 pass. Band 5 is recalibrated in a new output directory,
`Band5_Products_0323`, from the LTA file.

Done (2026-10-25 18:35 to 2026-10-26 01:23 SAST, 5 h 43 min of
calibration): the flux-model provenance has 3C48 by setjy with
Perley-Butler 2017 and 0323+055 by setjy-manual at 3.08 Jy at 1250 MHz
(setjy's 2.751 Jy at 1450 MHz is that model there); fluxscale
not_required with both as references; every target's gains from
0323+055. diagnose_flux_scale on the new data: 0323+055 at 3.20-3.31 Jy
on its six good scans (1.04-1.08 of its model in each sub-band), 3C48
as before (1.05-1.10): the two calibrators now agree. Each scan is 1.41
times what the old calibration gave (3.08/2.16 = 1.43), and the targets'
natural-weight noise from their data rose by 1.42-1.43 in every field
(G0222-0321 13.7 to 19.5 µJy/beam): the images' units were corrected,
their signal to noise is as it was.

### 156. Band 4 is low by the same factor

Band 4 (710 MHz MFS, made before findings 147-149), 15" apertures in the
primary-beam-corrected images, against each source's flux at 710 MHz
from RACS-low (887.5 MHz, 2019) and RACS-mid (1367.5 MHz, 2021) with the
spectral index the two give (sources with Fpeak/Ftot above 0.85 in
RACS-low and 0.8 in RACS-mid, at least 4 mJy, alpha between -2 and 1;
peeled sources left out): median 0.761 over 28 sources in 9 fields, 0.759
within 10', fields 0.72-0.93. Band 5 before finding 154 gave 0.767
against RACS-mid. Two observations calibrated from 3C48 coming out a
quarter low alike points to something systematic in the 3C48 scans or in
their path through conversion and calibration, not to chance; Band 4
wants the same recalibration, through its phase calibrator at its survey
flux.

The Band 4 model (2026-10-28). Band 4 (observed 2020-11-22, GWB 550-950
MHz, 4096 x 97.7 kHz) has 3C48 at the start, 3C147 at the end and
0323+055 between the targets. No survey measures 0323+055 inside 550-850
MHz; the nearest are MRC at 408 MHz (7.13 Jy, NED), RACS-low at 887.5 MHz
(4.05 Jy) and RACS-mid at 1367.5 MHz (2.88 Jy), with NVSS 2.79 at 1.4 GHz.
A log-parabola through those four gives 4.86 Jy at 700 MHz and a local
index of -0.72 at 550 MHz to -0.77 at 850 MHz; the power law 4.86 Jy at
700 MHz, alpha -0.75, departs from it by under 0.5% across the band, so
the curvature is left out: `calibration_calibrator_flux_models:
"0323+055=4.86@700,-0.75"`. The Band 5 model carried down to 700 MHz
gives 4.79 Jy (1.5% less), so the two bands' models agree. NVSS (1990s),
FIRST and RACS-mid (2021) agree at 1.4 GHz within 4%, so the source is
taken as steady over the epochs used. No code change: `ms`, `outdir` and
the flux-model line of config.yaml are set for Band 4 and the calibration
is run from Band4.fits into Band4_Products_0323. The observing log
records a source-catalogue error that spoiled scans 2 and 3 (the first
0323+055 scan and G0212-0335's first); the scan table's positions are
right and G0213-0251 was added last to it. Scan 2's solutions can reach
only scan 3, which select_scans already leaves out. Expected: 0323+055
at 1.04-1.08 of its model in diagnose_flux_scale, as in Band 5; the
targets up by about 1.31 (1/0.761); Band 4 robust-0 noise on the true
scale about 21-32 µJy/beam (16-24 before).

### 157. Images for extended emission: tapers, not natural weights

The recalibrated Band 5 targets (Band5_targets_0323, finished
2026-10-26): robust 0 at 26.7-31.7 µJy/beam for the two-scan targets
(44.0 G0229-0556, 46.2 G0231-0533), robust 0.5 at 22.8-27.3 (36.0,
37.7). The robust-2 images all diverged: 0.7 mJy/beam to 2.4e11 µJy/beam
of "noise", beams such as 14.4" x 1.93" (G0212-0335); untapered, natural
weights give the central square nearly all the weight and a beam with a
plateau that the clean cannot follow. `run_targets.py` now refuses
`--extra-robust` above 1 (`MOST_UNTAPERED_ROBUST`).

The sources at the pointing centres of G0216-0604, G0220-0318,
G0222-0321 and G0231-0533 were not seen. In NVSS they are 10.5, 4.0, 3.9
and 3.9 mJy over 65, 73, 96 and 136": in a 2.3" beam a few µJy each,
below any noise these data reach; Band 4 shows them with three times
the beam area and at a frequency where they are brighter. What they want
is surface brightness, a larger beam. `experiment_wsclean_image.py
--taper ARCSEC` weights the visibilities with WSClean's
`-taper-gaussian` (both I and V) and names the images for it
(`image_label`: robust_0_taper15). `run_targets.py --extra-taper ARCSEC
...` makes, for each target, images tapered to those beams at robust 0
on a grid of five pixels to the beam over the same field
(`taper_grid`: 15" gives 3" cells, 1024 pixels for the 6144 x 0.5"
field), so that the cleaning scales of 0, 5, 15 and 30 pixels reach six
beams; each is corrected for the primary beam and recorded like the
other extra images. The report now lists the extra images a run was
asked for. 6 new tests. The MeerKAT beam and WSClean test modules no
longer import on this machine (Anaconda's numpy became 1.26.4 on
2026-10-27 under astropy 8.0, which needs numpy 2), so 1186 tests ran;
all pass.

Run for all nine with `--extra-taper 6 15 30`. Expected: about 40-60
µJy/beam at 15" for the two-scan targets; G0216-0604's centre clear at
15", G0220-0318's and G0222-0321's at about 5-8 times the noise per beam,
G0231-0533's marginal at 15" and better at 30".

Made (2026-10-28, all nine by 11:58 SAST). The noise was four times what
I expected: at 6" (beams about 8" x 5") 51-116 µJy/beam, at 15" (16-25"
x 12-15") 159-295, at 30" (40-44" x 24-27") 325-578, and Stokes V 38-96,
85-202 and 132-406: the short baselines a taper favours are noisier
than the data's weights say. The centres, in primary-beam-corrected
images, against the same aperture on 24 empty places 4-8.5' out:
G0216-0604 found, 15 times the local noise at the peak at 15" and 10.1
+- 1.8 mJy within 54" (NVSS 10.5 mJy at 1.4 GHz); G0220-0318 found, 40
times the noise at 6" (60 at robust 0.5) and 4.4-5.3 mJy (NVSS 4.0);
G0222-0321 and G0231-0533 show a core at 4-10 times the noise, and their
extended emission (NVSS 3.9 mJy over 96" and 136") is not recovered: an
aperture of 72-97" is uncertain by 4-5 mJy here.

`diagnose_baseline_noise.py` (new) takes channel differences of Stokes I
and V, which keep the noise and lose the sky, scales them by the noise
the weights predict, and gives their scatter for ranges of baseline
length: one where the weights are right, more where a range carries
noise its weights do not know of. 1 new test; 1262 pass.

Its reading (2026-10-28): on G0222-0321 the scatter over the weights is
0.97-1.04 from 1 to 200 klambda and 1.09-1.29 below 1 klambda, on
G0231-0533 0.70-0.97 everywhere: the weights are right, the short
baselines are noisier per visibility (0.35 Jy against 0.16 for the
longest) and their weights say so. What limits the tapered images is how
few of them there are: 11% of the data lie within 4 klambda and 0.7%
within 1, so a 15" taper leans on about a fifth of the data, and its
Stokes V noise, 85-202 µJy/beam, is about what that gives. My 40-60
µJy/beam had assumed far more short baselines than 24 minutes give this
array. Neither reweighting nor flagging is a remedy. G0222-0321's and
G0231-0533's extended emission, about 0.15-0.3 mJy in a 15" beam, stays
below these data.

The tapered Stokes I images are, in their quietest tenth, 1.2-2 times as
noisy as Stokes V (15": G0216-0604 169 against 85 µJy/beam, G0222-0321
166 against 122): something besides the noise. The likeliest is the
sources in the primary beam's first sidelobe, 30-45' out, beyond the 51'
of the image, which a taper sees at their brightest and which are not
cleaned. `run_targets.py --taper-field-scale FACTOR` makes the tapered
images that many times as wide (at 3: 7680 x 1.2", 3072 x 3", 1536 x 6",
154'), in folders of their own (`final_robust_0_taper15_field3`); the
report records cell, size and field. 2 new tests; 1264 pass. Expected:
the 15" noise down by about a third, toward Stokes V.

Result (2026-10-28, all nine by 13:39 SAST). In the quietest tenth at 15",
G0216-0604 went from 169 to 101 µJy/beam and G0222-0321 from 166 to 127,
and Stokes I over Stokes V fell from 1.2-2 to 1.04-1.55 (6": 1.03-1.21;
30": 1.04-1.71, the most in G0212-0335 and G0216-0604). Image noise over
the nine: 6" 40-91, 15" 89-231, 30" 174-457 µJy/beam, against 51-116,
159-295 and 325-578. Stokes V also came out 10-20% lower, and the beams
differ a little (G0216-0604 at 15": 16.1 x 15.0" before, 16.6 x 12.6"
now), because Briggs weighting depends on the field's size through the
uv cell; so the ratio to V is the fairer measure of what the cleaning
gained. Around the four centres (3-8' out, corrected images) the noise
fell 20-47%. G0216-0604: peak 22 sigma at 15", 9.1 +/- 1.0 mJy at 30".
G0220-0318: 51 sigma at 6". G0222-0321: core only, 7.5 sigma at 6", its
72" aperture uncertain by 3.0-3.6 mJy. G0231-0533: no detection (peaks of
3-4.4 sigma; the 97" aperture is uncertain by 2.3-3.6 mJy). The aperture
fluxes moved by up to 2.7 mJy between the two sets, within their
uncertainties. To find G0222-0321's and G0231-0533's 4 mJy at 5 sigma, the
aperture uncertainty would have to fall to about 0.8 mJy, three to four
times lower than these data allow. Primary-beam correction still stops at
14.4' (beam 0.41): the wider field only gives the outer sources to the
clean.

### 158. A named secondary with two standards left no primary

Found while the Band 4 recalibration (finding 156) was flagging,
2026-10-28 16:05 SAST, before it reached calibration. Band 4 observed
both 3C48 (first scan) and 3C147 (last); inspection lists both as flux
and bandpass calibrators and records its strategy as 3C147 primary and
bandpass, 3C48 secondary. With `secondary_flux_calibrator: 0323+055` and
no `flux_calibrator`, `get_flux_calibrators` inferred the primary with
`require_unique_field` over both standards and would have stopped at
`print_summary`: "Multiple primary flux calibrators are available:
3C147, 3C48". Band 5 had 3C48 alone, so it never met this. Setting
`flux_calibrator: 3C147` in `config.yaml` would have answered it but cost
the flagging: `flux_calibrator` is among the flagging recipe's science
arguments, so finding 155's check would have refused the completed
flagging stage and asked for a new output directory. Now, with a named
secondary and more than one standard, the primary is inspection's
`primary_flux_calibrator` when that is still one of the vetted flux
calibrators, as `get_bandpass_calibrator` already does for the bandpass;
otherwise it fails as before. The named secondary replaces inspection's
own, so for Band 4 the references are 3C147 (Perley-Butler 2017, delay
and bandpass) and 0323+055 (setjy-manual), and 3C48 has no role. 4 new
tests; 1268 pass. `calibration.py` only; calibration reads it from disk
when the stage starts, and its recipe does not include the code, so the
running pipeline takes it up if it is in place before flagging ends (or
resumes with it after the stop: a code-only change to a failed attempt
with an unchanged recipe is accepted).

But 3C147 is the weaker choice for these data. The observing log:
"03:48:11 Started scan on 3C147", "03:53:00 W04 Tracking problem, BATCH
restarted", "04:00:00 ... observing session over". A bandpass from an
antenna off source takes on the primary beam's frequency dependence, and
the delay is solved on the same scan. 3C48 (scan 1, 10 min at the start)
is clean and is what Band 5 used. Choosing it cannot be done in
`Band4_Products_0323`: the manifest's `field_selection` must match on
every resume ("Output directory belongs to a different MS, field
selection, or pipeline configuration"), and on a resume without a
positional input omitted options are restored from it, so a changed
calibrator in `config.yaml` alone would be overridden there. In a new
output directory no code change is needed: `flux_calibrator: 3C48` and
`bandpass_calibrator: 3C48` with `secondary_flux_calibrator: 0323+055`
give references 3C48 and 0323+055 and bandpass 3C48 (the explicit
bandpass is needed: inspection's strategy, 3C147, otherwise wins the
bandpass tie).

### 159. A diverged clean was accepted as an image

Found in the Band 4 target run `Band4_targets_0323`, read 2026-10-30
09:10 SAST. G0213-0251's robust-0.5 image (one 22-minute scan, 5120
pixels, briggs 0.5, scales 0,5,15,30 px, sub-images of 2560) came out at
47.7 mJy/beam against 21.4 uJy/beam for the data: 2200 times the noise.
Its WSClean log shows the residual's noise at the start of each major
iteration going 118, 72.1, 58.6, 62.9, 109, 258, 769 uJy/beam, 1.92, 7.5
and 18.9 mJy/beam, with 22 "Multiscale clean diverged ... resetting"
warnings from the fourth iteration on, and the 30-pixel scale holding
components of -32 to +47 Jy per sub-image. WSClean's own reset of a
sub-image did not stop it. The run took 1 h 21 min of deconvolution and
WSClean exited 0, so `experiment_wsclean_image.py` recorded the attempt as
completed, `run_targets.py` primary-beam corrected it and listed it
beside the good images. The robust-0 image of the same data with the
same scales converged (34.3 uJy/beam). In 254 other cleans of more than
one major iteration in this run (the models of the self-calibration
rounds, the peel images and the final images, read from the cluster
2026-10-30) that noise never went above its lowest.

Now `run_command` takes a `watch`, given each line of output, and stops
the command (SIGTERM, then SIGKILL after 60 s) when the watch returns a
reason, raising `Diverged` (a `CommandFailed` that is not a crash). The
Stokes I clean is watched by `DivergenceWatch`: WSClean's "Estimated
standard deviation of background noise" (nJy, uJy, mJy or Jy) above
twice its lowest so far is divergence. On G0213 that stops the run at
258 uJy/beam, the sixth major iteration, about an hour earlier. The
clean is then made again with `calmer`: the largest scale dropped (0,5,15
here); with fewer than three scales left, without multi-scale, which is
the crash ladder's next rung and is not tried twice. When nothing is left
to drop the failure is raised and the target stops, so an image like this
is never accepted. The diverged attempt's log is kept as
`wsclean_<label>_I.diverged<N>.log`, the report's attempt says
`diverged: ...` and the reading says how many times and with what. Every
clean made through `experiment_wsclean_image.py` is covered: the
recipe's images, the peel images and the extra images. 12 new tests;
1280 pass. `experiment_wsclean_image.py` only.

Not explained: why robust 0.5 diverged where robust 0 did not. Robust
0.5 gives a broader beam with higher sidelobes for one scan's coverage,
and the components that ran away were all at 30 px. Removing the 30-pixel
scale should allow for that. It is a prediction, not yet tested on the
data.

On data (job 14533663, 2026-10-30 11:22-12:10 SAST): the first clean was
stopped at its fifth major iteration (118, 72.1, 58.6, 63.0, 124.6
uJy/beam), its log kept as `wsclean_robust_0p5_I.diverged1.log`; the clean
without the 30-pixel scale converged (118, 72.8, 49.1, 40.1, 37.3) to 37.5
uJy/beam, 31.1 in the quietest tenth, against 21.4 for the data.

### 160. What is left around the brightest sources of G0213-0251 and G0218-0245

Read from the cluster, 2026-10-30 13:30-14:30 SAST, after the Band 4
target run (robust 0: G0218 22.9 against 15.4 uJy/beam, G0213 34.3
against 22.5). No change to code.

**G0218-0245.** In the final residual the noise is 47 uJy/beam at the
317 mJy source (10.5 arcmin out, peeled), 42-44 within 100 pixels, 37 at
150-250 and 32 at 400-600, against 21-25 at the same distances from the
point mirrored through the centre. Displayed, it is spokes running out
from the source, the array's fan beams. Its mirror-symmetric and
antisymmetric parts about the source are alike (33 and 27-37
uJy/beam): amplitude and phase errors both. The source itself is taken
out exactly: +60 uJy/beam at its position in the joined image (0.02%),
and in the twelve sub-band images 0.02-0.74 mJy (at most 0.2%) except the
top one (-2.7 mJy, 0.8%). In the sub-band images the noise around it is
177-262 uJy/beam against 60-84 at the mirrored point from 581 to 791 MHz,
and 302-357 and 713-818 against 106-125 and 186-240 in the top two. The
peel's own image of the source (1024 pixels) is at 43.3 against 24.1 for
the data, and its rounds towards the source moved nothing (amplitudes
0.99-1.01, phases within 1 degree, 2 solutions per antenna; bandpass
5.5% rms at most). The field-wide part of the excess (quietest tenth
18.2 against 15.4) is a cross-hatched ripple over the whole image.

**G0213-0251.** The excess is in the inner field rather than about one
source: 39-49 uJy/beam within 600 pixels of the 584 mJy extended source
(9 arcmin out) but also 35-37 at the mirrored point, 40-52 about the
103 mJy source (24 arcmin) against 30-33 mirrored, and 26-28 at the
edges, against 22.5 for the data. Displayed: long straight streaks across
the whole image and spokes from the 103 mJy source. The peel's image of
the extended source is at 53.8 against 33.2 for the data and its rounds
moved nothing (amplitudes 0.97-1.02).

**Ruled out before.** A gain every minute for the whole field (finding
122) and twice the sub-bands (finding 123) left the rings about sources
where they were. The errors are not in the band-averaged or sub-band
flux of the sources, and per-scan gains and a gain every 30 channels
towards them find nothing.

**Tests set (no code change, in `peel_tests/`).** The extraction of each
source again as in the run, then its recipe with one change: for G0218's
source a last round of amplitude and phase every 2 minutes towards it
(`--time-solint 2min`); for G0213's, scales to 120 pixels
(`--scales 0 5 15 30 60 120`), and separately the 2-minute round. The
"before" image of each repeats the run's (43.1 and 54.0 uJy/beam) and is
the control. If neither moves its image by more than about 10%, what is
left is not an antenna's gain in that direction at any time scale the
data support, nor the clean of an extended source: baseline-based errors
are then the candidate.

**Results (jobs 14535309 and 14535310, 30 October, read 3 November).** The
controls repeated: 43.0 (G0218) and 54.1/53.5 (G0213) uJy/beam. The
2-minute round found amplitudes moving within the scan (up to 21% rms on
W04 towards G0218's source, 11% on C12 towards G0213's) and gave final
images of 40.0 against 22.6 for the data (G0218, was 43.0 against 24.7)
and 49.2 against 29.8 (G0213, was 54.1 against 34.0): the image and the
data's noise fell together, the ratio did not (1.74 to 1.77, 1.59 to
1.65). Scales to 120 pixels: 52.6 against 33.3 (1.57), no change; my
prediction of 45-50 was wrong. Neither moved the excess, so by the rule
above baseline-based errors are the next candidate for Band 4.

### 161. Band 5 rerun: manual flags for the start of 3C48 and 0323+055 scan 10

Read from the 2026-10-25 Band 5 calibration (`Band5_Products_0323`) on
3 November. Every table passed its check, but the log of the amplitude
solves (120 s) shows 44-46 of 52 solutions flagged for SNR < 3 in each of
the first five intervals of the 27-minute 3C48 scan (07:37-07:45 UTC),
then 1 of 52, then none: a 16-20 Jy source cannot give SNR < 3 on good
data. Only C00, C06, C08 and S04 solved throughout. The 10-second phase
solves over the same minutes lost 1-7 of about 50, so the phases held and
the amplitudes did not (an array still being set up at the start of the
observation is the likely cause). 32% of 3C48's amplitude solutions were
flagged against 9% of 0323+055's, and round-1 applycal took 3C48 from
13.6% to 58.9% flagged. The bandpass is solved over the whole scan, these
minutes included. The second half of 0323+055 scan 10 (09:54-10:01 UTC)
lost 15-35 of 52 amplitude solutions, phase solves there missed whole
antennas, and its calibrated 1300 MHz flux was half the other scans'
(1.45 against 3.1 Jy); scan 11, of 0323+055 again, follows at once.

`flagging.py` and `config.py`: `flagging_manual_commands` in config.yaml,
a list of flagdata commands in CASA's list syntax (`mode='manual'` only),
applied after the edge channels and before the automatic flaggers. Each
command is parsed before anything is flagged, counted before and after,
and reported in `flagging/manual_flag_report.json`; one that selects no
data or more than a quarter of the observation stops the stage instead.
Absent from FLAGGING_CONFIG unless given, so earlier runs' flagging
recipes are unchanged. `diagnose_flux_scale.py`: `--uv-bins` (each
calibrator's Stokes I by baseline length and sub-band, all scans
together), `--model` (a stated model for 0323+055), `--time-profile`
(each antenna's uncorrected amplitude minute by minute, every flag
ignored, against its own second half). Both calibrators came out above
their models on the central square (<5 kilolambda) by 5-10% (3C48) and
3-9% (0323+055), rising with frequency, against 0-1% in Band 4: the
profile tests whether they are resolved on the longer baselines, in which
case their point models bias the antenna gains. 16 tests
(`test_manual_flags.py`, `test_diagnose_flux_scale.py`); 1296 pass.

**Measured (3 November, 13:28 SAST, on the 21 October calibration
`Band5_Products/ingest/Band5.ms`, `flux_scale/Band5_profiles.json`).**
Time profiles, uncorrected data, flags ignored. 3C48 (scan 1, from 07:36:31
UTC): C00, C06, C08 and S04 at 0.01-0.07 of their later level, that is
without signal, for the first 11 minutes; at 07:47:30 every antenna stepped
up by 5-15% (C08 2.8, C04 2.5, C06 1.9) and settled by 07:52:30; C04 was at
2-3.5 times its later level until then; E03 lost its signal from about 07:41
until 07:54:30. So the first 16 of the scan's 27 minutes are an array still
being set up, and the bandpass of both earlier calibrations was solved over
them. 0323+055 scan 10 (from 09:49:40): every working antenna at 5-10 times
its later level for 5 minutes, then noise alone from 09:54:40 to the end of
the scan (C03, C05, C12 and W06 are at 1 throughout: they are dead, as their
quarantine said). Scan 11 (from 10:08:49) is steady. Manual flags for the
rerun: 3C48 07:36:00-07:54:00, E03 on 3C48 to 07:55:30, scan 10 from
09:54:30.

Against baseline length (all scans but 10): neither calibrator is resolved.
Beyond 5 kilolambda both are flat to 1-2% out to 130 kilolambda (3C48 1.00
of PB2017 at 1100-1300 MHz, 1.035 at 1400). The excess sits only on the
central square: under 5 kilolambda 3C48 is 1.064 and 0323+055 1.051 times
its own longer-baseline level, alike in every sub-band. Proportional to the
source (1.2 and 0.12 Jy), so not additive interference; confined to the
shortest baselines, so not antenna-based (an antenna's gain would carry it
to every baseline it is on): a baseline-based amplitude error of about 5-6%
on central-square baselines, which no antenna-based calibration or
calibrator model can remove. It is the kind of error finding 160 named as
the next candidate for Band 4. The calibration model is left as it was.

### 162. A dead antenna past the 5% bandpass limit stopped the Band 5 rerun

`Band5_Products_1103`, 3 November: flagging 13:44-14:20 SAST (the three
manual commands took 6.40, 0.04 and 2.40% of the observation, as planned),
then calibration failed at 14:28 SAST at the bandpass check: C03 90.0% and
C05 90.5% of their solutions flagged, over the per-antenna maximum of 80%.
W06 (delay) and C12 (pre-bandpass) had been excluded as before. Both are
dead (finding 161: noise alone through scans 1 and 10). On the 27-minute
3C48 scan of the earlier runs they kept 0 and 1% of their channels and fell
under the fixed exclusion limit of 5% usable (finding 142); on the scan's
9.7 settled minutes noise alone passed SNR 3 in 10% of their channels,
between that limit and the check's, and the run stopped. `calibration.py`:
the bandpass excludes an antenna whenever it alone fails the per-antenna
check, at or below `1 - max_caltable_group_flagged_fraction` usable (20% by
default), under the same limits on how many are excluded; the quarantine
record says what share counted. Not part of the calibration recipe, so the
failed stage continues from its delay and pre-bandpass checkpoints. 4 tests
(`test_bandpass_quarantine.py`); 1300 pass.

### 163. Band 5 rerun results, and what robust 1 does to Band 5

Calibration (attempt 2, 14:49-18:28 SAST, 218.5 min): every table passed;
amplitude solves with no SNR failures (2 failed solutions bridged for the
targets, against 294 in the 2026-10-25 calibration); amp/flux tables 18.9%
flagged (33.9%), the six dead or absent antennas being 18.75%. Targets in
`Band5_targets_1103` (jobs 14547621-34, 18:48 SAST to about 01:30 SAST),
with the same options as the 26 October run (G0216/G0229
`--phase-only-after-peel`, G0231 `--phase-only --time-phase-solint none`)
and `--extra-robust 0.5 1`. Robust 0 / its data noise, robust 0.5 / data and
beam (uJy/beam): G0212 28.9/27.3, 24.2/22.3 3.43x2.14"; G0213 30.1/28.1,
26.6/24.8; G0216 31.5/23.3, 32.5/19.2; G0218 29.6/27.4, 24.9/22.6; G0220
30.3/29.1, 25.6/24.1; G0222 31.3/29.5, 27.6/25.7; G0224 27.4/26.7,
23.3/22.4; G0229 41.9/39.8, 34.6/32.5; G0231 46.1/43.1, 37.4/35.1. Against
26 October the robust-0 images are 4-6% lower in six fields (G0212, G0220,
G0222 at the field stage, G0224 at the field stage, G0229) and alike in the
rest; the natural-weight estimates of the data (19.0-20.1) are unchanged.

Robust 1, untapered: the data's own noise is only 1-4% below robust 0.5
(G0212 21.8 against 22.3, G0220 23.3/24.1, G0224 22.1/22.4, G0231 34.2/35.1)
while the beam stretches (G0212 8.14x1.80", G0229 5.27x2.33", G0231
4.00x1.94"), and the clean is unstable at the 30-pixel scale: G0213 and
G0229 diverged and were stopped and retried without it (finding 159
working on Band 5), G0231 drifted 1.2 times above its lowest and was kept.
Images: better than robust 0.5 only on G0220 (24.2 against 25.6), alike on
G0229, worse on the rest (G0212 32.4, G0218 26.1, G0222 29.3, G0224 29.0,
G0231 71.5). My prediction of 20-23 uJy/beam was wrong. So robust 0.5 is
within a few per cent of the most sensitive weighting these data allow
and keeps a usable beam: 23-28 uJy/beam for the two-scan fields is the
floor, not ~20.

Observed, not yet tested: in every field that kept amplitude
self-calibration after peeling, the target stage raised the data's own
noise against the field stage (G0224 24.4 to 26.7, G0222 28.1 to 29.5,
G0213 27.4 to 28.1, G0220 28.6 to 29.1, G0218 27.2 to 27.4), while
G0229, phase only after the peel, kept 39.8. Once the bright sources are
out the field may be too faint for amplitude solutions, whose scatter then
reaches the weights. Test: the target stage of G0224 and G0222 again with
`--phase-only-after-peel` in a copy, comparing noise and the compact
sources' peaks (finding 152's measure).

Robust -0.3 (`--extra-robust 0.5 -0.3`, jobs 14552198-211, 4 November
04:16-06:28 SAST; image (quietest tenth) / data, beam): G0212 34.1 (33.4) /
32.5, 2.41x1.66"; G0213 34.8 (33.6) / 32.7, 2.33x1.66"; G0216 34.1 (30.0) /
27.4, 2.21x1.67"; G0218 35.3 (33.6) / 32.7, 2.32x1.62"; G0220 35.8 (35.0) /
34.3, 2.40x1.62"; G0222 37.0 (35.1) / 34.8, 2.39x1.65"; G0224 32.2 (31.5) /
31.4, 2.24x1.66"; G0229 50.6 (48.9) / 48.2, 2.33x1.80"; G0231 58.6 (54.8) /
53.2, 2.16x1.64". 1.16-1.21 times the robust-0 noise (G0231 1.27, G0216
1.08), as predicted (1.15-1.25); images 1.03-1.10 times their data noise,
and G0216 1.24 against 1.35 at robust 0 and 1.69 at robust 0.5, so less
weight on the central square leaves its bright source's residuals
relatively smaller. The beam narrowed on its minor axis (1.62-1.80"
against my 1.8") more than on its major (2.16-2.41" against my 2.0").
G0216 robust 1 diverged once and, retried, gave 34.9 uJy/beam against
18.3 for the data.

### 164. MeerKAT's gain calibrator was read as a second flux standard

First MeerKAT data through ORBIT: capture block 1708779077 (TDE J2344
obs 2, L band, 2024-02-24, 63 antennas; 63 min with J1939-6342 as flux
and bandpass calibrator, J2341-3506 as gain calibrator and J2344 as the
target), converted on ilifu with katdal's mvftoms at full resolution and
averaged to 1024 channels with CASA mstransform (`MeerKAT/1708779077.ms`,
29 GB; the kernsuite container's katdal cannot average channels itself,
finding of 2026-11-11). Inspection (16:33 SAST) found the telescope, band
(MEERKAT_L), 63 antennas, 7 scans and the roles from the intents, but
made J2341-3506 the phase calibrator and a second flux calibrator at
once. katdal tags a gain calibrator `CALIBRATE_PHASE,CALIBRATE_AMPLI`,
and `classify_fields` counted every `CALIBRATE_AMPLI` field as a flux
calibrator (the ALMA and VLA reading). The MeerKAT flux policy then
refuses the second one at setjy ("not one of the built-in
Stevens-Reynolds primary aliases"), so calibration would have stopped
after flagging. `inspect.py`: a field with both intents beside a distinct
`CALIBRATE_FLUX` field is the gain calibrator; `CALIBRATE_AMPLI` alone,
or with no flux field anywhere, still marks the flux scale. The
inspection report changes, so the run restarts in a fresh output folder.
4 tests (`test_inspection_intents.py`); 1304 pass.

**Rerun with the fix (16:49-17:43 SAST):** inspection names J1939-6342
(flux, bandpass), J2341-3506 (phase) and J2344 (target) with no warning.
Flagging in 11 minutes: 2.4% online flags, 10.1% outside 900-1670 MHz,
31.5% L-band masks, tfcrop 0.22% and extension 0.23%, 42.9% in all.
Calibration in 43 minutes: Stevens-Reynolds 2016 on J1939-6342; the
T phase solve failed its check on m000 and G was used; every table passed
with no flagged solutions (the bandpass's 41.6% are the masked channels);
fluxscale gives J2341-3506 2.010 +- 0.004 Jy at 1284 MHz; phase transfer
over 57-minute gaps a median 11 deg, 90th percentile 14 deg; post-cal
flagging 0.56%; statwt is off for MeerKAT by profile. First MeerKAT
calibration through ORBIT, with no change to the calibration code.

### 165. The MeerKAT MFS beam ignored WSClean's weights and the RFI masks

First MeerKAT image through ORBIT (1708779077, 18:47-20:02 SAST on
2026-11-11): WSClean 3.6 from the ilifu container, 6144 pixels of 1.5",
robust -0.5, eight joined sub-bands, 73 minutes; katbeam
MKAT-AA-L-JIM-2020 0.1 applied in 37 s from a venv named by `python:` in
config.yaml, because the Jupyter Python has no katbeam. Beam 6.13" x 4.89".
The apparent-flux image has 13.4 uJy/beam in the median 12.8' box and
11.4 in the quietest tenth. It has 22-24 within 20' of the centre, where
calibration errors around bright sources dominate before
self-calibration, and 0.1-1.5 mJy/beam within 12' of J2341-3506.

J2341-3506 lies 34.6' from the pointing centre. The corrected MFS image
gave it 1.63 Jy (peak) and 1.79 Jy (integrated), against fluxscale's
2.010 Jy. The seven corrected sub-band images fit a power law of
2.07 Jy at 1284 MHz with alpha -1.08 (3.5% rms), so katbeam and the flux
scale are right and the MFS beam is not. WSClean's MFS image is the mean
of the sub-band images weighted by their imaging weights (WSCIMGWG); that
weighted mean reproduces its J2341-3506 flux to 0.6%. Under Briggs -0.5
the weights rise steeply with frequency: the 1444 MHz sub-band carried
2.8 times the 1016 MHz one with the same visibilities. The 1177-1284 MHz
sub-band, inside the 1163-1299 MHz mask, carried none, and its image is
NaN. `correct_products` instead averaged the sub-band beams by their
coverage of 900-1670 MHz alone: equal weights, the empty sub-band
included, each at its nominal centre (1551 MHz for a sub-band whose data
lie at 1498-1524 MHz). At J2341-3506 that beam is 0.486 where WSClean's
weights give 0.426. The corrected MFS image was therefore 12% low there
for a flat spectrum, and by a Gaussian scaling of the measured beams 4%
at 20' and about 18% at the 0.3 cut-off. The 1498-1605 MHz sub-band's
own correction was 6% high at 35'.

`meerkat_beam.py`: the MFS beam is weighted by each sub-band's WSCIMGWG,
or by coverage when a header lacks one. Each sub-band is evaluated at the
centre of what survives the science range and the masks the flagging
stage applied, which `wsclean_imaging.py` reads from
`flagging/known_frequency_mask_report.json`. A sub-band with no weight,
no unflagged frequency or no finite pixel is skipped, and products an
earlier run left for it are removed. The record gives the weighting, the
masks and the MFS image's effective frequency, 1400 MHz here against the
1284 MHz in its header. Provenance schema 2. `pipeline.py` adds the scheme
and the mask report to the MeerKAT WSClean imaging recipe and refuses
schema-1 products. A completed WSClean imaging stage whose recipe changed
now reruns in place instead of demanding a new --outdir: its only
checkpoint binds WSClean's outputs to the command, MS and data column, so
only the correction is redone. CASA imaging still fails closed.
20 tests (`test_meerkat_beam.py`, `test_meerkat_wsclean.py`,
`test_pipeline_dependency_safety.py`); 1324 pass. The old weighting and
mask-blind frequencies each fail them.

Expected from the rerun: J2341-3506 at 1.86 Jy (peak) and 2.04 Jy
(integrated) in the corrected MFS image. What remains above its weighted
true flux of 1.96 Jy is its spectral index, which the MFS correction does
not remove, by design. The pass should take a minute or two.

**Rerun (12:20-12:23 SAST, 2026-11-12):** 2.6 minutes, WSClean not run
again. J2341-3506 1.86 Jy peak, 2.04 Jy integrated; effective frequency
1400.0 MHz; the 1177-1284 MHz sub-band skipped; sub-bands at 920.1,
1116.4, 1344.9, 1510.9 and 1650.0 MHz, as predicted.

### 166. The per-target chain, with peeling, made to run MeerKAT fields

The beam-corrected image of 1708779077 is 23 uJy/beam at the centre at
robust -0.5, twice the floor of the quietest parts (11.4). A residual
image rendered at 8x8 binning shows why: J2341-3506, 0.8 Jy apparent and
35' out, throws concentric ripples over the whole 2.56-degree field,
across the pointing centre too. Those errors are direction-dependent
(its beam response rotates with the sky, and 8 s and 836 kHz smear it by
up to 8% on the longest baselines there), so whole-field
self-calibration cannot remove them and peeling can. The chain that did
this for Band 5 (`run_targets.py`: scans, self-calibrated field, bright
off-axis sources peeled with their own gains, target, primary beam) was
uGMRT-only in five places, now telescope-aware:

- `run_targets.py` takes the telescope and band from the products'
  inspection and calibration reports. A MeerKAT field gets the band's
  grid (6144 x 1.5" at L), the channels inside the guaranteed range from
  the flagging stage's record (53-973 here), and reference antennas
  starting with the calibration's (m059) followed by MeerKAT's preferred
  ones, all overridable. A band without a katbeam model and products
  that contradict `--telescope` are refused. uGMRT keeps 5120 x 1", B4
  and its required `--channel-range`.
- The reference antennas reach every solve: the field and target
  recipes, and through `peel_sources.py --refant` each source's recipe.
  The recipe's own default names uGMRT antennas, which CASA would refuse.
- `correct_primary_beam.py --telescope meerkat` applies katbeam through
  `meerkat_beam.py` exactly as the imaging stage does (finding 165), with
  the masks from `<products>/flagging/known_frequency_mask_report.json`.
  It writes `-MFS-image-pb.fits`, which `run_targets.py` now looks for.
- `experiment_wsclean_image.py` chooses the dirty image that gives the
  data's noise from the correlations the inspection recorded: Stokes V for
  circular hands or all four linear correlations, and none for XX and YY
  alone. WSClean cannot form V from those, and Q holds the beams'
  off-axis leakage of J2341-3506. `--noise-stokes` overrides the choice.
  The report keeps `stokes_v` for V and adds `data_noise`, and the
  recipe's and the chain's report lines omit the data's noise when there
  is none.

29 tests (`test_run_targets.py`, `test_correct_primary_beam.py`,
`test_experiment_wsclean_image.py`, `test_selfcal_recipe.py`,
`test_peel_sources.py`); 1353 pass. Restoring the uGMRT antennas, the
uGMRT beam step, an unconditional Stokes V image, or a peel that drops
`--refant` each fails them. A dry run against copies of 1708779077's
records gives the expected commands at every step. None of the five
files is part of the pipeline's stage fingerprints, so its resume state
is unaffected.

Predicted for J2344 with `--peel-above 30 --keep-radius 20`, at the
chain's robust 0: J2341-3506 and the 37 mJy source 34.7' out are peeled,
while the sources inside 20' stay in the field for its self-calibration.
The centre should reach 8-9 uJy/beam, with a beam of about 7-8". A few
times the noise will probably remain within about 5' of J2341-3506,
where the smearing is in the data. The run should take about a day,
with up to about 100 GB of working copies at a time.

### 167. One rule sizes every band's images; the chain reads each band

Only uGMRT Bands 4 and 5 and MeerKAT L have been run on real data. The
other bands carried grids that no run had checked. Band 4's grid was set
by a rule (its comment in `config.py`): enough pixels across the beam,
and a field reaching out to where the measured beam has fallen to about
a tenth at the band's lowest frequency (0.12 at 550 MHz for 5120 x 1").
The same rule, with four or five pixels across lambda/B of the longest
baselines at the band's centre (25 km for the GMRT, 7.7 km for MeerKAT),
is met by the runs that were made. Band 5 at 6144 x 0.5" ends where the
beam at 1050 MHz is 0.12, and MeerKAT L at 6144 x 1.5" where it is 0.14
at 900 MHz. MeerKAT's beam is approximated by the cosine taper of Mauch
et al. (2020) for this. The untested grids fail it:

| Band | Grid held | Fault |
|---|---|---|
| uGMRT B2 | 8192 x 1.0" | 13 pixels to the 13" beam; beam 0.80 at the edge at 125 MHz |
| uGMRT B3 | 8192 x 0.8" | 8 pixels; beam 0.54 at 250 MHz |
| uGMRT B5 (pipeline's own) | 8192 x 0.3" | beam 0.26 at 1050 MHz |
| MeerKAT UHF | 6144 x 2.0" | beam 0.25 at 580 MHz |
| MeerKAT S0 | 6144 x 0.6" | beam 0.33 at 1750 MHz |

`config.py` now gives B2 8192 x 3.0", B3 8192 x 1.5", B5 6144 x 0.5" (as
its targets were imaged), UHF 7168 x 2.0" and S0-S4 8192 x 0.6", for
both the self-calibration and the final images. B4 and L are unchanged.
No products of these bands exist, so no recipe of a finished stage
changes.

`run_targets.py` now reads uGMRT products as it reads MeerKAT's. The
band comes from the inspection report, the grid from one table that
agrees with `config.py`, and the channel range from what the flagging
stage recorded masking at the band's edges, divided by the targets'
channel averaging (`calibration/target_channel_average.json`). For
Band5_Products_1103 that gives 122-3974, the range given by hand then,
and 31-993 once averaged by four. A band without a measured beam
(MeerKAT S) runs only with `--no-primary-beam`, and its report says the
image is uncorrected. A band the products contradict is refused. The
uGMRT recipe keeps its own reference antennas. `run_target_slurm.sh`
documents MeerKAT, whose Python must have katbeam (`PIPELINE_PYTHON`).
The two DD tests that sat outside `tests/`, which failed when run
together, are now in the suite.

`tests/test_band_coverage.py` checks every uGMRT and MeerKAT band. Each
must have every stage's settings and the chain's grid. Each grid, for
the chain and both of the pipeline's images, must pass the rule. The old
grids must fail it, and a primary beam must be measured wherever one is
applied. 12 tests (`test_band_coverage.py`, `test_run_targets.py`); 1365
pass. Bands 2 and 3, UHF and S still need a first real run, which
`SCIENCE_VALIDATION.md` counts as their validation.

## Orchestration checks that came back clean

Worth recording, since these are the claims most worth doubting:

- MPI is restricted to the imaging stage with the CASA backend
  (`stage_uses_mpi`). The README's claim that self-calibration deliberately
  avoids MPI because its `tclean` calls use `savemodel='modelcolumn'` is
  actually enforced, not merely stated.
- Automatic science-target selection is fail-closed and order-independent: it
  requires finite, internally consistent CASA counts for every candidate,
  rejects the whole selection if any are malformed, applies a flagged-fraction
  ceiling, and records the full ranking with a hash of the statistics it read.
- The imaging recipe fingerprint includes `FINAL_IMAGING_CONFIG[band]` and
  `IMAGING_TASK_CONFIG`, so the `pblimit` and restored `robust` values from
  findings 10 and 23 correctly invalidate a resume rather than being reused
  silently.
- `preflight()` gates telescope/input combinations, GMRT flag-file usage,
  and raw-LOFAR FITS before any stage runs.

## Behavioural changes to be aware of

- In-band sub-band imaging is now opt-in for every telescope and profile.
  Previously the orchestrator enabled it by default for uGMRT runs. Enable it
  with `--make-subbands` (plus `--subbands N` for the count) or
  `make_subbands: true` in `config.yaml`; a bare `--subbands N` prints a note
  and does nothing. Continuum imaging, primary-beam correction and source
  detection are unaffected.

- `calibration/phase_scan.cal` is a new required product. Existing output
  directories will recompute the calibration stage; they cannot be resumed
  into the new recipe, which is correct — the applied table list changed.
- `refantmode='strict'` can fail a run that previously completed. That is the
  intended fail-closed outcome: it means the chosen reference antenna is not
  available throughout, and a better one should be selected.
- LOFAR imaging is slower by roughly the channel count in the deconvolution.
- Self-calibration is slower and more memory-hungry because w-projection is
  now actually being performed.
- A lower `pblimit` enlarges the region auto-multithresh searches, so imaging
  and self-calibration do modestly more work. With `gridder='wproject'` and
  the default `normtype='flatnoise'` the restored image is not divided by the
  primary beam, so a smaller `pblimit` cannot amplify edge pixels.
- uGMRT B2 runs (and B3 at the bottom of its band) will now print a warning
  that the image does not reach `pbmin`. That is a pre-existing property of
  the configured `imsize`/`cell`, not a regression; it was simply never
  reported before.
- Amplitude self-calibration now needs DR >= 100 rather than DR >= 30-50. It
  is off by default, so this only affects runs that opted in -- where it was
  previously being attempted on fields too poorly modelled to support it.
- A DD run whose imaging stage did not publish a measured beam now fails
  closed instead of quietly using CASA's generic beam. Re-run imaging, or use
  `--skip-imaging` and correct the DDFacet output separately.
- MeerKAT and LOFAR imaging now use their intended negative Briggs robust
  instead of +0.5. This changes resolution, beam size and RMS for those
  telescopes -- it is the point of the fix, but any MeerKAT image made with
  the repository as previously shipped is not comparable to one made now.
- Flagging now aborts on a Measurement Set whose POLARIZATION table does not
  provide both parallel hands of one feed basis, where it previously
  proceeded across all correlations.
- Reference-antenna recommendations may differ from previous runs. The
  recorded ranking shows both the profile choice and the least-flagged one.
- The inspection report is now schema_version 3. The added keys are additive
  (`field_intents`, `calibrator_catalogues`, `time_coverage`), but a
  Measurement Set carrying CALIBRATE_AMPLI intents will classify differently
  than before: that field moves from the phase list to the flux list.
- Inspection may emit new warnings about calibrator/target bracketing. They
  describe the observation, not a pipeline failure, and do not stop the run.
- Source detection now aborts when the imaging report names a product for the
  requested role but its FITS export is missing, where it previously produced
  a catalogue from a different image. Re-run the imaging stage's FITS export.
- FITS ingest now rejects an imported Measurement Set with an unusable
  spectral or polarization axis at ingest rather than several stages later.
- uGMRT B4 and MeerKAT UHF images that span the full digitizer output now
  carry a different reference frequency, so their primary-beam correction and
  spectral-index reference both shift. This is the point of finding 38, but it
  means such images are not directly comparable to ones made before it.
- A CASA task failure now surfaces its own error message rather than a
  compatibility-layer message. Some previously reported "does not support
  required science keyword" failures were really data or selection errors.
- Output directories from before this audit should not be reused: their
  products were made under a different recipe. Use a new output directory
  rather than forcing reuse. (This bullet previously credited that to finding
  40; see its correction.)
- A killed self-calibration now continues on its own. Products from the
  unfinished step are kept as `<name>.incomplete.<N>` and listed in the
  self-calibration report; they are not reused and can be deleted once
  reviewed.
- Each stage now prints one memory line before it starts (job limit, job
  usage, node available) and may print a warning when the node is short of
  memory. A SIGKILL message now names its likely cause and points to a JSON
  file under `diagnostics/`.
- Self-calibration now refuses to start when the process's memory limit is
  below the minimum its imaging needs, instead of being killed partway.
- A target whose stage fails no longer stops the remaining targets. The run
  still fails at the end, naming each target that did not complete.
- Every calibrated science target is now self-calibrated and imaged, each
  under `Products/targets/<name>/`. Previously one target was selected and
  imaged. Use `--target`, `--targets` or `--target-selection auto` for less.
- When auto-multithresh opens no mask, self-calibration and final imaging now
  clean the compact emission inside the primary beam to 6 sigma with a capped
  iteration count, instead of publishing the dirty map. Turn it off with
  `automask_fallback_enabled`.
- Self-calibration measures the coherence of the calibrated target data
  before imaging and reports it, and the relaxed-mask retry now fires at
  5 sigma instead of 8, which recovers fields whose auto-multithresh mask
  never opens.
- Self-calibration reports plainly when CLEAN found nothing to deconvolve:
  such an image is the dirty map, and self-calibration had no model to solve
  against. The residual pass is skipped in that case and can never fail the
  stage.
- Self-calibration flags residual interference before the final image, and
  calibration grows its corrected-data flags over neighbouring channels and
  times. Both change the science recipe, so they apply to stages run from now
  on; both can be turned off (`residual_flagging`, `postcal_extend`).
- Primary-beam and PB-corrected images are bound by the content of their data
  files, not their metadata or timestamps (findings 53 and 55). A beam whose
  pixels no longer match its record is rebuilt rather than stopping imaging,
  and PB correction and FITS export rerun after it.
- The self-calibration recipe fingerprint is new (schema 2 after finding 52).
  A failed self-calibration recorded under an earlier recipe is set aside as
  `selfcal.superseded.<N>` and restarted rather than resumed. Completed
  self-calibration stages are adopted as before. Superseded directories hold
  full MS copies; delete them once reviewed if disk space is tight.

## Verification performed

- `python3 -m py_compile` over every module.
- The full regression suite: `tests/` (209 tests, run from inside that
  directory), plus `test_ddcal` (5) and `test_pipeline_ddcal` (2). **216 tests,
  all passing.** The same suite was run against the pre-audit sources first to
  establish a 209/209 baseline, so every failure could be attributed. See
  finding 41.
- After finding 65: `tests/` 358, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **365 tests, all passing.**
- After finding 64: `tests/` 354, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **361 tests, all passing.**
- After finding 63: `tests/` 353, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **360 tests, all passing.**
- After finding 62: `tests/` 350, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **357 tests, all passing.**
- After finding 61: `tests/` 348, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **355 tests, all passing.**
- After finding 60: `tests/` 346, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **353 tests, all passing.** A failed imaging
  attempt whose recipe changed is now set aside and restarted, as a failed
  self-calibration already was (finding 52); imaging does not write to its
  input either.
- After finding 59: `tests/` 336, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **343 tests, all passing.**
- After finding 58: `tests/` 333, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **340 tests, all passing.**
- After findings 56 and 57: `tests/` 324, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **331 tests, all passing.** A failed imaging
  attempt now resumes when its recorded science recipe matches, as calibration
  and flagging already did, instead of only through the legacy audit that
  treats any newer configuration file as a science change.
- After finding 55: `tests/` 297, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **304 tests, all passing.** CASA 6.7.5-18 was
  also installed here as a standalone package to replay the primary-beam flow
  against real CASA; see finding 55.
- After finding 54: `tests/` 295, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **302 tests, all passing.**
- After finding 53: `tests/` 292, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **299 tests, all passing.**
- After finding 52: `tests/` 284, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **291 tests, all passing.**
- After finding 51: `tests/` 281, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **288 tests, all passing.**
- After finding 50: `tests/` 274, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **281 tests, all passing.**
- After findings 42-49: `tests/` 263, plus `test_ddcal` (5) and
  `test_pipeline_ddcal` (2). **270 tests, all passing.** Findings 46-49 added
  25 of them. The 11 orchestrator tests were also run against the previous
  `pipeline.py`: 8 failed, and the 3 that passed are fail-closed checks that
  should hold in both versions. The two truncated-column fixtures were run
  through the previous `model_column_status()`, which accepted both.
- Configuration wiring checked by loading `config.py` and asserting the new
  calibration keys and the per-band `wprojplanes` synchronisation.
- The uGMRT aperture arithmetic was computed from the repository's own beam
  model and imaging defaults, and `effective_aperture()` was extracted and run
  against real B2/B4/B5 geometry; its verdicts match that arithmetic.
- The DD beam-selection bug was reproduced against the real
  `product_candidates()` ordering before the fix, and all four post-fix
  branches (recorded provider, suffix fallback, fail-closed, genuine CASA
  provider) were exercised.
- The RA-wrap fix and the uncovered-source recording were checked on
  constructed inputs; the direction-selection JSON and the generated DD
  primary-beam script were produced end to end from the existing test fixture
  and the script compiles.
- The QA grouping rewrite was benchmarked and checked for identical group
  statistics against the previous implementation.
- The config.yaml override effect was measured by loading `config.py` with and
  without the shipped YAML and diffing every per-band imaging value; the fix
  was then re-measured, and the override mechanism re-checked by deliberately
  uncommenting a flat key.
- `recommend_refant()` was extracted and exercised on three constructed
  antenna-statistics cases.
- The inspection scan-coverage and bracketing logic was extracted and run
  against mocked `msmetadata` for a healthy calibrator cycle and for a
  calibrator that stops before the target does.
- Stage reordering was checked for a reversed list, an already-ordered list,
  and a list containing a duplicate.
- The source-detection substitution was reproduced before the fix and all
  three post-fix branches (missing export, present export, restored
  preference) were exercised, including the recorded role provenance.
- The recipe fingerprint was checked for stage and configuration coverage
  before and after the fix, and the self-calibration payload hash confirmed to
  change when one self-calibration setting changes.
- The compatibility layer was exercised on a real-CASA-style error, a genuine
  `unexpected keyword argument`, an end-to-end task failure, and a true
  signature mismatch.
- Telescope inference was re-checked for the full uGMRT array, a south-arm
  sub-array, LOFAR station names including `SE607`, and MeerKAT.
- The reference-frequency helpers were exercised on six band/coverage
  combinations, confirming they are a no-op wherever the CASA default is
  already correct.
- `compare_fits_and_ms()` was checked on five constructed FITS/MS structures
  covering matching single-IF and multi-IF imports, a truncated import, a
  dropped-correlation import, and a non-UVFITS summary.
- `nearest_pointing_offset_arcmin()` was checked against identical directions,
  a 1 arcmin declination offset, a 1 arcmin great-circle offset in RA at
  dec=-60 (cos(dec) term), an RA=0 wrap, and nearest-of-several selection.

No CASA, no Measurement Set and no reference reduction were available in this
environment. These changes are recipe-level corrections verified by static
review and by the repository's own tests; they have not been exercised against
real telescope data.

## Not changed, but worth recording

- Only Taylor-0 is primary-beam corrected on the CASA path. For a 2:1
  fractional bandwidth the MFS image is a frequency-weighted average while the
  beam is evaluated at one reference frequency. This is already documented as
  a limitation and is the reason the independently corrected sub-band products
  exist.
- `gaintype='T'` for the uGMRT and MeerKAT gain solves averages the two
  parallel hands. That is defensible for Stokes I and is left as configured.
- The image sizes themselves are unchanged. Reaching PB=0.3 in uGMRT B2 would
  need roughly 18k pixels at 1.0", which is a resource decision for the
  observer, not something to change silently. The pipeline now reports the
  aperture it actually delivers so that decision can be made from evidence.
- Which beam model CASA uses internally for a given telescope was not
  verifiable here (no CASA in this environment). The fix removes the
  dependency rather than characterising it, which is the more robust outcome,
  but confirming the `.pb` extent on real data per telescope is still worth
  doing once.
- `--Facets-NFacets` is set to the number of DD solution directions plus one.
  Those are different quantities in DDFacet: the solution directions come from
  the cluster catalogue passed as `--NodesFile`, while `NFacets` controls the
  imaging facet tiling. The current coupling is not obviously harmful but it
  is not principled either. Left alone because DDFacet's exact `NFacets`
  semantics could not be verified here, and this is part of a command sequence
  the project treats as validated; `SCIENCE_VALIDATION.md` Gate 6 now requires
  the two numbers to be recorded and justified separately.
- `killms_nchan_solutions()` derives the number of frequency solution bands
  from `max(fmax) - min(fmin)` across spectral windows. For non-contiguous
  SPWs that spans the gaps and over-estimates the usable bandwidth, giving
  more solution bands than the data support.
- `antenna_stats()` runs one full `flagdata` summary per antenna, so it makes
  64 passes over a MeerKAT Measurement Set. Correct, but the dominant cost of
  the flagging stage's QA.
- `detect_bad_antennas()` skips antennas whose selected visibility count is
  below half the median, which is exactly the population of antennas that
  dropped out partway through an observation. A genuinely bad antenna present
  for only part of the run can therefore escape the threshold.
- The self-calibration acceptance metric is peak/RMS dynamic range, which
  rewards a rising peak as much as a falling noise floor. It is the
  conventional choice and is left as configured, but residual RMS reduction is
  the more direct measure of whether a loop helped.
