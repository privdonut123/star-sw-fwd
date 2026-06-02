# Forward STAR Alignment: Detailed End-of-Day Summary

## 1. Current State

Work was done on branch:

```bash
xihe-align
```

This branch is based on `dev`. Before the final end-of-day commit, it had four pushed alignment commits:

```bash
d926ef0d05 Add opt-in Forward alignment diagnostics scaffold
a0e5ba5833 Add planar residual metadata for FST alignment
222ce3cb0b Update forward afterburner alignment run defaults
c29236f3ad Add Forward alignment pull diagnostics
```

The final end-of-day commit adds these remaining source/documentation changes:

```bash
StRoot/StFwdTrackMaker/StFwdTrackMaker.cxx
StRoot/StFwdTrackMaker/StFwdTrackMaker.h
StRoot/StFwdTrackMaker/include/Tracker/TrackFitter.h
StRoot/StFwdTrackMaker/macro/mudst/fwd_afterburner.C
fwd_alignment_residual_qa.C
FORWARD_ALIGNMENT_SUMMARY.md
```

Those final changes add track-level post-selection metadata, fix the PV/FST sorting ambiguity, update the QA sorting check, keep the afterburner default at effectively all events, and add this Markdown summary.

Generated files such as `align_test.root`, `fGeom.root`, `fwd_align_qa.pdf`, logs, MuDst/PicoDst outputs, and temporary output directories are still untracked and should not be committed unless explicitly needed.

Local ignore rules were added in `.git/info/exclude` for:

```bash
.sl79_gcc11/
temp_gccflags.c
```

## 2. Files Modified Or Added

The alignment work touched the Forward tracking maker, the GenFit track fitter path, the afterburner macros, the geometry-cache macro, and a new QA macro.

Main files changed:

```bash
StRoot/StFwdTrackMaker/StFwdTrackMaker.cxx
StRoot/StFwdTrackMaker/StFwdTrackMaker.h
StRoot/StFwdTrackMaker/include/Tracker/TrackFitter.h
StRoot/StFwdTrackMaker/macro/mudst/fwd_afterburner.C
StRoot/StFwdTrackMaker/macro/build_geom.C
fwd_alignment_residual_qa.C
```

The most important new file is:

```bash
fwd_alignment_residual_qa.C
```

It is a standalone ROOT macro that reads `align_test.root` and produces alignment QA plots.

## 3. What Was Added To `StFwdTrackMaker`

The main new feature is an opt-in alignment diagnostics tree.

The maker now has controls:

```cpp
void setFillAlignment(bool fill = true);
void setAlignmentOutputFilename(std::string fn);
```

When enabled, `StFwdTrackMaker` creates:

```cpp
TTree *fwdAlign
```

inside an output ROOT file. In the afterburner workflow this file is configured as:

```cpp
align_test.root
```

The purpose of this tree is to save one row per fitted track measurement, with enough information to study alignment residuals offline without rerunning reconstruction for every plot.

The central new method is:

```cpp
void StFwdTrackMaker::FillAlignment();
```

It loops over:

```cpp
mForwardTracker->getTrackResults()
```

then over each GenFit track point and raw measurement, and fills detector ID, hit ID, residuals, pulls, FST sensor identifiers, track quality, and track kinematics.

The alignment filling is called after Forward track fitting, so the tree reflects the fitted GenFit track state and residual information.

## 4. Alignment Tree Schema Added

The `fwdAlign` tree now stores event-level information:

```cpp
run
event
nSeeds
nFitTracks
```

Track-level information, repeated on every hit row from that track:

```cpp
trackIndex
chi2
ndf
pval
fitConverged
fitConvergedFully
fitConvergedPartially
trackNHitsFit
trackPx
trackPy
trackPz
trackP
trackPt
trackEta
```

Measurement-level information:

```cpp
pointIndex
measurementIndex
detId
hitId
measurementDim
residualDim
hasResidual
sorting
meas0
meas1
meas2
```

FST geometry grouping information:

```cpp
fstGlobalSensor
fstDisk
fstWedge
fstSensor
```

Residual information from GenFit:

```cpp
resBiased0
resBiased1
resBiased2
resUnbiased0
resUnbiased1
resUnbiased2
```

Residual uncertainties and pulls:

```cpp
resBiasedSigma0
resBiasedSigma1
resBiasedSigma2
pullBiased0
pullBiased1
pullBiased2

resUnbiasedSigma0
resUnbiasedSigma1
resUnbiasedSigma2
pullUnbiased0
pullUnbiased1
pullUnbiased2
```

Invalid or unavailable quantities are stored as:

```cpp
-99999.0
```

This is represented in code as:

```cpp
constexpr float kInvalidAlignValue = -99999.0f;
```

Analysis cuts should reject these values, for example:

```cpp
resUnbiased0 > -90000
pullUnbiased0 > -90000
trackEta > -90000
```

## 5. Helper Functions Added

Several helper functions were added in `StFwdTrackMaker.cxx`.

### `setAlignmentVector`

Purpose: copy a `TVectorD` into three fixed tree branches.

Reason: GenFit measurements and residuals may have dimension 2 or 3. The tree uses fixed branches `x0`, `x1`, and `x2`, so missing components must be filled with the invalid sentinel.

Behavior:

```cpp
x0 = source[0] if available, else -99999
x1 = source[1] if available, else -99999
x2 = source[2] if available, else -99999
```

This is used for raw measurement coordinates and residual vectors.

### `setAlignmentPulls`

Purpose: compute residual uncertainties and pulls from GenFit residual objects.

It reads:

```cpp
residual.getState()
residual.getCov()
```

Then for each available residual dimension:

```cpp
sigma = sqrt(cov(i, i))
pull  = residual_i / sigma
```

If the covariance is missing, non-positive, or unusable, the pull and sigma are stored as `-99999`.

This lets the QA macro check whether a residual is large only in absolute terms or also large compared with the uncertainty expected by GenFit.

### `setAlignmentTrackKinematics`

Purpose: derive track-level momentum quantities for post-selection.

It reads the fitted momentum stored in:

```cpp
gtr.mMomentum
```

and fills:

```cpp
trackPx = momentum.X()
trackPy = momentum.Y()
trackPz = momentum.Z()
trackPt = momentum.Perp()
trackP  = momentum.Mag()
trackEta = 0.5 * log((p + pz) / (p - pz))
```

Rapidity was intentionally removed. We only keep pseudorapidity `trackEta`, because it does not require a particle mass hypothesis.

## 6. Where Momentum Comes From

The track momentum used in the alignment tree comes from `GenfitTrackResult`.

When a GenFit track is stored, `GenfitTrackResult` sets:

```cpp
mMomentum = mTrack->getCardinalRep()->getMom(
    mTrack->getFittedState(0, mTrack->getCardinalRep())
);
```

So the alignment tree momentum is not a seed momentum. It is the fitted GenFit momentum at fitted state `0`.

The alignment tree then copies and derives from this value:

```cpp
trackPx
trackPy
trackPz
trackP
trackPt
trackEta
```

These values are repeated on every hit/residual row from that track. That repetition is intentional because it allows simple residual selections without joining to a separate track tree.

Example:

```cpp
detId == 45 &&
hasResidual > 0 &&
trackNHitsFit >= 4 &&
trackPt > 0.5
```

## 7. Planar FST/FTT Measurement Change

For alignment work, Forward silicon/tracker measurements were changed so FST and FTT are represented as planar GenFit measurements instead of spacepoints.

The important setting is:

```cpp
static constexpr bool kUseSpacePoints = false;
```

Effect:

- PV still remains a spacepoint.
- FST and FTT hits become planar measurements.
- FST residuals become two-dimensional local-plane residuals.
- `resUnbiased0` and `resUnbiased1` are now meaningful local planar residual components.
- `meas0` and `meas1` are local planar measurement coordinates.
- `meas2`, `res*2`, and `pull*2` are usually invalid for planar FST rows.

This was needed because alignment should be based on local detector-plane residuals, not 3D spacepoint residuals.

The current work focuses on FST. FTT is intentionally not interpreted yet because its plane indexing and sorting behavior still need separate checking.

## 8. FST Sensor Metadata

FST rows are decoded into disk/wedge/sensor identifiers.

The conversion helper is:

```cpp
FwdHit::fstSensorWedgeDiskFromGlobalIndex(globalIndex, disk, wedge, sensor)
```

The mapping is:

```cpp
globalSensor = disk * 12 * 3 + wedge * 3 + sensor
```

with:

```cpp
disk   = 0..2
wedge  = 0..11
sensor = 0..2
globalSensor = 0..107
```

This means:

```cpp
globalSensor 0   -> disk 0, wedge 0,  sensor 0
globalSensor 36  -> disk 1, wedge 0,  sensor 0
globalSensor 107 -> disk 2, wedge 11, sensor 2
```

The alignment tree stores all four values so that downstream QA and alignment extraction do not have to redo the conversion.

## 9. PV/FST Sorting Ambiguity Fix

A real ambiguity was found: previously both PV and FST global sensor 0 could have:

```cpp
sorting = 0
```

That made `sorting` ambiguous if viewed alone.

The final sorting fix reserves sorting value `0` for PV:

```cpp
PV sorting  = 0
FST sorting = 1..108
FTT sorting = 109...
```

The FST global sensor decode now does:

```cpp
globalSensor = sorting - 1
```

So the tree keeps the desired FST convention:

```cpp
fstGlobalSensor = 0..107
```

while avoiding a collision with PV.

After this fix, the QA expectation is:

```cpp
sorting - fstGlobalSensor = 1
```

for FST rows.

The disk/wedge/sensor conversion remains correct because the conversion helper still receives the original `globalSensor` value in the range `0..107`.

## 10. Afterburner And Geometry Macro Changes

`build_geom.C` was updated so the default geometry tag is:

```cpp
y2022
```

and the default geometry cache output is:

```cpp
fGeom.root
```

The afterburner was updated to enable alignment output:

```cpp
fwdTrack->setFillAlignment(true);
fwdTrack->setAlignmentOutputFilename("align_test.root");
```

It also uses the geometry cache:

```cpp
fwdTrack->setGeoCache("fGeom.root");
```

and turns on the Forward tracking chain.

The final afterburner change in `fwd_afterburner.C` sets:

```cpp
nEvents default: 100 -> 999999999
```

This makes the default afterburner run over effectively all available events unless a smaller event count is passed explicitly.

## 11. QA Macro Added

A new ROOT macro was added at the top of the repo:

```cpp
fwd_alignment_residual_qa.C
```

Usage:

```bash
root4star -l -b -q 'fwd_alignment_residual_qa.C("align_test.root","fwd_align_qa")'
```

Outputs:

```bash
fwd_align_qa.root
fwd_align_qa.pdf
```

The macro focuses on FST only for now. FTT is intentionally not interpreted yet because its plane/sorting behavior still needs separate checking.

The main FST selection is:

```cpp
detId == 45
hasResidual > 0
residualDim == 2
fstGlobalSensor >= 0
fitConverged > 0
```

The QA macro checks that required branches exist and conditionally adds pull plots if the pull branches are present.

If pull branches are missing, it still produces residual-only QA and prints a message saying that the input file needs to be regenerated with the updated maker.

## 12. What Each QA Plot Is For

The summary page checks:

- total alignment rows
- FST rows
- FST residual rows after cuts
- whether pull branches are available
- whether the file schema matches the expected alignment tree

The detector/dimension page checks:

- `detId` distribution
- `measurementDim` vs `residualDim`
- `hasResidual`
- `sorting - fstGlobalSensor`

For the current fixed sorting scheme, FST should peak at:

```cpp
sorting - fstGlobalSensor = 1
```

Residual distribution plots show:

- unbiased residual 0
- unbiased residual 1
- biased residual 0
- biased residual 1

Unbiased residuals are the main alignment observable. Biased residuals are useful for comparison but are less appropriate for deriving corrections because the hit being tested participates in the fit.

Pull distribution plots show:

```cpp
pull = residual / residualSigma
```

Good behavior is approximately:

```cpp
mean near 0
RMS near 1
```

If pull RMS is much larger than 1, the uncertainties may be underestimated or there may be unmodeled misalignment/tails.

If pull RMS is much smaller than 1, uncertainties may be overestimated or correlations may be too strong.

Occupancy plots check whether each disk/wedge/sensor has enough statistics. These are necessary before trusting any alignment constants.

Mean residual by disk/wedge/global sensor is the first actual alignment signal.

Residual versus local coordinate plots are used to distinguish:

- constant offset: translation-like misalignment
- slope versus local coordinate: rotation-like misalignment
- structured patterns: possible geometry, hit model, or weak-mode effects

Disk-split sensor profiles are intended as the first crude per-sensor alignment diagnostic.

## 13. How To Run The Current Workflow

Build geometry cache:

```bash
cd /Users/xihehan/Alignment/star-sw-fwd
root4star -l -b -q 'StRoot/StFwdTrackMaker/macro/build_geom.C("y2022","fGeom.root")'
```

Rebuild STAR code after source changes:

```bash
cons
```

Run afterburner:

```bash
root4star -l -b -q 'StRoot/StFwdTrackMaker/macro/mudst/fwd_afterburner.C("pp500.MuDst.root",100)'
```

Run QA:

```bash
root4star -l -b -q 'fwd_alignment_residual_qa.C("align_test.root","fwd_align_qa")'
```

Example ROOT selection after the new track branches exist:

```cpp
detId==45 &&
hasResidual>0 &&
residualDim==2 &&
trackNHitsFit>=4 &&
trackPt>0.5 &&
abs(trackEta)<4
```

Example draw:

```cpp
fwdAlign->Draw(
  "resUnbiased0*10000",
  "detId==45 && hasResidual>0 && trackNHitsFit>=4 && trackPt>0.5"
);
```

## 14. Current Validation Status

Confirmed locally:

- The QA macro runs on the existing `align_test.root`.
- The existing old file does not have the newest track kinematic branches.
- The existing old file does not have the newest sorting fix.
- The existing old file may not have branches added after it was produced.
- `git diff --check` passed for the modified code during the work.

Not yet confirmed:

- Full STAR rebuild after the newest final changes.
- New `align_test.root` produced after the newest sorting and track-kinematic changes.
- QA output after the newest tree schema.
- Whether `MeasurementOnPlane::getCov()` compiles cleanly in the exact STAR/GenFit environment.

## 15. CMS Paper Alignment Ideas To Adopt

The CMS tracker alignment paper describes a track-based alignment strategy built around residual minimization. The core idea is to adjust detector geometry parameters so that reconstructed hits agree with fitted track predictions.

The paper emphasizes several concepts that map well to Forward STAR:

1. Use track-hit residuals as the basic observable.
2. Minimize normalized residuals, not just raw residuals.
3. Use unbiased residuals where possible.
4. Work hierarchically: large structures first, then smaller modules.
5. Watch for weak modes and systematic distortions.
6. Validate with independent track-quality and physics-quality observables.
7. Use track samples with different topologies to break degeneracies.
8. Iterate: align, rerun reconstruction, remeasure residuals, repeat.

For STAR Forward, the nearest equivalent is:

```cpp
alignment observable = GenFit unbiased planar residual
normalization        = GenFit residual uncertainty
grouping             = FST disk / wedge / sensor
track quality        = fit convergence, chi2/ndf, pval, nHitsFit, pt, eta
```

The STAR implementation should be smaller and more staged than CMS. CMS solves many alignment parameters at once. For Forward STAR, it is safer to start with residual QA and simple corrections, then increase the number of degrees of freedom only after the residuals are understood.

## 16. Proposed STAR Forward Alignment Strategy

### Stage 1: Build trustworthy residual samples

First, produce a clean `align_test.root` with:

- planar FST residuals
- unbiased residuals
- residual uncertainties
- pulls
- FST disk/wedge/sensor IDs
- track kinematic and quality branches

Use only good rows:

```cpp
detId == kFstId
hasResidual > 0
residualDim == 2
fitConverged > 0
trackNHitsFit >= chosen threshold
trackPt > chosen threshold
reasonable trackEta range
valid residual and pull values
```

This stage is about making sure the residuals themselves are meaningful before solving any constants.

### Stage 2: Validate FST geometry bookkeeping

Before extracting constants, verify:

```cpp
fstGlobalSensor = 0..107
fstDisk = 0..2
fstWedge = 0..11
fstSensor = 0..2
sorting - fstGlobalSensor = 1
measurementDim = residualDim = 2
```

If these checks fail, alignment constants would be assigned to the wrong physical sensors.

### Stage 3: Start with coarse alignment

The first correction should not be full sensor-level six-degree-of-freedom alignment.

Start with coarse translations:

```cpp
disk-level mean residual 0
disk-level mean residual 1
```

Then inspect whether entire disks show coherent offsets.

If one disk has a nonzero mean residual while the others are near zero, that suggests a disk-level shift or reference-frame mismatch.

### Stage 4: Move to sensor-level translations

After disk-level behavior is understood, derive per-sensor corrections from:

```cpp
mean(resUnbiased0) by fstGlobalSensor
mean(resUnbiased1) by fstGlobalSensor
```

Initial approximation:

```cpp
delta_u_sensor ~= -mean(resUnbiased0)
delta_v_sensor ~= -mean(resUnbiased1)
```

The sign should be verified with a controlled test: apply a known small artificial shift and confirm the residual response.

Do not apply sensor-level constants from low-stat sensors.

### Stage 5: Add rotations only when justified

Rotations should be considered only if residuals show slopes versus local coordinates.

Examples:

```cpp
resUnbiased0 vs meas1 slope -> possible in-plane rotation contribution
resUnbiased1 vs meas0 slope -> possible in-plane rotation contribution
residual trends versus local radius/phi -> possible FST geometry orientation issue
```

For the first alignment pass, translations are safer than rotations.

### Stage 6: Use pulls to judge significance

Raw residual means tell us the size of deviations.

Pull means tell us whether the deviations are significant relative to the fitted uncertainty.

Useful checks:

```cpp
mean(pullUnbiased0) by sensor
mean(pullUnbiased1) by sensor
pull RMS by disk
pull RMS by trackPt
pull RMS by trackEta
```

A sensor with a large residual but also large uncertainty may not be urgent.

A sensor with a consistent nonzero pull is a stronger alignment candidate.

### Stage 7: Control weak modes

CMS emphasizes weak modes: geometry distortions that leave track residuals deceptively good while biasing track parameters.

For STAR Forward, possible weak-mode-like problems include:

- coherent disk shifts
- coherent rotations around the beamline
- radial scale distortions
- charge-dependent curvature biases
- eta-dependent residual trends
- pt-dependent residual trends

Validation should split residuals by:

```cpp
track charge
trackPt
trackEta
disk
wedge
sensor
track type if available
```

If residuals improve globally but become charge-dependent or eta-dependent, the alignment may be absorbing a tracking/modeling bias rather than detector geometry.

### Stage 8: Iterate

The alignment loop should be:

1. Build geometry cache.
2. Run afterburner.
3. Produce `align_test.root`.
4. Run QA macro.
5. Extract residual means/slopes.
6. Produce trial correction constants.
7. Rerun reconstruction with corrections.
8. Compare before/after QA.
9. Keep only corrections that improve residuals and do not introduce weak-mode signatures.

## 17. Near-Term Next Work

The next concrete implementation step should be a residual extraction macro.

It should read:

```cpp
align_test.root
```

and produce tables of:

```cpp
mean residual 0 by disk
mean residual 1 by disk
mean residual 0 by global sensor
mean residual 1 by global sensor
mean pull 0 by global sensor
mean pull 1 by global sensor
entries per sensor
```

It should apply configurable cuts:

```cpp
trackNHitsFit
trackPt
trackEta
chi2/ndf
fitConverged or fitConvergedFully
```

The output should initially be diagnostic only, not automatically applied to geometry.

Recommended first product:

```cpp
fst_alignment_constants_test.C
```

with output:

```bash
fst_alignment_constants_test.root
fst_alignment_constants_test.txt
```

The text output should contain one row per sensor:

```text
globalSensor disk wedge sensor n meanRes0_um errRes0_um meanRes1_um errRes1_um meanPull0 meanPull1
```

## 18. Assumptions And Defaults

Current assumptions:

- FST alignment is the first priority.
- FTT alignment is postponed until its plane/sorting behavior is understood.
- PV is a constraint and should not be treated as an alignable FST hit.
- `trackEta` is sufficient; rapidity is not needed.
- Planar residuals are the correct residual type for FST alignment.
- Unbiased residuals are the main alignment observable.
- Pulls are diagnostic and should not be used alone to derive corrections.
- The first correction level should be disk/sensor translations, not full rotations.
- Geometry constants should not be updated automatically until the residual extraction is validated.
- The sign convention for applying residual-derived shifts must be verified with a controlled artificial-shift test.

## 19. Open Items

After this final commit, before using the newest tree schema for alignment, check:

- STAR rebuild succeeds.
- New `align_test.root` contains `trackNHitsFit`, `trackPt`, and `trackEta`.
- New `align_test.root` has FST `sorting - fstGlobalSensor = 1`.
- `fstDisk`, `fstWedge`, and `fstSensor` still decode correctly.
- QA macro runs on the new file and produces pull pages.
- Decide later whether the committed afterburner default of `999999999` events should remain long-term or be changed back to a smaller default for quick tests.

After that, the next alignment-specific task is to write the residual extraction macro and start producing candidate disk/sensor correction tables.
