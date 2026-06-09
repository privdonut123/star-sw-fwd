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
inspect_fst_ftus_geometry.C
toy_fst_inner_disk_xygamma.C
FORWARD_ALIGNMENT_SUMMARY.md
```

Those changes add track-level post-selection metadata, fix the PV/FST sorting ambiguity, update the QA sorting check, configure the afterburner alignment workflow, add geometry and in-plane alignment toy macros, and keep this Markdown summary as the running alignment log. The earlier delta-z toy macro was later removed because the available branches did not yet provide a trustworthy per-hit local slope model for a real z solve.

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

Additional standalone macros used during the current alignment debugging are:

```bash
inspect_fst_ftus_geometry.C
toy_fst_inner_disk_xygamma.C
```

`toy_fst_inner_disk_xygamma.C` is a toy in-plane disk-level solve for `deltaX`, `deltaY`, and `gammaZ` using only inner FST sensors. It intentionally does not solve `deltaZ`.

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

The afterburner alignment configuration sets:

```cpp
nEvents default: 1000
```

This keeps the default run long enough for alignment QA while still bounded for quick iteration. A larger event count can still be passed explicitly at the ROOT macro call site.

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
- Decide later whether the committed afterburner default of `1000` events should remain long-term or be changed for large production-style alignment tests.

After that, the next alignment-specific task is to write the residual extraction macro and start producing candidate disk/sensor correction tables.

## 20. June 6, 2026 Update: Sensor-Local Coordinates, Outer FST Debugging, And First Delta-Z Toy

This update records the follow-up alignment work done after the initial residual QA and alignment tree scaffolding. The focus was not to solve a full alignment yet, but to clarify the FST coordinate model, make the FST planar measurement more alignment-natural, and build the first toy macro for a very restricted `delta z` alignment test.

### 20.1 Main Conceptual Clarification

The key conceptual separation is now:

```text
raw-ish FST strip measurement -> detector-local measurement coordinates
geometry cache / later alignment constants -> detector placement O/U/V
GenFit residual -> local difference between measured hit and fitted track crossing
```

This matters because an alignment procedure should move sensor placement parameters, not rewrite the measured strip coordinate. For GenFit planar measurements, the detector-local hit should be stable:

```text
meas0, meas1 = where the cluster fired inside this sensor
```

while the plane placement should carry the geometry:

```text
global hit model = O + meas0 * U + meas1 * V
```

where:

- `O` is the sensor active-center origin in STAR global coordinates.
- `U` is the sensor local x/radial-like direction expressed as a global vector.
- `V` is the sensor local y/phi-like direction expressed as a global vector.

This is the model we want for alignment iterations:

```text
same measured meas0/meas1
new O/U/V after applying alignment
new residuals
```

### 20.2 Geometry Origin And FTUS Sensor Ordering

The FST plane origin and axes are still created in:

```cpp
StRoot/StFwdTrackMaker/include/Tracker/FwdGeomUtils.h
```

The function is:

```cpp
FwdGeomUtils::getFstSensorOrigin(int index, TVector3 &u, TVector3 &v)
```

The important correction is that STAR hit sensor ordering and AGML `FTUS` copy ordering are not the same:

```cpp
// STAR event/hit sensor order:
sensor 0 = inner
sensor 1 = outer
sensor 2 = outer

// AGML FTUS copy order:
FTUS_1 = outer
FTUS_2 = outer
FTUS_3 = inner
```

The code now maps:

```cpp
static const int kEventSensorToFtusCopy[3] = {3, 1, 2};
```

This means:

```text
event sensor 0 -> FTUS_3
event sensor 1 -> FTUS_1
event sensor 2 -> FTUS_2
```

The origin is also no longer the raw `FTUS` node translation. The `FTUS` node translation is near the wedge/disk origin; the active silicon is an offset `TGeoTubeSeg` inside the node. The code now inspects the active shape and computes its active center:

```cpp
rCenter = 0.5 * (tube->GetRmin() + tube->GetRmax());
phiCenter = 0.5 * (tube->GetPhi1() + tube->GetPhi2());
```

then transforms that local active center to global coordinates:

```cpp
_matrix->LocalToMaster(activeLocal, activeMaster);
origin.SetXYZ(activeMaster[0], activeMaster[1], activeMaster[2]);
```

This fixed the earlier problem where all sensor origins were effectively at the wedge/disk origin rather than the active silicon center.

### 20.3 U And V Axes

`U` and `V` still come from the current geometry cache:

```cpp
u = column 0 of the TGeo rotation matrix
v = column 1 of the TGeo rotation matrix
```

The code extracts them as:

```cpp
u.SetXYZ(rot[0], rot[3], rot[6]);
v.SetXYZ(rot[1], rot[4], rot[7]);
```

Then it normalizes `V` so the local `V` direction is consistently counterclockwise / positive phi-like:

```cpp
if (u.Cross(v).Z() < 0) v = -v;
```

Important interpretation:

```text
U and V are local detector axes expressed in global STAR coordinates.
They are not local coordinate values.
```

The local coordinate values are `meas0` and `meas1`.

### 20.4 Direct Sensor-Local FST Measurement Conversion

The FST planar measurement calculation is now in:

```cpp
StRoot/StFwdTrackMaker/include/Tracker/TrackFitter.h
```

inside:

```cpp
TrackFitter::createTrackPointFromPlanarMeasurement(...)
```

The current FST-specific branch starts from:

```cpp
const int globalSensor = static_cast<int>(fh->_genfit_plane_index);
const int disk = globalSensor / (kFstNumWedgePerDisk * kFstNumSensorsPerWedge);
const int electronicWedge = (globalSensor / kFstNumSensorsPerWedge) % kFstNumWedgePerDisk;
const int sensor = globalSensor % kFstNumSensorsPerWedge;

const double r = fh->_localPosition[0];
const double stripPhi = fh->_localPosition[1];
```

Here:

```text
r        = radial strip center in cm
stripPhi = meanPhiStrip * kFstStripPitchPhi
```

`stripPhi` is not global STAR phi. It is the FST phi-strip index expressed as an angle in the wedge coordinate.

The full FST wedge has:

```text
128 phi bins per radial row
8 radial rows per wedge
```

The sensors are:

```text
sensor 0 inner: full 30 degree wedge, 4 radial rows x 128 phi bins
sensor 1 outer: one outer half, 4 radial rows x 64 phi bins
sensor 2 outer: the other outer half, 4 radial rows x 64 phi bins
```

However, `meanPhiStrip` is effectively a wedge-level `0..127` coordinate, not a local `0..63` coordinate for each outer sensor. This was an important debugging point.

The code defines:

```cpp
stripSign = kFstzFilp[disk] * kFstzDirct[electronicWedge];
halfWedgePhi = 0.5 * kFstNumPhiSegPerWedge * kFstStripPitchPhi;
edgeToCenterPhi = halfWedgePhi - 0.5 * kFstStripPitchPhi;
```

Meanings:

```text
stripSign
    +1 or -1 depending on disk and wedge orientation.
    It tells whether increasing strip number moves toward positive local phi/V.

halfWedgePhi
    15 degrees, because each wedge is 30 degrees wide.

edgeToCenterPhi
    angular distance from strip-0 center to wedge center.
    It is 15 degrees minus half a strip pitch because strip 0 is a strip center,
    not a physical wedge edge.
```

The local angular coordinate is:

```cpp
dphi
```

This is:

```text
hit angular offset relative to the wedge-center radial axis
```

It is not global phi.

The current formula is:

```cpp
double dphi = stripSign * (stripPhi - edgeToCenterPhi);
if (sensor == 1) {
    dphi = stripSign * (edgeToCenterPhi - stripPhi + 0.5 * kFstStripGapPhi);
} else if (sensor == 2) {
    dphi = stripSign * (edgeToCenterPhi - stripPhi - 0.5 * kFstStripGapPhi);
}
```

For the sensor center:

```cpp
const double sensorRSpan = 0.5 * kFstNumRStripsPerWedge * kFstStripPitchR;
const double centerR = (sensor == 0)
    ? kFstrStart[0] + 0.5 * sensorRSpan
    : kFstrStart[kFstNumRStripsPerWedge / 2] + 0.5 * sensorRSpan;
const double outerCenterDphi = 0.5 * (halfWedgePhi + kFstStripGapPhi);
```

The active center angles are:

```text
sensor 0: 0 degrees relative to wedge center
sensor 1: +8 degrees times stripSign
sensor 2: -8 degrees times stripSign
```

The local Cartesian measurement is then:

```cpp
hitOnPlane[0] = r * cos(dphi) - centerR * cos(centerDphi);
hitOnPlane[1] = r * sin(dphi) - centerR * sin(centerDphi);
```

Interpretation:

```text
hitOnPlane[0] = radial-like local coordinate relative to the active sensor center
hitOnPlane[1] = phi-like local coordinate relative to the active sensor center
```

This removed dependence on `plane->getO()` from the FST measurement-value calculation. That is more natural for alignment because moving the plane later should not change the strip-measured local coordinate.

### 20.5 Outer Sensor Sign Debugging

There was a short but important debugging loop on the outer sensors.

The geometry cache inspection showed, for disk 0 / wedge 0:

```text
sensor 0 active center: 75 degrees
sensor 1 active center: 83 degrees
sensor 2 active center: 67 degrees
```

Relative to the wedge center at 75 degrees:

```text
sensor 0:  0 degrees
sensor 1: +8 degrees
sensor 2: -8 degrees
```

For disk 1, `kFstzFilp` flips this:

```text
sensor 1: -8 degrees
sensor 2: +8 degrees
```

A temporary sign change wrongly treated sensor 2 as if its strip center were around strip `31.5`, like sensor 1. That was incorrect because `meanPhiStrip` is a wedge-level coordinate. Sensor 2's center is around strip `95.5`, not `31.5`.

The bad temporary behavior was effectively:

```text
stripSign = +1:
sensor 0 center: 0 degrees
sensor 1 hit center: +8 degrees
sensor 2 hit center: +7 degrees  <-- wrong side
```

The corrected behavior is:

```text
stripSign = +1:
sensor 0 center: 0 degrees
sensor 1 hit center: +8 degrees
sensor 2 hit center: -8 degrees

stripSign = -1:
sensor 0 center: 0 degrees
sensor 1 hit center: -8 degrees
sensor 2 hit center: +8 degrees
```

This is consistent with the `fGeom.root` active sensor centers and with the final code state.

The user's QA plot after the bad temporary patch showed very large alternating `pullUnbiased1` values. That plot should be discarded because it was produced with the wrong sign convention. The subsequent correction should be used for the next afterburner run.

### 20.6 Geometry Inspection Macro

Added:

```cpp
inspect_fst_ftus_geometry.C
```

Purpose:

```text
Print FST disk, wedge, event sensor id, GEANT FSTW copy,
FTUS copy, node origin, active center, radial range, phi range, and path.
```

Example usage:

```bash
root -l -b -q 'inspect_fst_ftus_geometry.C("fGeom.root")'
root -l -b -q 'inspect_fst_ftus_geometry.C("fGeom.root",0,0)'
root -l -b -q 'inspect_fst_ftus_geometry.C("fGeom.root",-1,-1,"fst_ftus_geometry.csv")'
```

This macro was used to verify that:

```text
FTUS node translation is not the active silicon center.
FTUS active shapes contain the correct R and phi ranges.
STAR sensor order and AGML FTUS copy order differ.
Outer sensor centers are at approximately +/-8 degrees from wedge center.
```

### 20.7 QA Macro Updates

Updated:

```cpp
fwd_alignment_residual_qa.C
```

Important current QA behavior:

```text
Measurement maps:
    Use all valid FST measurement rows.
    Do not apply the track momentum/eta cut.

Residual and pull plots:
    Require trackEta >= 2.5
    Require trackEta <= 4.0
    Require trackP > 0.5
    Require fit convergence
    Require residualDim == 2
```

The local measurement-map binning is now:

```text
meas0: -6 cm to +6 cm, 1 cm bins
meas1: -12 cm to +12 cm, 1 cm bins
```

The macro now includes:

```text
meas1 vs meas0 for all FST hits
meas1 vs meas0 by disk
meas1 vs meas0 by wedge
meas1 vs meas0 by sensor-in-wedge
meas1 vs meas0 by global sensor
residual vs local meas0/meas1
pull vs local meas0/meas1
mean residual / pull by sensor
occupancy diagnostics
```

The local measurement maps were crucial for seeing that the inner sensors centered naturally first, while the outer sensors were sensitive to the manual half-wedge conversion.

### 20.8 Delta-Z Toy Removed

The earlier file:

```cpp
toy_fst_inner_delta_z.C
```

was removed from the working tree.

The reason is conceptual, not just cleanup: the available `fwdAlign` branches currently contain track-level momentum, but not the fitted per-hit local track direction or the exact local derivatives at each measurement plane. A real z solve needs the response of the local residual to a displacement of the plane along global z:

```text
d(residual0) / dz
d(residual1) / dz
```

Using only global track momentum as a proxy produced toy numbers that were too easy to misread as alignment constants. Those numbers are no longer part of the recommended workflow.

Generated delta-z artifacts were also removed:

```text
toy_fst_delta_z_solver.*
toy_fst_delta_z_inner_res0.*
toy_fst_delta_z_inner_res1.*
toy_fst_inner_delta_z.root
toy_fst_inner_delta_z.pdf
```

### 20.9 In-Plane Inner-Sensor Toy Solver

The remaining toy solver is:

```cpp
toy_fst_inner_disk_xygamma.C
```

Purpose:

```text
Use only fstSensor == 0 rows to solve one in-plane rigid correction per FST disk.
```

The solved parameters are:

```text
deltaX
deltaY
gammaZ
```

The first-order model is:

```text
residual0 ~= -U dot [(deltaX, deltaY, 0) + gammaZ * (zhat x P)]
residual1 ~= -V dot [(deltaX, deltaY, 0) + gammaZ * (zhat x P)]
```

where:

```text
P = O + meas0 * U + meas1 * V
```

This is deliberately an in-plane toy. It does not solve `deltaZ`.

### 20.10 June 9 Closure Status

The newest useful closure diagnostic is:

```text
measGlobal = O + meas0 * U + meas1 * V
closure    = measGlobal - fstHitGlobal
```

After fixing the outer-sensor half-gap sign convention, the in-plane closure components are the trusted check:

```text
closureU = closure dot U
closureV = closure dot V
```

The current interpretation is:

```text
closureU/V passing -> local strip-to-sensor-coordinate mapping is probably consistent.
closureZ failing   -> sensor plane z/source-geometry convention still needs investigation.
```

This is important because FST `U` and `V` are nearly transverse:

```text
U_z ~= 0
V_z ~= 0
```

Therefore:

```text
measGlobal.z ~= O.z
closureZ ~= O.z - fstHitGlobalZ
```

So `closureZ` mainly tests whether the z coordinate used by `FwdGeomUtils::getFstSensorOrigin()` matches the z coordinate used when the original `StFwdHit` global position was built. It is not primarily a test of the `r/phi` local measurement conversion.

The failed z closure is not evidence by itself that `meas0/meas1` are wrong. It is evidence that z placement, shape-origin convention, or hit global z creation must be reconciled before any delta-z alignment attempt should be trusted.

### 20.11 Current Recommended Next Step

Before applying any correction constants:

1. Keep using `fwd_alignment_residual_qa.C` to monitor residuals, pulls, local measurement maps, and closure.
2. Treat `closureU/V` as the primary validation for the FST local measurement conversion.
3. Investigate `fstPlaneOriginZ` versus `fstHitGlobalZ` directly from `align_test.root`.
4. Trace the source of `fstHitGlobalZ` in the FST hit-making path.
5. Compare that source to the `FwdGeomUtils::getFstSensorOrigin()` active-center z calculation.
6. Do not revive a delta-z solve until the z convention is understood and per-hit local track slopes are available.

The conservative alignment path remains:

```text
first:  inner-sensor-only in-plane disk corrections
next:   inner-sensor-only per-sensor in-plane corrections
later:  z and outer-sensor studies after the z convention and outer geometry are stable
```
