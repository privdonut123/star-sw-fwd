# StFwdTrackMaker
The `StFrwdTrackMaker` performs tracking with the forward silicon tracker (FST) and the forward sTGC tracker (FTT). The tracking consists of several stages:
1. Load & sort space points from FST and FTT
2. Find seed tracks using FTT
3. Fit seed tracks using FTT + PV if available
4. Refit tracks with FST hits if  available
5. Fill StEvent structures

# Chain Options
- `fwdTrack` : the forward track maker
- Simulation:
    - `fttFastSim` : Ftt fast simulator, produces space points directly
    - `fstFastSim` : Fst fast simulator, produces space points directly
    - `fcsSim` : Fcs simulator, produces raw data to be processed by standard offline chain

# Development setup
1. Checkout the code 
```
git clone --no-checkout git@github.com:jdbrice/star-sw-1.git
cd star-sw-1
git config core.sparseCheckout true

touch .git/info/sparse-checkout

echo "StRoot/RTS/" >> .git/info/sparse-checkout
echo "StRoot/StBFChain/" >> .git/info/sparse-checkout
echo "StRoot/StEvent/" >> .git/info/sparse-checkout
echo "StRoot/StFstClusterMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFstDbMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFstHitMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFstRawHitMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFstSimMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFstUtil/" >> .git/info/sparse-checkout
echo "StRoot/StFttClusterMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFttDbMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFttHitCalibMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFttPointMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFttQAMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFttRawHitMaker/" >> .git/info/sparse-checkout
echo "StRoot/StFwdTrackMaker/" >> .git/info/sparse-checkout

echo "StDb" >> .git/info/sparse-checkout
echo "StarDb" >> .git/info/sparse-checkout

git checkout fwd-tracking

```
2. initialize the environment:
```
source StRoot/StFwdTrackMaker/macro/env.sh
cons # build the code if not done already
```
Note: this environment doesnt have the normal macro paths, so you might need to add (for DEV):
```
gSystem->Load( "libStarRoot.so" );
gROOT->SetMacroPath(".:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
```
to your macros if they load others (e.g. 'bfc.C').

3. Build the geometry cache
```
./StRoot/StFwdTrackMaker/macro/build_geom.C
```

Symlink macros for use. From project root
```
$ ls 
StarDb
StDb
StRoot
```

## For simulation:
```
ln -s StRoot/StFwdTrackMaker/macro/sim sim
```

### Test track seed finder
```
./sim/run_batch_seed <JOB_ID>
``` 

### Test track fitting (MC track finding)
```
./sim/run_batch_fast <JOB_ID>
``` 

### Test track finding + fitting 
```
./sim/run_batch_full <JOB_ID>    # WIP
``` 


## For DAQ data files:
```
ln -s StRoot/StFwdTrackMaker/macro/daq daq
```



# XML Configuration
The forward tracker can be configured via an XML configuration file.
The basic structure is:
```xml
<?xml version="1.0" encoding="UTF-8"?>
<config>
    <Output url="full_track.root" />
    <Geometry>fGeom.root</Geometry>
    <Source ... />
    
    <TrackFinder ... >
        ...
    </TrackFinder>
    
    <TrackFitter ... >
        ...
    </TrackFitter>
</config>
```

where the `<Output url="..." />` attribute determines the name of the 
