import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
from pylorentz import Momentum4

kFst = 24
kStgc = 26

recoFields = ['reco.mDidFitConverge',
 'reco.mDidFitConvergeFully',
 'reco.mNumberOfFailedPoints',
 'reco.mNumberOfSeedPoints',
 'reco.mNumberOfFitPoints',
 'reco.mChi2',
 'reco.mNDF',
 'reco.mPval',
 'reco.mCharge',
 'reco.mPrimaryMomentum.fUniqueID',
 'reco.mPrimaryMomentum.fBits',
 'reco.mPrimaryMomentum.fX',
 'reco.mPrimaryMomentum.fY',
 'reco.mPrimaryMomentum.fZ',
 'reco.mIdTruth',
 'reco.mQATruth',
 'reco.mDCA[3]',
 'reco.mVtxIndex']

class Single:
    """
    A class to represent a single nested object (track, McTrack etc.) from the ROOT file.
    """
    def __init__(self, event_record, index, prefix=""):
        self.event_record = event_record
        self.index = index
        self.prefix = prefix
        fields = event_record.fields
        self.record = ak.Record({field: event_record[field][index] for field in fields})        

    def __getitem__(self, key):
        # Construct the full field name by prepending "mcTracks."
        full_key = f"{self.prefix}{key}"
        # Return the value from the record if the field exists
        if full_key in self.record.fields:
            return self.record[full_key]
        raise KeyError(f"Field '{full_key}' not found in record")
    
    @property
    def lv(self, E=1):
        # if "mPxyz.mX1" in self:
        return Momentum4(self["mE"], self["mPxyz.mX1"], self["mPxyz.mX2"], self["mPxyz.mX3"])
        # elif "mPrimaryMomentum.fY" in self:
            # return Momentum4(E, self["mPrimaryMomentum.fX"], self["mPrimaryMomentum.fY"], self["mPrimaryMomentum.fZ"])

f = uproot.open("single_particle_gun_muminus_500Events_4PerEvent_Pt_0.10to5.00_Eta_2.50to4.00_Phi_0.00to6.28.FwdTree.root")
# print(f.keys())
FwdTree = f["fwd"]
# Get the branches
branches = FwdTree.keys()
# Print the branches
# print("Branches in FwdTree:", branches)


# get the mcTrack info
mcTrackBranch = FwdTree["mcTracks"]
mcTracks = mcTrackBranch.arrays( ak_add_doc=True, library="ak" )
print( "Tree has %d events" % (f["fwd"].num_entries) )

# mask_eta = (mcTracks["mcTracks.mEta"] > 2.5) & (mcTracks["mcTracks.mEta"] < 4.0)


muon_eta = np.zeros(1)
numFST = np.zeros(f["fwd"].num_entries)
fields = mcTracks.fields
print("fields:", fields)

recoBranch = FwdTree["reco"]
recoTracks = recoBranch.arrays( recoFields, ak_add_doc=True, library="ak" )

#['mcTracks.mGePid', 'mcTracks.mEta']
for iEvent in np.arange( 0, f["fwd"].num_entries ):

    numMcTracksThisEvent = mcTrackBranch.array()[iEvent]
    McTracks = [ Single( mcTracks[fields][iEvent], iMcTrack, prefix="mcTracks." ) for iMcTrack in np.arange(0, numMcTracksThisEvent) ]
    for iMcTrack in np.arange(0, numMcTracksThisEvent):

        mct = Single( mcTracks[fields ][iEvent], iMcTrack, prefix="mcTracks." )
        if mct.lv.eta > 2.5 and mct.lv.eta < 4.0:
            numFST[iEvent] = mct["mHits[30]"][kFst]
        else:
            numFST[iEvent] = 0

        # Hits in both FST and sTGC
        if mct["mHits[30]"][kFst] > 0 and mct["mHits[30]"][kStgc] > 0:
            muon_eta = np.append(muon_eta, mct.lv.eta)

    for iRecoTrack in np.arange(0, f["fwd"]["reco"].array()[iEvent]):
        reco = Single( recoTracks[ recoFields ][iEvent], iRecoTrack, prefix="reco." )
        # print( "reco.mIdTruth = ", reco["mIdTruth"] )
        mct = McTracks[ reco["mIdTruth"] - 1 ]
        # print( "mct.mId = ", mct["mId"] )


# plot hist of muon_eta
plt.hist( muon_eta.tolist(), bins=100, range=(-1, 4.0), histtype='step', label='Muon Eta')
plt.show()
# print( muon_eta )

plt.hist( numFST, bins=10, range=(0, 10), histtype='step', label='Num FST Hits')
plt.xlabel("Num FST Hits")
plt.ylabel("Num Tracks")
plt.title("Num FST Hits{ 2.5 < eta < 4.0 }")
plt.legend()
plt.show()

plt.hist( ak.flatten(mcTracks["mcTracks.mCharge"]), bins=10, range=(-1, 1), histtype='stepfilled', label='Charge')
plt.xlabel("Charge")
plt.ylabel("Num Tracks")
plt.title("Charge{ 2.5 < eta < 4.0 }")  
plt.legend()
plt.show()

['header',
 'header/TObject',
 'header/TObject/fUniqueID',
 'header/TObject/fBits',
 'header/pv',
 'header/pv/pv.TObject',
 'header/pv/pv.TObject/pv.fUniqueID',
 'header/pv/pv.TObject/pv.fBits',
 'header/pv/pv.fX',
 'header/pv/pv.fY',
 'header/pv/pv.fZ',
 'header/run',
 'header/event',
 'header/tofmult',
 'header/vpdVz',
 'mcTracks',
 'mcTracks/mcTracks.mGePid',
 'mcTracks/mcTracks.mId',
 'mcTracks/mcTracks.mIsShower',
 'mcTracks/mcTracks.mHits[30]',
 'mcTracks/mcTracks.mItrmdVertex',
 'mcTracks/mcTracks.mIdVx',
 'mcTracks/mcTracks.mIdVxEnd',
 'mcTracks/mcTracks.mCharge',
 'mcTracks/mcTracks.mE',
 'mcTracks/mcTracks.mEta',
 'mcTracks/mcTracks.mPxyz.mX1',
 'mcTracks/mcTracks.mPxyz.mX2',
 'mcTracks/mcTracks.mPxyz.mX3',
 'mcTracks/mcTracks.mpT',
 'mcTracks/mcTracks.mPtot',
 'mcTracks/mcTracks.mRapidity',
 'nSeedTracks',
 'fstHits',
 'fstHits/fstHits.fUniqueID',
 'fstHits/fstHits.fBits',
 'fstHits/fstHits.mId',
 'fstHits/fstHits.mIdTruth',
 'fstHits/fstHits.mApv',
 'fstHits/fstHits.mMaxTimeBin',
 'fstHits/fstHits.mMeanRStrip',
 'fstHits/fstHits.mMeanPhiStrip',
 'fstHits/fstHits.mCharge',
 'fstHits/fstHits.mChargeErr',
 'fstHits/fstHits.mNRawHits',
 'fstHits/fstHits.mNRawHitsR',
 'fstHits/fstHits.mNRawHitsPhi',
 'fstHits/fstHits.mLocalPosition[3]',
 'fstHits/fstHits.mHardwarePosition',
 'fstHits/fstHits.mXYZ.fUniqueID',
 'fstHits/fstHits.mXYZ.fBits',
 'fstHits/fstHits.mXYZ.fX',
 'fstHits/fstHits.mXYZ.fY',
 'fstHits/fstHits.mXYZ.fZ',
 'fttPoints',
 'fttPoints/fttPoints.fUniqueID',
 'fttPoints/fttPoints.fBits',
 'fttPoints/fttPoints.mPlane',
 'fttPoints/fttPoints.mQuadrant',
 'fttPoints/fttPoints.mX',
 'fttPoints/fttPoints.mY',
 'fttPoints/fttPoints.mClusters',
 'fttPoints/fttPoints.mXYZ.fUniqueID',
 'fttPoints/fttPoints.mXYZ.fBits',
 'fttPoints/fttPoints.mXYZ.fX',
 'fttPoints/fttPoints.mXYZ.fY',
 'fttPoints/fttPoints.mXYZ.fZ',
 'fttPoints/fttPoints.mIdTruth',
 'fttPoints/fttPoints.mQaTruth',
 'fttClusters',
 'fttClusters/fttClusters.fUniqueID',
 'fttClusters/fttClusters.fBits',
 'fttClusters/fttClusters.mId',
 'fttClusters/fttClusters.mPlane',
 'fttClusters/fttClusters.mQuadrant',
 'fttClusters/fttClusters.mOrientation',
 'fttClusters/fttClusters.mNStrips',
 'fttClusters/fttClusters.mSumAdc',
 'fttClusters/fttClusters.mX',
 'fttClusters/fttClusters.mSigma',
 'fttClusters/fttClusters.mRawHits',
 'fttClusters/fttClusters.mNeighbors',
 'fttClusters/fttClusters.mPoints',
 'fttClusters/fttClusters.mIdTruth',
 'fttClusters/fttClusters.mQaTruth',
 'wcalClusters',
 'wcalClusters/wcalClusters.fUniqueID',
 'wcalClusters/wcalClusters.fBits',
 'wcalClusters/wcalClusters.mXYZ.fUniqueID',
 'wcalClusters/wcalClusters.mXYZ.fBits',
 'wcalClusters/wcalClusters.mXYZ.fX',
 'wcalClusters/wcalClusters.mXYZ.fY',
 'wcalClusters/wcalClusters.mXYZ.fZ',
 'wcalClusters/wcalClusters.mClu',
 'hcalClusters',
 'hcalClusters/hcalClusters.fUniqueID',
 'hcalClusters/hcalClusters.fBits',
 'hcalClusters/hcalClusters.mXYZ.fUniqueID',
 'hcalClusters/hcalClusters.mXYZ.fBits',
 'hcalClusters/hcalClusters.mXYZ.fX',
 'hcalClusters/hcalClusters.mXYZ.fY',
 'hcalClusters/hcalClusters.mXYZ.fZ',
 'hcalClusters/hcalClusters.mClu',
 'wcalHits',
 'wcalHits/wcalHits.fUniqueID',
 'wcalHits/wcalHits.fBits',
 'wcalHits/wcalHits.mXYZ.fUniqueID',
 'wcalHits/wcalHits.mXYZ.fBits',
 'wcalHits/wcalHits.mXYZ.fX',
 'wcalHits/wcalHits.mXYZ.fY',
 'wcalHits/wcalHits.mXYZ.fZ',
 'wcalHits/wcalHits.mHit',
 'hcalHits',
 'hcalHits/hcalHits.fUniqueID',
 'hcalHits/hcalHits.fBits',
 'hcalHits/hcalHits.mXYZ.fUniqueID',
 'hcalHits/hcalHits.mXYZ.fBits',
 'hcalHits/hcalHits.mXYZ.fX',
 'hcalHits/hcalHits.mXYZ.fY',
 'hcalHits/hcalHits.mXYZ.fZ',
 'hcalHits/hcalHits.mHit',
 'epdHits',
 'epdHits/epdHits.fUniqueID',
 'epdHits/epdHits.fBits',
 'epdHits/epdHits.mXYZ.fUniqueID',
 'epdHits/epdHits.mXYZ.fBits',
 'epdHits/epdHits.mXYZ.fX',
 'epdHits/epdHits.mXYZ.fY',
 'epdHits/epdHits.mXYZ.fZ',
 'epdHits/epdHits.mHit',
 'reco',
 'reco/reco.fUniqueID',
 'reco/reco.fBits',
 'reco/reco.mProjections',
 'reco/reco.mFTTPoints',
 'reco/reco.mFSTPoints',
 'reco/reco.mEcalClusters',
 'reco/reco.mHcalClusters',
 'reco/reco.mDidFitConverge',
 'reco/reco.mDidFitConvergeFully',
 'reco/reco.mNumberOfFailedPoints',
 'reco/reco.mNumberOfSeedPoints',
 'reco/reco.mNumberOfFitPoints',
 'reco/reco.mChi2',
 'reco/reco.mNDF',
 'reco/reco.mPval',
 'reco/reco.mCharge',
 'reco/reco.mPrimaryMomentum.fUniqueID',
 'reco/reco.mPrimaryMomentum.fBits',
 'reco/reco.mPrimaryMomentum.fX',
 'reco/reco.mPrimaryMomentum.fY',
 'reco/reco.mPrimaryMomentum.fZ',
 'reco/reco.mIdTruth',
 'reco/reco.mQATruth',
 'reco/reco.mDCA[3]',
 'reco/reco.mVtxIndex',
 'seeds',
 'seeds/seeds.fUniqueID',
 'seeds/seeds.fBits',
 'seeds/seeds.mXYZ.fUniqueID',
 'seeds/seeds.mXYZ.fBits',
 'seeds/seeds.mXYZ.fX',
 'seeds/seeds.mXYZ.fY',
 'seeds/seeds.mXYZ.fZ',
 'seeds/seeds.mTrackId',
 'seeds/seeds.mSector',
 'seeds/seeds.mCov[9]']