#ifndef TRACK_FITTER_H
#define TRACK_FITTER_H

#include "GenFit/ConstField.h"
#include "GenFit/EventDisplay.h"
#include "GenFit/Exception.h"
#include "GenFit/FieldManager.h"
#include "GenFit/KalmanFitStatus.h"
#include "GenFit/GblFitter.h"

#include "TDatabasePDG.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"

#include <vector>
#include <memory>

#include "StFwdTrackMaker/Common.h"

#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/STARField.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"

#include "StarGenerator/UTIL/StarRandom.h"
#include "FitterUtils.h"

/* Class for interfacing with GenFit for fitting tracks 
 *
 */
class TrackFitter {
// Accessors and options
  public:
    std::shared_ptr<genfit::Track> getTrack() { return mFitTrack; }

    // this is used rarely for debugging purposes, especially to check/compare plane misalignment
    static const bool kUseSpacePoints = true; // use spacepoints instead of planar measurements
    static const int kVerbose = 1; // verbosity level for debugging

  public:

    /**
     * @brief Construct a new Track Fitter object
     *
     * @param _mConfig : Config object
     * @param geoCache : Geometry cache filename
     */
    TrackFitter(FwdTrackerConfig _mConfig, TString geoCache) : mConfig(_mConfig), mGeoCache(geoCache), mFitTrack(nullptr) {}

    /**
     * @brief Setup the tracker object
     * Load geometry
     * Setup Material Effects
     * Setup the magnetic field
     * Setup the fitter
     * Setup the fit planes
     */
    void setup() {

        // the geometry manager that GenFit will use
        TGeoManager * gMan = nullptr;

        // Setup the Geometry used by GENFIT
        LOG_INFO << "StFwdTrackMaker is loading the geometry cache: " << mConfig.get<string>("Geometry", mGeoCache.Data()).c_str() << endm;
        TGeoManager::Import(mConfig.get<string>("Geometry", mGeoCache.Data()).c_str());
        gMan = gGeoManager;
        // Set up the material interface and set material effects on/off from the config
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

        // Set Material Stepper debug level
        genfit::MaterialEffects::getInstance()->setDebugLvl( mConfig.get<int>("TrackFitter.MaterialEffects:DebugLvl", 0) );

        genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch( mConfig.get<int>("TrackFitter.MaterialEffects.EnergyLossBetheBloch", true) );
        genfit::MaterialEffects::getInstance()->setNoiseBetheBloch( mConfig.get<int>("TrackFitter.MaterialEffects.NoiseBetheBloch", true) );
        genfit::MaterialEffects::getInstance()->setNoiseCoulomb( mConfig.get<int>("TrackFitter.MaterialEffects.NoiseCoulomb", true) );
        genfit::MaterialEffects::getInstance()->setEnergyLossBrems( mConfig.get<int>("TrackFitter.MaterialEffects.EnergyLossBrems", true) );
        genfit::MaterialEffects::getInstance()->setNoiseBrems( mConfig.get<int>("TrackFitter.MaterialEffects.NoiseBrems", true) );
        genfit::MaterialEffects::getInstance()->ignoreBoundariesBetweenEqualMaterials( mConfig.get<int>("TrackFitter.MaterialEffects.ignoreBoundariesBetweenEqualMaterials", true) );

        // do this last to override
        genfit::MaterialEffects::getInstance()->setNoEffects( !mConfig.get<bool>("TrackFitter:MaterialEffects", false)); // negated, true means defaul is all effects on (noEffects off)
        if (!mConfig.get<bool>("TrackFitter:MaterialEffects", false)){
            LOG_INFO << "Turning OFF GenFit Material Effects in stepper" << endm;
        }

        // Determine which Magnetic field to use
        // Either constant field or real field from StarFieldAdaptor
        if (mConfig.get<bool>("TrackFitter:constB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 5.0)); // 0.5 T Bz
            LOG_INFO << "StFwdTrackMaker: Tracking with constant magnetic field" << endm;
        } else if (mConfig.get<bool>("TrackFitter:zeroB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 0.)); // ZERO FIELD
            LOG_INFO << "StFwdTrackMaker: Tracking with ZERO magnetic field" << endm;
        } else {
            mBField = std::unique_ptr<genfit::AbsBField>(new StarFieldAdaptor());
            LOG_INFO << "StFwdTrackMaker: Tracking with StarFieldAdapter" << endm;
        }
        // we must have one of the two available fields at this point
        // note, the pointer is still bound to the lifetime of the TackFitter
        genfit::FieldManager::getInstance()->init(mBField.get());

        // initialize the main mFitter using a KalmanFitter with reference tracks
        mFitter = std::unique_ptr<genfit::AbsKalmanFitter>(new genfit::KalmanFitter());

        // Here we load several options from the config,
        // to customize the mFitter behavior
        mFitter->setMaxFailedHits(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxFailedHits", -1)); // default -1, no limit
        mFitter->setDebugLvl(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:DebugLvl", 0)); // default 0, no output
        mFitter->setMaxIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxIterations", 4000)); // default 4 iterations
        mFitter->setMinIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MinIterations", 10)); // default 0 iterations

        // Set the fit convergence paramters
        mFitter->setRelChi2Change( mConfig.get<double>("TrackFitter.KalmanFitterRefTrack:RelChi2Change", 1e-3) );
        // mFitter->setAbsChi2Change( mConfig.get<double>("TrackFitter.KalmanFitterRefTrack:AbsChi2Change", 1e-3) );
        mFitter->setDeltaPval( mConfig.get<double>("TrackFitter.KalmanFitterRefTrack:DeltaPval", 1e-3) );
        mFitter->setBlowUpFactor( mConfig.get<double>("TrackFitter.KalmanFitterRefTrack:BlowUpFactor", 1e9) );

        // FwdGeomUtils looks into the loaded geometry and gets detector z locations if present
        FwdGeomUtils fwdGeoUtils( gMan );

        // Create the genfit Planes for the FST sensors
        createAllFstPlanes( fwdGeoUtils );
        createAllFttPlanes( fwdGeoUtils );
        LOG_DEBUG << "Created all FST and FTT planes" << endm;
    }

    /**
     * @brief Convert the 3x3 covmat to 2x2 by dropping z
     *
     * @param h : hit with cov matrix
     * @return TMatrixDSym : cov matrix 2x2
     */
    TMatrixDSym CovMatPlane(KiTrack::IHit *h){
        TMatrixDSym cm(2);
        cm(0, 0) = static_cast<FwdHit*>(h)->_covmat(0, 0);
        cm(1, 1) = static_cast<FwdHit*>(h)->_covmat(1, 1);
        cm(0, 1) = static_cast<FwdHit*>(h)->_covmat(0, 1);
        return cm;
    }

    /**
     * @brief Get projection to a given plane
     *
     * @param fstPlane : plane index
     * @param fitTrack : track to project
     * @return genfit::MeasuredStateOnPlane
     */
    genfit::MeasuredStateOnPlane projectToPlane(genfit::SharedPlanePtr plane, std::shared_ptr<genfit::Track> fitTrack, int iState = 0) {
        genfit::MeasuredStateOnPlane nil;
        if (plane == nullptr) {
            LOG_ERROR << "Plane is null, cannot project" << endm;
            return nil;
        }
        if (fitTrack == nullptr) {
            LOG_ERROR << "Track is null, cannot project" << endm;
            return nil;
        }

        genfit::MeasuredStateOnPlane tst = fitTrack->getFittedState(iState);
        fitTrack->getCardinalRep()->extrapolateToPlane(tst, plane);
        return tst;   
    }


    /**
     * @brief Get projection to given FST plane
     *
     * @param fstPlane : plane index
     * @param fitTrack : track to project
     * @return genfit::MeasuredStateOnPlane
     */
    genfit::MeasuredStateOnPlane projectToFst(size_t fstSensorPlaneIndex, std::shared_ptr<genfit::Track> fitTrack) {
        if (fstSensorPlaneIndex > mFstSensorPlanes.size()) {
            genfit::MeasuredStateOnPlane nil;
            LOG_ERROR << "FST plane index out of range: " << fstSensorPlaneIndex << endm;
            return nil;
        }
        return projectToPlane(mFstSensorPlanes[fstSensorPlaneIndex], fitTrack);
    }

    /**
     * @brief Get projection to given FTT plane
     *
     * @param iFttPlane : plane index
     * @param fitTrack : track to project
     * @return genfit::MeasuredStateOnPlane
     */
    genfit::MeasuredStateOnPlane projectToFtt(size_t iFttPlane, std::shared_ptr<genfit::Track> fitTrack) {
        if (iFttPlane > mFttPlanes.size()) {
            LOG_ERROR << "FTT plane index out of range: " << iFttPlane << endm;
            genfit::MeasuredStateOnPlane nil;
            return nil;
        }
        return projectToPlane(mFttPlanes[iFttPlane], fitTrack);
    }

    /**
     * @brief setup the track from the given seed and optional primary vertex
     * @param trackSeed : seed points
     * @param seedMom : seed momentum
     * @param seedPos : seed position
     * @param Vertex : primary vertex
     */
    bool setupTrack(Seed_t trackSeed ) {
        
        // setup the track fit seed parameters
        GenericFitSeeder gfs;
        int seedQ = 1;
        TVector3 seedPos(0, 0, 0);
        TVector3 seedMom(0, 0, 10); // this default seed actually works better than a half-bad guess
        gfs.makeSeed( trackSeed, seedPos, seedMom, seedQ );
        LOG_DEBUG << "Setting track fit seed position = " << TString::Format( "(%f, %f, %f)", seedPos.X(), seedPos.Y(), seedPos.Z() ) << endm; 
        LOG_DEBUG << "Setting track fit seed momentum = " << TString::Format( "(%f, %f, %f)", seedMom.X(), seedMom.Y(), seedMom.Z() ) << endm;
        LOG_DEBUG << "Setting track fit seed charge = " << seedQ << endm;

        if ( seedQ == 0 ) {
            LOG_ERROR << "Seed charge is zero, skipping track -> usually means collinear points" << endm;
            return false;
        }

        // create the track representations
        // Note that multiple track reps differing only by charge results in a silent failure of GenFit
        auto theTrackRep = new genfit::RKTrackRep(mPdgMuon * -1 * seedQ); // bc pos PDG codes are for neg particles

        // Create the track
        mFitTrack = std::make_shared<genfit::Track>(theTrackRep, seedPos, seedMom);
        // now add the points to the track

        int hitId(0);       // hit ID
        size_t planeId(0);     // detector plane ID

        // initialize the hit coords on plane and spacepoint for PV
        TVectorD hitOnPlane(2);
        TVectorD spacepoint(3);

        /******************************************************************************************************************
		 * loop over the hits, add them to the track
		 ******************************************************************************************************************/
        // use these to enforce our sorting parameters
        size_t idxFst = 0; // index of the FST hit
        size_t idxFtt = 0; // index of the FTT hit
        for (auto h : trackSeed) {
            auto fh = dynamic_cast<FwdHit*>(h);
            if (fh == nullptr) {
                LOG_ERROR << "Hit is not a FwdHit, cannot add to track" << endm;
                continue;
            }
            

            /******************************************************************************************************************
            * If the Primary vertex is included
            ******************************************************************************************************************/
            if ( kUseSpacePoints || fh->isPV() ) {
                LOG_DEBUG << "Treating hit as a spacepoint" << endm;
                if ( fh->isPV() ){
                    LOG_DEBUG << "Including primary vertex in fit" << endm;
                }
                TVectorD pv(3);
                pv[0] = h->getX();
                pv[1] = h->getY();
                pv[2] = h->getZ();
                LOG_DEBUG << "x = " << pv[0] << "+/- " << fh->_covmat(0,0) << ", y = " << pv[1] << " +/- " << fh->_covmat(1,1) << ", z = " << pv[2] << " +/- " << fh->_covmat(2,2) << endm;

                auto tp = new genfit::TrackPoint();
                genfit::SpacepointMeasurement *measurement = new genfit::SpacepointMeasurement(pv, fh->_covmat, fh->_detid, ++hitId, tp);
                tp->addRawMeasurement(measurement);
                tp->setTrack(mFitTrack.get());

                // Set the sorting parameter
                if ( fh->isPV() ){
                    tp->setSortingParameter(0);
                }
                // These below are only used if kUseSpacePoints is true
                if ( fh->isFtt() ){
                    tp->setSortingParameter(4 + idxFtt);
                    idxFtt++;
                }
                if ( fh->isFst() ){
                    tp->setSortingParameter(1 + idxFst);
                    idxFst++;
                }
                // add the spacepoint to the track
                mFitTrack->insertPoint( tp );
                continue;
            }

            // Otherwise we treat the measurement as a planar measurement
            hitOnPlane[0] = h->getX();
            hitOnPlane[1] = h->getY();
            auto tp = new genfit::TrackPoint();
            genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitOnPlane, CovMatPlane(h), fh->_detid, ++hitId, tp);
            genfit::SharedPlanePtr plane = getPlaneFor( fh );
            planeId = fh->_genfit_plane_index;
            // I do this to make the planeId unique between FST and FTT
            // But I do not know if it is needed
            if (fh->isFtt()) {
                planeId = kFstNumSensors + fh->_genfit_plane_index;
            }          
            measurement->setPlane(plane, planeId);
            tp->addRawMeasurement(measurement);
            tp->setTrack(mFitTrack.get());
            tp->setSortingParameter(planeId); // or use the hitId?
            mFitTrack->insertPoint( tp );
        } // loop on trackSeed
        return true;
    } // setupTrack

    /** @brief performs the fit on a track
     *  @param t : track to fit
    */
    void performFit( std::shared_ptr<genfit::Track> t ){
        /******************************************************************************************************************
		 * Do the fit
		 ******************************************************************************************************************/
        try {

            // prepare the track for fitting
            // int nFailedPoints = 0;
            // bool changed = false;
            // changed = dynamic_cast<genfit::KalmanFitterRefTrack*>( mFitter.get() )->prepareTrack( mFitTrack.get(), mFitTrack->getCardinalRep(), false, nFailedPoints);
            // LOG_DEBUG << "Track prepared for fit with " << nFailedPoints << " failed points, changed? = " << changed << endm;

            // check the track for consistency
            mFitTrack->checkConsistency();
            // do the fit
            mFitter->processTrack(mFitTrack.get());

            // check the track for consistency
            mFitTrack->checkConsistency();

            // find track rep with smallest chi2
            mFitTrack->determineCardinalRep();
            // update the seed
            // t->udpateSeed();

            auto status = mFitTrack->getFitStatus();
            if (status == nullptr) {
                LOG_ERROR << "Fit status is null, cannot get fit status" << endm;
                return;
            }
            if ( kVerbose > 0 ) {
                LOG_INFO << "Fit status:  " << status->isFitConverged() << endm;
                LOG_INFO << "-Fit pvalue: " << status->getPVal() << endm;
                LOG_INFO << "-Fit Chi2:   " << status->getChi2() << endm;
            }

            if ( status->isFitConverged() ){
             
                auto cr = mFitTrack->getCardinalRep();
                auto p = cr->getMom( mFitTrack->getFittedState( 0, cr ));
                int rcQ = status->getCharge(); 
                if ( kVerbose > 0 ) { 
                    LOG_INFO << "Fit momentum: " << p.X() << ", " << p.Y() << ", " << p.Z() << endm;
                    LOG_INFO << "\tFit Pt: " << p.Pt() << ", eta: " << p.Eta() << ", phi: " << p.Phi() << endm;
                }
            }


        } catch (genfit::Exception &e) {
            LOG_ERROR << "Exception on fit update" << e.what() << endm;
        }
        if ( kVerbose > 0 ) {
            LOG_INFO << "Track fit update complete!" << endm;
        }
    }

    /*
     * @brief Get all FST planes
     *
     * @return std::vector<genfit::SharedPlanePtr>
     */
    void createAllFstPlanes( FwdGeomUtils &fwdGeoUtils )
    {
        // create FWD GeomUtils to get the plane locations
        for (int globalSensorIndex = 0; globalSensorIndex < kFstNumSensors; globalSensorIndex++)
        {
            TVector3 u(1, 0, 0); 
            TVector3 v(0, 1, 0);
            TVector3 o = fwdGeoUtils.getFstSensorOrigin(globalSensorIndex, u, v);
            if (kVerbose > 0) {
                LOG_INFO << "Adding FST Sensor " << globalSensorIndex << " at " << o.X() << ", " << o.Y() << ", " << o.Z() << endm;
                LOG_INFO << "\tSensor " << globalSensorIndex << " U = " << u.X() << ", " << u.Y() << ", " << u.Z() << endm;
                LOG_INFO << "\tSensor " << globalSensorIndex << " V = " << v.X() << ", " << v.Y() << ", " << v.Z() << endm;
            }
            mFstSensorPlanes.push_back(
                genfit::SharedPlanePtr(
                    new genfit::DetPlane(o, u, v)));
        }
    }

    /*
     * @brief Get all FTT planes
     *
     * @return std::vector<genfit::SharedPlanePtr>
     */
    void createAllFttPlanes( FwdGeomUtils &fwdGeoUtils ) {
        // create FWD GeomUtils to get the plane locations
        // 4 planes, 4 quadrants, 2 planes per quadrant
        for (int globalPlaneIndex = 0; globalPlaneIndex < 32; globalPlaneIndex++)
        {
            TVector3 u(1, 0, 0); 
            TVector3 v(0, 1, 0);
            TVector3 o = fwdGeoUtils.getFttQuadrant(globalPlaneIndex, u, v);
            if (kVerbose > 0) { 
                LOG_INFO << "Adding FTT Plane " << globalPlaneIndex << " at " << o.X() << ", " << o.Y() << ", " << o.Z() << endm;
                LOG_INFO << "\tPlane " << globalPlaneIndex << " U = " << u.X() << ", " << u.Y() << ", " << u.Z() << endm;
                LOG_INFO << "\tPlane " << globalPlaneIndex << " V = " << v.X() << ", " << v.Y() << ", " << v.Z() << endm;
            }
            mFttPlanes.push_back(
                genfit::SharedPlanePtr(
                    new genfit::DetPlane(o, u, v)));
        }
    }
    

    /**
     * @brief Primary track fitting routine
     *
     * @param trackSeed :
     * @param Vertex : Primary Vertex
     * @param seedMomentum : seed momentum (can be from MC)
     * @return void : the results can be accessed via the getTrack() method
     */
    long long fitTrack(Seed_t trackSeed, TVector3 *seedMomentum = 0) {
        long long itStart = FwdTrackerUtils::nowNanoSecond();
        LOG_DEBUG << "Fitting track with " << trackSeed.size() << " FWD Measurements" << endm;

        /******************************************************************************************************************
		 * First sort the seed, bc GENFIT seemingly cannot handle out of order points
		 ******************************************************************************************************************/
        std::sort(trackSeed.begin(), trackSeed.end(), 
            [](KiTrack::IHit *a, KiTrack::IHit *b) 
                { return a->getZ() < b->getZ(); }
        );

        /******************************************************************************************************************
		 * Setup the track fit seed parameters and objects
		 ******************************************************************************************************************/
        bool valid = setupTrack(trackSeed);
        if ( !valid ){
            LOG_ERROR << "Failed to setup track for fit" << endm;
            return -1;
        }
        LOG_DEBUG << "Ready to fit with " << mFitTrack->getNumPoints() << " track points" << endm;

        /******************************************************************************************************************
		 * Do the fit
		 ******************************************************************************************************************/
        performFit( mFitTrack );
        long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
        return duration;
    } // fitTrack

    genfit::SharedPlanePtr getPlaneFor( FwdHit * fh ){
        
        // sTGC
        if ( fh->isFtt() ){
            if ( fh->_genfit_plane_index > mFttPlanes.size() ) {
                LOG_ERROR << "Invalid FTT genfit plane index: " << fh->_genfit_plane_index << endm;
                return nullptr;
            }
            return mFttPlanes[fh->_genfit_plane_index];
        }

        // FST
        if ( fh->isFst() ){
            if ( fh->_genfit_plane_index > mFstSensorPlanes.size() ) {
                LOG_ERROR << "Invalid FST sensor genfit plane index: " << fh->_genfit_plane_index << endm;
                return nullptr;
            }
            return mFstSensorPlanes[ fh->_genfit_plane_index ];
        }
        LOG_ERROR << "Unknown FwdHit type, cannot get plane - if this is a PV then use an abs measurement" << endm;
        return nullptr;
    }

    // Store the planes for FTT and FST
    vector<genfit::SharedPlanePtr> mFttPlanes;
    vector<genfit::SharedPlanePtr> mFstSensorPlanes; // 108 planes, one for each sensor

  protected:
    std::unique_ptr<genfit::AbsBField> mBField;

    FwdTrackerConfig mConfig; // main config object
    TString mGeoCache; // the name of the geometry cache file

    // Main GenFit fitter instance
    std::unique_ptr<genfit::AbsKalmanFitter> mFitter = nullptr;

    // PDG codes for the default plc type for fits
    const int mPdgPiPlus = 211;
    const int mPdgPiMinus = -211;
    const int mPdgPositron = 11;
    const int mPdgElectron = -11;
    const int mPdgMuon = 13;
    const int mPdgAntiMuon = -13;

    // GenFit state - resused
    std::shared_ptr<genfit::Track> mFitTrack;
};


#endif
