#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/StFwdHitMaker.h"


void StFwdHitMaker::loadFttHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    LOG_DEBUG << "Loading FTT Hits" << endm;
    LOG_DEBUG << "Ftt From Source: " << mFttDataSource << endm;
    if ( StFwdHitLoader::DataSource::GEANT == mFttDataSource ){
        LOG_DEBUG << "Loading FTT Hits from GEANT" << endm;
        loadFttHitsFromGEANT( mcTrackMap, hitMap, count );
    } else if ( StFwdHitLoader::DataSource::StEvent == mFttDataSource ){
        LOG_DEBUG << "Loading FTT Hits from StEvent" << endm;
        loadFttHitsFromStEvent( mcTrackMap, hitMap, count );
    } else if ( StFwdHitLoader::DataSource::MuDst == mFttDataSource ){
        LOG_DEBUG << "Loading FTT Hits from MuDst (Not Yet Implemented)" << endm;
        // loadFttHitsFromMuDst( mcTrackMap, hitMap, count );
    }
    return;
} // loadFttHits

void StFwdHitMaker::loadFttHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    LOG_DEBUG << "Loading FTT Hits from Data" << endm;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    StFttCollection *col = event->fttCollection();
    size_t numFwdHitsPrior = mFwdHitsFtt.size();

    if ( col && col->numberOfPoints() > 0 ){
        LOG_DEBUG << "The Ftt Collection has " << col->numberOfPoints() << " points" << endm;
        TMatrixDSym hitCov3(3);
        const double sigXY = 0.2; //
        hitCov3(0, 0) = sigXY * sigXY;
        hitCov3(1, 1) = sigXY * sigXY;
        hitCov3(2, 2) = 4; // unused since they are loaded as points on plane
        static const double mm_to_cm = 0.1;
        for ( auto point : col->points() ){

            float xcm = point->xyz().x()*mm_to_cm;
            float ycm = point->xyz().y()*mm_to_cm;
            float zcm = point->xyz().z();

            // get the track id
            int track_id = point->idTruth();
            shared_ptr<McTrack> mcTrack = nullptr;
            if ( mcTrackMap.count(track_id) ) {
                mcTrack = mcTrackMap[track_id];
                LOG_DEBUG << "Adding McTrack to FTT hit: " << track_id << endm;
            }

            mFwdHitsFtt.push_back(FwdHit(count++, // id
                xcm, ycm, zcm,
                -point->plane(), // volume id
                kFttId, // detid
                track_id, // track id
                hitCov3, // covariance matrix
                mcTrack) // mcTrack
                );
            mFttHits.push_back( TVector3( xcm, ycm, zcm)  );
        } // end of loop over points
    } else {
        LOG_DEBUG << "The Ftt Collection is EMPTY points" << endm;
    }

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFtt.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFtt[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
        // add to MC track map
        if ( hit->getMcTrack() ){
            hit->getMcTrack()->addFttHit(hit);
        }
    }

    if ( numFwdHitsPost != numFwdHitsPrior ){
        LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FTT hits from StEvent" << endm;
    }
}

void StFwdHitMaker::loadFttHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    /************************************************************/
    // STGC Hits
    St_g2t_fts_hit *g2t_stg_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");

    size_t numFwdHitsPrior = mFwdHitsFtt.size();
    if (!g2t_stg_hits){
        LOG_WARN << "geant/g2t_stg_hit is empty" << endm;
        return;
    }

    // make the Covariance Matrix once and then reuse
    TMatrixDSym hitCov3(3);
    const double sigXY = 0.02;
    hitCov3(0, 0) = sigXY * sigXY;
    hitCov3(1, 1) = sigXY * sigXY;
    hitCov3(2, 2) = 0.1; // unused since they are loaded as points on plane

    int nstg = g2t_stg_hits->GetNRows();

    LOG_DEBUG << "This event has " << nstg << " stg hits in geant/g2t_stg_hit " << endm;

    bool filterGEANT = mFwdConfig.get<bool>( "Source:fttFilter", false );
    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_stg_hits->At(i);
        if (0 == git)
            continue; // geant hit
        int track_id = git->track_p;
        int volume_id = git->volume_id;
        int plane_id = (volume_id - 1) / 100;           // from 1 - 16. four chambers per station

        // only use the hits on the front modules
        if ( volume_id % 2 ==0 )
            continue;

        float x = git->x[0];// + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float y = git->x[1];// + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float z = git->x[2];

        if (plane_id < 0 || plane_id >= 4) {
            continue;
        }
        mFwdHitsFtt.push_back(
            FwdHit(
                count++, // id
                x, y, z, // position
                -plane_id, // volume id
                kFttId, // detid
                track_id, // track id
                hitCov3, // covariance matrix
                mcTrackMap[track_id] // mcTrack
                )
            );
        mFttHits.push_back( TVector3( x, y, z )  );
    } // loop on hits

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFtt.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFtt[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);

        if ( dynamic_cast<FwdHit*>(hit)->_mcTrack ){
            dynamic_cast<FwdHit*>(hit)->_mcTrack->addFttHit(hit);
        }
    }

    if ( numFwdHitsPost != numFwdHitsPrior ){
        LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from MuDst" << endm;
    }

} // loadFttHits



/**
 * @brief Loads FST hits from various sources into the hitmap and McTrackMap (if availabale)
 *
 * Order of precedence:
 * MuDst StMuFstCollection (Data)
 * StEvent StFstHitCollection (Data or slowsim)
 * StEvent StRndHitCollection (fast sim)
 * GEANT St_g2t_fts_hit (starsim only) - note if rasterizer is active this takes priority over FastSim
 *
 * @param mcTrackMap : MC track map if running sim
 * @param hitMap : FST hitmap to populate
 * @param count  : number of hits loaded
 */
int StFwdHitMaker::loadFstHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){

    int count = loadFstHitsFromMuDst(mcTrackMap, hitMap);
    if ( count > 0 ) return count; // only load from one source at a time

    count += loadFstHitsFromStEvent(mcTrackMap, hitMap);
    if ( count > 0 ) return count; // only load from one source at a time

    bool siRasterizer = mFwdConfig.get<bool>( "SiRasterizer:active", false );

    if ( !siRasterizer ) count += loadFstHitsFromStRnDHits( mcTrackMap, hitMap );
    if ( count > 0 ) return count; // only load from one source at a time

    return loadFstHitsFromGEANT( mcTrackMap, hitMap );
} // loadFstHits

int StFwdHitMaker::loadFstHitsFromMuDst( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    int count = 0;
    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker) {
        LOG_WARN << " No MuDstMaker ... bye-bye" << endm;
        return 0;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        LOG_WARN << " No MuDst ... bye-bye" << endm;
        return 0;
    }

    StMuFstCollection * fst = mMuDst->muFstCollection();
    if (!fst) {
        LOG_WARN << "No StMuFstCollection ... bye-bye" << endm;
        return 0;
    }

    size_t numFwdHitsPrior = mFwdHitsFst.size();
    LOG_INFO << "Loading " << fst->numberOfHits() << " StMuFstHits" << endm;
    TMatrixDSym hitCov3(3);
    for ( unsigned int index = 0; index < fst->numberOfHits(); index++){
        StMuFstHit * muFstHit = fst->getHit( index );

        float vR = muFstHit->localPosition(0);
        float vPhi = muFstHit->localPosition(1);
        float vZ = muFstHit->localPosition(2);

        int wedgeIndex  = muFstHit->getWedge();
        int sensorIndex = muFstHit->getSensor();
        int diskIndex   = muFstHit->getDisk();
        int globalIndex = FwdHit::fstGlobalSensorIndex( diskIndex, wedgeIndex, sensorIndex );

        float x0 = vR * cos( vPhi );
        float y0 = vR * sin( vPhi );
        hitCov3 = makeSiCovMat( TVector3( x0, y0, vZ ), mFwdConfig );

        
        mFstHits.push_back( TVector3( x0, y0, vZ)  );

        // we use d+4 so that both FTT and FST start at 4
        mFwdHitsFst.push_back(
            FwdHit(
                count++, // id
                x0, y0, vZ, // position
                diskIndex+4, // volume id
                kFstId, // detid
                0, // track id
                hitCov3, // covariance matrix
                nullptr // mcTrack
            )
        );
        mFwdHitsFst.back()._genfit_plane_index = globalIndex;
    } // index

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
    }

    if ( numFwdHitsPost != numFwdHitsPrior ){
        LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from MuDst" << endm;
    }

    // TODO add to hitmap
    return count;
} // loadFstHitsFromMuDst

int StFwdHitMaker::loadFstHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    int count = 0;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (!event) {
        LOG_WARN << "No StEvent, cannot load FST hits from StEvent StFstHitCollection" << endm;
        return 0;
    }
    LOG_DEBUG << "Got StEvent, loading Fst Hits" << endm;
    StFstHitCollection *fstHitCollection = event->fstHitCollection();
    size_t numFwdHitsPrior = mFwdHitsFst.size();

    if ( fstHitCollection && fstHitCollection->numberOfHits() > 0){
        // reuse this to store cov mat
        TMatrixDSym hitCov3(3);
        LOG_DEBUG << "StFstHitCollection is NOT NULL, loading hits" << endm;
        for ( unsigned int iw = 0; iw < kFstNumWedges; iw++ ){
            StFstWedgeHitCollection * wc = fstHitCollection->wedge( iw );
            if ( !wc ) continue;
            for ( unsigned int is = 0; is < kFstNumSensorsPerWedge; is++ ){
                StFstSensorHitCollection * sc = wc->sensor( is );
                if ( !sc ) continue;
                StSPtrVecFstHit fsthits = sc->hits();
                for ( unsigned int ih = 0; ih < fsthits.size(); ih++ ){
                    float vR   = fsthits[ih]->localPosition(0);
                    float vPhi = fsthits[ih]->localPosition(1);
                    float vZ   = fsthits[ih]->localPosition(2);
                    LOG_INFO << TString::Format("FST local position: %f %f %f", vR, vPhi, vZ) << endm;
                    
                    int wedgeIndex  = iw % kFstNumWedgePerDisk;
                    int sensorIndex = is % kFstNumSensorsPerWedge;
                    int diskIndex   = iw / kFstNumWedgePerDisk;
                    int globalIndex = FwdHit::fstGlobalSensorIndex( diskIndex, wedgeIndex, sensorIndex );

                    LOG_DEBUG << "diskIndex = " << diskIndex << ", wedgeIndex = " << wedgeIndex << ", sensorIndex = " << sensorIndex << ", globalIndex = " << globalIndex << endm;
                    float x0 = vR * cos( vPhi );
                    float y0 = vR * sin( vPhi );
                    hitCov3 = makeSiCovMat( TVector3( x0, y0, vZ ), mFwdConfig );

                       
                    mFstHits.push_back( TVector3( x0, y0, vZ)  );
                    int track_id = fsthits[ih]->idTruth();
                    shared_ptr<McTrack> mcTrack = nullptr;
                    if ( mcTrackMap.count(track_id) ) {
                        mcTrack = mcTrackMap[track_id];
                        LOG_DEBUG << "FST Hit: idTruth = " << track_id << endm;
                        LOG_DEBUG << "Adding McTrack to FST hit" << endm;
                    }

                    // we use d+4 so that both FTT and FST start at 4 for historical reasons
                    mFwdHitsFst.push_back(
                        FwdHit(
                            count++, // id
                            x0, y0, vZ, // position
                            diskIndex+4, // volume id
                            kFstId, // detid
                            track_id, // mc track id
                            hitCov3, // covariance matrix
                            mcTrack // mcTrack
                        )
                    );
                    // store a pointer to the original StFstHit
                    mFwdHitsFst.back()._hit = fsthits[ih];
                    mFwdHitsFst.back()._genfit_plane_index = globalIndex;
                }
            } // loop is
        } // loop iw
        LOG_DEBUG << " FOUND " << mFstHits.size() << " FST HITS in StFstHitCollection" << endm;
    } // fstHitCollection
    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
        // add to MC track map
        if ( hit->getMcTrack() ){
            hit->getMcTrack()->addFstHit(hit);
        }
    }
    if ( numFwdHitsPost != numFwdHitsPrior ){
        LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from StEvent" << endm;
    }
    return count;
} //loadFstHitsFromStEvent

int StFwdHitMaker::loadFstHitsFromStRnDHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    LOG_DEBUG << "Looking for FST hits in StEvent StRnDHit Collection" << endm;
    int count = 0;
    // Get the StEvent handle
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (!event) {
        LOG_DEBUG << "No StEvent, cannot load FST FastSim hits from StEvent StRndHitCollection" << endm;
        return 0;
    }

    size_t numFwdHitsPrior = mFwdHitsFst.size();
    StRnDHitCollection *rndCollection = event->rndHitCollection();
    if (!rndCollection) return 0;

    const StSPtrVecRnDHit &hits = rndCollection->hits();

    // we will reuse this to hold the cov mat
    TMatrixDSym hitCov3(3);

    for (unsigned int fsthit_index = 0; fsthit_index < hits.size(); fsthit_index++) {
        StRnDHit *hit = hits[fsthit_index];

        if ( hit->layer() > 6 ){
            // skip sTGC hits here
            continue;
        }

        const StThreeVectorF pos = hit->position();

        StMatrixF covmat = hit->covariantMatrix();

        // copy covariance matrix element by element from StMatrixF
        hitCov3(0,0) = covmat[0][0]; hitCov3(0,1) = covmat[0][1]; hitCov3(0,2) = covmat[0][2];
        hitCov3(1,0) = covmat[1][0]; hitCov3(1,1) = covmat[1][1]; hitCov3(1,2) = covmat[1][2];
        hitCov3(2,0) = covmat[2][0]; hitCov3(2,1) = covmat[2][1]; hitCov3(2,2) = covmat[2][2];

        mFwdHitsFst.push_back(
            FwdHit(
                count++, // id
                hit->position().x(), hit->position().y(), hit->position().z(), // position
                hit->layer(), // volume id
                kFstId, // detid
                hit->idTruth(), // mc track id
                hitCov3, // covariance matrix
                mcTrackMap[hit->idTruth()] // mcTrack
            )
        );
        mFstHits.push_back( TVector3( hit->position().x(), hit->position().y(), hit->position().z())  );
    }

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
    }
    if ( numFwdHitsPost != numFwdHitsPrior ){
        LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from StEvent FastSim" << endm;
    }

    return count;
} //loadFstHitsFromStEvent

int StFwdHitMaker::loadFstHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){
    int count = 0;
    LOG_DEBUG << "Looking for FST hits in geant struct" << endm;
    /************************************************************/
    // Load FSI Hits from GEANT
    St_g2t_fts_hit *g2t_fsi_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");

    if ( !g2t_fsi_hits ){
        LOG_DEBUG << "No g2t_fts_hits, cannot load FST hits from GEANT" << endm;
        return 0;
    }

    int nfsi = g2t_fsi_hits->GetNRows();
    size_t numFwdHitsPrior = mFwdHitsFst.size();

    // reuse this to store cov mat
    TMatrixDSym hitCov3(3);

    for (int i = 0; i < nfsi; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_fsi_hits->At(i);

        if (0 == git)
            continue; // geant hit

        int track_id = git->track_p;
        int volume_id = git->volume_id;  // 4, 5, 6
        int d = volume_id / 1000;        // disk id

        int plane_id = d - 4;
        float x = git->x[0];
        float y = git->x[1];
        float z = git->x[2];

        if (mSiRasterizer->active()) {
            TVector3 rastered = mSiRasterizer->raster(TVector3(git->x[0], git->x[1], git->x[2]));
            LOG_INFO << TString::Format("Rastered: %f %f %f -> %f %f %f", git->x[0], git->x[1], git->x[2], rastered.X(), rastered.Y(), rastered.Z()) << endm;
            x = rastered.X();
            y = rastered.Y();
        } else {
            LOG_INFO << "Using GEANT FST hit positions without rasterization" << endm;
        }

        if (plane_id > 3 || plane_id < 0) {
            continue;
        }

        hitCov3 = makeSiCovMat( TVector3( x, y, z ), mFwdConfig );
        mFwdHitsFst.push_back(
            FwdHit(
                count++, // id
                x, y, z, // position
                d, // volume id
                kFstId, // detid
                track_id, // mc track id
                hitCov3, // covariance matrix
                mcTrackMap[track_id] // mcTrack
            )
        );
        mFstHits.push_back( TVector3( x, y, z )  );
    }

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);

        // add to MC track map
        if ( hit->getMcTrack() )
            hit->getMcTrack()->addFstHit(hit);
    }
    if ( numFwdHitsPost != numFwdHitsPrior ){
        LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from GEANT" << endm;
    }

    return count;
} // loadFstHitsFromGEANT
