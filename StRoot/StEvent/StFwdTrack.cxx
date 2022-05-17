/***************************************************************************
 *
 * $Id: StFwdTrack.cxx
 *
 * Author: jdb, Feb 2022
 ***************************************************************************
 *
 * Description: StFwdTrack stores the Forward tracks built from Fst and Ftt
 *
 ***************************************************************************/

#include "StEvent/StFwdTrack.h"
#include "GenFit/Track.h"
// #include "GenFit/Exception.h"
#include "StMessMgr.h"

ClassImp( StFwdTrack )

StFwdTrack::StFwdTrack( genfit::Track *t ) : mGenfitTrack(t), mProjections(0) {

}

/* momentum
 * get the track momentum at the first point (PV if included)
 */
const StThreeVectorD StFwdTrack::momentum() const{
    auto cr = mGenfitTrack->getCardinalRep();
    TVector3 p = cr->getMom(mGenfitTrack->getFittedState(0, cr));
    return StThreeVectorD( p.X(), p.Y(), p.Z() );
}

/* momentumAt
 * get the track momentum at the nthh point (if available)
 */
const StThreeVectorD StFwdTrack::momentumAt(int _id) const{
    if (!mGenfitTrack)
        return StThreeVectorD( 0, 0, 0 );
    
    int id = _id;
    if (id >= mGenfitTrack->getNumPoints()){
        id = mGenfitTrack->getNumPoints() - 1;
    }

    while ( !mGenfitTrack->getPoint( id ) || !mGenfitTrack->getPoint( id )->getFitterInfo() ){
        id = id -1;
        if ( id < 0 ){
            LOG_WARN << " Cannot get momentum at point: " << _id << endm;
            return StThreeVectorD( 0, 0, 0 );
        }
    } // while valid point

    auto cr = mGenfitTrack->getCardinalRep();
    TVector3 p = cr->getMom(mGenfitTrack->getFittedState(id, cr));
    LOG_INFO << "Momentum at point: " << id << endm;
    return StThreeVectorD( p.X(), p.Y(), p.Z() );
}
const char StFwdTrack::charge() const{
    auto cr = mGenfitTrack->getCardinalRep();
    return cr->getCharge(mGenfitTrack->getFittedState(0, cr)); // at the first fit point
}

const unsigned int StFwdTrack::numberOfFitPoints() const {
    if ( !mGenfitTrack )
        return 0;
    return mGenfitTrack->getNumPoints();
}