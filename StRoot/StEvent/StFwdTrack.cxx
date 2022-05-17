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

const StThreeVectorF StFwdTrack::momentum() const{
    auto cr = mGenfitTrack->getCardinalRep();
    TVector3 p = cr->getMom(mGenfitTrack->getFittedState(0, cr));
    return StThreeVectorF( p.X(), p.Y(), p.Z() );
}
const StThreeVectorF StFwdTrack::momentumAt(int _id) const{
    if (!mGenfitTrack)
        return StThreeVectorF( 0, 0, 0 );
    
    int id = _id;
    if (id > mGenfitTrack->getNumPoints ()){
        id = mGenfitTrack->getNumPoints();
    }
    {
        auto cr = mGenfitTrack->getCardinalRep();
        TVector3 p = cr->getMom(mGenfitTrack->getFittedState(id, cr));
        return StThreeVectorF( p.X(), p.Y(), p.Z() );
    // } catch ( int e ) {
        LOG_WARN << " Cannot get momentum at point: " << id << endm;
    }
    return StThreeVectorF( 0, 0, 0 );
}
const char StFwdTrack::charge() const{
    {
        auto cr = mGenfitTrack->getCardinalRep();
        return cr->getCharge(mGenfitTrack->getFittedState(0, cr)); // at the first fit point
        LOG_WARN << "Cannot get track charge " << endm;
    }
    return -99;
}