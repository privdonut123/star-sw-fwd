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

ClassImp( StFwdTrack )

StFwdTrack::StFwdTrack( genfit::Track *t ) : mGenfitTrack(t), mProjections(0) {

}

