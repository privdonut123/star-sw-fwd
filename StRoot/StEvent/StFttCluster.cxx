/***************************************************************************
 *
 * $Id: StFttCluster.cxx,v 1.0 2021/11/18 14:53:48 jdb Exp $
 *
 * Author: jdb, Nov 2021
 ***************************************************************************
 *
 * Description: 
 *
 ***************************************************************************/ 
#include "StEvent/StFttCluster.h"


StFttCluster::StFttCluster() :
mId(-1),
mOrientation(kFttUnknownOrientation),
mNStrips(0),
mSumAdc(0.0),
mX(0.0),
mSigma(0.0),
mHits(0),            // Tower hits of the current cluster
mNeighbors(0)
{

}


StFttCluster::~StFttCluster(){}