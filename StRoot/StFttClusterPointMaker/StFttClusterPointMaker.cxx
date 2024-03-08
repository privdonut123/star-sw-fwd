/***************************************************************************
 *
 * StFttClusterPointMaker.cxx
 *
 * Author: Zhen 2023
 ***************************************************************************
 *
 * Description: StFttClusterPointMaker - class to fill the StFttPoint in StEvent
 * this version save all the clusters as point, with prescise resolution in the direction perp to the strip
 * 
 * 
 *
 ***************************************************************************/
#include <vector>
#include <map>
#include <array>
#include <algorithm>


#include "StEvent.h"
#include "StEnumerations.h"

#include "StFttClusterPointMaker.h"


#include "StEvent/StFttRawHit.h"
#include "StEvent/StFttCluster.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"

#include "StFttDbMaker/StFttDb.h"


//_____________________________________________________________
StFttClusterPointMaker::StFttClusterPointMaker( const char* name )
: StMaker( name ),
  mEvent( 0 ),          /// pointer to StEvent
  mDebug( false ),       /// print out of all full messages for debugging
//   mDebug( true ),       /// print out of all full messages for debugging
  mUseTestData( false ),
  mFttDb( nullptr )
{
    LOG_DEBUG << "StFttClusterPointMaker::ctor"  << endm;
    LOG_INFO << "******** StFttClusterPointMaker::StFttClusterPointMaker = "<<name<<endm;
}

//_____________________________________________________________
StFttClusterPointMaker::~StFttClusterPointMaker()
{  /* no op */ }

//_____________________________________________________________
Int_t
StFttClusterPointMaker::Init()
{
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterPointMaker::InitRun( Int_t runnumber )
{ 
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterPointMaker::FinishRun( Int_t runnumber )
{ 
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterPointMaker::Finish()
{ 
    return kStOk;
}


//_____________________________________________________________
Int_t
StFttClusterPointMaker::Make()
{ 
    LOG_DEBUG << "StFttClusterPointMaker::Make()" << endm;

    mEvent = (StEvent*)GetInputDS("StEvent");
    if(mEvent) {
        LOG_DEBUG<<"Found StEvent"<<endm;
    } else {
        return kStOk;
    }
    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        return kStOk;
    }

    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb"));

    assert( mFttDb );

    if ( mUseTestData )
        InjectTestData();

    // next we will need them in even more detail
    // per strip group, but start here
    // Zhen: remove the defination to the .h files
    for (int i = 0; i<16; i++)
    {
        for (int j = 0; j<4; j++)
        {
            clustersPerRob[i][j].clear();
        }
    }

    for ( StFttCluster* clu : mFttCollection->clusters() ) 
    {
        // group clusters by quadrant, hor, vert, hdiag, vdiag
       

        UChar_t rob = mFttDb->rob( clu );
        // if ( clu->nStrips() < 2 ) continue;// add cluster width limit
         if(mDebug)
        {
            LOG_INFO << "rob = " << (int)rob << endm;
            LOG_INFO << "direction = " << (int)clu->orientation() << endm;
            LOG_INFO << "cluster x = " << clu->x() << endm;
        }
        //group clusters in differnt rob and orientation
        clustersPerRob[ (int)rob ][ clu->orientation() ].push_back( clu );// clustersPerRob[ rob ] [ orientation ]
        // clustersPerRob[ clu->orientation() ] [ rob ].push_back( clu ); 
    } // loop on cluster

    for (int i = 0; i < 16; i++) 
    {
        if (mDebug) 
        {
            LOG_INFO << "Now at ROB " << i << endm;
            LOG_INFO << "nCluster kFttVertical = " << clustersPerRob[ i ][ kFttVertical ].size() << endm;
            LOG_INFO << "nCluster kFttHorizontal = " << clustersPerRob[ i ][ kFttHorizontal ].size() << endm;
            LOG_INFO << "nCluster kFttDiagonalV = " << clustersPerRob[ i ][ kFttDiagonalV ].size() << endm;
            LOG_INFO << "nCluster kFttDiagonalH = " << clustersPerRob[ i ][ kFttDiagonalH ].size() << endm;
        }

        MakeLocalPoints((UChar_t)i); // make local points for each Quadrand
    }
    MakeGlobalPoints();
    LOG_INFO << "StFttClusterPointMaker made " << mFttCollection->numberOfPoints() << " points this event" << endm;

    return kStOk;
}

void StFttClusterPointMaker::InjectTestData(){
    mFttCollection->rawHits().clear();

    // TODO: inject clean strip hits to test cluster finder
    // StFttRawHit *hit = new StFttRawHit( sec, rdo, feb, vm, vmm[0].ch, vmm[0].adc, vmm[0].bcid, vmm[0].tb );
    // hit->setMapping( plane, quadrant, row, strip ) 
}

void StFttClusterPointMaker::MakeGlobalPoints() {
    for ( StFttPoint * p : mFttCollection->points() ){

        float x = p->x();
        float y = p->y();
        float z = 0.0;
        StThreeVectorD global;

        // dx is a local shift
        float dx = 0, dy = 0, dz = 0;
        // sx is only {1,-1} -> reflected or normal
        float sx = 0, sy = 0, sz = 0;
        mFttDb->getGloablOffset( p->plane(), p->quadrant(), dx, sx, dy, sy, dz, sz );
        global.set( (x) * sx +dx, (y) * sy +dy, (z + dz) * sz );
        p->setXYZ( global );
    }
}

//--------------------------------------------------------------
//for  the loacl coordinate, if using the the center of pin hole as (0,0)
//center of first strip of V&H strips is 15.95mm
//center of first strip of dia strips is 19.42mm
//in the cluster point maker, the thing I need to do is just loop all the clusters 
void StFttClusterPointMaker::MakeLocalPoints(UChar_t Rob)
{
    //init the point 
    StFttPoint* point;
    double x = -999.;
    double y = -999.;
    //init the number of cluster
    size_t nClusters_X = 0;size_t nClusters_Y = 0;size_t nClusters_DX = 0;size_t nClusters_DY = 0;
    nClusters_X = clustersPerRob[(UChar_t)Rob][kFttVertical].size();
    nClusters_Y = clustersPerRob[(UChar_t)Rob][kFttHorizontal].size();

    if(mDebug)
    {
        LOG_INFO << "rob = " << (int)Rob << endm;
        LOG_INFO << "nClusterX = " << nClusters_X << " nClusterY = " << nClusters_Y << endm;
        LOG_INFO << "nCluster dV = " << clustersPerRob[(UChar_t)Rob][kFttDiagonalV].size()<< endm;
        LOG_INFO << "nCluster dH = " << clustersPerRob[(UChar_t)Rob][kFttDiagonalH].size()<< endm;
    }

    //loop all the clusters and save the cluster information
    for (int iClu_X = 0; iClu_X < nClusters_X; iClu_X++)
    {
        point = new StFttPoint();
        //get the x information of cluster
        StFttCluster* clu_x = clustersPerRob[(UChar_t)Rob][kFttVertical][iClu_X];
        if (clu_x->row() > 2)
        {
            LOG_ERROR << "Wrong cluster row information for cluster at Rob " << Rob << " direction : " << mFttDb->Direction_name[kFttVertical] <<  "!!!!!!" << endm;
            continue;
        }
        
        //x from cluster reconstruction, y at the center of strip with max ADC
        point->setX(clu_x->x());
        point->setY( mFttDb->YX_StripGroupEdge[clu_x->row()]+clu_x->maxStripLength()/2. );
        //x sigma from the cluster reconstruction, y sigma from uniform distribution
        point->setSigmaX(clu_x->sigma());
        point->setSigmaY(clu_x->maxStripLength()/sqrt(12));
        point->setPlane(clu_x->plane());
        point->setQuadrant(clu_x->quadrant());
        point->addCluster(clu_x,kFttVertical);
        mFttPoint.push_back(point);
        mFttCollection->addPoint(point);
    }

    //loop all the clusters and save tht cluster information
    for (int iClu_Y = 0; iClu_Y < nClusters_Y; iClu_Y++)
    {
        point = new StFttPoint();
        StFttCluster* clu_y = clustersPerRob[(UChar_t)Rob][kFttHorizontal][iClu_Y];
        if (clu_y->row() > 2)
        {
            LOG_ERROR << "Wrong cluster row information for cluster at Rob " << Rob << " direction : " << mFttDb->Direction_name[kFttHorizontal] <<  "!!!!!!" << endm;
            continue;
        }

        //y from reconstruction, x at the center of strip with max ADC
        point->setX( mFttDb->YX_StripGroupEdge[clu_y->row()]+clu_y->maxStripLength()/2. );
        point->setY( clu_y->x() );
        //y sigma from cluster reconstruction, x sigma from uniform distribution
        point->setSigmaX(clu_y->maxStripLength()/sqrt(12.));
        point->setSigmaY(clu_y->sigma());
        point->setPlane(clu_y->plane());
        point->setQuadrant(clu_y->quadrant());
        point->addCluster(clu_y,kFttHorizontal);
        mFttPoint.push_back(point);
        mFttCollection->addPoint(point);
    }
    
    
    // need to discuss with Daniel about how to sace the diagonal cluster 
    //loop all the diagonalH clusters and save the cluster information
    for (int iClu_DH = 0; iClu_DH < clustersPerRob[(UChar_t)Rob][kFttDiagonalH].size(); iClu_DH++)
    {
        point = new StFttPoint();
        StFttCluster* clu_DH = clustersPerRob[(UChar_t)Rob][kFttDiagonalH][iClu_DH];
        if (clu_DH->row() <= 2)
        {
            LOG_ERROR << "Wrong cluster row information for cluster at Rob " << Rob << " direction : " << mFttDb->Direction_name[kFttDiagonalH] <<  "!!!!!!" << endm;
            continue;
        }

        //TODO: Now, for clusters in the diagonal direction, the X-axis is along the diagonal direction in the first quadrant, not the horizontal direction. 
        //TODO: Do I need to do the transfer of coordinate now or using it?
        //TODO: now the y all set to positive, how to set the region x>y and the region x<y
        point->setX( clu_DH->x() );
        point->setY( clu_DH->maxStripLength()/2. );
        //y sigma from cluster reconstruction, x sigma from uniform distribution
        point->setSigmaX(clu_DH->sigma());
        point->setSigmaY(clu_DH->maxStripLength()/sqrt(12));
        point->setPlane(clu_DH->plane());
        point->setQuadrant(clu_DH->quadrant());
        point->addCluster(clu_DH,kFttDiagonalH);
        mFttPoint.push_back(point);
        mFttCollection->addPoint(point);
    }

    // need to discuss with Daniel about how to sace the diagonal cluster 
    //loop all the diagonalH clusters and save the cluster information
    for (int iClu_DV = 0; iClu_DV < clustersPerRob[(UChar_t)Rob][kFttDiagonalV].size(); iClu_DV++)
    {
        point = new StFttPoint();
        StFttCluster* clu_DV = clustersPerRob[(UChar_t)Rob][kFttDiagonalV][iClu_DV];
        if (clu_DV->row() <= 2)
        {
            LOG_ERROR << "Wrong cluster row information for cluster at Rob " << Rob << " direction : " << mFttDb->Direction_name[kFttDiagonalV] <<  "!!!!!!" << endm;
            continue;
        }

        //TODO: Now, for clusters in the diagonal direction, the X-axis is along the diagonal direction in the first quadrant, not the horizontal direction. 
        //TODO: Do I need to do the transfer of coordinate now or using it?
        //TODO: now the y all set to positive, how to set the region x>y and the region x<y
        point->setX( clu_DV->x() );
        point->setY( clu_DV->maxStripLength()/2. );
        //y sigma from cluster reconstruction, x sigma from uniform distribution
        point->setSigmaX(clu_DV->sigma());
        point->setSigmaY(clu_DV->maxStripLength()/sqrt(12));
        point->setPlane(clu_DV->plane());
        point->setQuadrant(clu_DV->quadrant());
        point->addCluster(clu_DV,kFttDiagonalV);
        mFttPoint.push_back(point);
        mFttCollection->addPoint(point);
    }
    
}


