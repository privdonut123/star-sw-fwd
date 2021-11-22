#ifndef STFTTCLUSTER_H
#define STFTTCLUSTER_H

#include "StObject.h"

#include "StContainers.h"  // For StPtrVecFttRawHit
#include "StEnumerations.h"

class StFttRawHit;
class StFttDb;

class StFttCluster : public StObject {
public:
    StFttCluster();
    ~StFttCluster();
    
    int id() const; // Cluster ID
    UChar_t orientation() const;
    int nStrips() const;
    int nNeighbors() const;
    float sumAdc() const;
    float x() const;  // Mean x ("center of gravity") in local grid coordinate (1st moment).
    float sigma() const; // Maximum 2nd moment (along major axis).

    void setId(int cluid);
    void setOrientation( UChar_t );
    void setNStrips(int numStrips);
    void setSumAdc(int theSumAdc);
    void setX(float x0);
    void setSigma(float sigma);

    StPtrVecFttRawHit& hits();
    const StPtrVecFttRawHit& hits() const;
    void addNeighbor(StFttCluster* neighbor);
    StPtrVecFttCluster& neighbor();
    const StPtrVecFttCluster& neighbor() const;    
    // StPtrVecFcsPoint& points();
    // const StPtrVecFcsPoint& points() const;
    // void addPoint(StFcsPoint* p);
    // void addPoint(StFcsPoint* p1, StFcsPoint* p2);
    // void print(Option_t *option="") const;

private:
    Int_t mId=-1;             // Eventwise cluster ID
    UChar_t mOrientation = kFttUnknownOrientation;        // Orientation of cluster
    Int_t mNStrips=0;         // Number of strips
    Float_t mSumAdc=0.0;      // Total ADC (0th moment)
    Float_t mX=0.0;             // Mean x ("center of gravity") in local grid coordinate (1st moment)
    Float_t mSigma=0.0;        // 2nd moment
    StPtrVecFttRawHit mHits;            // Tower hits of the current cluster
    StPtrVecFttCluster mNeighbors;    // Neighbor clusters
    // StPtrVecFcsPoint mPoints;        // Fitted points (photons) in the cluster

    ClassDef(StFttCluster, 1)
};


inline int StFttCluster::id() const { return mId; } // Cluster ID
inline UChar_t StFttCluster::orientation() const { return mOrientation; }
inline int StFttCluster::nStrips() const { return mNStrips; }
inline int StFttCluster::nNeighbors() const { return mNeighbors.size(); }
// inline int StFttCluster::nPoints() const { return mPoints.size(); }
inline float StFttCluster::sumAdc() const { return mSumAdc; }
inline float StFttCluster::x() const { return mX; } // Mean x ("center of gravity") in local grid coordinate (1st moment).
inline float StFttCluster::sigma() const { return mSigma; } // 2nd moment

inline void StFttCluster::setOrientation( UChar_t so ) { mOrientation = so; }
inline void StFttCluster::setNStrips(int numStrips) { mNStrips = numStrips; }
inline void StFttCluster::setSumAdc(int theSumAdc) { mSumAdc = theSumAdc; }
inline void StFttCluster::setX(float x0) { mX = x0; }
inline void StFttCluster::setSigma(float sigma) { mSigma = sigma; }

inline void StFttCluster::setId(int cluid) { mId = cluid; }

inline StPtrVecFttRawHit& StFttCluster::hits() { return mHits; }
inline const StPtrVecFttRawHit& StFttCluster::hits() const { return mHits; }
inline StPtrVecFttCluster& StFttCluster::neighbor() { return mNeighbors; }
inline const StPtrVecFttCluster& StFttCluster::neighbor() const { return mNeighbors; }
// inline StPtrVecFcsPoint& StFttCluster::points() { return mPoints; }
// inline const StPtrVecFcsPoint& StFttCluster::points() const { return mPoints; }

#endif  // STFTTCLUSTER_H
