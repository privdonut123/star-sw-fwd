/***************************************************************************
 *
 * $Id: StFttRawHit.cxx,v 2.1 2021/11/18 14:53:48 jdb Exp $
 *
 * Author: jdb, Nov 2021
 ***************************************************************************
 *
 * Description: 
 *
 ***************************************************************************/ 
#include "StFttRawHit.h"
#include <cmath>


StFttRawHit::StFttRawHit() 
: mSector(0),
mRDO(0),
mFEB(0),
mVMM(0),
mChannel(0),
mADC(0),
mBCID(0),
mTB(-32000),
mPlane(255),
mQuadrant(255),
mRow(255),
mStrip(255)
{ /*noop*/ }

StFttRawHit::StFttRawHit(   UChar_t mSector, UChar_t mRDO, UChar_t mFEB, 
                            UChar_t mVMM, UChar_t mChannel, UShort_t mADC, 
                            UShort_t mBCID, Short_t mTB){
    setRaw( mSector, mRDO, mFEB, mVMM, mChannel, mADC, mBCID, mTB);
}

void StFttRawHit::setRaw(   UChar_t mSector, UChar_t mRDO, UChar_t mFEB, 
                            UChar_t mVMM, UChar_t mChannel, UShort_t mADC, 
                            UShort_t mBCID, Short_t mTB){
    this->mSector  = mSector;
    this->mRDO     = mRDO;
    this->mFEB     = mFEB;
    this->mVMM     = mVMM;
    this->mChannel = mChannel;
    this->mADC     = mADC;
    this->mBCID    = mBCID;
    this->mTB      = mTB;
}

void StFttRawHit::setMapping(   UChar_t mPlane, UChar_t mQuadrant, 
                                UChar_t mRow, UChar_t mStrip, UChar_t mOrientation ){
    this->mPlane       = mPlane;
    this->mQuadrant    = mQuadrant;
    this->mRow         = mRow;
    this->mStrip       = mStrip;
    this->mOrientation = mOrientation;
}



ostream&
operator<<( ostream &os, const StFttRawHit& rh )
{
    os << " StFttRawHit( " << endl;
    os << "\tmSector = "   << rh.sector()   << endl;
    os << "\tmRDO = "      << rh.rdo()      << endl;
    os << "\tmFEB = "      << rh.feb()      << endl;
    os << "\tmVMM = "      << rh.vmm()      << endl;
    os << "\tmChannel = "  << rh.channel()  << endl;
    os << "\tmADC = "      << rh.adc()      << endl;
    os << "\tmBCID = "     << rh.bcid()     << endl;
    os << "\tmTB = "       << rh.tb()       << endl;
    os << "\tmPlane = "    << rh.plane()    << endl;
    os << "\tmQuadrant = " << rh.quadrant() << endl;
    os << "\tmRow = "      << rh.row()      << endl;
    os << "\tmStrip = "    << rh.strip()    << " ) " << endl;


    return os;
}