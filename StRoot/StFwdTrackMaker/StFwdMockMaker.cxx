#include "StFwdTrackMaker/StFwdMockMaker.h"
#include "StFwdMockMaker.h"

#include "StEvent/StEvent.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StEvent/StFwdTrack.h"

ClassImp(StFwdMockMaker);

#include <fstream>

#include <iomanip>
#include "StMemStat.h"  // Ensure this is available and correctly included

#ifndef __CINT__
#include <chrono>
class MemoryLogger {
    public:
        MemoryLogger(const std::string& filename)
            : start_time_(std::chrono::steady_clock::now()), log_file_(filename, std::ios::app) {
            if (log_file_.tellp() == 0) {
                log_file_ << "# Time(s)\tProgSize(MB)\n";
            }
        }
    
        void log(const std::string& comment = "") {
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - start_time_;
            double progSizeMB = StMemStat::ProgSize();  // assuming kB
    
            log_file_ << std::fixed << std::setprecision(6)
                      << elapsed.count() << "\t" << progSizeMB;
    
            if (!comment.empty()) {
                log_file_ << "\t# " << comment;
            }
    
            log_file_ << "\n";
        }
    
    private:
        std::chrono::steady_clock::time_point start_time_;
        std::ofstream log_file_;
    };

    MemoryLogger MemLogger2("mem_log2.txt");

    #endif // CINT


int StFwdMockMaker::Make() {
    MemLogger2.log("StFwdMockMaker::Make() [");
    LOG_INFO << "StFwdMockMaker::Make()" << endm;

    mFwdHitsFtt.clear();
    for ( int i = 0; i < 1000; i++ ) {
        mFwdHitsFtt.push_back( TVector3( i, i, i ) );
    }

    FillEvent();

    MemLogger2.log("] StFwdMockMaker::Make()");
    return kStOk;
}

StFwdTrack * StFwdMockMaker::makeStFwdTrack( size_t indexTrack ){
    LOG_DEBUG << "StFwdMockMaker::makeStFwdTrack()" << endm;
    StFwdTrack *fwdTrack = new StFwdTrack( );

    TVector3 tv3, mom;
    float cov[9] = { 0.01, 0, 0,
                     0, 0.01, 0,
                     0, 0, 0.01 };
    for ( int i = 0; i < 100; i++){
        tv3.SetXYZ( i, i, i );
        mom.SetXYZ( i, i, i );
        fwdTrack->mProjections.push_back( StFwdTrackProjection( 0, StThreeVectorF( tv3.X(), tv3.Y(), tv3.Z() ), StThreeVectorF( mom.X(), mom.Y(), mom.Z() ), cov) );
    }

    for ( int i = 0; i < 100; i++ ){
        StFwdTrackSeedPoint p(
            StThreeVectorD( i, i, i ),
            0, // detectorId * 10 + sector
            indexTrack,
            cov
        );
        fwdTrack->mFTTPoints.push_back( p );
    }

    return fwdTrack;
} // makeStFwdTrack


void StFwdMockMaker::FillEvent() {
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if ( !ftc ){
        LOG_INFO << "Creating the StFwdTrackCollection" << endm;
        ftc = new StFwdTrackCollection();
        stEvent->setFwdTrackCollection( ftc );
    }

    size_t indexTrack = 0;
    for ( int i = 0; i < 1000; i++ ) {
            StFwdTrack* fwdTrack = makeStFwdTrack( indexTrack );
            indexTrack++;
            if (nullptr == fwdTrack)
                continue;
            ftc->addTrack( fwdTrack );
    }

    LOG_INFO << "StFwdTrackCollection has " << ftc->numberOfTracks() << " tracks now" << endm;

}