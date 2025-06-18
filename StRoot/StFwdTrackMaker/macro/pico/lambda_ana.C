#include "StPicoMcTrack.h"
#include "StPicoMcVertex.h"
#include "StPicoFwdTrack.h"
#include "StPicoFwdVertex.h"
#include "StPicoFcsHit.h"
#include "StPicoFcsCluster.h"

bool inAcc(StPicoMcTrack *track) {
    // Example acceptance criteria
    return (track->nHitsFts() > 2 || track->nHitsStg() > 4);
}

int testTrackType = 1; // 0=Global, 1=BLC, 2=Primary, 3=FwdVertex

// make a TClonesArray for McTrack
TClonesArray *mcTracks = new TClonesArray("StPicoMcTrack", 1000);
TClonesArray *mcVertices = new TClonesArray("StPicoMcVertex", 1000);
TClonesArray *fwdTracks = new TClonesArray("StPicoFwdTrack", 1000);
TClonesArray *fwdVertices = new TClonesArray("StPicoFwdVertex", 1000);
TClonesArray *fcsHits = new TClonesArray("StPicoFcsHit", 1000);
TClonesArray *fcsClusters = new TClonesArray("StPicoFcsCluster", 1000);

void lambda_ana(){

    TFile *fOutput = new TFile("lambda_analysis.root", "RECREATE");
    TH2F * hNumFstStgc = new TH2F("hNumFstStgc", "Number of FST and STGC hits; nFST; nSTGC", 15, 0, 15, 15, 0, 15);
    TH1F * hNumFST = new TH1F("hNumFST", "Number of FST hits", 100, 0, 100);
    TH1F * hNumHCAL = new TH1F("hNumHCAL", "Number of HCAL hits", 100, 0, 100);
    TH1F * hNumWCAL = new TH1F("hNumWCAL", "Number of WCAL hits", 100, 0, 100);
    TH1F * hNumSTGC = new TH1F("hNumSTGC", "Number of STGC hits", 100, 0, 100);
    TH1F * hM0 = new TH1F("hM0", "Mass of all tracks", 100, 1, 2);
    TH1F * hFwdM0 = new TH1F("hFwdM0", "Mass of all tracks", 100, 1, 2);
    TH1F * hFwdM0Correct = new TH1F("hFwdM0Correct", "Mass of all tracks", 100, 1, 2);
    TH1F * hFwdM0Incorrect = new TH1F("hFwdM0Incorrect", "Mass of all tracks", 100, 1, 2);
    TH1F * hEta = new TH1F("hEta", "Pseudorapidity of tracks", 100, -5, 5);
    TH1F * hEtaAcc = new TH1F("hEtaAcc", "Pseudorapidity of tracks that hit FWD", 100, -5, 5);
    TH1F * hIdTruthEtaAcc = new TH1F("hIdTruthEtaAcc", "Pseudorapidity of tracks that hit FWD", 100, -5, 5);

    TH1F * hPt = new TH1F("hPt", "Transverse momentum of tracks", 100, 0, 10);
    TH1F * hPtAcc = new TH1F("hPtAcc", "Transverse momentum of tracks that hit FWD", 100, 0, 10);
    TH1F * hIdTruthPtAcc = new TH1F("hIdTruthPtAcc", "Transverse momentum of tracks that hit FWD", 100, 0, 10);


    TChain *chain = new TChain("PicoDst");
    chain->Add("*.picoDst.root");

    chain->SetBranchAddress("McTrack", &mcTracks);
    chain->SetBranchAddress("McVertex", &mcVertices);
    chain->SetBranchAddress("FwdTracks", &fwdTracks);
    chain->SetBranchAddress("FwdVertices", &fwdVertices);
    chain->SetBranchAddress("FcsHits", &fcsHits);
    chain->SetBranchAddress("FcsClusters", &fcsClusters);

    size_t nEntries = chain->GetEntries();
    cout << "Number of events: " << nEntries << endl;
    for (size_t i = 0; i < nEntries; ++i) {
        mcTracks->Clear(); // Clear the TClonesArray for each entry
        mcVertices->Clear(); // Clear the TClonesArray for each entry
        fwdTracks->Clear(); // Clear the TClonesArray for each entry
        fwdVertices->Clear(); // Clear the TClonesArray for each entry
        fcsHits->Clear(); // Clear the TClonesArray for each entry
        fcsClusters->Clear(); // Clear the TClonesArray for each entry

        chain->GetEntry(i);
        if (i % 1000 == 0) {
            cout << "Processing entry: " << i << endl;
        }

        // Loop over the McTrack objects in the TClonesArray
        for (int j = 0; j < mcTracks->GetEntriesFast(); ++j) {
            StPicoMcTrack *t1 = static_cast<StPicoMcTrack*>(mcTracks->At(j));
            if (!t1) continue;
            if ( t1->id() > 2 ) continue; // Skip tracks with id > 2
            TLorentzVector lv1 = t1->fourMomentum();
            hNumFstStgc->Fill(t1->nHitsFts(), t1->nHitsStg());
            StPicoMcVertex *vtx1 = static_cast<StPicoMcVertex*>(mcVertices->At(t1->idVtxStart() - 1));
            if (!vtx1) continue;

            // printf("Vertex ID: %d, Position: (%.2f, %.2f, %.2f)\n", 
            //        vtx1->id(), vtx1->position().X(), vtx1->position().Y(), vtx1->position().Z());

            for (int k = j; k < mcTracks->GetEntriesFast(); ++k) {
                if ( j == k ) continue; // Skip self-comparison
                StPicoMcTrack *t2 = static_cast<StPicoMcTrack*>(mcTracks->At(k));
                if (!t2) continue;
                StPicoMcVertex *vtx2 = static_cast<StPicoMcVertex*>(mcVertices->At(t2->idVtxStart() - 1));
                if (!vtx2) continue;
                
                if ( t2->id() > 2 ) continue; // Skip tracks with id > 2
                TLorentzVector lv2 = t2->fourMomentum();
    
                TLorentzVector lv = lv1 + lv2;
                hM0->Fill(lv.M());
                hEta->Fill(lv.PseudoRapidity());
                hPt->Fill(lv.Pt());
                if (inAcc(t1) && inAcc(t2)) {
                    hEtaAcc->Fill(lv.PseudoRapidity());
                    hPtAcc->Fill(lv.Pt());
                }
            }
        }

        // Loop over the FwdTracks
        for (int j = 0; j < fwdTracks->GetEntriesFast(); ++j) {
            StPicoFwdTrack *fwdTrack = static_cast<StPicoFwdTrack*>(fwdTracks->At(j));
            if ( fwdTrack->trackType() != testTrackType ) continue; // Only consider global tracks
            if (!fwdTrack) continue;

            // get matched Mc Track
            if ( fwdTrack->idTruth() < 1 || fwdTrack->idTruth() > 2 ) continue; // Skip tracks without MC match
            auto mcTrack = static_cast<StPicoMcTrack*>(mcTracks->At(fwdTrack->idTruth() - 1));

            for (int k = j; k < fwdTracks->GetEntriesFast(); k++) {
                if ( j == k ) continue; // Skip self-comparison
                StPicoFwdTrack *fwdTrack2 = static_cast<StPicoFwdTrack*>(fwdTracks->At(k));
                if (!fwdTrack2) continue;
                if ( fwdTrack2->trackType() != testTrackType ) continue; // Only consider global tracks

                // get matched Mc Track
                if ( fwdTrack2->idTruth() < 1 || fwdTrack2->idTruth() > 2 ) continue; // Skip tracks without MC match
                auto mcTrack2 = static_cast<StPicoMcTrack*>(mcTracks->At(fwdTrack2->idTruth() - 1));

                TLorentzVector mclv1 = mcTrack->fourMomentum();
                TLorentzVector mclv2 = mcTrack2->fourMomentum();
                TLorentzVector mclv = mclv1 + mclv2;
                
                if (inAcc(mcTrack) && inAcc(mcTrack2)) {
                    hIdTruthEtaAcc->Fill(mclv.PseudoRapidity());
                    hIdTruthPtAcc->Fill(mclv.Pt());
                }

                // printf( "Matched mcTrack GeantId(): mcTrack1=%d, mcTrack2=%d\n", mcTrack->geantId(), mcTrack2->geantId() );

                
                TLorentzVector fwdlv1 = fwdTrack->fourMomentum( 0.13957); // Assuming pion mass for lv calculation);
                TLorentzVector fwdlv2 = fwdTrack2->fourMomentum( 0.93827); // Assuming proton mass for lv calculation);
                TLorentzVector fwdlv = fwdlv1 + fwdlv2;
                hFwdM0->Fill(fwdlv.M());
                

                TLorentzVector fwdlv1B = fwdTrack2->fourMomentum( 0.13957); // Assuming pion mass for lv calculation);
                TLorentzVector fwdlv2B = fwdTrack->fourMomentum( 0.93827); // Assuming proton mass for lv calculation);
                TLorentzVector fwdlvB = fwdlv1B + fwdlv2B;
                hFwdM0->Fill(fwdlvB.M());

                if ( mcTrack->geantId() == 9 && mcTrack2->geantId() == 14 ){
                    hFwdM0Correct->Fill(fwdlv.M());
                    hFwdM0Incorrect->Fill(fwdlvB.M());
                } else if ( mcTrack2->geantId() == 9 && mcTrack->geantId() == 14 ) {
                    hFwdM0Incorrect->Fill(fwdlv.M());
                    hFwdM0Correct->Fill(fwdlvB.M());
                }
            }

        }

    } // end of loop over entries


    return;
}