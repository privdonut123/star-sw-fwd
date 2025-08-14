#include "StPicoMcTrack.hpp"
#include "StPicoMcVertex.hpp"
#include "StPicoFwdTrack.hpp"
#include "StPicoFwdVertex.hpp"
#include "StPicoFcsHit.hpp"
#include "StPicoFcsCluster.hpp"

bool inAcc(StPicoMcTrack *track) {
    // Example acceptance criteria
    return (track->nHitsFts() > 2);
}

int testTrackType = 2; // 0=Global, 1=BLC, 2=Primary, 3=FwdVertex

// make a TClonesArray for McTrack
TClonesArray *mcTracks = new TClonesArray("StPicoMcTrack", 1000);
TClonesArray *mcVertices = new TClonesArray("StPicoMcVertex", 1000);
TClonesArray *fwdTracks = new TClonesArray("StPicoFwdTrack", 1000);
TClonesArray *fwdVertices = new TClonesArray("StPicoFwdVertex", 1000);
TClonesArray *fcsHits = new TClonesArray("StPicoFcsHit", 1000);
TClonesArray *fcsClusters = new TClonesArray("StPicoFcsCluster", 1000);


void format_eff_plot(TH1* hNum, TH1 * hDen, TString name, TString title, int color = kBlack) {
    TH1 * hEff = (TH1*)hNum->Clone(name);
    hEff->Divide(hNum,hDen, 1, 1, "B");
    hEff->SetTitle(title);
    hEff->SetLineColor(color);
    hEff->SetLineWidth(4);
}

void qa(    TString dataDir = "/Users/brandenburg.89/star/ssw/data/muon_minus/", 
            int tt = -1, // track type
            TString label = "",
            int maxEntries = 50000  // limit the number of entries to process, -1 for all
         ){
    
    TChain *chain = new TChain("PicoDst");
    chain->Add(TString::Format("%s*.picoDst.root", dataDir.Data()));
    if (chain->GetEntries() == 0) {
        printf("No entries found in chain for %s\n", dataDir.Data());
        return;
    }
    chain->SetBranchAddress("McTrack", &mcTracks);
    chain->SetBranchAddress("McVertex", &mcVertices);
    chain->SetBranchAddress("FwdTracks", &fwdTracks);
    chain->SetBranchAddress("FwdVertices", &fwdVertices);
    chain->SetBranchAddress("FcsHits", &fcsHits);
    chain->SetBranchAddress("FcsClusters", &fcsClusters);

    if (tt >= 0) {
        testTrackType = tt; // Set the track type if specified
    }

    TFile *fOutput = new TFile( TString::Format( "efficiency_trackType%d_%s.root", testTrackType, label.Data() ) , "RECREATE");
    fOutput->cd();

    int nEtaBins = 200;
    float etaMin = -5.0;
    float etaMax = 5.0;

    int nPtBins = 40;
    float ptMin = 0.0;
    float ptMax = 2.0;

    int nPhiBins = 100;
    float phiMin = -3.24;
    float phiMax = 3.24;

    int nMultBins = 100;
    float multMin = 0.0;
    float multMax = 100.0;

    TH1 * hMcPrimaryEta = new TH1F("hMcPrimaryEta", "MC Primary Eta", nEtaBins, etaMin, etaMax);
    TH1 * hMcPrimaryEtaInAcc = new TH1F("hMcPrimaryEtaInAcc", "MC Primary (in acc) Eta", nEtaBins, etaMin, etaMax);

    TH1 * hMcPrimaryPt = new TH1F("hMcPrimaryPt", "MC Primary Pt", nPtBins, ptMin, ptMax);
    TH1 * hMcPrimaryPtInAcc = new TH1F("hMcPrimaryPtInAcc", "MC Primary (in acc) Pt", nPtBins, ptMin, ptMax);

    TH1 * hMcPrimaryPhi = new TH1F("hMcPrimaryPhi", "MC Primary Phi", nPhiBins, phiMin, phiMax);
    TH1 * hMcPrimaryPhiInAcc = new TH1F("hMcPrimaryPhiInAcc", "MC Primary (in acc) Phi", nPhiBins, phiMin, phiMax);
    
    TH1 * hMcPrimaryMult = new TH1F("hMcPrimaryMult", "MC Primary Mult", nMultBins, multMin, multMax);
    TH1 * hMcPrimaryMultInAcc = new TH1F("hMcPrimaryMultInAcc", "MC Primary (in acc) Mult", nMultBins, multMin, multMax);

    TH1 * hMatchedMcEta = new TH1F("hMatchedMcEta", "Matched MC Eta", nEtaBins, etaMin, etaMax);
    TH1 * hMatchedMcPt = new TH1F("hMatchedMcPt", "Matched MC Pt", nPtBins, ptMin, ptMax);
    TH1 * hMatchedMcPhi = new TH1F("hMatchedMcPhi", "Matched MC Phi", nPhiBins, phiMin, phiMax);

    TH1 * hMatchedMcPrimaryEta = new TH1F("hMatchedMcPrimaryEta", "Matched MC Primary Eta", nEtaBins, etaMin, etaMax);
    TH1 * hMatchedMcPrimaryPt = new TH1F("hMatchedMcPrimaryPt", "Matched MC Primary Pt", nPtBins, ptMin, ptMax);
    TH1 * hMatchedMcPrimaryPhi = new TH1F("hMatchedMcPrimaryPhi", "Matched MC Primary Phi", nPhiBins, phiMin, phiMax);

    TH1 * hMatchedMcPrimaryMult = new TH1F("hMatchedMcPrimaryMult", "Matched MC Primary Mult", nMultBins, multMin, multMax);

    TH2 * hMatchedMcPrimaryEtaCharge = new TH2F("hMatchedMcPrimaryEtaCharge", "Matched MC Primary Eta vs Charge", nEtaBins, etaMin, etaMax, 3, -0.5, 2.5);
    TH2 * hMatchedMcPrimaryPtCharge = new TH2F("hMatchedMcPrimaryPtCharge", "Matched MC Primary Pt vs Charge", nPtBins, ptMin, ptMax, 3, -0.5, 2.5);
    TH2 * hMatchedMcPrimaryPhiCharge = new TH2F("hMatchedMcPrimaryPhiCharge", "Matched MC Primary Phi vs Charge", nPhiBins, phiMin, phiMax, 3, -0.5, 2.5);

    TH1 * hCurveResolution = new TH1F("hCurveResolution", "Curve Resolution", 100, -5.0, 5.0);
    TH2 * hCurveResolutionEta = new TH2F("hCurveResolutionEta", "Curve Resolution vs Eta", nEtaBins, etaMin, etaMax, 100, -5.0, 5.0);
    TH2 * hCurveResolutionPt = new TH2F("hCurveResolutionPt", "Curve Resolution vs Pt", nPtBins, ptMin, ptMax, 100, -5.0, 5.0);
    TH2 * hCurveResolutionPhi = new TH2F("hCurveResolutionPhi", "Curve Resolution vs Phi", nPhiBins, phiMin, phiMax, 100, -5.0, 5.0);


    TH1 * hQATruth = new TH1F("hQATruth", "QA Truth", 102, -1.0, 101.0);
    TH1 * hQATruthGoodChi2 = new TH1F("hQATruthGoodChi2", "QA Truth", 102, -1.0, 101.0);

    TH1 * hFwdTrackType = new TH1F("hFwdTrackType", "Forward Track Type", 4, -0.5, 3.5);

    size_t fwdPerfect = 0;
    size_t fwdGoodFits = 0;

    size_t nEntries = chain->GetEntries();
    // nEntries = 500000; // Limit to 10,000 entries for testing
    cout << "Number of events: " << nEntries << endl;
    if (maxEntries > 0 && maxEntries < nEntries) {
        nEntries = maxEntries; // Limit the number of entries to process
    }
    cout << "Processing " << nEntries << " entries." << endl;
    for (size_t i = 0; i < nEntries; ++i) {
        mcTracks->Clear(); // Clear the TClonesArray for each entry
        mcVertices->Clear(); // Clear the TClonesArray for each entry
        fwdTracks->Clear(); // Clear the TClonesArray for each entry
        fwdVertices->Clear(); // Clear the TClonesArray for each entry
        fcsHits->Clear(); // Clear the TClonesArray for each entry
        fcsClusters->Clear(); // Clear the TClonesArray for each entry

        chain->GetEntry(i);
        if (i % 10000 == 0) {
            cout << "Processing entry: " << i << endl;
        }

        size_t fwdMult = 0; // Forward multiplicity counter
        size_t mcMult = 0; // MC multiplicity counter
        // Loop over the McTrack objects in the TClonesArray
        for (int j = 0; j < mcTracks->GetEntriesFast(); ++j) {
            StPicoMcTrack *t1 = static_cast<StPicoMcTrack*>(mcTracks->At(j));
            if (!t1) continue;
            // printf("Track with idVtxStart = %d\n", t1->idVtxStart());
            
            TLorentzVector lv1 = t1->fourMomentum();
            StPicoMcVertex *vtx1 = static_cast<StPicoMcVertex*>(mcVertices->At(t1->idVtxStart() - 1));
            if (!vtx1) {
                // printf("Vertex ID: %d not found for track ID: %d\n", t1->idVtxStart(), t1->id());
                continue;
            }

            if ( t1->idVtxStart() != 1 ){
                // printf("Skipping track with idVtxStart != 1: %d\n", t1->idVtxStart());
                continue; // Only consider primary tracks
            }
            mcMult++; // Increment MC multiplicity counter
            hMcPrimaryEta->Fill(lv1.PseudoRapidity());
            hMcPrimaryPt->Fill(lv1.Pt());
            hMcPrimaryPhi->Fill(lv1.Phi());
            if (inAcc(t1)) {
                hMcPrimaryEtaInAcc->Fill(lv1.PseudoRapidity());
                hMcPrimaryPtInAcc->Fill(lv1.Pt());
                hMcPrimaryPhiInAcc->Fill(lv1.Phi());
                fwdMult++;
            }
        }

        size_t fwdMultReco = 0; // Reconstructed forward multiplicity counter
        
        // Loop over the FwdTracks
        for (int j = 0; j < fwdTracks->GetEntriesFast(); ++j) {
            StPicoFwdTrack *fwdTrack = static_cast<StPicoFwdTrack*>(fwdTracks->At(j));
            
            if (!fwdTrack) continue;
            hFwdTrackType->Fill(fwdTrack->trackType());

            if ( fwdTrack->trackType() != testTrackType ) continue; // Only tracks of chosen type


            // get matched Mc Track
            if ( fwdTrack->idTruth() < 1  || fwdTrack->idTruth() > mcTracks->GetEntriesFast() ) {
                // printf("FwdTrack with idTruth < 1 || idTruth > len(mcTracks): %d -> status = %d, chi2 = %f \n", fwdTrack->idTruth(), (int)fwdTrack->status(), fwdTrack->chi2());
                continue;
            }

            hQATruth->Fill(fwdTrack->qaTruth());

            // if ( fwdTrack->chi2() < 0.001 ) continue; // Skip tracks with chi2 < 0.01
            hQATruthGoodChi2->Fill(fwdTrack->qaTruth());

            auto mcTrack = static_cast<StPicoMcTrack*>(mcTracks->At(fwdTrack->idTruth() - 1));
            if (!mcTrack) { 
                printf("No matched MC track for FwdTrack with idTruth: %d\n", fwdTrack->idTruth());
                continue;
            }

            TLorentzVector lvMc = mcTrack->fourMomentum();
            hMatchedMcEta->Fill(lvMc.PseudoRapidity());
            hMatchedMcPt->Fill(lvMc.Pt());
            hMatchedMcPhi->Fill(lvMc.Phi());

            if ( mcTrack->idVtxStart() == 1 && fwdTrack->chi2() > 0.0001 && fabs(fwdTrack->numberOfFitPoints()) > 1 ){
                hMatchedMcPrimaryEta->Fill(lvMc.PseudoRapidity());
                hMatchedMcPrimaryPt->Fill(lvMc.Pt());
                hMatchedMcPrimaryPhi->Fill(lvMc.Phi());

                int weightCharge = 1.0 - (int)(fwdTrack->charge() == mcTrack->charge());
                
                hMatchedMcPrimaryEtaCharge->Fill(lvMc.PseudoRapidity(), weightCharge);
                hMatchedMcPrimaryPtCharge->Fill(lvMc.Pt(), weightCharge);
                hMatchedMcPrimaryPhiCharge->Fill(lvMc.Phi(), weightCharge);

                float fwdCurve = 1.0 / (fwdTrack->momentum().Pt());
                float mcCurve = 1.0 / (lvMc.Pt());
                float curveRes = (fwdCurve - mcCurve) / mcCurve;
                hCurveResolution->Fill(curveRes);
                hCurveResolutionEta->Fill(lvMc.PseudoRapidity(), curveRes);
                hCurveResolutionPt->Fill(lvMc.Pt(), curveRes);
                hCurveResolutionPhi->Fill(lvMc.Phi(), curveRes);

                fwdMultReco++;

                if(fwdTrack->qaTruth()>= 99){
                    fwdPerfect++;
                }
                if ( fwdTrack->chi2() > 0.0001 ){
                    fwdGoodFits++;
                }
            }
        } // end of FwdTracks loop

        if ( fwdMult > 0){
            // compute the average efficiency vs multiplicity
            hMcPrimaryMult->Fill(mcMult);
            hMcPrimaryMultInAcc->Fill(fwdMult);
            hMatchedMcPrimaryMult->Fill(fwdMult, (float)(fwdMultReco/(float)fwdMult));
        }
    } // end of event loop

    printf( "Processed %zu entries, found %0.1f MC Primary tracks (%0.1f in FWD Acc), %0.1f Fwd tracks (%zu perfect seeds, %lu good fits) of type %d\n", 
            nEntries, 
            hMcPrimaryEta->GetEntries(), 
            hMcPrimaryEtaInAcc->GetEntries(),
            hMatchedMcPrimaryEta->GetEntries(), 
            fwdPerfect,
            fwdGoodFits,
            testTrackType );

    // Matched over McPrimary
    format_eff_plot(hMatchedMcEta, hMcPrimaryEta, "Eta_Matched_McPrimary", "Efficiency vs #eta; #eta; Efficiency", kRed);
    format_eff_plot(hMatchedMcPt, hMcPrimaryPt, "Pt_Matched_McPrimary", "Efficiency vs p_{T}; p_{T} (GeV/c); Efficiency", kRed);
    format_eff_plot(hMatchedMcPhi, hMcPrimaryPhi, "Phi_Matched_McPrimary", "Efficiency vs #phi; #phi (rad); Efficiency", kRed);

    // Matched McPrimary over McPrimary 
    format_eff_plot(hMatchedMcPrimaryEta, hMcPrimaryEtaInAcc, "Eta_MatchedPrimary_McPrimary", "Efficiency vs #eta; #eta; Efficiency", kRed);
    format_eff_plot(hMatchedMcPrimaryPt, hMcPrimaryPtInAcc, "Pt_MatchedPrimary_McPrimary", "Efficiency vs p_{T}; p_{T} (GeV/c); Efficiency", kRed);
    format_eff_plot(hMatchedMcPrimaryPhi, hMcPrimaryPhiInAcc, "Phi_MatchedPrimary_McPrimary", "Efficiency vs #phi; #phi (rad); Efficiency", kRed);

    format_eff_plot(hMatchedMcPrimaryMult, hMcPrimaryMultInAcc, "Mult_MatchedPrimary_McPrimary", "Efficiency vs Multiplicity; Multiplicity; Efficiency", kRed);

    // Make the QATruth histograms
    hQATruth->SetTitle("QA Truth; QA Truth; Count");
    hQATruth->SetLineColor(kRed);
    hQATruth->SetFillColorAlpha(kRed, 0.3);
    hQATruth->SetLineWidth(2);
    hQATruthGoodChi2->SetTitle("QA Truth Good Chi2; QA Truth; Count");
    hQATruthGoodChi2->SetLineColor(kRed);
    hQATruthGoodChi2->SetFillColorAlpha(kRed, 0.3);
    hQATruthGoodChi2->SetLineWidth(2);

    // Create cumulative histograms
    hQATruth->Scale(1.0 / hQATruth->Integral()); // Normalize to 1
    hQATruthGoodChi2->Scale(1.0 / hQATruthGoodChi2->Integral()); 
    TH1 * hQaTruthC = hQATruth->GetCumulative(false);
    TH1 * hQaTruthGoodChi2C = hQATruthGoodChi2->GetCumulative(false);
    
    // Write histograms to the output file
    fOutput->Write();
}