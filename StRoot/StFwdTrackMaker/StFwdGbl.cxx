//-*-mode: C++; c-basic-offset: 2; -*-
/* Copyright 2013
 *  Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *  This is an interface to General Broken Lines
 *
 *  Version: 5 (Tadeas)
 *  - several bug-fixes:
 *    - Scatterers at bad points
 *    - Jacobians at a point before they should be (code reorganized)
 *    - Change of sign of residuals
 *  Version: 4 (Tadeas)
 *  Fixed calculation of equvivalent scatterers (solution by C. Kleinwort)
 *  Now a scatterer is inserted at each measurement (except last) and between each two measurements.
 *  TrueHits/Clusters. Ghost (1D) hits ignored. With or without magnetic field.
 *  Version: 3 (Tadeas)
 *  This version now supports both TrueHits and Clusters for VXD.
 *  It can be used for arbitrary material distribution between
 *  measurements. Moments of scattering distribution are computed
 *  and translated into two equivalent thin GBL scatterers placed
 *  at computed positions between measurement points.
 *  Version: 2 ... never published (Tadeas)
 *  Scatterer at each boundary (tooo many scatterers). TrueHits/Clusters. Without global der.&MP2 output.
 *  Version: 1 (Sergey & Tadeas)
 *  Scatterers at measurement planes. TrueHits
 *  Version 0: (Sergey)
 *  Without scatterers. Genfit 1.
 *
 *  This file is part of GENFIT.
 *
 *  GENFIT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GENFIT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */

/* July 8th, 2022
 * Gavin Wilks
 * Modify GFGbl.h/cc class to work with Forward Tracking System at STAR.
 *
 * This is used for the alignment of the forward system, which includes
 * the Forward Silicon Tracker (FST) and Forward Thickgem Tracker (FTT).
 *
 * This class is used to create the mille.dat file, which is required
 * for alignment with Millepede II.
 */

#include "StRoot/StFwdTrackMaker/StFwdGbl.h"
#include "GenFit/GblTrajectory.h"
#include "GenFit/GblPoint.h"
#include "GenFit/MilleBinary.h"
//#include "MyDebugTools.h"

#include "GenFit/AbsMeasurement.h"
#include "GenFit/PlanarMeasurement.h"
#include "GenFit/SpacepointMeasurement.h"
#include "GenFit/MeasurementOnPlane.h"
#include "GenFit/KalmanFitterInfo.h"
#include "GenFit/AbsTrackRep.h"
//#include "GenFit/MaterialProperties.h"
#include "GenFit/Track.h"
#include "GenFit/FieldManager.h"
#include "GenFit/HMatrixU.h"
#include "GenFit/HMatrixV.h"
#include "GenFit/GblFitter.h"

#include <list>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <string>
#include <Math/SMatrix.h>
#include <TMatrixD.h>
//#include <TVectorDfwd.h>
#include <TMatrixT.h>

#include <TVector3.h>

//#define DEBUG
//#define OUTPUT


#ifdef DEBUG
//ofstream debug("gbl.debug");
#endif

//#ifdef OUTPUT





//bool writeHistoDataForLabel(double label, TVectorD res, TVectorD measErr, TVectorD resErr, TVectorD downWeights, TVectorD localPar, TMatrixDSym localCov)
//{
//  if (label < 1.) return false;
//  
//  unsigned int id = floor(label);
//  // skip segment (5 bits)
//  id = id >> 5;
//  unsigned int sensor = id & 7;
//  id = id >> 3;
//  unsigned int ladder = id & 31;
//  id = id >> 5;
//  unsigned int layer = id & 7;
//  if (layer == 7 && ladder == 2) {
//    label = sensor;
//  } else if (layer == 7 && ladder == 3) {
//    label = sensor + 9 - 3;
//  } else {
//    label = layer + 3;
//  }
//  
//  if (label > 12.) return false;
//  
//  int i = int(label);
//  
//  #ifdef OUTPUT
//  resHistosU[i - 1]->Fill(res[0]);
//  resHistosV[i - 1]->Fill(res[1]);
//  mhistosU[i - 1]->Fill(res[0] / measErr[0]);
//  mhistosV[i - 1]->Fill(res[1] / measErr[1]);
//  ghistosU[i - 1]->Fill(res[0] / resErr[0]);
//  ghistosV[i - 1]->Fill(res[1] / resErr[1]);
//  downWeightsHistosU[i - 1]->Fill(downWeights[0]);
//  downWeightsHistosV[i - 1]->Fill(downWeights[1]);
//  localPar1[i - 1]->Fill(localPar(0));
//  localPar2[i - 1]->Fill(localPar(1));
//  localPar3[i - 1]->Fill(localPar(2));
//  localPar4[i - 1]->Fill(localPar(3));
//  localPar5[i - 1]->Fill(localPar(4));
//  #endif
//  
//  
//  return true;
//}
//#endif



using namespace gbl;
using namespace std;
using namespace genfit;

StFwdGbl::StFwdGbl() :
genfit::AbsFitter(), m_milleFileName("millefile.dat"), m_gblInternalIterations(""), m_pValueCut(0.), m_minNdf(1),
m_chi2Cut(0.),
m_enableScatterers(true),
m_enableIntermediateScatterer(false)
{
  
}

void StFwdGbl::beginRun()
{
  milleFile = new gbl::MilleBinary(m_milleFileName);
  
  //#ifdef OUTPUT
  //diag = new TFile(rootFileName.c_str(), "RECREATE");
  //char name[20];
  //
  //for (int i = 0; i < 12; i++) {
  //  sprintf(name, "res_u_%i", i + 1);
  //  resHistosU[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
  //  sprintf(name, "res_v_%i", i + 1);
  //  resHistosV[i] = new TH1F(name, "Residual (V)", 1000, -0.1, 0.1);
  //  sprintf(name, "meas_pull_u_%i", i + 1);
  //  mhistosU[i] = new TH1F(name, "Res/Meas.Err. (U)", 1000, -20., 20.);
  //  sprintf(name, "meas_pull_v_%i", i + 1);
  //  mhistosV[i] = new TH1F(name, "Res/Meas.Err. (V)", 1000, -20., 20.);
  //  sprintf(name, "pull_u_%i", i + 1);
  //  ghistosU[i] = new TH1F(name, "Res/Res.Err. (U)", 1000, -20., 20.);
  //  sprintf(name, "pull_v_%i", i + 1);
  //  ghistosV[i] = new TH1F(name, "Res/Res.Err. (V)", 1000, -20., 20.);
  //  sprintf(name, "downWeights_u_%i", i + 1);
  //  downWeightsHistosU[i] = new TH1F(name, "Down-weights (U)", 1000, 0., 1.);
  //  sprintf(name, "downWeights_v_%i", i + 1);
  //  downWeightsHistosV[i] = new TH1F(name, "Down-weights (V)", 1000, 0., 1.);
  //  sprintf(name, "localPar1_%i", i + 1);
  //  
  //  localPar1[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
  //  sprintf(name, "localPar2_%i", i + 1);
  //  localPar2[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
  //  sprintf(name, "localPar3_%i", i + 1);
  //  localPar3[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
  //  sprintf(name, "localPar4_%i", i + 1);
  //  localPar4[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
  //  sprintf(name, "localPar5_%i", i + 1);
  //  localPar5[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
  //}
  //fitResHisto = new TH1I("fit_result", "GBL Fit Result", 21, -1, 20);
  //ndfHisto = new TH1I("ndf", "GBL Track NDF", 41, -1, 40);
  /*for(int i = 0; i < 7; i++){ 
    residualsUkf[i]  = new TH1F(Form("residualsUkf_%d", i) , "u (cm)", 100, -3.,3.);
    residualsVkf[i]  = new TH1F(Form("residualsVkf_%d", i) , "v (cm)", 100, -.5,.5);
    residualsUgbl[i] = new TH1F(Form("residualsUgbl_%d",i), "u (cm)", 100, -.2,.2);
    residualsVgbl[i] = new TH1F(Form("residualsVgbl_%d",i), "v (cm)", 100, -2.,2.);
       merrUgbl[i] = new TH1F(Form("merrUgbl_%d",   i),"",100,-.1,.1);
       merrVgbl[i] = new TH1F(Form("merrVgbl_%d",   i),"",100,-.1,.1);
       rerrUgbl[i] = new TH1F(Form("rerrUgbl_%d",   i),"",100,-.1,.1);
       rerrVgbl[i] = new TH1F(Form("rerrVgbl_%d",   i),"",100,-.1,.1);
    dweightUgbl[i] = new TH1F(Form("dweightUgbl_%d",i),"",100,-2.,2.);
    dweightVgbl[i] = new TH1F(Form("dweightVgbl_%d",i),"",100,-2.,2.); 
  }*/
  chi2OndfHistoGBL = new TH1F("chi2_ndfGBL", "Track Chi2/NDF", 100, 0., 10.);
  pValueHistoGBL = new TH1F("p_valueGBL", "Track P-value", 100, 0., 1.);
  chi2OndfHisto = new TH1F("chi2_ndf", "Track Chi2/NDF", 100, 0., 10.);
  pValueHisto = new TH1F("p_value", "Track P-value", 100, 0., 1.);
  
  //stats = new TH1I("stats", "0: NDF>0 | 1: fTel&VXD | 2: all | 3: VXD | 4: SVD | 5: all - PXD | 6: fTel&SVD | 7: bTel", 10, 0, 10);
  
  //#endif
}

void StFwdGbl::endRun()
{
  //#ifdef OUTPUT
  //diag->cd();
  //diag->Write();
  //diag->Close();
  //#endif
  // This is needed to close the file before alignment starts
  cout << "Before Millefile" << endl;
  if (milleFile)
    delete milleFile;
  cout << "Deleted Millefile" << endl;
}

/**
 * @brief Evaluates moments of radiation length distribution from list of
 * material steps and computes parameters describing a corresponding thick scatterer.
 *
 * Based on input from Claus Kleinwort (DESY),
 * adapted for continuous material distribution represented by
 * a sum of step functions. Returned thick scatterer can be represented by two GBL scattering points
 * at (s - ds) and (s + ds) with variance of theta equal to theta/sqrt(2) for both points.
 * Calculates variance of theta from total sum of radiation lengths
 * instead of summimg squares of individual deflection angle variances.
 *
 * @param length returned: Length of the track
 * @param theta returned: Variation of distribution of deflection angle
 * @param s returned: First moment of material scattering distribution
 * @param ds returned: Second moment (variance) of material scattering distribution
 * @param p Particle momentum magnitude (GeV/c)
 * @param mass Mass of particle (GeV/c/c)
 * @param steps Vector of material steps from (RKTrackRep) extrapolation
 * @return void
 */
void getScattererFromMatList(double& length,
                             double& theta, 
                             double& s, 
                             double& ds, 
                             const double p, 
                             const double mass, 
                             const double charge, 
                             const std::vector<genfit::MatStep>& steps)
{
  theta = 0.; s = 0.; ds = 0.;
  
  //cout << "if (steps.empty()) return;" << endl;
  if (steps.empty()) return;
  //cout << "passed test" << endl;
  
  // sum of step lengths
  double len = 0.;
  // normalization
  double sumxx = 0.;
  // first moment (non-normalized)
  double sumx2x2 = 0.;
  // (part of) second moment / variance (non-normalized)
  double sumx3x3 = 0.;
  
  // cppcheck-suppress unreadVariable
  double xmin = 0.;
  double xmax = 0.;
  
  for (unsigned int i = 0; i < steps.size(); i++) {
    const genfit::MatStep step = steps.at(i);
    // inverse of material radiation length ... (in 1/cm) ... "density of scattering"
    double rho = 1. / step.materialProperties_.getRadLen();
    //double rho = 1. / 1.; //step.materialProperties_.getRadLen();
    len += fabs(step.stepSize_);
    xmin = xmax;
    xmax = xmin + fabs(step.stepSize_);
    // Compute integrals
    
    // integral of rho(x)
    sumxx   += rho * (xmax - xmin);
    // integral of x*rho(x)
    sumx2x2 += rho * (xmax * xmax - xmin * xmin) / 2.;
    // integral of x*x*rho(x)
    sumx3x3 += rho * (xmax * xmax * xmax - xmin * xmin * xmin) / 3.;
  }
  // This ensures PDG formula still gives positive results (but sumxx should be > 1e-4 for it to hold)
  //cout << "if (sumxx < 1.0e-10) return;" << endl;
  if (sumxx < 1.0e-10) return;
  //cout << "passed" << endl;
  // Calculate theta from total sum of radiation length
  // instead of summimg squares of individual deflection angle variances
  // PDG formula:
  double beta = p / sqrt(p * p + mass * mass);
  theta = (0.0136 / p / beta) * fabs(charge) * sqrt(sumxx) * (1. + 0.038 * log(sumxx));
  //theta = (0.015 / p / beta) * fabs(charge) * sqrt(sumxx);
  
  // track length
  length = len;
  // Normalization factor
  double N = 1. / sumxx;
  // First moment
  s  = N * sumx2x2;
  // Square of second moment (variance)
  // integral of (x - s)*(x - s)*rho(x)
  double ds_2 = N * (sumx3x3 - 2. * sumx2x2 * s + sumxx * s * s);
  ds = sqrt(ds_2);
  
  //#ifdef DEBUG
  cout << "Thick scatterer parameters:" << endl;
  cout << "Variance of theta: " << theta << endl;
  cout << "Mean s           : " << s << endl;
  cout << "Variance of s    : " << ds << endl;
  
  //#endif
}


void StFwdGbl::processTrackWithRep(genfit::Track* trk, const genfit::AbsTrackRep* rep, bool resortHits)
{
  // Flag used to mark error in raw measurement combination
  // measurement won't be considered, but scattering yes
  bool skipMeasurement = false;
  // Chi2 of Reference Track
  double trkChi2 = 0.;
  double trkNdf = 0;
 
  auto cardinalRep = trk->getCardinalRep();
  auto cardinalStatus = trk->getFitStatus(cardinalRep);
  trkChi2 = cardinalStatus->getChi2();
  trkNdf = cardinalStatus->getNdf();
  // This flag enables/disables fitting of q/p parameter in GBL
  // It is switched off automatically if no B-field at (0,0,0) is detected.
  bool fitQoverP = true;
  //TODO: Use clever way to determine zero B-field
  double Bfield = genfit::FieldManager::getInstance()->getFieldVal(TVector3(0., 0., 0.)).Mag();
  //cout << "Bfield @ (0,0,0) = " << Bfield << endl;

  if (!(Bfield > 0.))
  {
    fitQoverP = false;
    //cout << "DO NOT FIT q/p!" << endl;
  }
 
  // Dimesion of repository/state
  int dim = rep->getDim();
  // current measurement point
  genfit::TrackPoint* point_meas;
  // current raw measurement
  genfit::AbsMeasurement* raw_meas;
  
  // We assume no initial kinks, this will be reused several times
  TVectorD scatResidual(2);
  scatResidual.Zero();
  
  // All measurement points in ref. track
  int npoints_meas = trk->getNumPointsWithMeasurement();
  
  #ifdef DEBUG
  int npoints_all = trk->getNumPoints();
  
  if (resortHits)
    //cout << "WARNING: Hits resorting in GBL interface not supported." << endl;
  
  cout << "-------------------------------------------------------" << endl;
  cout << "               GBL processing genfit::Track            " << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << " # Ref. Track Points  :  " << npoints_all  << endl;
  cout << " # Meas. Points       :  " << npoints_meas << endl;
  
  #endif
  // List of prepared GBL points for GBL trajectory construction
  std::vector<gbl::GblPoint> listOfPoints;
  
  std::vector<double> listOfLayers;
  // index of point with seed information (0 for none)
  unsigned int seedLabel = 0;
  // Seed covariance
  // TMatrixDSym clCov(dim);
  // Seed state
  TMatrixDSym clSeed(dim);
  
  // propagation Jacobian to next point from current measurement point
  TMatrixD jacPointToPoint(dim, dim);
  jacPointToPoint.UnitMatrix();
  //cout << "jacPointToPoint" << endl;
  //jacPointToPoint.Print();  

  int n_gbl_points = 0;
  int n_gbl_meas_points = 0;
  int Ndf = 0;
  double Chi2 = 0.;
  double lostWeight = 0.;
  
  // Momentum of track at current plane
  double trackMomMag = 0.;
  // Charge of particle at current plane :-)
  double particleCharge = 1.;
  
  bool si1 = false;
  bool si2 = false;
  bool si3 = false;
  bool siO1_1 = false;
  bool siO1_2 = false;
  bool siO2_1 = false;
  bool siO2_2 = false;

  bool stgcP1Q0 = false;

  bool layer1 = false;
  bool layer2 = false;
  bool layer3 = false;
 
  bool innerSi1 = false;
  bool innerSi2 = false;
  bool innerSi3 = false;

  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    genfit::TrackPoint* point_meas_temp = trk->getPointWithMeasurement(ipoint_meas);
    genfit::PlanarMeasurement* measPlanar = dynamic_cast<genfit::PlanarMeasurement*>(point_meas_temp->getRawMeasurement(0));
    int sid;
    if (measPlanar) sid = measPlanar->getPlaneId();
  
    if((sid-1000)/36 == 0) layer1 = true;
    if((sid-1000)/36 == 1) layer2 = true;
    if((sid-1000)/36 == 2) layer3 = true;
    if (sid == 1036) si2 = true; 
    if (sid == 1037) siO2_1 = true; 
    if (sid == 1038) siO2_2 = true;
    if (sid == 1000) si1 = true; 
    if (sid == 1001) siO1_1 = true; 
    if (sid == 1002) siO1_2 = true;
    if (sid == 4) stgcP1Q0 = true;
    if((sid-1000)/36 == 0 && (sid-1000)%3 == 0) innerSi1 = true;
    if((sid-1000)/36 == 1 && (sid-1000)%3 == 0) innerSi1 = true;
    if((sid-1000)/36 == 2 && (sid-1000)%3 == 0) innerSi1 = true;
  }

  //if(m_allFST)
  //{
  //  if(!layer1 || !layer2 || !layer3 || !stgcP1Q0 || !innerSi1 || !innerSi2 || !innerSi3 ) 
  //  {
  //    cout << "WARNING!!! THE REQUIRED SENSOR IS NOT PRESENT, SKIP TRACK!" << endl;
  //    return;
  //  }
  //}
  if(m_allFST)
  {
    if(!layer1 || !layer2 || !layer3) 
    {
      cout << "WARNING!!! THE REQUIRED SENSOR IS NOT PRESENT, SKIP TRACK!" << endl;
      return;
    }
  }
  /*else
  {
    if(!siO1_1) 
    {
      cout << "WARNING!!! THE REQUIRED SENSOR IS NOT PRESENT, SKIP TRACK!" << endl;
      return;
    }
  }*/
  //if(!siO2_2) 
  //{
  //  cout << "WARNING!!! THE REQUIRED SENSOR IS NOT PRESENT, SKIP TRACK!" << endl;
  //  return;
  //}
  bool includeVertex = false;
   
  //cout << "Before point loop" << endl;
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    //cout << "Before grabbing PointWithMeasurement" << endl;
    //if(m_includeVertex && ipoint_meas == 0) continue;

    point_meas = trk->getPointWithMeasurement(ipoint_meas);
    
    //cout << "Before checking fitter info" << endl;
    if (!point_meas->hasFitterInfo(rep)) {
      cout << " ERROR: Measurement point does not have a fitter info. Track will be skipped." << endl;
      return;
    }
    // Fitter info which contains Reference state and plane
    //cout << "Before  grabbing fitter info" << endl;
    genfit::KalmanFitterInfo* fi = dynamic_cast<genfit::KalmanFitterInfo*>(point_meas->getFitterInfo(rep));
    if (!fi) {
      cout << " ERROR: KalmanFitterInfo (with reference state) for measurement does not exist. Track will be skipped." << endl;
      return;
    }
    // Current detector plane
    //cout << "Before grabbing plane" << endl;
    genfit::SharedPlanePtr plane = fi->getPlane();
    if (!fi->hasReferenceState()) {
      cout << " ERROR: Fitter info does not contain reference state. Track will be skipped." << endl;
      return;
    }
    // Reference StateOnPlane for extrapolation
    //cout << "Before grabbing referencestatonplane" << endl;
    genfit::ReferenceStateOnPlane* reference = new genfit::ReferenceStateOnPlane(*fi->getReferenceState());//(dynamic_cast<const genfit::ReferenceStateOnPlane&>(*fi->getReferenceState()));
    // Representation state at plane
    //cout << "Before grabbing state" << endl;
    TVectorD state = reference->getState();
    // track direction at plane (in global coords)
    //cout << "Before grabbing trackdir" << endl;
    TVector3 trackDir = rep->getDir(*reference);
    // track momentum vector at plane (in global coords)
    //cout << "Before grabbing mag" << endl;
    trackMomMag = rep->getMomMag(*reference);
    // charge of particle
    particleCharge = rep->getCharge(*reference);
    // mass of particle
    double particleMass = rep->getMass(*reference);
    
    // Parameters of a thick scatterer between measurements
    double trackLen = 0.;
    double scatTheta = 0.;
    double scatSMean = 0.;
    double scatDeltaS = 0.;
    // Parameters of two equivalent thin scatterers
    double theta1 = 0.;
    double theta2 = 0.;
    double s1 = 0.;
    double s2 = 0.;
    
    TMatrixDSym noise;
    TVectorD deltaState;
    // jacobian from s2 to M2
    TMatrixD jacMeas2Scat(dim, dim);
    jacMeas2Scat.UnitMatrix();
    
    
    // Now get measurement. First have a look if we need to combine SVD clusters...
    // Load Jacobian from previous extrapolation
    //rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
    // Try to get VxdId of current plane
    int sensorId = 0;
    //cout << "before grabbing measurement" << endl;
    ///genfit::MeasurementOnPlane* ;


    genfit::PlanarMeasurement* measPlanar;
    genfit::SpacepointMeasurement* measSpace;
    if(point_meas->getRawMeasurement(0)->getDim() <= 2) measPlanar = dynamic_cast<genfit::PlanarMeasurement*>(point_meas->getRawMeasurement(0));
    if(point_meas->getRawMeasurement(0)->getDim() == 3) measSpace = dynamic_cast<genfit::SpacepointMeasurement*>(point_meas->getRawMeasurement(0));
    if (measPlanar && point_meas->getRawMeasurement(0)->getDim() <= 2) sensorId = measPlanar->getPlaneId();
    if (measSpace  && point_meas->getRawMeasurement(0)->getDim() == 3) {
      sensorId = measSpace->getDetId();
      includeVertex = true;
    } 
    //if(sensorId == 1000) si1 = true;
    //if(sensorId == 1036) si2 = true;
    //if(sensorId == 1072) si3 = true;
    //WARNING: Now we only support 2D measurements. If 2 raw measurements are stored at the track
    // point, these are checked if they correspond to "u" and "v" measurement (for SVD clusters) and these
    // measurements are combined. SLANTED SENSORS NOT YET SUPPORTED!!!
    //cout << "before determining dimension of the measurement" << endl;
    if (point_meas->getRawMeasurement(0)->getDim() != 2
      && trk->getPointWithMeasurement(ipoint_meas)->getNumRawMeasurements() == 2
      && point_meas->getRawMeasurement(0)->getDim() == 1
      && point_meas->getRawMeasurement(1)->getDim() == 1) {
      genfit::AbsMeasurement* raw_measU = 0;
      genfit::AbsMeasurement* raw_measV = 0;
    
      // cout << " Two 1D Measurements encountered. " << endl;
    
      int sensorId2 = -1;
      genfit::PlanarMeasurement* measPlanar2 = dynamic_cast<genfit::PlanarMeasurement*>(point_meas->getRawMeasurement(0));
      if (measPlanar2) sensorId2 = measPlanar->getPlaneId();
    
      // We only try to combine if at same sensor id (should be always, but who knows)
      // otherwise ignore this point
      if (sensorId != sensorId2) {
	skipMeasurement = true;
	cout << " ERROR: Incompatible sensorIDs at measurement point " << ipoint_meas << ". Measurement will be skipped." << endl;
      }
    
      // We have to combine two SVD 1D Clusters at the same plane into one 2D recohit
      genfit::AbsMeasurement* raw_meas1 = point_meas->getRawMeasurement(0);
      genfit::AbsMeasurement* raw_meas2 = point_meas->getRawMeasurement(1);
      // Decide which cluster is u and which v based on H-matrix
      if (raw_meas1->constructHMatrix(rep)->isEqual(genfit::HMatrixU())
	  && raw_meas2->constructHMatrix(rep)->isEqual(genfit::HMatrixV())) {
	// right order U, V
	raw_measU = raw_meas1;
	raw_measV = raw_meas2;
      } else if (raw_meas2->constructHMatrix(rep)->isEqual(genfit::HMatrixU())
		 && raw_meas1->constructHMatrix(rep)->isEqual(genfit::HMatrixV())) {
        // inversed order V, U
        raw_measU = raw_meas2;
	raw_measV = raw_meas1;
      } else {
	// Incompatible measurements ... skip track ... I saw this once and just before a segfault ...
	cout << " ERROR: Incompatible 1D measurements at meas. point " << ipoint_meas << ". Track will be skipped." << endl;
	return;
      }
      // Combine raw measurements
      TVectorD _raw_coor(2);
      _raw_coor(0) = raw_measU->getRawHitCoords()(0);
      _raw_coor(1) = raw_measV->getRawHitCoords()(0);
      // Combine covariance matrix
      TMatrixDSym _raw_cov(2);
      _raw_cov.Zero();
      _raw_cov(0, 0) = raw_measU->getRawHitCov()(0, 0);
      _raw_cov(1, 1) = raw_measV->getRawHitCov()(0, 0);
      // Create new combined measurement
      raw_meas = new genfit::PlanarMeasurement(_raw_coor, _raw_cov, raw_measU->getDetId(), raw_measU->getHitId(), point_meas);
    } else if(point_meas->getRawMeasurement(0)->getDim() == 2) {
      // Default behavior
      raw_meas = point_meas->getRawMeasurement(0);
      TVectorD raw_coor = raw_meas->getRawHitCoords();
      float x = raw_coor(0);
      float y = raw_coor(1);
      if(sensorId >= 1000 && sensorId < 9000) {
        if((sensorId-1000)%3 == 0) x+=10.75;
        if((sensorId-1000)%3 != 0) x+=22.25;
        float r = TMath::Sqrt(x*x+y*y);
        if(sensorId >= 1000 && sensorId < 1036 && r >= 22.250) skipMeasurement = true;
        if(sensorId >= 1036 && sensorId < 1072 && r <=  7.875) skipMeasurement = true;
        if(sensorId >= 1036 && sensorId < 1072 && r >= 25.125) skipMeasurement = true;
        if(sensorId >= 1072 && sensorId < 1108 && r <=  7.875) skipMeasurement = true;
      }
    } else {
      //measSpace->constructPlane(state);
      raw_meas = point_meas->getRawMeasurement(0);
      //cout << "Vertex information" << endl;
      //measSpace->constructMeasurementsOnPlane(*dynamic_cast<genfit::StateOnPlane*>(reference))[0]->getCov().Print();
      //measSpace->constructMeasurementsOnPlane(*dynamic_cast<genfit::StateOnPlane*>(reference))[0]->getState().Print();
    }
    //TODO: We only support 2D measurements in GBL (ot two 1D combined above)
    /*if (raw_meas->getRawHitCoords().GetNoElements() != 2 && raw_meas_vertex->getRawHitCoords().GetNoElements() != 3) {
      skipMeasurement = true;
      #ifdef DEBUG
      cout << " WARNING: Measurement " << (ipoint_meas + 1) << " is not 2D or 3D. Measurement Will be skipped. " << endl;
      #endifi
    }*/
      
    // Now we have all necessary information, so lets insert current measurement point
    // if we don't want to skip it (e.g. ghost SVD hit ... just 1D information)
    //cout << "decided the dimension of the measurement" << endl;
    if (!skipMeasurement) {
      //cout << "inside if statment for !skipMeasurment" << endl;
      
      
      // 2D hit coordinates
      TVectorD raw_coor = raw_meas->getRawHitCoords();
      TVectorD raw_coorNew(2); 
      if(point_meas->getRawMeasurement(0)->getDim() == 2)
      {
        raw_coorNew(0) = raw_coor(0); 
        raw_coorNew(1) = raw_coor(1); 
      }
      if(point_meas->getRawMeasurement(0)->getDim() == 3)
      {
        //cout << "3D measurement" << endl;
        raw_coorNew = raw_meas->constructMeasurementsOnPlane(*dynamic_cast<genfit::StateOnPlane*>(reference))[0]->getState();
      }
      //cout << "Raw Coordinates " <<endl;
      //raw_coorNew.Print();
      // Covariance matrix of measurement
      TMatrixDSym raw_cov(2); 
      if(point_meas->getRawMeasurement(0)->getDim() == 2)
      {
        raw_cov = raw_meas->getRawHitCov();
      }
      if(point_meas->getRawMeasurement(0)->getDim() == 3)
      {
        //cout << "3D measurement" << endl;
        raw_cov = raw_meas->constructMeasurementsOnPlane(*dynamic_cast<genfit::StateOnPlane*>(reference))[0]->getCov();
        auto planeSpacePoint = measSpace->constructPlane(*dynamic_cast<genfit::StateOnPlane*>(reference));
        //cout << "Plane information" << endl;
        //planeSpacePoint->Print();
        TVectorD measured = raw_meas->constructMeasurementsOnPlane(*dynamic_cast<genfit::StateOnPlane*>(reference))[0]->getState();
        TVector2 rcpoint(measured(0),measured(1));
        //cout << "Vertex in lab" << endl;
        //planeSpacePoint->toLab(rcpoint).Print();
        std::unique_ptr<const genfit::AbsHMatrix> HitHMatrix3d(raw_meas->constructHMatrix(rep));
        TVectorD planestate = HitHMatrix3d->Hv(state);
        //cout << "Track Vertex in lab" << endl;
        TVector2 trackpoint(planestate(0),planestate(1));
        //planeSpacePoint->toLab(trackpoint).Print();
        //state.Print();
      }
      //cout << "Raw Covariance" << endl;
      //raw_cov.Print();
      // Projection matrix from repository state to measurement coords
      std::unique_ptr<const genfit::AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(rep));
      // Residual between measured position and reference track position
      TVectorD residual = -1.*(raw_coorNew - HitHMatrix->Hv(state));
      //cout << "state" << endl;
      //state.Print();

      //cout << "Track Projection" << endl;
      //HitHMatrix->Hv(state).Print();
      //cout << "right before filling the residuals" << endl;
      //cout << "residuals = " << residual(0) << ", " << residual(1) << endl;
      //residualsUkf[ipoint_meas]->Fill(residual(0));
      //residualsVkf[ipoint_meas]->Fill(residual(1));

      //cout << "grabbed information and filled residuals" << endl;
      //TVectorD residual = (raw_coor - HitHMatrix->Hv(state));

      //trkChi2 += residual(0) * residual(0) / raw_cov(0, 0) + residual(1) * residual(1) / raw_cov(1, 1);
        
      // Measurement point
      //cout << "jacPointToPoint" << endl;
      //jacPointToPoint.Print();
      TMatrixDSym raw_covNew(2);
      raw_covNew(0,0) = raw_cov(0,0);
      raw_covNew(0,1) = raw_cov(0,1);
      raw_covNew(1,0) = raw_cov(1,0);
      raw_covNew(1,1) = raw_cov(1,1);
      //cout << "2x2 Covariant" << endl;
      //raw_covNew.Print();
      //cout << "Residuals = " << endl;
      //residual.Print();
      gbl::GblPoint measPoint(jacPointToPoint);
      // Projection from local (state) coordinates .to measurement coordinates (inverted)
      // 2x2 matrix ... last block of H matrix (2 rows x 5 columns)
      //TMatrixD proL2m = HitHMatrix->getMatrix().GetSub(0, 1, 3, 4);
      TMatrixD proL2m(2, 2);
      proL2m.Zero();
      proL2m(0, 0) = HitHMatrix->getMatrix()(0, 3);
      proL2m(0, 1) = HitHMatrix->getMatrix()(0, 4);
      proL2m(1, 1) = HitHMatrix->getMatrix()(1, 4);
      proL2m(1, 0) = HitHMatrix->getMatrix()(1, 3);
      //cout << "HitHMatrix" << endl;
      //HitHMatrix->getMatrix().Print();
      proL2m.Invert();
      //raw_cov *= 100.;
      //cout << "proL2m" << endl;
      //proL2m.Print();
      //cout << "raw_cov" << endl;
      //raw_cov.Print();
      measPoint.addMeasurement(proL2m, residual, raw_covNew.Invert());
      //measPoint.addMeasurement(HitHMatrix->getMatrix(), residual, raw_cov.Invert());
      //cout << "inverted raw_cov" << endl;
      //raw_cov.Print();
      //cout << "residual" << endl;
      //residual.Print();

      //Add global derivatives to the point
        
      // sensor label = sensorID * 10, then last digit is label for global derivative for the sensor
      int label;
      int wedgelabel = -1;
      int halflabel = -1;
      double lu, lv, s;

      if(sensorId < 1000) 
      {
        label = sensorId * 10;
        //lu = lu_pent;
        //lv = lv_pent;
        //s  = (lu + lv)/2.0;
      }
      if(sensorId >= 1000)
      {
        if((sensorId - 1000)%3 == 0)
        {
          label = 1000 + (sensorId-1000) * 10;
          //lu = lu_iSi;
          //lv = lv_iSi;
          //s  = (lu + lv)/2.0;
        }
        else if((sensorId - 1000)%3 != 0) 
        {
          label = 1000 + (sensorId-1000) * 10;
          //lu = lu_oSi;
          //lv = lv_oSi;
          //s  = (lu + lv)/2.0;
        }
        //else if((sensorId - 1000)%3 == 2) label = 1000 + (sensorId-1001) * 10; // joins the outer sensors together with one set of alignment parameters
        wedgelabel = 3000 + ((sensorId-1000)/3)*10;
        halflabel  = 4000 + (((sensorId-1000)/18)%2)*10;
      }
      if(sensorId > 9000)
      {
        label = 9000;
      }

      //cout << "SensorID = " << sensorId << "     label = " << label << "      wedgelabel = " << wedgelabel << "      halflabel = " << halflabel << endl;
      // values for global derivatives
      TMatrixD derGlobalFTT(2, 6);
      TMatrixD derGlobalFST(2, 18);
        
      // labels for global derivatives
      std::vector<int> labGlobal;
        
      // track direction in global coords
      TVector3 tDir = trackDir;
      // sensor u direction in global coords
      TVector3 uDir = plane->getU();
      // sensor v direction in global coords
      TVector3 vDir = plane->getV();
      // sensor normal direction in global coords
      TVector3 nDir = plane->getNormal();
      //file << sensorId << endl;
      //outputVector(uDir, "U");
      //outputVector(vDir, "V");
      //outputVector(nDir, "Normal");
      // track direction in local sensor system
      //tDir.Print();
      //uDir.Print();
      //vDir.Print();
      //nDir.Print();

      //TVector3 tLoc = TVector3(uDir.Dot(tDir), vDir.Dot(tDir), nDir.Dot(tDir));
      TVector3 tLoc = TVector3(tDir.Dot(uDir), tDir.Dot(vDir), tDir.Dot(nDir));
        
      // track u-slope in local sensor system
      double uSlope = tLoc[0] / tLoc[2];
      // track v-slope in local sensor system
      double vSlope = tLoc[1] / tLoc[2];
        
      // Measured track u-position in local sensor system
      double uPos = raw_coor[0];
      // Measured track v-position in local sensor system
      double vPos = raw_coor[1];

      //cout << "uSlope = " << uSlope << endl;
      //cout << "vSlope = " << vSlope << endl;
        
     
      //This follows from CMS derivatives, excluding sensor surface parameters
      //Global derivatives for alignment in sensor local coordinates
      derGlobalFST(0, 0) = 1.0;                    // du/du
      derGlobalFST(0, 1) = 0.0;                    // du/dv
      derGlobalFST(0, 2) = -uSlope;               // du/dw
      derGlobalFST(0, 3) = -vPos * uSlope;          // du/d(alpha)
      derGlobalFST(0, 4) = uPos * uSlope;         // du/d(beta)
      derGlobalFST(0, 5) = -vPos;                   // du/d(gamma)
         
      derGlobalFST(1, 0) = 0.0;                    // dv/du
      derGlobalFST(1, 1) = 1.0;                    // dv/dv
      derGlobalFST(1, 2) = -vSlope;               // dv/dw
      derGlobalFST(1, 3) = -vPos * vSlope;          // dv/d(alpha)
      derGlobalFST(1, 4) = uPos * vSlope;         // dv/d(beta)
      derGlobalFST(1, 5) = uPos;                  // dv/d(gamma)

      derGlobalFTT(0, 0) = 1.0;                    // du/du
      derGlobalFTT(0, 1) = 0.0;                    // du/dv
      derGlobalFTT(0, 2) = -uSlope;               // du/dw
      derGlobalFTT(0, 3) = -vPos * uSlope;          // du/d(alpha)
      derGlobalFTT(0, 4) = uPos * uSlope;         // du/d(beta)
      derGlobalFTT(0, 5) = -vPos;                   // du/d(gamma)
         
      derGlobalFTT(1, 0) = 0.0;                    // dv/du
      derGlobalFTT(1, 1) = 1.0;                    // dv/dv
      derGlobalFTT(1, 2) = -vSlope;               // dv/dw
      derGlobalFTT(1, 3) = -vPos * vSlope;          // dv/d(alpha)
      derGlobalFTT(1, 4) = uPos * vSlope;         // dv/d(beta)
      derGlobalFTT(1, 5) = uPos;                  // dv/d(gamma)

      labGlobal.push_back(label + 1); // u
      labGlobal.push_back(label + 2); // v
      labGlobal.push_back(label + 3); // w
      labGlobal.push_back(label + 4); // alpha
      labGlobal.push_back(label + 5); // beta
      labGlobal.push_back(label + 6); // gamma

      //cout << "derGlobal" << label; 
      //derGlobal.Print();

      /*
      // Global derivatives for movement of whole detector system in global coordinates
      //TODO: Usage of this requires Hierarchy Constraints to be provided to MP2
  
      // sensor centre position in global system
      TVector3 detPos = plane->getO();
      //cout << "detPos" << endl;
      //detPos.Print();

      // global prediction from raw measurement
      TVector3 pred = detPos + uPos * uDir + vPos * vDir;
      //cout << "pred" << endl;
      //pred.Print();

      double xPred = pred[0];
      double yPred = pred[1];
      double zPred = pred[2];

      // scalar product of sensor normal and track direction
      double tn = tDir.Dot(nDir);
      //cout << "tn" << endl;
      //cout << tn << endl;

      // derivatives of local residuals versus measurements
      TMatrixD drdm(3, 3);
      drdm.UnitMatrix();
      for (int row = 0; row < 3; row++)
	for (int col = 0; col < 3; col++)
	  drdm(row, col) -= tDir[row] * nDir[col] / tn;

      //cout << "drdm" << endl;
      //drdm.Print();

      // derivatives of measurements versus global alignment parameters
      TMatrixD dmdg(3, 6);
      dmdg.Zero();
      dmdg(0, 0) = 1.; dmdg(0, 4) = -zPred; dmdg(0, 5) = yPred;
      dmdg(1, 1) = 1.; dmdg(1, 3) = zPred;  dmdg(1, 5) = -xPred;
      dmdg(2, 2) = 1.; dmdg(2, 3) = -yPred; dmdg(2, 4) = xPred;

      //cout << "dmdg" << endl;
      //dmdg.Print();

      // derivatives of local residuals versus global alignment parameters
      TMatrixD drldrg(3, 3);
      drldrg.Zero();
      drldrg(0, 0) = uDir[0]; drldrg(0, 1) = uDir[1]; drldrg(0, 2) = uDir[2];
      drldrg(1, 0) = vDir[0]; drldrg(1, 1) = vDir[1]; drldrg(1, 2) = vDir[2];

      //cout << "drldrg" << endl;
      //drldrg.Print();

      //cout << "drdm * dmdg" << endl;
      //(drdm * dmdg).Print();

      // derivatives of local residuals versus rigid body parameters
      TMatrixD drldg(3, 6);
      drldg = drldrg * (drdm * dmdg);

      //cout << "drldg" << endl;
      //drldg.Print();

      // offset to determine labels for sensor sets or individual layers
      // 0: PXD, TODO 1: SVD, or individual layers
      // offset 0 is for alignment of whole setup
      int offset = 0;
      //if (sensorId > 16704) offset = 20; // svd, but needs to introduce new 6 constraints: sum rot&shifts of pxd&svd = 0
      */
      if(sensorId >= 1000 && sensorId < 9000)
      {
        TMatrixD drdps(6, 2); // d(residuals)/d(sensor alignment parameters)

        drdps(0, 0) = 1.0;                    // du/du
        drdps(1, 0) = 0.0;                    // du/dv
        drdps(2, 0) = -uSlope;                // du/dw
        drdps(3, 0) = -vPos * uSlope;         // du/d(alpha)
        drdps(4, 0) = uPos * uSlope;          // du/d(beta)
        drdps(5, 0) = -vPos;                  // du/d(gamma)
       
        drdps(0, 1) = 0.0;                    // dv/du
        drdps(1, 1) = 1.0;                    // dv/dv
        drdps(2, 1) = -vSlope;                // dv/dw
        drdps(3, 1) = -vPos * vSlope;         // dv/d(alpha)
        drdps(4, 1) = uPos * vSlope;          // dv/d(beta)
        drdps(5, 1) = uPos;                   // du/d(gamma)

        // This depends on which sensor we are on
        TMatrixD dpsdpw(6, 6); // d(sensor alignment parameters)/d(wedge alignment parameters)

        int is = sensorId - 1000;

        int half = (is / 18) % 2; // 0 (left +x half), 1 (right -x half) 

        // FST disk (integer division rounds down)
        int disk = is / 36; // 0-2

        // FST wedge
        int wedge = is / 3; // 0-35
        int wedgeId = wedge%12; // 0-12

        // FST sensor
        int sid = is % 3; // 0 (inner), 1 (outer), 2 (outer)

        bool isPosY = false;            
        if(disk == 0 || disk == 2)
        {
          if(wedge%2 == 0)
          {
            if(sid == 1)      isPosY = true; //phiAxisShift = +8.0 * TMath::Pi() / 180.0;
            else if(sid == 2) isPosY = false; //phiAxisShift = -8.0 * TMath::Pi() / 180.0;
          }
          else if(wedge%2 == 1)
          {
            if(sid == 1)      isPosY = false; //phiAxisShift = -8.0 * TMath::Pi() / 180.0;
            else if(sid == 2) isPosY = true; //phiAxisShift = +8.0 * TMath::Pi() / 180.0;
          }
        }
        else if(disk == 1)
        {
          if(wedge%2 == 0)
          {
            if(sid == 1)      isPosY = false; //phiAxisShift = -8.0 * TMath::Pi() / 180.0;
            else if(sid == 2) isPosY = true; //phiAxisShift = +8.0 * TMath::Pi() / 180.0;
          }
          else if(wedge%2 == 1)
          {
            if(sid == 1)      isPosY = true; //phiAxisShift = +8.0 * TMath::Pi() / 180.0;
            else if(sid == 2) isPosY = false; //phiAxisShift = -8.0 * TMath::Pi() / 180.0;
          }
        }

        if(sid == 0)
        {
          dpsdpw(0, 0) = 1.0; dpsdpw(0, 1) = 0.0;   dpsdpw(0, 2) = 0.0;   dpsdpw(0, 3) = 0.0; dpsdpw(0, 4) = 0.0; dpsdpw(0, 5) = 0.0; //dps/duw
          dpsdpw(1, 0) = 0.0; dpsdpw(1, 1) = 1.0;   dpsdpw(1, 2) = 0.0;   dpsdpw(1, 3) = 0.0; dpsdpw(1, 4) = 0.0; dpsdpw(1, 5) = 0.0; //dps/dvw
          dpsdpw(2, 0) = 0.0; dpsdpw(2, 1) = 0.0;   dpsdpw(2, 2) = 1.0;   dpsdpw(2, 3) = 0.0; dpsdpw(2, 4) = 0.0; dpsdpw(2, 5) = 0.0; //dps/dww
          dpsdpw(3, 0) = 0.0; dpsdpw(3, 1) = 0.0;   dpsdpw(3, 2) = 0.0;   dpsdpw(3, 3) = 1.0; dpsdpw(3, 4) = 0.0; dpsdpw(3, 5) = 0.0; //dps/dalphaw
          dpsdpw(4, 0) = 0.0; dpsdpw(4, 1) = 0.0;   dpsdpw(4, 2) = 5.75;  dpsdpw(4, 3) = 0.0; dpsdpw(4, 4) = 1.0; dpsdpw(4, 5) = 0.0; //dps/dbetaw
          dpsdpw(5, 0) = 0.0; dpsdpw(5, 1) = -5.75; dpsdpw(5, 2) = 0.0;   dpsdpw(5, 3) = 0.0; dpsdpw(5, 4) = 0.0; dpsdpw(5, 5) = 1.0; //dps/dgammaw
        }
        if(sid != 0)
        {
          double cos8 = TMath::Cos(8.*TMath::Pi()/180.);
          double sin8 = TMath::Sin(8.*TMath::Pi()/180.);
          double r = 6.34099129;
          double x = 5.53346453;
          double y = 3.0966015;
          double costm8 = TMath::Cos(21.2319658*TMath::Pi()/180.);
          double sintm8 = TMath::Sin(21.2319658*TMath::Pi()/180.);
          if(isPosY)
          {
            dpsdpw(0, 0) = cos8;     dpsdpw(0, 1) = -sin8;     dpsdpw(0, 2) = 0.0; dpsdpw(0, 3) = 0.0;  dpsdpw(0, 4) = 0.0;   dpsdpw(0, 5) = 0.0;//dps/duw
            dpsdpw(1, 0) = sin8;     dpsdpw(1, 1) =  cos8;     dpsdpw(1, 2) = 0.0; dpsdpw(1, 3) = 0.0;  dpsdpw(1, 4) = 0.0;   dpsdpw(1, 5) = 0.0;//dps/dvw
            dpsdpw(2, 0) = 0.0;      dpsdpw(2, 1) = 0.0;       dpsdpw(2, 2) = 1.0; dpsdpw(2, 3) = 0.0;  dpsdpw(2, 4) = 0.0;   dpsdpw(2, 5) = 0.0;//dps/dww
            dpsdpw(3, 0) = 0.0;      dpsdpw(3, 1) = 0.0;       dpsdpw(3, 2) = y  ; dpsdpw(3, 3) = 1.0;  dpsdpw(3, 4) = -sin8; dpsdpw(3, 5) = 0.0;//dps/dalphaw
            dpsdpw(4, 0) = 0.0;      dpsdpw(4, 1) = 0.0;       dpsdpw(4, 2) = -x ; dpsdpw(4, 3) = sin8; dpsdpw(4, 4) = 1.0;   dpsdpw(4, 5) = 0.0;//dps/dbetaw
            dpsdpw(5, 0) = r*costm8; dpsdpw(5, 1) = -r*sintm8; dpsdpw(5, 2) = 0.0; dpsdpw(5, 3) = 0.0;  dpsdpw(5, 4) = 0.0;   dpsdpw(5, 5) = 1.0;//dps/dgammaw
          }
          else
          {
            dpsdpw(0, 0) = cos8;     dpsdpw(0, 1) = sin8;      dpsdpw(0, 2) = 0.0; dpsdpw(0, 3) = 0.0;   dpsdpw(0, 4) = 0.0;   dpsdpw(0, 5) = 0.0;//dps/duw
            dpsdpw(1, 0) = -sin8;    dpsdpw(1, 1) = cos8;      dpsdpw(1, 2) = 0.0; dpsdpw(1, 3) = 0.0;   dpsdpw(1, 4) = 0.0;   dpsdpw(1, 5) = 0.0;//dps/dvw
            dpsdpw(2, 0) = 0.0;      dpsdpw(2, 1) = 0.0;       dpsdpw(2, 2) = 1.0; dpsdpw(2, 3) = 0.0;   dpsdpw(2, 4) = 0.0;   dpsdpw(2, 5) = 0.0;//dps/dww
            dpsdpw(3, 0) = 0.0;      dpsdpw(3, 1) = 0.0;       dpsdpw(3, 2) = -y ; dpsdpw(3, 3) = 1.0;   dpsdpw(3, 4) = sin8;  dpsdpw(3, 5) = 0.0;//dps/dalphaw
            dpsdpw(4, 0) = 0.0;      dpsdpw(4, 1) = 0.0;       dpsdpw(4, 2) = -x ; dpsdpw(4, 3) = -sin8; dpsdpw(4, 4) = 1.0;   dpsdpw(4, 5) = 0.0;//dps/dbetaw
            dpsdpw(5, 0) = r*costm8; dpsdpw(5, 1) = r*sintm8;  dpsdpw(5, 2) = 0.0; dpsdpw(5, 3) = 0.0;   dpsdpw(5, 4) = 0.0;   dpsdpw(5, 5) = 1.0;//dps/dgammaw
          }
        }


        // Both halves share the same default coordinate system
        TMatrixD dpwdph(6, 6); // d(wedge alignment parameters)/d(half alignment parameters)

        double angle = 5./12.*TMath::Pi() - 1./6.*TMath::Pi()*wedgeId;
        double cosangle = TMath::Cos(angle);
        double cos90pangle = TMath::Cos(TMath::Pi()/2.+angle);
        double sinangle = TMath::Sin(angle);
        dpwdph(0, 0) = cosangle; dpwdph(0, 1) = -sinangle; dpwdph(0, 2) = 0.0;            dpwdph(0, 3) = 0.0;          dpwdph(0, 4) = 0.0;       dpwdph(0, 5) = 0.0;
        dpwdph(1, 0) = sinangle; dpwdph(1, 1) =  cosangle; dpwdph(1, 2) = 0.0;            dpwdph(1, 3) = 0.0;          dpwdph(1, 4) = 0.0;       dpwdph(1, 5) = 0.0;
        dpwdph(2, 0) = 0.0;      dpwdph(2, 1) = 0.0;       dpwdph(2, 2) = 1.0;            dpwdph(2, 3) = 0.0;          dpwdph(2, 4) = 0.0;       dpwdph(2, 5) = 0.0;
        dpwdph(3, 0) = 0.0;      dpwdph(3, 1) = 0.0;       dpwdph(3, 2) = 16.5*sinangle;  dpwdph(3, 3) = 1.0;          dpwdph(3, 4) = -sinangle; dpwdph(3, 5) = 0.0;
        dpwdph(4, 0) = 0.0;      dpwdph(4, 1) = 0.0;       dpwdph(4, 2) = -16.5*cosangle; dpwdph(4, 3) = -cos90pangle; dpwdph(4, 4) = 1.0;       dpwdph(4, 5) = 0.0;
        dpwdph(5, 0) = 0.0;      dpwdph(5, 1) = 16.5;      dpwdph(5, 2) = 0.0;            dpwdph(5, 3) = 0.0;          dpwdph(5, 4) = 0.0;       dpwdph(5, 5) = 1.0;
 

        TMatrixD drdpw = dpsdpw * drdps;         //d(residuals)/d(wedge alignment parameters)   6x2
        TMatrixD drdph = dpwdph * dpsdpw *drdps; //d(residuals)/d(half  alignment parameters)   6x2

        derGlobalFST(0, 0)  = 1.0;               derGlobalFST(1, 0)  = 0.0;           
        derGlobalFST(0, 1)  = 0.0;               derGlobalFST(1, 1)  = 1.0;           
        derGlobalFST(0, 2)  = -uSlope;           derGlobalFST(1, 2)  = -vSlope;       
        derGlobalFST(0, 3)  = -vPos * uSlope;    derGlobalFST(1, 3)  = -vPos * vSlope;
        derGlobalFST(0, 4)  = uPos * uSlope;     derGlobalFST(1, 4)  = uPos * vSlope; 
        derGlobalFST(0, 5)  = -vPos;             derGlobalFST(1, 5)  = uPos;          

        derGlobalFST(0, 6)  = drdpw(0, 0);       derGlobalFST(1, 6)  = drdpw(0, 1); 
        derGlobalFST(0, 7)  = drdpw(1, 0);       derGlobalFST(1, 7)  = drdpw(1, 1); 
        derGlobalFST(0, 8)  = drdpw(2, 0);       derGlobalFST(1, 8)  = drdpw(2, 1); 
        derGlobalFST(0, 9)  = drdpw(3, 0);       derGlobalFST(1, 9)  = drdpw(3, 1); 
        derGlobalFST(0, 10) = drdpw(4, 0);       derGlobalFST(1, 10) = drdpw(4, 1); 
        derGlobalFST(0, 11) = drdpw(5, 0);       derGlobalFST(1, 11) = drdpw(5, 1); 
         
        derGlobalFST(0, 12) = drdph(0, 0);       derGlobalFST(1, 12) = drdph(0, 1); 
        derGlobalFST(0, 13) = drdph(1, 0);       derGlobalFST(1, 13) = drdph(1, 1); 
        derGlobalFST(0, 14) = drdph(2, 0);       derGlobalFST(1, 14) = drdph(2, 1); 
        derGlobalFST(0, 15) = drdph(3, 0);       derGlobalFST(1, 15) = drdph(3, 1); 
        derGlobalFST(0, 16) = drdph(4, 0);       derGlobalFST(1, 16) = drdph(4, 1); 
        derGlobalFST(0, 17) = drdph(5, 0);       derGlobalFST(1, 17) = drdph(5, 1); 
          
        labGlobal.push_back(wedgelabel + 1); // u
        labGlobal.push_back(wedgelabel + 2); // v
        labGlobal.push_back(wedgelabel + 3); // w
        labGlobal.push_back(wedgelabel + 4); // alpha
        labGlobal.push_back(wedgelabel + 5); // beta
        labGlobal.push_back(wedgelabel + 6); // gamma

        labGlobal.push_back(halflabel + 1); // u
        labGlobal.push_back(halflabel + 2); // v
        labGlobal.push_back(halflabel + 3); // w
        labGlobal.push_back(halflabel + 4); // alpha
        labGlobal.push_back(halflabel + 5); // beta
        labGlobal.push_back(halflabel + 6); // gamma

        measPoint.addGlobals(labGlobal, derGlobalFST);
        
        //cout << "Global Derivatives" << endl;
        //derGlobalFST.Print();        

      }
      if(sensorId < 1000 || sensorId > 9000)
      {
        measPoint.addGlobals(labGlobal, derGlobalFTT);
      }
     

      listOfPoints.push_back(measPoint);
      listOfLayers.push_back((unsigned int) sensorId);
      n_gbl_points++;
      n_gbl_meas_points++;
    } else {
      // Incompatible measurement, store point without measurement
      gbl::GblPoint dummyPoint(jacPointToPoint);
      listOfPoints.push_back(dummyPoint);
      listOfLayers.push_back((unsigned int) sensorId);
      n_gbl_points++;
      skipMeasurement = false;
      #ifdef DEBUG
      cout << " Dummy point inserted. " << endl;
      #endif
    }
    //cout << "end of adding a point" << endl;


    //cout << " Starting extrapolation..." << endl;
    try {

      // Extrapolate to next measurement to get material distribution
      if (ipoint_meas < npoints_meas - 1) {
      //if ((!includeVertex && ipoint_meas < npoints_meas - 1) || (includeVertex && ipoint_meas > 0 && ipoint_meas < npoints_meas - 1)) {
	// Check if fitter info is in place
	if (!trk->getPoint(ipoint_meas + 1)->hasFitterInfo(rep)) {
	  cout << " ERROR: Measurement point does not have a fitter info. Track will be skipped." << endl;
	  return;
	}
	// Fitter of next point info which is only used now to get the plane
	genfit::KalmanFitterInfo* fi_i_plus_1 = dynamic_cast<genfit::KalmanFitterInfo*>(trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep));
	if (!fi_i_plus_1) {
	  cout << " ERROR: KalmanFitterInfo (with reference state) for measurement does not exist. Track will be skipped." << endl;
	  return;
	}
	genfit::StateOnPlane refCopy(*reference);
	// Extrap to point + 1, do NOT stop at boundary
	rep->extrapolateToPlane(refCopy, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, false);
	getScattererFromMatList(trackLen,
				scatTheta,
				scatSMean,
				scatDeltaS,
				trackMomMag,
				particleMass,
				particleCharge,
				rep->getSteps());
	// Now calculate positions and scattering variance for equivalent scatterers
	// (Solution from Claus Kleinwort (DESY))
	s1 = 0.;
	s2 = scatSMean + scatDeltaS * scatDeltaS / (scatSMean - s1);
        //cout << "scatTheta = " << scatTheta << endl;
        //cout << "scatDeltaS = " << scatDeltaS << endl;
        //cout << "scatSMean = " << scatSMean << endl;
	theta1 = sqrt(scatTheta * scatTheta * scatDeltaS * scatDeltaS / (scatDeltaS * scatDeltaS + (scatSMean - s1) * (scatSMean - s1)));
	theta2 = sqrt(scatTheta * scatTheta * (scatSMean - s1) * (scatSMean - s1) / (scatDeltaS * scatDeltaS + (scatSMean - s1) * (scatSMean - s1)));

	if (m_enableScatterers && !m_enableIntermediateScatterer) {
	  theta1 = scatTheta;
	  theta2 = 0;
	} else if (!m_enableScatterers) {
	  theta1 = 0.;
	  theta2 = 0.;
	}

	if (s2 < 1.e-4 || s2 >= trackLen - 1.e-4 || s2 <= 1.e-4) {
	  cout << " WARNING: GBL points will be too close. GBLTrajectory construction might fail. Let's try it..." << endl;
	}

      }
      // Return back to state on current plane
      delete reference;
      reference = new ReferenceStateOnPlane(*fi->getReferenceState());

      // If not last measurement, extrapolate and get jacobians for scattering points between this and next measurement
      //if ((!includeVertex && ipoint_meas < npoints_meas - 1) || (includeVertex && ipoint_meas > 0 && ipoint_meas < npoints_meas - 1)) {
      if (ipoint_meas < npoints_meas - 1) {
	if (theta2 > scatEpsilon) {
	  // First scatterer will be placed at current measurement point (see bellow)

	  // theta2 > 0 ... we want second scatterer:
	  // Extrapolate to s2 (remember we have s1 = 0)
	  rep->extrapolateBy(*reference, s2, false, true);
	  rep->getForwardJacobianAndNoise(jacMeas2Scat, noise, deltaState);
	  // Finish extrapolation to next measurement
	  double nextStep = rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, true);
	  //if (0. > nextStep) {
	  //cout << " ERROR: The extrapolation to measurement point " << (ipoint_meas + 2) << " stepped back by " << nextStep << "cm !!! Track will be skipped." << endl;
	  //  return;
	  //}
	  rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);

	} else {
	  // No scattering: extrapolate whole distance between measurements
	  rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, true);
	  //NOTE: we will load the jacobian from this extrapolation in next loop into measurement point
	  //jacPointToPoint.Print();
	  rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
	  //jacPointToPoint.Print();
	}
      }
    } catch (...) {
      cout << " ERROR: Extrapolation failed. Track will be skipped." << endl;
      return;
    }
    //cout << " Extrapolation finished." << endl;


    // Now store scatterers if not last measurement and if we decided
    // there should be scatteres, otherwise the jacobian in measurement
    // stored above is already correct
    if (ipoint_meas < npoints_meas - 1) {
    //if ((!includeVertex && ipoint_meas < npoints_meas - 1) || (includeVertex && ipoint_meas > 0 && ipoint_meas < npoints_meas - 1)) {

      if (theta1 > scatEpsilon) {
	// We have to insert first scatterer at measurement point
	// Therefore (because state is perpendicular to plane, NOT track)
	// we have non-diaonal matrix of multiple scattering covariance
	// We have to project scattering into plane coordinates
	double c1 = trackDir.Dot(plane->getU());
	double c2 = trackDir.Dot(plane->getV());
	TMatrixDSym scatCov(2);
	scatCov(0, 0) = 1. - c2 * c2;
	scatCov(1, 1) = 1. - c1 * c1;
	scatCov(0, 1) = c1 * c2;
	scatCov(1, 0) = c1 * c2;
	scatCov *= theta1 * theta1 / (1. - c1 * c1 - c2 * c2) / (1. - c1 * c1 - c2 * c2) ;

        //cout << "theta1 = " << theta1 << endl;
        //cout << "c1 = " << c1 << endl;
        //cout << "c2 = " << c2 << endl;
        //cout << "1. - c1 * c1 - c2 * c2 = " << 1. - c1 * c1 - c2 * c2 << endl;  

        //cout << "scatCov" << endl;
        //scatCov.Print();
	// last point is the just inserted measurement (or dummy point)
	gbl::GblPoint& lastPoint = listOfPoints.at(n_gbl_points - 1);
	lastPoint.addScatterer(scatResidual, scatCov.Invert());

      }

      if (theta2 > scatEpsilon) {
	// second scatterer is somewhere between measurements
	// TrackRep state is perpendicular to track direction if using extrapolateBy (I asked Johannes Rauch),
	// therefore scattering covariance is diagonal (and both elements are equal)
	TMatrixDSym scatCov(2);
	scatCov.Zero();
	scatCov(0, 0) = theta2 * theta2;
	scatCov(1, 1) = theta2 * theta2;

	gbl::GblPoint scatPoint(jacMeas2Scat);
	scatPoint.addScatterer(scatResidual, scatCov.Invert());
	listOfPoints.push_back(scatPoint);
	listOfLayers.push_back(((unsigned int) sensorId) + 0.5);
	n_gbl_points++;
      }


    }
    // Free memory on the heap
    delete reference;
  }
  // We should have at least two measurement points to fit anything
  if (n_gbl_meas_points > 1) {
    int fitRes = -1;
    double pvalue = 0.;
    gbl::GblTrajectory* traj = 0;
    try {
      // Construct the GBL trajectory, seed not used
      traj = new gbl::GblTrajectory(listOfPoints, seedLabel, clSeed, fitQoverP);
      // Fit the trajectory
      fitRes = traj->fit(Chi2, Ndf, lostWeight, m_gblInternalIterations);
      if (fitRes != 0) {
        #ifdef DEBUG
        //traj->printTrajectory(100);
        //traj->printData();
        //traj->printPoints(100);
        #endif
      }
    } catch (...) {
      cout << "GBL FAILED" << endl;
      // Gbl failed critically (usually GblPoint::getDerivatives ... singular matrix inversion)
      return;
    }
    
    pvalue = TMath::Prob(Chi2, Ndf);
    
    //traj->printTrajectory(100);
    //traj->printData();
    //traj->printPoints(100);
    
    #ifdef OUTPUT
    // Fill histogram with fit result
    fitResHisto->Fill(fitRes);
    ndfHisto->Fill(Ndf);
    #endif
    
    //#ifdef DEBUG
    cout << " Ref. Track Chi2      :  " << trkChi2 << endl;
    cout << " Ref. end momentum    :  " << trackMomMag << " GeV/c ";
///    if (abs(trk->getCardinalRep()->getPDG()) == 11) {
///      if (particleCharge < 0.)
///        cout << "(electron)";
///      else
///        cout << "(positron)";
///    }
    cout << endl;
    cout << "------------------ GBL Fit Results --------------------" << endl;
    cout << " Fit q/p parameter    :  " << ((fitQoverP) ? ("True") : ("False")) << endl;
    cout << " Valid trajectory     :  " << ((traj->isValid()) ? ("True") : ("False")) << endl;
    cout << " Fit result           :  " << fitRes << "    (0 for success)" << endl;
    cout << " # GBL meas. points   :  " << n_gbl_meas_points << endl;
    cout << " # GBL all points     :  " << n_gbl_points << endl;
    cout << " GBL track NDF        :  " << Ndf << "    (-1 for failure)" << endl;
    cout << " GBL track Chi2       :  " << Chi2 << endl;
    cout << " GBL track P-value    :  " << pvalue;
    if (pvalue < m_pValueCut)
      cout << " < p-value cut " << m_pValueCut;
    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    //#endif
    
    #ifdef OUTPUT
    bool hittedLayers[12];
    for (int hl = 0; hl < 12; hl++) {
      hittedLayers[hl] = false;
    }
    #endif
    
    // GBL fit succeded if Ndf >= 0, but Ndf = 0 is useless
    //TODO: Here should be some track quality check
    //    if (Ndf > 0 && fitRes == 0) {
    if (traj->isValid() && fitRes == 0/*&& pvalue >= m_pValueCut && Ndf >= m_minNdf*/) {
      
      // In case someone forgot to use beginRun and dind't reset mille file name to ""
      if (!milleFile && m_milleFileName != "")
        milleFile = new gbl::MilleBinary(m_milleFileName);
      
      // Loop over all GBL points
      int imeas = 0;
      for (unsigned int p = 0; p < listOfPoints.size(); p++) {
        unsigned int label = p + 1;
        //cout << "label = " << label << endl;
        unsigned int numRes;
        TVectorD residuals(5);
        TVectorD measErrors(2);
        TVectorD resErrors(2);
        TVectorD downWeights(2);
        //TODO: now we only provide info about measurements, not kinks
        if (!listOfPoints.at(p).hasMeasurement())
          continue;
        
        #ifdef OUTPUT
        // Decode VxdId and get layer in TB setup
        unsigned int l = 0;
        unsigned int id = listOfLayers.at(p);
        // skip segment (5 bits)
        id = id >> 5;
        unsigned int sensor = id & 7;
        id = id >> 3;
        unsigned int ladder = id & 31;
        id = id >> 5;
        unsigned int layer = id & 7;
        
        if (layer == 7 && ladder == 2) {
          l = sensor;
        } else if (layer == 7 && ladder == 3) {
          l = sensor + 9 - 3;
        } else {
          l = layer + 3;
        }
        
        hittedLayers[l - 1] = true;
        #endif
        TVectorD localPar(5);
        TMatrixDSym localCov(5);
        traj->getResults(label, localPar, localCov);
        // Get GBL fit results
        traj->getMeasResults(label, numRes, residuals, measErrors, resErrors, downWeights);
        //residuals.Print();
        //residualsUgbl[imeas]->Fill(residuals[0]);
        //residualsVgbl[imeas]->Fill(residuals[1]);
        //merrUgbl[imeas]->Fill(measErrors[0]);
        //merrVgbl[imeas]->Fill(measErrors[1]);
        //rerrUgbl[imeas]->Fill(resErrors[0]);
        //rerrVgbl[imeas]->Fill(resErrors[1]);
        //dweightUgbl[imeas]->Fill(downWeights[0]);
        //dweightVgbl[imeas]->Fill(downWeights[1]);
        //cout <<"imeas = " << imeas << endl;
        //cout << "residuals = " << residuals[0] << "    "       << residuals[0]   << endl;;
        //cout << "measErrors = " << measErrors[0] << "    "     << measErrors[0]  << endl;;
        //cout << "resErrors = " <<  resErrors[0] << "    "      << resErrors[0]   << endl;;
        //cout << "downWeights = " <<  downWeights[0] << "    "  << downWeights[0] << endl;;






        imeas++;
        if (m_chi2Cut && (fabs(residuals[0] / resErrors[0]) > m_chi2Cut || fabs(residuals[1] / resErrors[1]) > m_chi2Cut))
          return;
        // Write layer-wise data
        #ifdef OUTPUT
        if (!writeHistoDataForLabel(listOfLayers.at(p), residuals, measErrors, resErrors, downWeights, localPar, localCov))
          return;
        #endif
        
      } // end for points
      
      // Write binary data to mille binary
      #ifndef DEBUG
      ///cout << "pvalue = " << pvalue << "   pvaluecut = " << m_pValueCut << endl;
      ///cout << "Ndf = " << Ndf << "    m_minNdf = " << m_minNdf << endl;
      #endif
      
      if (milleFile && m_milleFileName != "" && fitRes == 0 /*&& pvalue >= m_pValueCut && Ndf >= m_minNdf*/ && traj->isValid()) {
        //cout << "Trajectory is valid? " << traj->isValid() << endl;
        setSuccessfulFitFlag(true);
        traj->milleOut(*milleFile);  
        chi2OndfHistoGBL->Fill(Chi2 / Ndf);
        pValueHistoGBL->Fill(TMath::Prob(Chi2, Ndf));
        chi2OndfHisto->Fill(trkChi2 / trkNdf);
        pValueHisto->Fill(TMath::Prob(trkChi2, trkNdf));
        //cout << "traj->printTrajectory(1)" << endl;;
        //traj->printTrajectory(1);
        //cout << "traj->printPoints(1)" << endl;;
        //traj->printPoints(1);
        //cout << "traj->printData();" << endl;
        //traj->printData();
        #ifdef DEBUG
        cout << " GBL Track written to Millepede II binary file." << endl;
        cout << "-------------------------------------------------------" << endl;
        #endif
      }
      
      #ifdef OUTPUT
      // Fill histograms
      chi2OndfHisto->Fill(Chi2 / Ndf);
      pValueHisto->Fill(TMath::Prob(Chi2, Ndf));
      // track counting
      stats->Fill(0);
      // hitted sensors statistics
      if (
        hittedLayers[0] &&
        hittedLayers[1] &&
        hittedLayers[2] &&
        hittedLayers[3] &&
        hittedLayers[4] &&
        hittedLayers[5] &&
        hittedLayers[6] &&
        hittedLayers[7] &&
        hittedLayers[8]
      ) {
        // front tel + pxd + svd
        stats->Fill(1);
      }
      
      if (
        hittedLayers[0] &&
        hittedLayers[1] &&
        hittedLayers[2] &&
        hittedLayers[3] &&
        hittedLayers[4] &&
        hittedLayers[5] &&
        hittedLayers[6] &&
        hittedLayers[7] &&
        hittedLayers[8] &&
        hittedLayers[9] &&
        hittedLayers[10] &&
        hittedLayers[11]
      ) {
        // all layers
        stats->Fill(2);
      }
      if (
        hittedLayers[3] &&
        hittedLayers[4] &&
        hittedLayers[5] &&
        hittedLayers[6] &&
        hittedLayers[7] &&
        hittedLayers[8]
      ) {
        // vxd
        stats->Fill(3);
      }
      if (
        hittedLayers[5] &&
        hittedLayers[6] &&
        hittedLayers[7] &&
        hittedLayers[8]
      ) {
        // svd
        stats->Fill(4);
      }
      if (
        hittedLayers[0] &&
        hittedLayers[1] &&
        hittedLayers[2] &&
        
        hittedLayers[5] &&
        hittedLayers[6] &&
        hittedLayers[7] &&
        hittedLayers[8] &&
        hittedLayers[9] &&
        hittedLayers[10] &&
        hittedLayers[11]
      ) {
        // all except pxd
        stats->Fill(5);
      }
      if (
        hittedLayers[0] &&
        hittedLayers[1] &&
        hittedLayers[2] &&
        
        hittedLayers[5] &&
        hittedLayers[6] &&
        hittedLayers[7] &&
        hittedLayers[8]
      ) {
        // front tel + svd
        stats->Fill(6);
      }
      if (
        hittedLayers[9] &&
        hittedLayers[10] &&
        hittedLayers[11]
      ) {
        // backward tel
        stats->Fill(7);
      }
      #endif
    }
    
    // Free memory
    delete traj;
  }
}

