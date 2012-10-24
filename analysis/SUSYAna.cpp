#include "SUSYAna.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include "TSystem.h"

#include "SUSYDump.h"

using std::cout;
using std::endl;
using namespace std;

const int met_index = 1; // 0 = RAW PFMET, 1 = TYPE1 corrected MET, 2 = TYPE 0 corrected MET

SUSYAna::SUSYAna(TTree *tree) : TreeContent(tree) 
{
  if (tree ==  0) {
    cerr << "SUSYAna: ERROR - Tree empty !!!" << endl;
    exit(0);
  }
}

// main event loop
void SUSYAna::Loop(TString fout, bool debug, TString type) {

  DEBUG = debug;

  ProcInfo_t info;

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  
  cout << "SUSYAna: Running over " << nentries << " events" << endl;

  if (DEBUG) 
    cout << "SUSYAna: Debug mode" << endl;

  if (type == "none") 
    cout << "SUSYAna: You are running with JES flag 'none'" << endl;
  else if (type == "up")
    cout << "SUSYAna: You are running with JES flag 'up'" << endl;
  else if (type == "down")
    cout << "SUSYAna: You are running with JES flag 'down'" << endl;
  else {
    cout << "SUSYAna: No JES flag specified ! Assuming JES flag 'none'" << endl;
    type = "none";
  }

  init(1);

  Long64_t nbytes = 0, nb = 0;

  // main event loop
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (fmod((double)jentry,(double)1000)==0) 
      cout << "Event " << jentry << endl;

    if (fmod((double)jentry,(double)10000)==0) {
      gSystem->GetProcInfo(&info);
      cout << " -> Memory in MB : " << info.fMemResident/1000. << " (resident) " 
	   << info.fMemVirtual/1000. << " (virtual) " << endl;
    }

    // Dump event/object information, see SUSYDump.h
    if (DEBUG) {
      BasicDump(jentry);
      TriggerDump("*"); // select e.g all muon trigger with "Mu"
      TruthDump();
      VertexDump();
      MuonDump(0);      // [less] detailed information 1 [0]
      TruthJetDump();
      PFJetDump();
      METDump();
      SCDump();
      EleDump(0);       // [less] detailed information 1 [0]
      PFEleDump(0);       // [less] detailed information 1 [0]
    }
    TriggerDump("*"); // select e.g all muon trigger with "Mu"
  }
  write(fout, 2);
}

// duplicates finder
bool SUSYAna::find_duplicate(int run, int evt, double x1, double x2) {

  pair<int,int> temp1(run,evt);
  pair<double,double> temp2(x1,x2);
  Key key (temp1, temp2);
  
  KeyIter pos = _keys.find (key);
  
  if (pos == _keys.end()) {
    _keys.insert (key);
    return false;
  }
  else {
    if (DEBUG) cout << "SUSYAna: duplicate run " << run << " , evt " << evt << " , scale " << x1 << " , vtx_z " << x2 << endl;
    return true;
  }
}

// initialize histograms
void SUSYAna::init(int nstages) {
  
  Char_t h_name[200];
  Char_t h_title[200];
  
  TH1F* h_temp;
  
  for (Int_t jj=0; jj<nstages; jj++) {
    
    strcpy(h_title,"\0");
    sprintf(h_title,"Stage %i", jj);
    
    strcpy(h_name,"\0");
    sprintf(h_name,"h1_%i_mu_pt", jj);
    h_temp = new TH1F(h_name, h_title, 350, 0., 7000.);
    h_temp->Sumw2();
    h_temp->SetXTitle("p_{T}^{#mu}  [GeV]");
    h1_mu_pt.push_back(h_temp);
    
    strcpy(h_name,"\0");
    sprintf(h_name,"h1_%i_mu_mu_Minv", jj);
    h_temp = new TH1F(h_name, h_title, 350, 0., 7000.);
    h_temp->Sumw2();
    h_temp->SetXTitle("M_{#mu #mu} [GeV]");
    h1_mu_Minv.push_back(h_temp);
  }

  h1_Z_mass = new TH1F("h1_Z_mass", "Z mass (hopefully!)", 200, 0, 200);
}

// write histograms to file
void SUSYAna::write(TString fout, int nstages) {
  
  TFile *f = new TFile(fout,"recreate");
 
  // for (Int_t jj=0; jj<nstages; jj++) {
    
  //   h1_mu_pt[jj]     -> Write();
  //   h1_mu_TrkIso[jj] -> Write();

  // }
  h1_Z_mass->Write();
  
  f->Close();
  delete f;
}
// helper
double SUSYAna::M2(double E1, double px1, double py1, double pz1, double E2, double px2, double py2, double pz2) {
  double InvMass2 =   (E1 + E2)*(E1 + E2)
    - (px1 + px2)*(px1 + px2)
    - (py1 + py2)*(py1 + py2)
    - (pz1 + pz2)*(pz1 + pz2);
  if (InvMass2 < 0) cout << "Mass Square negative: InvariantMass" << InvMass2 << endl; 
  return std::sqrt(InvMass2);
}

double SUSYAna::M2v2(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2){
  TLorentzVector lorentz1;
  TLorentzVector lorentz2;
  lorentz1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
  lorentz2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
  TLorentzVector lorentz1plus2 = lorentz1 + lorentz2;
  return lorentz1plus2.M();
}

double SUSYAna::AngleMuons(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2){
  TLorentzVector lorentz1;
  TLorentzVector lorentz2;
  lorentz1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
  lorentz2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
   
  double part1_px = lorentz1.Px();
  double part1_py = lorentz1.Py();
  double part1_pz = lorentz1.Pz();
  double part1_p = lorentz1.P();

  double part2_px = lorentz2.Px();
  double part2_py = lorentz2.Py();
  double part2_pz = lorentz2.Pz();
  double part2_p = lorentz2.P();

  //see exotica muon twiki for documentation
  double angle = acos(((-1.*part1_px * part2_px)+(-1.*part1_py * part2_py)+(-1.*part1_pz * part2_pz))/(part1_p * part2_p));
  return angle;

}

Double_t SUSYAna::DeltaPhi(double a, double b) {

  double temp = fabs(a-b);
  if (temp <= TMath::Pi())
    return temp;
  else
    return  2.*TMath::Pi() - temp;
}

Double_t SUSYAna::mT(double et1, double phi1, double et2, double phi2) {

  double mm = 2 * et1 * et2 * ( 1. - cos(phi1 - phi2) );
  return sqrt(mm);

}

Double_t SUSYAna::MinDEt(const std::vector<TLorentzVector> & objects, 
			 std::vector<UInt_t> * lista, 
			 std::vector<UInt_t> * listb) {
  
  // Find the combination with the lowest DEt
  UInt_t n = objects.size();
  if (n==0) return 0.;
  if (n==1) return objects[0].Et();
  if (n>10)  {
    cout << "MinDEt: n too big : " << n << endl;
    return -1;
  }
  if (lista!=0 && listb!=0) { lista->clear(); listb->clear(); }
  
  double mindiff = 1000000000., diff = 0.;
  
  // Determine the combination that minimises the difference
  std::vector< std::vector<UInt_t> > combinationset1;
  std::vector< std::vector<UInt_t> > combinationset2;
  Combinations::mycombinations(n, combinationset1, combinationset2);
  
  if (combinationset1.size() != combinationset2.size() ) {
    cout << "MinDEt: Combination set sizes to not match - something has gone wrong..." << endl;
  }
  
  for (UInt_t set = 0; set<combinationset1.size(); set++) {
    
    std::vector<UInt_t> la = combinationset1[set]; //!< Temporary list a for calculating best combo
    std::vector<UInt_t> lb = combinationset2[set]; //!< Temporary list b for calculating best combo
    
    Double_t aEt = 0., bEt = 0.;
    for (std::vector<UInt_t>::iterator ia=la.begin();ia!=la.end();++ia) {
      //      cout << (*ia) << " ";
      aEt += objects[ (*ia) ].Et();
    }
    //    cout << ", ";
    for (std::vector<UInt_t>::iterator ib=lb.begin();ib!=lb.end();++ib) {
      bEt += objects[ (*ib) ].Et();
      //      cout << (*ib) << " ";
    }
    //    cout << endl;
    diff = fabs(aEt - bEt);
    //    cout << "Difference in Et is " << diff << endl;
    if (diff < mindiff) {
      mindiff = diff;
      if (lista!=0 && listb!=0) { *lista = la; *listb = lb; }
    }
    la.clear(); lb.clear();
    
  } // end of loop over combination sets
  
  //  cout << "Minimum difference is " << mindiff << endl << endl << endl;
  //  cout << "===========================================" << endl;
  
  return mindiff;  
}

Double_t SUSYAna::AlphaT(const std::vector<TLorentzVector> & objects) {

  return 0.5*((SumET(objects) - MinDEt(objects))/(MT(objects)));

}

Double_t SUSYAna::SumET(const std::vector<TLorentzVector> & objects) {

  Double_t sEt = 0;    
  for (std::vector<TLorentzVector>::const_iterator o=objects.begin();o!=objects.end();++o) { 
    sEt+=o->Et(); 
  }
  return sEt;

}

Double_t SUSYAna::MT(const std::vector<TLorentzVector> & objects) {

  Double_t sEt = 0, sPx = 0, sPy = 0;
  for (std::vector<TLorentzVector>::const_iterator o=objects.begin();o!=objects.end();++o) { 
    sEt+=o->Et();sPx+=o->Px();sPy+=o->Py(); 
  }
  Double_t MTsq = sEt*sEt - sPx*sPx - sPy*sPy;
  return MTsq >= 0. ? sqrt(MTsq) : -sqrt(-MTsq);

}

// JES variation: recalculate Jets and propagate into MET
void SUSYAna::doJESandRecalculateMET(TString corr) {

  // WARNING : Not propagated into TCMET since this need track corrected Jets
  Double_t pfscale;

  Double_t dMEx;
  Double_t dMEy;
  Double_t dSumEt;

  // JME-10-010-pas-v8.pdf
  if (corr == "none")
    return;
  else if (corr == "up") {
    pfscale   = 1.04;
  }  
  else if (corr == "down") {
    pfscale   = 0.96;
  }
  else {
    cout << "SUSYAna::doJESandRecalculateMET: Option '" << corr << "' not known." << endl;
    cout << "                                 No correction applied." << endl;
    return;
  }
  if (DEBUG) cout << "SUSYAna::doJESandRecalculateMET called with option '" << corr << "'" << endl;
  
  // PF Jets
  dMEx   = 0.;
  dMEy   = 0.;
  dSumEt = 0.;
  
  for (int i=0; i<pfjet_n; i++) {

    // cleaned PF jets, have ptraw > 10 GeV
    dMEx   += (pfscale - 1.) * pfjet_px[i];
    dMEy   += (pfscale - 1.) * pfjet_py[i];
    dSumEt += (pfscale - 1.) * pfjet_Et[i];
    
    pfjet_E[i]  *= pfscale;
    pfjet_Et[i] *= pfscale;  
    pfjet_p[i]  *= pfscale;   
    pfjet_pt[i] *= pfscale;  
    pfjet_px[i] *= pfscale;
    pfjet_py[i] *= pfscale;
    pfjet_pz[i] *= pfscale;
  }
  
  met_sumet[met_index] += dSumEt;
  met_ex[met_index]    -= dMEx;
  met_ey[met_index]    -= dMEy;
  
  TVector3 mpf(met_ex[met_index], met_ey[met_index], 0.);

  met_phi[met_index] = mpf.Phi();
  met_et[met_index]  = mpf.Perp();

}
