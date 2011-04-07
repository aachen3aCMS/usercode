#define SUSYAna_cxx
#include "SUSYAna.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TSystem.h"

#include "SUSYDump.h"

using std::cout;
using std::endl;
using namespace std;

SUSYAna::SUSYAna(TTree *tree) : TreeContent(tree) {

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

  if (DEBUG) cout << "SUSYAna: Debug mode" << endl;

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

  TString trigger;

  init(2);

  Long64_t nbytes = 0, nb = 0;

  // main event loop
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //    if (DEBUG && jentry == 10000) break;
    
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
      CaloJetDump();
      PFJetDump();
      METDump();
      SCDump();
      EleDump(0);       // [less] detailed information 1 [0]
      PFEleDump(0);       // [less] detailed information 1 [0]
    }

    doJESandRecalculateMET(type);

    if (DEBUG) {
      if (type == "up" || type == "down") {
	CaloJetDump();
	PFJetDump();
	METDump();
      }
      continue;
    }

    double v=0.;
    if (vtx_n>0) v = vtx_z[0]; 

    if (find_duplicate(global_run, global_event, pdf_scale, v)) {
      cout << "SUSYAna: found duplicate " << endl;
    }

    // do your analysis
    cout << endl;
    cout << "Analysis example output " << endl;

    // list Trigger Objects
    if (DEBUG) {
      for (int i=0; i<trig_n; i++) {
	cout << "  " << unpack(trig_name[i]) << "  "  << trig_L1prescale[i] << "  " << trig_HLTprescale[i] 
	     << "  " << trig_pt[i] << "  " << trig_eta[i] << "  " << trig_phi[i] << endl;
      }
    }

    for (int j=0; j<muo_n; j++) {

      cout << "pt eta phi = " << muo_pt[j] << " " << muo_eta[j] << " " << muo_phi[j] << "    ID [ ";
      for (int k=0; k<24; k++)
	cout << muo_ID[j][k] << " ";
      cout << "] " << endl;
      if (muo_ID[j][6] == 1)
	cout << "  -> global prompt tight"<< endl;
      
      // list matched trigger objects
      for (int i=0; i<muo_trign[j]; i++) {
	int id = muo_trig[j][i];
	TString name = unpack((int*)trig_name[id]);
	TString filt = unpack(trig_filter[id]);
	cout << "     #" << i << "  " << name << "  " << filt << "  " << trig_pt[id] 
	     << "  " << trig_phi[id] << "  " << trig_pt[id] <<endl;
      }


      h1_mu_pt[0]     -> Fill(muo_pt[j]);
      h1_mu_TrkIso[0] -> Fill(muo_TrkIso[j]);

      if (muo_TrkIso[j] < 6. && muo_TrkChiNormCm[j]<10) {
	h1_mu_pt[1]     -> Fill(muo_pt[j]);
	h1_mu_TrkIso[1] -> Fill(muo_TrkIso[j]);

      }
    }
    
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
    h_temp = new TH1F(h_name, h_title, 200, 0., 200.);
    h_temp->Sumw2();
    h_temp->SetXTitle("p_{T}^{#mu}  [GeV]");
    h1_mu_pt.push_back(h_temp);
    
    strcpy(h_name,"\0");
    sprintf(h_name,"h1_%i_mu_TrkIso", jj);
    h_temp = new TH1F(h_name, h_title, 200, 0., 200.);
    h_temp->Sumw2();
    h_temp->SetXTitle("TrkIso");
    h1_mu_TrkIso.push_back(h_temp);
   }
}

// write histograms to file
void SUSYAna::write(TString fout, int nstages) {
  
  TFile *f = new TFile(fout,"recreate");
 
  for (Int_t jj=0; jj<nstages; jj++) {
    
    h1_mu_pt[jj]     -> Write();
    h1_mu_TrkIso[jj] -> Write();
    
  }
  
  f->Close();
   
}
// helper
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
  Double_t caloscale;
  Double_t pfscale;

  Double_t dMEx;
  Double_t dMEy;
  Double_t dSumEt;

  // JME-10-010-pas-v8.pdf
  if (corr == "none")
    return;
  else if (corr == "up") {
    caloscale = 1.04;
    pfscale   = 1.04;
  }  
  else if (corr == "down") {
    caloscale = 0.96;
    pfscale   = 0.96;
  }
  else {
    cout << "SUSYAna::doJESandRecalculateMET: Option '" << corr << "' not known." << endl;
    cout << "                                 No correction applied." << endl;
    return;
  }
  if (DEBUG) cout << "SUSYAna::doJESandRecalculateMET called with option '" << corr << "'" << endl;
  
  dMEx   = 0.;
  dMEy   = 0.;
  dSumEt = 0.;
  
  // Calo Jets
  for (int i=0; i<calojet_n; i++) {

    // not final, best guess
    if (calojet_pt_raw[i]>20. && calojet_fem[i]<0.9) {
      dMEx   += (caloscale - 1.) * calojet_px[i];
      dMEy   += (caloscale - 1.) * calojet_py[i];
      dSumEt += (caloscale - 1.) * calojet_Et[i];
    }
    calojet_E[i]  *= caloscale;
    calojet_Et[i] *= caloscale;  
    calojet_p[i]  *= caloscale;   
    calojet_pt[i] *= caloscale;  
    calojet_px[i] *= caloscale;
    calojet_py[i] *= caloscale;
    calojet_pz[i] *= caloscale;
  }

  met_sumet[0] += dSumEt;
  met_ex[0]    -= dMEx;
  met_ey[0]    -= dMEy;

  TVector3 mcalo(met_ex[0], met_ey[0], 0.);

  met_phi[0] = mcalo.Phi();
  met_et[0]  = mcalo.Perp();

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
  
  met_sumet[3] += dSumEt;
  met_ex[3]    -= dMEx;
  met_ey[3]    -= dMEy;
  
  TVector3 mpf(met_ex[3], met_ey[3], 0.);

  met_phi[3] = mpf.Phi();
  met_et[3]  = mpf.Perp();

}
