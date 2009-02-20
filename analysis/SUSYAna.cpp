#define SUSYAna_cxx
#include "SUSYAna.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TSystem.h"

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

  if (DEBUG) cout << "SUSYAna: DEBUG modus " << endl;

  if (type == "data") 
    cout << "SUSYAna: You are running with flag 'data'" << endl;
  else if (type == "mc")
    cout << "SUSYAna: You are running with flag 'mc'" << endl;
  else 
    cout << "SUSYAna: No flag specified !" << endl;

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

    trigger = global_HLT;
    if (trigger.Contains(":HLT_Mu7:")) {
      //      cout << "SUSYAna: " << trigger.SubString("HLT_Mu7") << endl;
    }

    double v=0.;
    if (vtx_n>0) v = vtx_z[0]; 

    if (find_duplicate(global_run, global_event, pdf_scale, v)) {
      cout << "SUSYAna: found duplicate " << endl;
    }

    // do your analysis
    for (int j=0; j<muo_n; j++) {

      h1_mu_pt[0]     -> Fill(muo_pt[j]);
      h1_mu_TrkIso[0] -> Fill(muo_TrkIso[j]);

      if (muo_TrkIso[j] < 6.) {
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
