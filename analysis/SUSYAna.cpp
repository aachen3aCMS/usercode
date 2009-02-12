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
void SUSYAna::Loop() {

  ProcInfo_t info;

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  
  cout << "Running over " << nentries << " events" << endl;

  if (DEBUG) cout << " -> in debug mode" << endl;

  init(2);

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
      cout << " -> Memory in kB : " << info.fMemResident << " (resident) " 
	   << info.fMemVirtual << " (virtual) " << endl;
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

  write(2);
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
void SUSYAna::write(int nstages) {
  
  TFile *f = new TFile(_fname,"recreate");
 
  for (Int_t jj=0; jj<nstages; jj++) {
    
    h1_mu_pt[jj]     -> Write();
    h1_mu_TrkIso[jj] -> Write();
    
  }
  
  f->Close();
   
}
