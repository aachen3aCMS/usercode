#include "Analysis.h"

// ROOT includes
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>
#include <TList.h>
#include <TKey.h>

// C / C++ includes
#include <iostream>
#include <sstream>
#include <cmath>

#include "Utilities.h"

using namespace std;

Analysis::Analysis(TTree & inputTree, TTree & outputTree) 
  : TreeContent(& inputTree) , fInputTree(inputTree), fOutputTree(outputTree)
{
  h_counters = 0;
}

// main event loop
void Analysis::Loop() 
{
  // process and memory info
  ProcInfo_t info;

  // set branch addresses in output tree if necessary
  SetBranchAddresses();
  CreateHistograms();

  // main event loop
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  INFO("Analysis: Chain contains " << nentries << " events");
  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) 
      break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    // some progress / process indicators
    if (fmod((double)jentry,(double)1000)==0) 
      INFO("Event " << jentry);
    if (fmod((double)jentry,(double)10000)==0) {
      gSystem->GetProcInfo(&info);
      INFO("Memory in MB : " << info.fMemResident/1000. << " (resident) " 
	   << info.fMemVirtual/1000. << " (virtual) ");
    }

    //////////////////////////////////////////////////////////////////////
    // Start here implementing your cuts. Use histograms to check later
    // if everything went as expected.

    // Fill histograms before applying a cut
    Fill("n_muo_n", muo_n);
    if (muo_n > 0) {
      Fill("n_muo_pt0", muo_pt[0]);
    }
    if (muo_n > 1) {
      Fill("n_muo_pt1", muo_pt[1]);
    }

    // select muons with pt cut    
    bool rejection = true;
    if (muo_n >= 2 && muo_pt[0] > 15. && muo_pt[1] > 7.) {
      rejection = false;
    }

    // apply cut
    if (rejection)
      continue;

    // fill histograms after applying a cut
    Fill("muo_n", muo_n);
    Fill("muo_pt0", muo_pt[0]);
    Fill("muo_pt1", muo_pt[1]);

    // write accepted events to output file
    fOutputTree.Fill();
  }
}

void Analysis::CreateHistograms()
{
  //////////////////////////////////////////////////////////////////////
  // Create histograms here that you are going to fill in the loop
  CreateHisto("n_muo_n", "Number of muons", 10, -0.5, 9.5);
  CreateHisto("n_muo_pt0", "pt of leading muon", 500, 0, 500);
  CreateHisto("n_muo_pt1", "pt of second muon", 500, 0, 500);
  CreateHisto("muo_n", "Number of muons", 10, -0.5, 9.5);
  CreateHisto("muo_pt0", "pt of leading muon", 500, 0, 500);
  CreateHisto("muo_pt1", "pt of second muon (skimmer cuts)", 500, 0, 500);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup)
{
  histo[name] = new TH1D(Form("h1_%d_%s", 0, name), title, nbinsx, xlow, xup);
}

void Analysis::CreateHisto(const char * name, const char * title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  histo2[name] = new TH2D(Form("h2_%d_%s", 0, name), title, 
			  nbinsx, xlow, xup, nbinsy, ylow, yup);
}

void Analysis::Fill(const char * name, double value)
{
  histo[name]->Fill(value);
}

void Analysis::Fill(const char * name, double x, double y)
{
  histo2[name]->Fill(x, y);
}

void Analysis::SetBranchAddresses()
{
  //////////////////////////////////////////////////////////////////////
  // initialize output tree branches
  INFO("Setting output branch addresses");
  TObjArray * branchlist = fInputTree.GetListOfBranches();
  TIter next(branchlist);
  while (TBranch * branch = (TBranch *) next()) {
    fOutputTree.SetBranchAddress(branch->GetName(), branch->GetAddress());
    LOG(4, "branch " << branch->GetName() << " has address " << (void *) branch->GetAddress());
  }
  // tree size
  const UInt_t fMaxTreeSize = 100000000;
  INFO("Setting maximum tree size to " << fMaxTreeSize);
  fOutputTree.SetMaxTreeSize(fMaxTreeSize);
}

Bool_t Analysis::Notify()
{
  INFO("A new file has been opened: " << fChain->GetTreeNumber());

  // list all histograms in current file
  TFile * f = fChain->GetCurrentFile();
  TH1F * h1 = (TH1F *) f->Get("ACSkimAnalysis/h_counters");
  if (h1 == 0) {
    ERROR("Could not read h_counters from file");
    return kFALSE;
  }
  if (h_counters != 0) {
    INFO("Found existing histogram " << h1->GetName());
    h_counters->Add(h1);
  }
  else {
    INFO("Found new histogram " << h1->GetName());
    h_counters = h1;
    h_counters->SetDirectory(fOutputTree.GetCurrentFile()->GetDirectory("ACSkimAnalysis"));
  }

  return kTRUE;
}
