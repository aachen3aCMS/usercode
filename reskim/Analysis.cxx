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
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TEntryList.h>

// C / C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/net/scratch_cms/institut_3a/padeken/lhapdfLib/include/LHAPDF/LHAPDF.h"
#include "Utilities.h"


using namespace std;
using namespace LHAPDF;

Analysis::Analysis(TTree & inputTree, TTree & outputTree) 
  : TreeContent(& inputTree) , fInputTree(inputTree), fOutputTree(outputTree)
{
  h_counters = 0;
}

// main event loop
void Analysis::Loop(TString pdfset) 
{
  // process and memory info
  ProcInfo_t info;

  bool data = false;
  bool cteq6ll = false;
  bool CT10 = false;
  if(pdfset == "data")data = true;
  else if(pdfset == "cteq6ll")cteq6ll = true;
  else if(pdfset == "CT10")CT10 = true;
  else{
    cout << "pdfset has to be 'cteq6ll','CT10' or 'data'" << endl;
    return;
  }

  // set branch addresses in output tree if necessary
  SetBranchAddresses();
  CreateHistograms();
  
  TString BadLaserFile= "/user/padeken/TauAna/AllBadHCALLaser.txt.gz";
  HcalLaser = new HCALLaserFilter(BadLaserFile);

  ///######################################################################################################
  ///   Initializing the neccessary stuff for LHAPDF
  ///   !!always!! do:
  ///   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/.automount/net_rw/net__scratch_cms/institut_3a/erdweg/LHAPDF/lib
  ///
  TString PDF_path = "/net/scratch_cms/institut_3a/padeken/lhapdfLib/share/lhapdf/PDFsets/";
  //TString PDF_path = "/.automount/net_rw/net__scratch_cms/institut_3a/erdweg/LHAPDF/share/lhapdf/PDFsets/";

  ///   PDF sets to use, beware that the maximum number of pdf sets to use is 5 at the moment
  TString PDF_1 = "NNPDF22_nlo_100.LHgrid";
  //TString PDF_1 = "MSTW2008lo68cl.LHgrid";
  TString PDF_2 = "MSTW2008nlo68cl.LHgrid";
  TString PDF_3 = "CT10.LHgrid";

  ///   Initialize PDF sets
  TString PDFset_1 = PDF_path+PDF_1;
  TString PDFset_2 = PDF_path+PDF_2;
  TString PDFset_3 = PDF_path+PDF_3;
  
  /// declare variables for saving weights, initialize
  Double_t weights_1[200];
  Double_t weights_2[200];
  Double_t weights_3[200];
  
  Double_t ori;
  TBranch *nweights_1;
  TBranch *bweights_1;
  TBranch *nweights_2;
  TBranch *bweights_2;
  TBranch *nweights_3;
  TBranch *bweights_3;
  if(!data){
    LHAPDF::initPDFSet(1,PDFset_1.Data());
    LHAPDF::initPDFSet(2,PDFset_2.Data());
    LHAPDF::initPDFSet(3,PDFset_3.Data());
    ///   reference PDF set which is used to produce the Monte Carlos
    if(cteq6ll)LHAPDF::initPDFSet(4,(PDF_path+"cteq6ll.LHpdf").Data());
    else if(CT10)LHAPDF::initPDFSet(4,PDFset_3.Data());

    NPDF_1 = 0;
    if (LHAPDF::numberPDF(1) == 1)
      NPDF_1 = 1;
    else
      NPDF_1 = LHAPDF::numberPDF(1) + 1;                    /// usually 0: mean, 1...2*N: N eigenvectors
    if ( NPDF_1 > 200 ) cout << "  Your weights array isn't big enough! Please increase it's size " << PDF_1 << endl;
    cout << PDF_1 << " size: " << NPDF_1 << endl;

    NPDF_2 = 0;
    if (LHAPDF::numberPDF(2) == 1)
      NPDF_2 = 1;
    else
      NPDF_2 = LHAPDF::numberPDF(2) + 1;                    /// usually 0: mean, 1...2*N: N eigenvectors
    if ( NPDF_2 > 200 ) cout << "  Your weights array isn't big enough! Please increase it's size " << PDF_2 << endl;
    cout << PDF_2 << " size: " << NPDF_2 << endl;

    NPDF_3 = 0;
    if (LHAPDF::numberPDF(3) == 1)
      NPDF_3 = 1;
    else
      NPDF_3 = LHAPDF::numberPDF(3) + 1;                    /// usually 0: mean, 1...2*N: N eigenvectors
    if ( NPDF_3 > 200 ) cout << "  Your weights array isn't big enough! Please increase it's size " << PDF_3 << endl;
    cout << PDF_3 << " size: " << NPDF_3 << endl;

    for (int i=0; i<NPDF_1; i++) {
      weights_1[i] = 0.;
    }
    for (int i=0; i<NPDF_2; i++) {
      weights_2[i] = 0.;
    }
    for (int i=0; i<NPDF_3; i++) {
      weights_3[i] = 0.;
    }

    /// new branches with weights for each PDF set
    PDF_1.ReplaceAll(".LHgrid","");
    TString temp = "";
    temp = Form("pdf_weights_%s[%d]/double",PDF_1.Data(),NPDF_1);
    bweights_1 = fOutputTree.Branch("pdf_weights_"+PDF_1, weights_1, temp.Data());
  
    //PDF_2.Remove(PDF_2.Length()-7,7);
    PDF_2.ReplaceAll(".LHgrid","");
    temp = Form("pdf_weights_%s[%d]/double",PDF_2.Data(),NPDF_2);
    bweights_2 = fOutputTree.Branch("pdf_weights_"+PDF_2, weights_2, temp.Data());
  
    //PDF_3.Remove(PDF_3.Length()-7,7);
    PDF_3.ReplaceAll(".LHgrid","");
    temp = Form("pdf_weights_%s[%d]/double",PDF_3.Data(),NPDF_3);
    bweights_3 = fOutputTree.Branch("pdf_weights_"+PDF_3, weights_3, temp.Data());
    ///
    ///   End of Initializing PDF stuff
    ///######################################################################################################
  }
  /// main event loop
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "entries: " << nentries << endl;
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
    if(!HcalLaser->filter(global_run, lumi_section, global_event) && global_isdata){
        continue;
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
    ///   some high-pt Muon ID cuts
    bool rejection = true;
    for(int i = 0; i < muo_n; i++){
      if (muo_n >= 1 && muo_pt[i] > 40. && muo_ID[i][1] == 1 && muo_ID[i][3] == 1 && muo_TrackerLayersMeasTk[i] > 5 && muo_StationsMatched[i] >= 2 && muo_ValidMuonHitsCm[i] >= 1 && muo_ValidPixelHitsCm[i] >= 1 && fabs(muo_d0Tk[i]) < 0.02 && fabs(muo_dzTk[i]) < 0.5 && muo_TevReco_ptError[i][1]/muo_TevReco_pt[i][1] < 0.3) {
        rejection = false;
        break;
      }
    }
    // apply cut
    if (rejection)
      continue;

    if(!data){
      ///######################################################################################################
      ///   Calculating the PDF weights
      ///
      LHAPDF::usePDFMember(4,0);         /// evaluate original PDF set's weight
      ori = LHAPDF::xfx(4, pdf_x1, pdf_scale, pdf_id1) * LHAPDF::xfx(4, pdf_x2, pdf_scale, pdf_id2);

      for (int i=0; i<NPDF_1; i++) {
        LHAPDF::usePDFMember(1,i);
        weights_1[i] =                 /// for every event: calculate a weight for every PDF subset -> array of size NPDF
          LHAPDF::xfx(1, pdf_x1, pdf_scale, pdf_id1) * LHAPDF::xfx(1, pdf_x2, pdf_scale, pdf_id2)/ori;
      }
      bweights_1->Fill();

      for (int i=0; i<NPDF_2; i++) {
        LHAPDF::usePDFMember(2,i);
        weights_2[i] =                 /// for every event: calculate a weight for every PDF subset -> array of size NPDF
          LHAPDF::xfx(2, pdf_x1, pdf_scale, pdf_id1) * LHAPDF::xfx(2, pdf_x2, pdf_scale, pdf_id2)/ori;
      }
      bweights_2->Fill();

      for (int i=0; i<NPDF_3; i++) {
        LHAPDF::usePDFMember(3,i);
        weights_3[i] =               /// for every event: calculate a weight for every PDF subset -> array of size NPDF
          LHAPDF::xfx(3, pdf_x1, pdf_scale, pdf_id1) * LHAPDF::xfx(3, pdf_x2, pdf_scale, pdf_id2)/ori;
      }

      bweights_3->Fill();
      ///
      ///   End of calculating the PDF weights
      ///######################################################################################################
    }
    // fill histograms after applying a cut
    Fill("muo_n", muo_n);
    Fill("muo_pt0", muo_pt[0]);
    Fill("muo_pt1", muo_pt[1]);

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
    if(fInputTree.GetBranchStatus(branch->GetName())){
        LOG(4, "branch " << branch->GetName() << " has been deactivated at adress " << (void *) branch->GetAddress());
    }else{
        fOutputTree.SetBranchAddress(branch->GetName(), branch->GetAddress());
        LOG(4, "branch " << branch->GetName() << " has address " << (void *) branch->GetAddress());
    }
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
