//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 26 16:25:32 2009 by ROOT version 5.18/00a
// from TTree allData/data after cuts
// found on file: out.root
//////////////////////////////////////////////////////////

#ifndef TreeContent_h
#define TreeContent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define INTSIZE 4

class TreeContent {
public :

  union u64 {
    char c[INTSIZE];
    int i;
  };
  
  int get_size(const int* a) {
    int size=0;
    u64 tmp;
    for(int i=0; ; ++i ) {
      tmp.i = a[i];
      size++;
      for(int j=0; j < INTSIZE ; ++j )
	if(tmp.c[j] == '\x0')
	  return size;
    }
    return -1;
  }
  
  std::string unpack(const int* a) {
    int size = get_size(a);
    u64 tmp;
    char* c = new char[size*INTSIZE];
    for(int i=0; i < size ; ++i) {
      tmp.i = a[i];
      for(int j=0; j < INTSIZE; ++j )
	c[i*INTSIZE+j] = tmp.c[j];
    }
    return c;
  }
  
  int* pack(const char* c) {
    u64 tmp;
    int j=0,count=0;
    int size = strlen(c)/INTSIZE+1;
    int* ii=new int[size];
    for(int i=0 ; ; i++) {
      tmp.c[j++]=c[i];
      ii[count]=tmp.i;
      if(j==INTSIZE) {
	j=0;
	count++;
      }
      if(c[i] == '\x0')
	break;
    }
    return ii;
  }


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        global_weight;
   Int_t           global_procID;
   Double_t        global_pthat;
   Int_t           global_store;
   Int_t           global_run;
   Int_t           global_event;
   Int_t           global_bx;
   Int_t           global_lumiblk;
   Int_t           global_orbit;
   Int_t           global_exp;
   Int_t           global_isdata;
   Char_t          global_HLT[100000];
   Int_t           trig_HLTName[50];
   Int_t           trig_n;
   Int_t           trig_prescale[200];
   Int_t           trig_name[200][100];   //[trig_n]
   Int_t           trig_filter[200][100];   //[trig_n]
   Double_t        trig_pt[200];   //[trig_n]
   Double_t        trig_eta[200];   //[trig_n]
   Double_t        trig_phi[200];   //[trig_n]
   Int_t           truth_n;
   Int_t           truth_pdgid[100];   //[truth_n]
   Int_t           truth_bvtxid[100];   //[truth_n]
   Int_t           truth_evtxid[100];   //[truth_n]
   Double_t        truth_E[100];   //[truth_n]
   Double_t        truth_Et[100];   //[truth_n]
   Double_t        truth_p[100];   //[truth_n]
   Double_t        truth_pt[100];   //[truth_n]
   Double_t        truth_px[100];   //[truth_n]
   Double_t        truth_py[100];   //[truth_n]
   Double_t        truth_pz[100];   //[truth_n]
   Double_t        truth_eta[100];   //[truth_n]
   Double_t        truth_phi[100];   //[truth_n]
   Double_t        truth_m[100];   //[truth_n]
   Int_t           truthl_n;
   Int_t           truthl_ori[100];   //[truthl_n]
   Int_t           truthl_pdgid[100];   //[truthl_n]
   Double_t        truthl_E[100];   //[truthl_n]
   Double_t        truthl_Et[100];   //[truthl_n]
   Double_t        truthl_p[100];   //[truthl_n]
   Double_t        truthl_pt[100];   //[truthl_n]
   Double_t        truthl_px[100];   //[truthl_n]
   Double_t        truthl_py[100];   //[truthl_n]
   Double_t        truthl_pz[100];   //[truthl_n]
   Double_t        truthl_eta[100];   //[truthl_n]
   Double_t        truthl_phi[100];   //[truthl_n]
   Int_t           pdf_id1;
   Int_t           pdf_id2;
   Float_t         pdf_x1;
   Float_t         pdf_x2;
   Float_t         pdf_f1;
   Float_t         pdf_f2;
   Float_t         pdf_scale;
   UInt_t          vtx_n;
   UInt_t          vtx_ntr[100];   //[vtx_n]
   Double_t        vtx_x[100];   //[vtx_n]
   Double_t        vtx_y[100];   //[vtx_n]
   Double_t        vtx_z[100];   //[vtx_n]
   Double_t        vtx_chi[100];   //[vtx_n]
   Double_t        met[3];
   Double_t        mex[3];
   Double_t        mey[3];
   Double_t        meteta[3];
   Double_t        metphi[3];
   Double_t        sumet[3];
   Double_t        sumetsig[3];
   UInt_t          jet_n;
   Double_t        jet_E[100];   //[jet_n]
   Double_t        jet_Et[100];   //[jet_n]
   Double_t        jet_p[100];   //[jet_n]
   Double_t        jet_pt[100];   //[jet_n]
   Double_t        jet_px[100];   //[jet_n]
   Double_t        jet_py[100];   //[jet_n]
   Double_t        jet_pz[100];   //[jet_n]
   Double_t        jet_eta[100];   //[jet_n]
   Double_t        jet_phi[100];   //[jet_n]
   Double_t        jet_fem[100];   //[jet_n]
   Int_t           jet_truth[100];   //[jet_n]
   UInt_t          truthjet_n;
   Double_t        truthjet_E[100];   //[truthjet_n]
   Double_t        truthjet_Et[100];   //[truthjet_n]
   Double_t        truthjet_p[100];   //[truthjet_n]
   Double_t        truthjet_pt[100];   //[truthjet_n]
   Double_t        truthjet_px[100];   //[truthjet_n]
   Double_t        truthjet_py[100];   //[truthjet_n]
   Double_t        truthjet_pz[100];   //[truthjet_n]
   Double_t        truthjet_eta[100];   //[truthjet_n]
   Double_t        truthjet_phi[100];   //[truthjet_n]
   UInt_t          ele_n;
   Double_t        ele_E[100];   //[ele_n]
   Double_t        ele_Et[100];   //[ele_n]
   Double_t        ele_p[100];   //[ele_n]
   Double_t        ele_pt[100];   //[ele_n]
   Double_t        ele_px[100];   //[ele_n]
   Double_t        ele_py[100];   //[ele_n]
   Double_t        ele_pz[100];   //[ele_n]
   Double_t        ele_eta[100];   //[ele_n]
   Double_t        ele_phi[100];   //[ele_n]
   Double_t        ele_charge[100];   //[ele_n]
   Double_t        ele_RelTrkIso[100];   //[ele_n]
   Double_t        ele_TrkIso[100];   //[ele_n]
   Double_t        ele_ECalIso[100];   //[ele_n]
   Double_t        ele_HCalIso[100];   //[ele_n]
   Double_t        ele_TrkIsoDep[100];   //[ele_n]
   Double_t        ele_ECalIsoDep[100];   //[ele_n]
   Double_t        ele_HCalIsoDep[100];   //[ele_n]
   Double_t        ele_AllIso[100];   //[ele_n]
   Double_t        ele_TrkChiNorm[100];   //[ele_n]
   Double_t        ele_d0[100];   //[ele_n]
   Double_t        ele_sd0[100];   //[ele_n]
   Int_t           ele_hits[100];   //[ele_n]
   Int_t           ele_truth[100];   //[ele_n]
   UInt_t          ele_ID[100][5];   //[ele_n]
   UInt_t          muo_n;
   Double_t        muo_E[100];   //[muo_n]
   Double_t        muo_Et[100];   //[muo_n]
   Double_t        muo_p[100];   //[muo_n]
   Double_t        muo_pt[100];   //[muo_n]
   Double_t        muo_px[100];   //[muo_n]
   Double_t        muo_py[100];   //[muo_n]
   Double_t        muo_pz[100];   //[muo_n]
   Double_t        muo_eta[100];   //[muo_n]
   Double_t        muo_phi[100];   //[muo_n]
   Double_t        muo_charge[100];   //[muo_n]
   Double_t        muo_RelTrkIso[100];   //[muo_n]
   Double_t        muo_TrkIso[100];   //[muo_n]
   Double_t        muo_ECalIso[100];   //[muo_n]
   Double_t        muo_HCalIso[100];   //[muo_n]
   Double_t        muo_TrkIsoDep[100];   //[muo_n]
   Double_t        muo_ECalIsoDep[100];   //[muo_n]
   Double_t        muo_HCalIsoDep[100];   //[muo_n]
   Double_t        muo_AllIso[100];   //[muo_n]
   Double_t        muo_TrkChiNorm[100];   //[muo_n]
   Double_t        muo_d0[100];   //[muo_n]
   Double_t        muo_sd0[100];   //[muo_n]
   Int_t           muo_prompttight[100];   //[muo_n]
   Int_t           muo_hits[100];   //[muo_n]
   Int_t           muo_truth[100];   //[muo_n]
   UInt_t          muo_trign[100];   //[muo_n]
   UInt_t          muo_trig[100][100];   //[muo_n]

   // List of branches
   TBranch        *b_global_weight;   //!
   TBranch        *b_global_procID;   //!
   TBranch        *b_global_pthat;   //!
   TBranch        *b_global_store;   //!
   TBranch        *b_global_run;   //!
   TBranch        *b_global_event;   //!
   TBranch        *b_global_bx;   //!
   TBranch        *b_global_lumiblk;   //!
   TBranch        *b_global_orbit;   //!
   TBranch        *b_global_exp;   //!
   TBranch        *b_global_isdata;   //!
   TBranch        *b_global_HLT;   //!
   TBranch        *b_trig_HLTName;   //!
   TBranch        *b_trig_n;   //!
   TBranch        *b_trig_prescale;   //!
   TBranch        *b_trig_name;   //!
   TBranch        *b_trig_filter;   //!
   TBranch        *b_trig_pt;   //!
   TBranch        *b_trig_eta;   //!
   TBranch        *b_trig_phi;   //!
   TBranch        *b_truth_n;   //!
   TBranch        *b_truth_pdgid;   //!
   TBranch        *b_truth_bvtxid;   //!
   TBranch        *b_truth_evtxid;   //!
   TBranch        *b_truth_E;   //!
   TBranch        *b_truth_Et;   //!
   TBranch        *b_truth_p;   //!
   TBranch        *b_truth_pt;   //!
   TBranch        *b_truth_px;   //!
   TBranch        *b_truth_py;   //!
   TBranch        *b_truth_pz;   //!
   TBranch        *b_truth_eta;   //!
   TBranch        *b_truth_phi;   //!
   TBranch        *b_truth_m;   //!
   TBranch        *b_truthl_n;   //!
   TBranch        *b_truthl_ori;   //!
   TBranch        *b_truthl_pdgid;   //!
   TBranch        *b_truthl_E;   //!
   TBranch        *b_truthl_Et;   //!
   TBranch        *b_truthl_p;   //!
   TBranch        *b_truthl_pt;   //!
   TBranch        *b_truthl_px;   //!
   TBranch        *b_truthl_py;   //!
   TBranch        *b_truthl_pz;   //!
   TBranch        *b_truthl_eta;   //!
   TBranch        *b_truthl_phi;   //!
   TBranch        *b_pdf_id1;   //!
   TBranch        *b_pdf_id2;   //!
   TBranch        *b_pdf_x1;   //!
   TBranch        *b_pdf_x2;   //!
   TBranch        *b_pdf_f1;   //!
   TBranch        *b_pdf_f2;   //!
   TBranch        *b_pdf_scale;   //!
   TBranch        *b_vtx_n;   //!
   TBranch        *b_vtx_ntr;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_chi;   //!
   TBranch        *b_met;   //!
   TBranch        *b_mex;   //!
   TBranch        *b_mey;   //!
   TBranch        *b_meteta;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_sumet;   //!
   TBranch        *b_sumetsig;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_Et;   //!
   TBranch        *b_jet_p;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_fem;   //!
   TBranch        *b_jet_truth;   //!
   TBranch        *b_truthjet_n;   //!
   TBranch        *b_truthjet_E;   //!
   TBranch        *b_truthjet_Et;   //!
   TBranch        *b_truthjet_p;   //!
   TBranch        *b_truthjet_pt;   //!
   TBranch        *b_truthjet_px;   //!
   TBranch        *b_truthjet_py;   //!
   TBranch        *b_truthjet_pz;   //!
   TBranch        *b_truthjet_eta;   //!
   TBranch        *b_truthjet_phi;   //!
   TBranch        *b_ele_n;   //!
   TBranch        *b_ele_E;   //!
   TBranch        *b_ele_Et;   //!
   TBranch        *b_ele_p;   //!
   TBranch        *b_ele_pt;   //!
   TBranch        *b_ele_px;   //!
   TBranch        *b_ele_py;   //!
   TBranch        *b_ele_pz;   //!
   TBranch        *b_ele_eta;   //!
   TBranch        *b_ele_phi;   //!
   TBranch        *b_ele_charge;   //!
   TBranch        *b_ele_TrkIso;   //!
   TBranch        *b_ele_RelTrkIso;   //!
   TBranch        *b_ele_ECalIso;   //!
   TBranch        *b_ele_HCalIso;   //!
   TBranch        *b_ele_TrkIsoDep;   //!
   TBranch        *b_ele_ECalIsoDep;   //!
   TBranch        *b_ele_HCalIsoDep;   //!
   TBranch        *b_ele_AllIso;   //!
   TBranch        *b_ele_TrkChiNorm;   //!
   TBranch        *b_ele_d0;   //!
   TBranch        *b_ele_sd0;   //!
   TBranch        *b_ele_hits;   //!
   TBranch        *b_ele_truth;   //!
   TBranch        *b_ele_ID;   //!
   TBranch        *b_muo_n;   //!
   TBranch        *b_muo_E;   //!
   TBranch        *b_muo_Et;   //!
   TBranch        *b_muo_p;   //!
   TBranch        *b_muo_pt;   //!
   TBranch        *b_muo_px;   //!
   TBranch        *b_muo_py;   //!
   TBranch        *b_muo_pz;   //!
   TBranch        *b_muo_eta;   //!
   TBranch        *b_muo_phi;   //!
   TBranch        *b_muo_charge;   //!
   TBranch        *b_muo_TrkIso;   //!
   TBranch        *b_muo_RelTrkIso;   //!
   TBranch        *b_muo_ECalIso;   //!
   TBranch        *b_muo_HCalIso;   //!
   TBranch        *b_muo_TrkIsoDep;   //!
   TBranch        *b_muo_ECalIsoDep;   //!
   TBranch        *b_muo_HCalIsoDep;   //!
   TBranch        *b_muo_AllIso;   //!
   TBranch        *b_muo_TrkChiNorm;   //!
   TBranch        *b_muo_d0;   //!
   TBranch        *b_muo_sd0;   //!
   TBranch        *b_muo_prompttight;   //!
   TBranch        *b_muo_hits;   //!
   TBranch        *b_muo_truth;   //!
   TBranch        *b_muo_trign;   //!
   TBranch        *b_muo_trig;   //!

   TreeContent(TTree *tree=0);
   virtual ~TreeContent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeContent_cxx
TreeContent::TreeContent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out.root");
      if (!f) {
         f = new TFile("out.root");
         f->cd("out.root:/ACSkimAnalysis");
      }
      tree = (TTree*)gDirectory->Get("allData");

   }
   Init(tree);
}

TreeContent::~TreeContent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeContent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeContent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeContent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("global_weight", &global_weight, &b_global_weight);
   fChain->SetBranchAddress("global_procID", &global_procID, &b_global_procID);
   fChain->SetBranchAddress("global_pthat", &global_pthat, &b_global_pthat);
   fChain->SetBranchAddress("global_store", &global_store, &b_global_store);
   fChain->SetBranchAddress("global_run", &global_run, &b_global_run);
   fChain->SetBranchAddress("global_event", &global_event, &b_global_event);
   fChain->SetBranchAddress("global_bx", &global_bx, &b_global_bx);
   fChain->SetBranchAddress("global_lumiblk", &global_lumiblk, &b_global_lumiblk);
   fChain->SetBranchAddress("global_orbit", &global_orbit, &b_global_orbit);
   fChain->SetBranchAddress("global_exp", &global_exp, &b_global_exp);
   fChain->SetBranchAddress("global_isdata", &global_isdata, &b_global_isdata);
   fChain->SetBranchAddress("global_HLT", global_HLT, &b_global_HLT);
   fChain->SetBranchAddress("trig_HLTName", trig_HLTName, &b_trig_HLTName);
   fChain->SetBranchAddress("trig_n", &trig_n, &b_trig_n);
   fChain->SetBranchAddress("trig_prescale", trig_prescale, &b_trig_prescale);
   fChain->SetBranchAddress("trig_name", trig_name, &b_trig_name);
   fChain->SetBranchAddress("trig_filter", trig_filter, &b_trig_filter);
   fChain->SetBranchAddress("trig_pt", trig_pt, &b_trig_pt);
   fChain->SetBranchAddress("trig_eta", trig_eta, &b_trig_eta);
   fChain->SetBranchAddress("trig_phi", trig_phi, &b_trig_phi);
   fChain->SetBranchAddress("truth_n", &truth_n, &b_truth_n);
   fChain->SetBranchAddress("truth_pdgid", truth_pdgid, &b_truth_pdgid);
   fChain->SetBranchAddress("truth_bvtxid", truth_bvtxid, &b_truth_bvtxid);
   fChain->SetBranchAddress("truth_evtxid", truth_evtxid, &b_truth_evtxid);
   fChain->SetBranchAddress("truth_E", truth_E, &b_truth_E);
   fChain->SetBranchAddress("truth_Et", truth_Et, &b_truth_Et);
   fChain->SetBranchAddress("truth_p", truth_p, &b_truth_p);
   fChain->SetBranchAddress("truth_pt", truth_pt, &b_truth_pt);
   fChain->SetBranchAddress("truth_px", truth_px, &b_truth_px);
   fChain->SetBranchAddress("truth_py", truth_py, &b_truth_py);
   fChain->SetBranchAddress("truth_pz", truth_pz, &b_truth_pz);
   fChain->SetBranchAddress("truth_eta", truth_eta, &b_truth_eta);
   fChain->SetBranchAddress("truth_phi", truth_phi, &b_truth_phi);
   fChain->SetBranchAddress("truth_m", truth_m, &b_truth_m);
   fChain->SetBranchAddress("truthl_n", &truthl_n, &b_truthl_n);
   fChain->SetBranchAddress("truthl_ori", truthl_ori, &b_truthl_ori);
   fChain->SetBranchAddress("truthl_pdgid", truthl_pdgid, &b_truthl_pdgid);
   fChain->SetBranchAddress("truthl_E", truthl_E, &b_truthl_E);
   fChain->SetBranchAddress("truthl_Et", truthl_Et, &b_truthl_Et);
   fChain->SetBranchAddress("truthl_p", truthl_p, &b_truthl_p);
   fChain->SetBranchAddress("truthl_pt", truthl_pt, &b_truthl_pt);
   fChain->SetBranchAddress("truthl_px", truthl_px, &b_truthl_px);
   fChain->SetBranchAddress("truthl_py", truthl_py, &b_truthl_py);
   fChain->SetBranchAddress("truthl_pz", truthl_pz, &b_truthl_pz);
   fChain->SetBranchAddress("truthl_eta", truthl_eta, &b_truthl_eta);
   fChain->SetBranchAddress("truthl_phi", truthl_phi, &b_truthl_phi);
   fChain->SetBranchAddress("pdf_id1", &pdf_id1, &b_pdf_id1);
   fChain->SetBranchAddress("pdf_id2", &pdf_id2, &b_pdf_id2);
   fChain->SetBranchAddress("pdf_x1", &pdf_x1, &b_pdf_x1);
   fChain->SetBranchAddress("pdf_x2", &pdf_x2, &b_pdf_x2);
   fChain->SetBranchAddress("pdf_f1", &pdf_f1, &b_pdf_f1);
   fChain->SetBranchAddress("pdf_f2", &pdf_f2, &b_pdf_f2);
   fChain->SetBranchAddress("pdf_scale", &pdf_scale, &b_pdf_scale);
   fChain->SetBranchAddress("vtx_n", &vtx_n, &b_vtx_n);
   fChain->SetBranchAddress("vtx_ntr", vtx_ntr, &b_vtx_ntr);
   fChain->SetBranchAddress("vtx_x", vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_chi", vtx_chi, &b_vtx_chi);
   fChain->SetBranchAddress("met", met, &b_met);
   fChain->SetBranchAddress("mex", mex, &b_mex);
   fChain->SetBranchAddress("mey", mey, &b_mey);
   fChain->SetBranchAddress("meteta", meteta, &b_meteta);
   fChain->SetBranchAddress("metphi", metphi, &b_metphi);
   fChain->SetBranchAddress("sumet", sumet, &b_sumet);
   fChain->SetBranchAddress("sumetsig", sumetsig, &b_sumetsig);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_Et", jet_Et, &b_jet_Et);
   fChain->SetBranchAddress("jet_p", jet_p, &b_jet_p);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_px", jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_fem", jet_fem, &b_jet_fem);
   fChain->SetBranchAddress("jet_truth", jet_truth, &b_jet_truth);
   fChain->SetBranchAddress("truthjet_n", &truthjet_n, &b_truthjet_n);
   fChain->SetBranchAddress("truthjet_E", truthjet_E, &b_truthjet_E);
   fChain->SetBranchAddress("truthjet_Et", truthjet_Et, &b_truthjet_Et);
   fChain->SetBranchAddress("truthjet_p", truthjet_p, &b_truthjet_p);
   fChain->SetBranchAddress("truthjet_pt", truthjet_pt, &b_truthjet_pt);
   fChain->SetBranchAddress("truthjet_px", truthjet_px, &b_truthjet_px);
   fChain->SetBranchAddress("truthjet_py", truthjet_py, &b_truthjet_py);
   fChain->SetBranchAddress("truthjet_pz", truthjet_pz, &b_truthjet_pz);
   fChain->SetBranchAddress("truthjet_eta", truthjet_eta, &b_truthjet_eta);
   fChain->SetBranchAddress("truthjet_phi", truthjet_phi, &b_truthjet_phi);
   fChain->SetBranchAddress("ele_n", &ele_n, &b_ele_n);
   fChain->SetBranchAddress("ele_E", ele_E, &b_ele_E);
   fChain->SetBranchAddress("ele_Et", ele_Et, &b_ele_Et);
   fChain->SetBranchAddress("ele_p", ele_p, &b_ele_p);
   fChain->SetBranchAddress("ele_pt", ele_pt, &b_ele_pt);
   fChain->SetBranchAddress("ele_px", ele_px, &b_ele_px);
   fChain->SetBranchAddress("ele_py", ele_py, &b_ele_py);
   fChain->SetBranchAddress("ele_pz", ele_pz, &b_ele_pz);
   fChain->SetBranchAddress("ele_eta", ele_eta, &b_ele_eta);
   fChain->SetBranchAddress("ele_phi", ele_phi, &b_ele_phi);
   fChain->SetBranchAddress("ele_charge", ele_charge, &b_ele_charge);
   fChain->SetBranchAddress("ele_RelTrkIso", ele_RelTrkIso, &b_ele_RelTrkIso);
   fChain->SetBranchAddress("ele_TrkIso", ele_TrkIso, &b_ele_TrkIso);
   fChain->SetBranchAddress("ele_ECalIso", ele_ECalIso, &b_ele_ECalIso);
   fChain->SetBranchAddress("ele_HCalIso", ele_HCalIso, &b_ele_HCalIso);
   fChain->SetBranchAddress("ele_TrkIsoDep", ele_TrkIsoDep, &b_ele_TrkIsoDep);
   fChain->SetBranchAddress("ele_ECalIsoDep", ele_ECalIsoDep, &b_ele_ECalIsoDep);
   fChain->SetBranchAddress("ele_HCalIsoDep", ele_HCalIsoDep, &b_ele_HCalIsoDep);
   fChain->SetBranchAddress("ele_AllIso", ele_AllIso, &b_ele_AllIso);
   fChain->SetBranchAddress("ele_TrkChiNorm", ele_TrkChiNorm, &b_ele_TrkChiNorm);
   fChain->SetBranchAddress("ele_d0", ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_sd0", ele_sd0, &b_ele_sd0);
   fChain->SetBranchAddress("ele_hits", ele_hits, &b_ele_hits);
   fChain->SetBranchAddress("ele_truth", ele_truth, &b_ele_truth);
   fChain->SetBranchAddress("ele_ID", ele_ID, &b_ele_ID);
   fChain->SetBranchAddress("muo_n", &muo_n, &b_muo_n);
   fChain->SetBranchAddress("muo_E", muo_E, &b_muo_E);
   fChain->SetBranchAddress("muo_Et", muo_Et, &b_muo_Et);
   fChain->SetBranchAddress("muo_p", muo_p, &b_muo_p);
   fChain->SetBranchAddress("muo_pt", muo_pt, &b_muo_pt);
   fChain->SetBranchAddress("muo_px", muo_px, &b_muo_px);
   fChain->SetBranchAddress("muo_py", muo_py, &b_muo_py);
   fChain->SetBranchAddress("muo_pz", muo_pz, &b_muo_pz);
   fChain->SetBranchAddress("muo_eta", muo_eta, &b_muo_eta);
   fChain->SetBranchAddress("muo_phi", muo_phi, &b_muo_phi);
   fChain->SetBranchAddress("muo_charge", muo_charge, &b_muo_charge);
   fChain->SetBranchAddress("muo_RelTrkIso", muo_RelTrkIso, &b_muo_RelTrkIso);
   fChain->SetBranchAddress("muo_TrkIso", muo_TrkIso, &b_muo_TrkIso);
   fChain->SetBranchAddress("muo_ECalIso", muo_ECalIso, &b_muo_ECalIso);
   fChain->SetBranchAddress("muo_HCalIso", muo_HCalIso, &b_muo_HCalIso);
   fChain->SetBranchAddress("muo_TrkIsoDep", muo_TrkIsoDep, &b_muo_TrkIsoDep);
   fChain->SetBranchAddress("muo_ECalIsoDep", muo_ECalIsoDep, &b_muo_ECalIsoDep);
   fChain->SetBranchAddress("muo_HCalIsoDep", muo_HCalIsoDep, &b_muo_HCalIsoDep);
   fChain->SetBranchAddress("muo_AllIso", muo_AllIso, &b_muo_AllIso);
   fChain->SetBranchAddress("muo_TrkChiNorm", muo_TrkChiNorm, &b_muo_TrkChiNorm);
   fChain->SetBranchAddress("muo_d0", muo_d0, &b_muo_d0);
   fChain->SetBranchAddress("muo_sd0", muo_sd0, &b_muo_sd0);
   fChain->SetBranchAddress("muo_prompttight", muo_prompttight, &b_muo_prompttight);
   fChain->SetBranchAddress("muo_hits", muo_hits, &b_muo_hits);
   fChain->SetBranchAddress("muo_truth", muo_truth, &b_muo_truth);
   fChain->SetBranchAddress("muo_trig", muo_trig, &b_muo_trig);
   fChain->SetBranchAddress("muo_trign", muo_trign, &b_muo_trign);
   Notify();
}

Bool_t TreeContent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeContent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeContent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeContent_cxx
