//
// Package:    SusyACSkimAnalysis
// Class:      SusyACSkimAnalysis
// 
// Description: Skeleton analysis for SUSY search with Lepton + Jets + MET
//
// Original Author:  Carsten Magass
//         Created:  November 2008
//

// System include files
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

// ROOT includes
#include <TNtuple.h>
#include <TH1F.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// SUSY include files
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"
#include "SusyAnalysis/EventSelector/interface/uncorrectionTypeMET.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

#include "PhysicsTools/Utilities/interface/PtComparator.h"

using namespace std;
using namespace pat;

////////////////////////////////
//
// Class declaration
//
class SusyACSkimAnalysis : public edm::EDFilter {
public:
  explicit SusyACSkimAnalysis(const edm::ParameterSet&);
  ~SusyACSkimAnalysis();
  
private:
  //*** CMSSW interface
  /// Called once per job, at start
  virtual void beginJob(const edm::EventSetup&) ;
  /// Called for each event
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  /// Called once per job, at end
  virtual void endJob();
  
  /// Helpers
  virtual bool isStable(int pdgid);
  virtual bool isDecaying(int pdgid);
  virtual bool isSUSY(int pdgid);

  /// Print a summary of counts for all selectors
  virtual void printSummary(void);

  //*** Plotting
  /// Define all plots
  virtual void initPlots();

private:

  // Plots
  TTree * mAllData; // Will contain the SUSY-AC specific data

  TH1F* h_counters;

  // Data tags
  edm::InputTag jetTag_;
  edm::InputTag metTag_;
  edm::InputTag elecTag_;
  edm::InputTag muonTag_;
  edm::InputTag genTag_;
  edm::InputTag trigTag_;
  edm::InputTag genJetTag_;
  edm::InputTag vertexTag_;

  std::string gen;

  pat::JetCorrFactors::CorrStep correction_;

  GreaterByPt<pat::Muon>      ptcomp_muo;
  GreaterByPt<pat::Electron>  ptcomp_ele;
  GreaterByPt<pat::Jet>       ptcomp_jet;
  GreaterByPt<reco::GenJet>   ptcomp_genjet;

  typedef std::pair<std::string,float> IdPair;

  int nele_;
  int nmuo_;
  int njet_;
  double muopt_;
  double muoeta_;
  double elept_;
  double eleeta_;
  double jetpt_;
  double jetfem_;
  double jeteta_;
  double met_;

  // Counters
  unsigned int nrEventTotalRaw_;
  unsigned int nrEventPassedRaw_;

  double localPi;

  // Tree variables
  int mTreerun;
  int mTreeevent;
  int mTreelumiblk;
  int mTreestore;
  int mTreebx;
  int mTreeorbit;
  int mTreeexp;
  int mTreedata;

  Char_t mTreeHLT[10000];

  double mTreeMET[3];
  double mTreeMEX[3];
  double mTreeMEY[3];
  double mTreeMETphi[3];
  double mTreeMETeta[3];
  double mTreeSumET[3];
  double mTreeSumETSignif[3];

  int mTreeNtruth;
  int mTreetruthpdgid[100];
  int mTreetruthbvtxid[100];
  int mTreetruthevtxid[100];
  double mTreetruthE[100];
  double mTreetruthEt[100];
  double mTreetruthp[100];
  double mTreetruthpt[100];
  double mTreetruthpx[100];
  double mTreetruthpy[100];
  double mTreetruthpz[100];
  double mTreetrutheta[100];
  double mTreetruthphi[100];
  double mTreetruthm[100];

  int mTreeNtruthl;
  int mTreetruthlpdgid[50];
  int mTreetruthlori[50];
  double mTreetruthlE[50];
  double mTreetruthlEt[50];
  double mTreetruthlp[50];
  double mTreetruthlpt[50];
  double mTreetruthlpx[50];
  double mTreetruthlpy[50];
  double mTreetruthlpz[50];
  double mTreetruthleta[50];
  double mTreetruthlphi[50];

  int mTreepdfid1;
  int mTreepdfid2;
  float mTreepdfx1;
  float mTreepdfx2;
  float mTreepdff1;
  float mTreepdff2;
  float mTreepdfscale;

  int    mTreeNjet;
  int    mTreeJetTruth[50];
  double mTreeJetEt[50];
  double mTreeJetPt[50];
  double mTreeJetP[50];
  double mTreeJetPx[50];
  double mTreeJetPy[50];
  double mTreeJetPz[50];
  double mTreeJetE[50];
  double mTreeJetEta[50];
  double mTreeJetPhi[50];
  double mTreeJetFem[50];

  int    mTreeNtruthjet;
  double mTreetruthJetEt[50];
  double mTreetruthJetPt[50];
  double mTreetruthJetP[50];
  double mTreetruthJetPx[50];
  double mTreetruthJetPy[50];
  double mTreetruthJetPz[50];
  double mTreetruthJetE[50];
  double mTreetruthJetEta[50];
  double mTreetruthJetPhi[50];

  int    mTreeNele;
  int    mTreeEleID[50][5];
  int    mTreeEleTruth[50];
  double mTreeEleEt[50];
  double mTreeEleP[50];
  double mTreeElePt[50];
  double mTreeElePx[50];
  double mTreeElePy[50];
  double mTreeElePz[50];
  double mTreeEleE[50];
  double mTreeEleEta[50];
  double mTreeElePhi[50];
  double mTreeEleTrkIso[50];
  double mTreeEleRelTrkIso[50];
  double mTreeEleECalIso[50];
  double mTreeEleHCalIso[50];
  double mTreeEleAllIso[50];
  double mTreeEleTrkChiNorm[50];
  double mTreeEleCharge[50];

  int    mTreeNmuo;
  int    mTreeMuoTruth[50];
  int    mTreeMuoHits[50];
  double mTreeMuoEt[50];
  double mTreeMuoP[50];
  double mTreeMuoPt[50];
  double mTreeMuoPx[50];
  double mTreeMuoPy[50];
  double mTreeMuoPz[50];
  double mTreeMuoE[50];
  double mTreeMuoEta[50];
  double mTreeMuoPhi[50];
  double mTreeMuoTrkIso[50];
  double mTreeMuoRelTrkIso[50];
  double mTreeMuoECalIso[50];
  double mTreeMuoHCalIso[50];
  double mTreeMuoAllIso[50];
  double mTreeMuoTrkChiNorm[50];
  double mTreeMuoCharge[50];
  double mTreeMuod0[50];
  double mTreeMuosd0[50];

  int    mTreeNvtx;
  int    mTreeVtxntr[50];
  double mTreeVtxx[50];
  double mTreeVtxy[50];
  double mTreeVtxz[50];
  double mTreeVtxchi[50];

  double mTreeEventWeight;
  int    mTreeProcID;
  double mTreePthat;

};

