//
// Package:    UserCode/aachen3a/ACSusyAnalysis
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
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "aachen3a/ACSusyAnalysis/interface/TriggerTools.h"

using namespace std;
using namespace pat;
using namespace ACSusyAnalysis;

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
  virtual void beginJob();
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

  double DeltaPhi(double a, double b);

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
  edm::InputTag genJetTag_;
  edm::InputTag vertexTag_;

  bool is_MC;
  bool is_SHERPA;
  bool do_fatjets;

  std::string cor_;
  std::string flav_;
  std::string btag_;
  pat::JetCorrFactors::CorrStep correction_;

  GreaterByPt<pat::Muon>      ptcomp_muo;
  GreaterByPt<pat::Electron>  ptcomp_ele;
  GreaterByPt<pat::Jet>       ptcomp_jet;
  GreaterByPt<reco::GenJet>   ptcomp_genjet;

  typedef std::pair<std::string,float> IdPair;

  double _jeta[100];
  double _jphi[100];

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

  double pthat_low_;
  double pthat_high_;

  // Counters
  unsigned int nrEventTotalRaw_;
  unsigned int nrEventPassedPthatRaw_;
  unsigned int nrEventPassedRaw_;

  double localPi;
  
  // Tree variables
  int mTreerun;
  int mTreeevent;
  int mTreestore;
  int mTreebx;
  int mTreeorbit;
  int mTreeexp;
  int mTreedata;

  int mTreelumiblk;
  double mTreedellumi;
  double mTreereclumi;
  double mTreedellumierr;
  double mTreereclumierr;

  int mTreenoisel;
  int mTreenoiset;
  int mTreenoiseh;

  double mTreeecalr9;
  double mTreeecale;
  double mTreeecalpt;
  double mTreeecalpx;
  double mTreeecalpy;
  double mTreeecalpz;
  double mTreeecaleta;
  double mTreeecalphi;
  double mTreeecaltime;
  double mTreeecalchi;
  int    mTreeecalieta;
  int    mTreeecaliphi;
  int    mTreeecalflag;

  Char_t mTreeHLT[100000];

  int mTreetrighltname[50];
  
  int mTreeNtrig;
  int mTreetrigpre[500];
  int mTreetrigname[500][100];
  int mTreefiltname[500][100];
  double mTreetrigpt[500];
  double mTreetrigeta[500];
  double mTreetrigphi[500];
  
  double mTreeMET[3];
  double mTreeMEX[3];
  double mTreeMEY[3];
  double mTreeMETphi[3];
  double mTreeMETeta[3];
  double mTreeSumET[3];
  double mTreeSumETSignif[3];

  int    mTreeNtracks;
  double mTreetrackshqf;

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
  int mTreetruthlpdgid[100];
  int mTreetruthlori[100];
  double mTreetruthlE[100];
  double mTreetruthlEt[100];
  double mTreetruthlp[100];
  double mTreetruthlpt[100];
  double mTreetruthlpx[100];
  double mTreetruthlpy[100];
  double mTreetruthlpz[100];
  double mTreetruthleta[100];
  double mTreetruthlphi[100];

  int mTreepdfid1;
  int mTreepdfid2;
  float mTreepdfx1;
  float mTreepdfx2;
  float mTreepdff1;
  float mTreepdff2;
  float mTreepdfscale;

  int    mTreeNjet;
  int    mTreeJetTruth[100];
  int    mTreeJetPart[100];
  double mTreeJetEt[100];
  double mTreeJetPt[100];
  double mTreeJetP[100];
  double mTreeJetPx[100];
  double mTreeJetPy[100];
  double mTreeJetPz[100];
  double mTreeJetE[100];
  double mTreeJetEta[100];
  double mTreeJetPhi[100];
  double mTreeJetFem[100];
  double mTreeJetFhad[100];
  double mTreeJetBtag[100];
  double mTreeJetCharge[100];
  double mTreeJetn90[100];
  double mTreeJetfhpd[100];
  double mTreeJetfrbx[100];

  int    mTreeNtruthjet;
  double mTreetruthJetEt[100];
  double mTreetruthJetPt[100];
  double mTreetruthJetP[100];
  double mTreetruthJetPx[100];
  double mTreetruthJetPy[100];
  double mTreetruthJetPz[100];
  double mTreetruthJetE[100];
  double mTreetruthJetEta[100];
  double mTreetruthJetPhi[100];

  int    mTreeNfatjet;
  int    mTreefatjetnsub[100];
  double mTreefatjetpt[100];
  double mTreefatjetpx[100];
  double mTreefatjetpy[100];
  double mTreefatjetpz[100];
  double mTreefatjete[100];
  double mTreefatjeteta[100];
  double mTreefatjetphi[100];
  double mTreefatjetsubpt[100][10];
  double mTreefatjetsubpx[100][10];
  double mTreefatjetsubpy[100][10];
  double mTreefatjetsubpz[100][10];
  double mTreefatjetsube[100][10];
  double mTreefatjetsubeta[100][10];
  double mTreefatjetsubphi[100][10];
  double mTreefatjetsubfem[100][10];
  double mTreefatjetsubfhad[100][10];
  double mTreefatjetsubbtag[100][10];
  double mTreefatjetsubn90[100][10];
  double mTreefatjetsubfhpd[100][10];
  double mTreefatjetsubfrbx[100][10];

  int    mTreeNele;
  int    mTreeEleID[100][5];
  int    mTreeEleTruth[100];
  int    mTreeEleHits[100];
  double mTreeEleEt[100];
  double mTreeEleP[100];
  double mTreeElePt[100];
  double mTreeElePx[100];
  double mTreeElePy[100];
  double mTreeElePz[100];
  double mTreeEleE[100];
  double mTreeEleEta[100];
  double mTreeElePhi[100];
  double mTreeEleTrkIso[100];
  double mTreeEleRelTrkIso[100];
  double mTreeEleECalIso[100];
  double mTreeEleHCalIso[100];
  double mTreeEleAllIso[100];
  double mTreeEleECalIsoDep[100];
  double mTreeEleHCalIsoDep[100];
  double mTreeEleTrkIsoDep[100];
  double mTreeEleTrkChiNorm[100];
  double mTreeEleCharge[100];
  double mTreeEled0[100];
  double mTreeElesd0[100];

  int    mTreeNmuo;
  int    mTreeNmuotrign[100];
  int    mTreeMuotrig[100][100];
  int    mTreeMuoTruth[100];
  int    mTreeMuoHitsCm[100];
  int    mTreeMuoHitsTk[100];
  int    mTreeMuoGood[100];
  double mTreeMuoEt[100];
  double mTreeMuoP[100];
  double mTreeMuoPt[100];
  double mTreeMuoPx[100];
  double mTreeMuoPy[100];
  double mTreeMuoPz[100];
  double mTreeMuoE[100];
  double mTreeMuoEta[100];
  double mTreeMuoPhi[100];
  double mTreeMuoTrkIso[100];
  double mTreeMuoRelTrkIso[100];
  double mTreeMuoECalIso[100];
  double mTreeMuoHCalIso[100];
  double mTreeMuoAllIso[100];
  double mTreeMuoECalIsoDep[100];
  double mTreeMuoHCalIsoDep[100];
  double mTreeMuoTrkIsoDep[100];
  double mTreeMuoTrkChiNormCm[100];
  double mTreeMuoTrkChiNormTk[100];
  double mTreeMuoCharge[100];
  double mTreeMuod0Cm[100];
  double mTreeMuod0Tk[100];
  double mTreeMuosd0Cm[100];
  double mTreeMuosd0Tk[100];
  double mTreeMuocalocomp[100];
  double mTreeMuocaltowe[100];

  int    mTreeNvtx;
  int    mTreeVtxntr[100];
  int    mTreeVtxfake[100];
  double mTreeVtxndf[100];
  double mTreeVtxx[100];
  double mTreeVtxy[100];
  double mTreeVtxz[100];
  double mTreeVtxchi[100];

  double mTreeEventWeight;
  int    mTreeProcID;
  double mTreePthat;

};

