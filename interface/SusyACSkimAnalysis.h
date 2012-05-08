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
#include <algorithm>

// ROOT includes
#include <TNtuple.h>
#include <TH1F.h>
#include <TVector3.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

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
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"

#include "aachen3a/ACSusyAnalysis/interface/TriggerTools.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"


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
  /// Called once per run, at start
  virtual void beginMyRun(const edm::Run&, const edm::EventSetup&);
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
  double DeltaR(double a, double b, double c, double d);
  //~ bool SortEleRefs(pat::ElectronRef& a, pat::ElectronRef& b); 
  //~ bool SortInt(int a, int b); 
  

private:

  // Plots
  TTree * mAllData; // Will contain the SUSY-AC specific data

  TH1F* h_counters;

  HLTConfigProvider hltConfig_;
  bool hltConfigInit_;
  L1GtUtils l1GtUtils_;
  int l1error;

	//structs for getting filter informations inside the Skimmer. Shamelessly stolen from MUSiC

   struct FilterDef {
      std::string name;
      unsigned int ID;
      bool active;
   };


   struct FilterResult {
      std::string   name;
      std::string   process;
      edm::InputTag results;
      edm::InputTag event;
      HLTConfigProvider config;
      std::vector< std::string > filter_names;
      std::vector< FilterDef > filter_infos;
   };
   std::vector< FilterResult > filters;

  // Data tags
  edm::InputTag calojetTag_;
  edm::InputTag pfjetTag_;
  edm::InputTag metTag_;
  edm::InputTag metTagPF_;
  edm::InputTag metTagTC_;
  edm::InputTag photonTag_;
  edm::InputTag elecTag_;
  edm::InputTag gsfelecTag_;
  edm::InputTag PFelecTag_;
  edm::InputTag muonTag_;
  edm::InputTag PFmuonTag_;
  edm::InputTag tauSrc_;
  edm::InputTag metTagPFnoPU_;
  edm::InputTag metTagJPFnoPUType1_ ;
  edm::InputTag metTagcorMetGlobalMuons_;
  edm::InputTag metTagHO_;
  edm::InputTag metTagNoHF_;
  edm::InputTag genTag_;
  edm::InputTag genJetTag_;
  edm::InputTag vertexTag_;
  edm::InputTag ebhitsTag_;
  edm::InputTag freducedBarrelRecHitCollection_;
  edm::InputTag freducedEndcapRecHitCollection_;
  edm::InputTag IsoDepElectron;
  edm::InputTag IsoValElectronPF;
  edm::InputTag IsoDepPhoton;
  edm::InputTag IsoValPhotonPF;
  
  
  std::vector<edm::InputTag> inputTagIsoDepElectrons_;
  std::vector<edm::InputTag> inputTagIsoDepPhotons_;
  //  std::vector<edm::InputTag> inputTagIsoValElectronsNoPFId_;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;   
  //~ std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;   
  std::vector<edm::InputTag> inputTagIsoValPhotonsPFId_;     

  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;


  bool is_PYTHIA8;
  bool is_MC;
  bool is_SHERPA;
  bool do_fatjets;
  bool matchAll_;
  bool susyPar_;
  bool doCaloJet_;
  bool doTaus_;
  //bool beamScaping_;

  std::string btag_;
  std::string processName_; 

  GreaterByPt<pat::Muon>      ptcomp_muo;
  GreaterByPt<pat::Tau>       ptcomp_tau;
  GreaterByPt<pat::Photon>    ptcomp_photon;
  GreaterByPt<pat::Electron>  ptcomp_ele;
  GreaterByPt<pat::Jet>       ptcomp_jet;
  GreaterByPt<reco::GenJet>   ptcomp_genjet;
  GreaterByPtRef<pat::ElectronRef>   ptcomp_EleRef;
  GreaterByPtRef<pat::PhotonRef>   ptcomp_PhotonRef;

  typedef std::pair<std::string,float> IdPair;

  JetIDSelectionFunctor jetId;
  std::string vers_;
  std::string qual_;
  JetIDSelectionFunctor::Version_t vers;
  JetIDSelectionFunctor::Quality_t qual;
  
  //needed for Combined Isolation
  EgammaTowerIsolation *hadDepth1Isolation03_, *hadDepth2Isolation03_;
  typedef std::vector< edm::Handle<edm::ValueMap<double> > > isoContainer;
  isoContainer *eIsoFromDepsValueMap_;
  isoContainer *muIsoFromDepsValueMap_;

  TString ACmuonID[24];
  TString ACtauID[18];

  double _jeta[100];
  double _jphi[100];

  double bfield;

  int npho_;
  int nele_;
  int npfele_;
  int nmuo_;
  int ncalojet_;
  int npfjet_;
  int ntau_;
  double muoptfirst_;
  double muoptother_;
  double muoeta_;
  double phopt_;
  double phoeta_;
  double elept_;
  double eleeta_;
  double pfelept_;
  double pfeleeta_;
  double calojetpt_;
  double calojeteta_;
  double pfjetpt_;
  double pfjeteta_;
  double metcalo_;
  double metpf_;
  double mettc_;
  double htc_;
  double PFhtc_;
  double taupt_;
  double taueta_;
  double qscale_low_;
  double qscale_high_;
  double muominv_;
  double muoDminv_;
  std::string trigger_;

  // Counters
  unsigned int nrEventTotalRaw_;
  unsigned int nrEventPassedQscaleRaw_;
  unsigned int nrEventPassedRaw_;

  int cpho_;
  int cele_;
  int cpfele_;
  int cmuo_;
  int ccalojet_;
  int cpfjet_;
  int ctau_;

  double localPi;
  
  int pre1;
  int pre2;
  
  double mumet_x;
  double mumet_y;
  double mumet_sumpt;
  
  double elemet_x;
  double elemet_y;
  double elemet_sumpt;
  

  //for ht cut
  double htcutsum_;
  double PFhtcutsum_;

  // Tree variables
  unsigned int mTreerun;
  unsigned int mTreeevent;
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
  
  
  

  

  Bool_t mTreenoiseHBHEFilterResult;
  Bool_t IsPythiaFiltered;
  int mTreeNFilter; 
  Bool_t mTreeFilterResults[100];
  int mTreeFilterName[100][20];

    
  int nTreePileUp;
  int mTreePUbx[3];
  int mTreePUNumInteractions[3];
  int mTreePUN[3];
  
  double mTreePUInstLumi[3][100];
  double mTreePUzPosi[3][100];
  double mTreePUsumPthi[3][100];
  double mTreePUsumPtlo[3][100];

  int mTreetrighltname[20];
  
  enum { nMaxTrigger = 7000 };
  int mTreeNtrig;
  int mTreetrigL1pre[nMaxTrigger];
  int mTreetrigHLTpre[nMaxTrigger];
  int mTreetrigname[nMaxTrigger][20];
  int mTreefiltname[nMaxTrigger][20];
  double mTreetrigpt[nMaxTrigger];
  double mTreetrigeta[nMaxTrigger];
  double mTreetrigphi[nMaxTrigger];
  
  double mTreeMET[10];
  double mTreeMEX[10];
  double mTreeMEY[10];
  double mTreeMETphi[10];
  double mTreeSumET[10];
  double mTreeSumETSignif[8];
  double mTreeMETSignif[8];
  double mTreeMETCaloMETInmHF[8];
  double mTreeMETCaloMETInpHF[8];
  double mTreeMETCaloMETPhiInmHF[8];
  double mTreeMETCaloMETPhiInpHF[8];
  double mTreeMETCaloSETInmHF[8];
  double mTreeMETCaloSETInpHF[8];
  double mTreeMETemEtFraction[8];
  double mTreeMETetFractionHadronic[8];
  double mTreeMETmaxEtInEmTowers[8];
  double mTreeMETmaxEtInHadTowers[8];
  double mTreeMETemEtInHF[8];
  double mTreeMETemEtInEE[8];
  double mTreeMETemEtInEB[8];
  double mTreeMEThadEtInHF[8];
  double mTreeMEThadEtInHE[8];
  double mTreeMEThadEtInHO[8];
  double mTreeMEThadEtInHB[8];
  
  double mTreeMETChargedEMEtFraction[10];
  double mTreeMETChargedHadEtFraction[10];
  double mTreeMETMuonEtFraction[10];
  double mTreeMETNeutralEMFraction[10];
  double mTreeMETNeutralHadEtFraction[10];
  double mTreeMETType6EtFraction[10];
  double mTreeMETType7EtFraction[10];

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

  int    mTreeNCalojet;
  int    mTreeCaloJetTruth[100];
  int    mTreeCaloJetPart[100];
  int    mTreeCaloJetConst[100];
  int    mTreeCaloJetn90[100];
  int    mTreeCaloJetn90hits[100];
  int    mTreeCaloJetID[100];
  double mTreeCaloJetEt[100];
  double mTreeCaloJetPt[100];
  double mTreeCaloJetPtRaw[100];
  double mTreeCaloJetP[100];
  double mTreeCaloJetPx[100];
  double mTreeCaloJetPy[100];
  double mTreeCaloJetPz[100];
  double mTreeCaloJetE[100];
  double mTreeCaloJetEta[100];
  double mTreeCaloJetPhi[100];
  double mTreeCaloJetFem[100];
  double mTreeCaloJetFhad[100];
  double mTreeCaloJetBtag[100];
  double mTreeCaloJetCharge[100];
  double mTreeCaloJetfhpd[100];
  double mTreeCaloJetfrbx[100];

  int    mTreeNPFjet;
  int    mTreePFJetTruth[100];
  int    mTreePFJetPart[100];
  int    mTreePFJetConst[100];
  int    mTreePFJetN[100][7];
  int    mTreePFJetn90[100];
  double mTreePFJetEt[100];
  double mTreePFJetPt[100];
  double mTreePFJetPtRaw[100];
  double mTreePFJetP[100];
  double mTreePFJetPx[100];
  double mTreePFJetPy[100];
  double mTreePFJetPz[100];
  double mTreePFJetE[100];
  double mTreePFJetEta[100];
  double mTreePFJetPhi[100];
  double mTreePFJetBtag[100];
  double mTreePFJetCharge[100];
  double mTreePFJetF[100][7];

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

  int    mTreeNSC;
  int    mTreeNSCtrign[100];
  int    mTreeSCTruth[200];
  double mTreeSCE[200];
  double mTreeSCPhi[200];
  double mTreeSCEta[200];
  double mTreeSCtrig[100][500];

  int    mTreeNpho;
  int    mTreePhoTruth[100];
  int    mTreeNphotrign[100];
  int    mTreePhotrig[100][500];

  double mTreePhoPFCandPx[100];
  double mTreePhoPFCandPy[100];
  double mTreePhoPFCandPz[100];
  double mTreePhoPFCandE[100];
  double mTreePhoPFCandeta[100];
  double mTreePhoPFCandphi[100];
  int mTreePhoPFCandpfid[100];
  double mTreePhoPFCandpfDeltaR[100];

  double mTreePhoEt[100];
  double mTreePhoP[100];
  double mTreePhoPt[100];
  double mTreePhoE[100];
  double mTreePhoRawE[100];
  double mTreePhoEta[100];
  double mTreePhoPhi[100];
  double mTreePhoPx[100];
  double mTreePhoPy[100];
  double mTreePhoPz[100];
  double mTreePhoSCEta[100];
  double mTreePhoDr03TkSumPt[100];
  double mTreePhoDr04TkSumPt[100];
  double mTreePhoSigmaIetaIeta[100];
  double mTreePhoEcaloIso[100];
  double mTreePhoHcaloIso[100];
  double mTreePhoTrackIso[100];
  double mTreePhoHCalOverEm[100];
  double mTreePhoe5x5[100];
  double mTreePhoe1x5[100];
  double mTreePhoe3x3[100];
  double mTreePhoSwissCross[100];
  bool mTreePhohasMatchedPromptElectron[100];
  double mTreePhoPFisoEG[100][3];


  double mTreePhoRoundness[100];
  double mTreePhoAngle[100];
  double mTreePhoR9[100];
  double mTreePhoSCEtaWidth[100];

  int    mTreePhoHasPixelSeed[100];
  int    mTreePhoHasConvTracks[100];
  int    mTreePhoisPF[100];


  int    mTreeNele;
  int    mTreeNeletrign[100];
  int    mTreeEletrig[100][500];
  int    mTreeEleTruth[100];
  int    mTreeEleSC[100];
  int    mTreeEleHits[100];
  int    mTreeEleValidHitFirstPxlB[100];
  int    mTreeEleTrkExpHitsInner[100];
  int    mTreeEleisECal[100];
  int    mTreeEleisTracker[100];
  int    mTreeEleClassification[100];
 
  double mTreeEleEt[100];
  double mTreeEleP[100];
  double mTreeEleTrackpt[100];
  double mTreeEleTrackptError[100];
  double mTreeElePt[100];
  double mTreeElePx[100];
  double mTreeElePy[100];
  double mTreeElePz[100];
  double mTreeEleE[100];
  double mTreeEleEta[100];
  double mTreeElePhi[100];
  double mTreeEleHCalOverEm[100];
  double mTreeEleHCalOverEmBc[100];
  double mTreeEleDr03TkSumPt[100];
  double mTreeEleDr04HCalTowerSumEt[100];
  double mTreeEleDr03HCalTowerSumEt[100];
  double mTreeEleDr04HCalTowerSumEtBc[100];
  double mTreeEleDr03HCalTowerSumEtBc[100];
  double mTreeEleDr04ECalRecHitSumEt[100];
  double mTreeEleDr03ECalRecHitSumEt[100];
  double mTreeEleSigmaIetaIeta[100];
  double mTreeEleDeltaEtaSuperClusterTrackAtVtx[100];
  double mTreeEleDeltaPhiSuperClusterTrackAtVtx[100];
  double mTreeEleTrkChiNorm[100];
  double mTreeEleCharge[100];
  double mTreeEled0vtx[100];
  double mTreeEled0bs[100];
  double mTreeElesd0[100];
  double mTreeEledzvtx[100];
  double mTreeEleConvdist[100];
  double mTreeEleConvdcot[100];
  double mTreeEleConvr[100];
  double mTreeElefbrem[100];
  double mTreeElePFisoEG[100][3];
  double mTreeElePFiso[100][3];
  double mTreeEledr03HcalDepth1[100];
  double mTreeEledr03HcalDepth2[100];
  double mTreeElee2x5Max[100];
  double mTreeElee5x5[100];
  double mTreeElee1x5[100];
  double mTreeEleCaloEt[100];
  double mTreeEleSCEta[100];
  double mTreeElehcalDepth1TowerSumEt03[100];
  double mTreeElehcalDepth2TowerSumEt03[100];
  double mTreeEleSwissCross[100];
  double mTreeEleEoverP[100];
  double mTreeEleECalEnergy[100];
  double mTreeEleTrackMomentumAtVtx[100];
  bool mTreeElehasMatchedConversion[100];
  double mTreeEleSCRawEt[100];

  
  double mTreeElePFCandPx[100];
  double mTreeElePFCandPy[100];
  double mTreeElePFCandPz[100];
  double mTreeElePFCandE[100];
  double mTreeElePFCandeta[100];
  double mTreeElePFCandphi[100];
  int mTreeElePFCandpfid[100];
  double mTreeElePFCandpfDeltaR[100];
  double mTreeElerho;
  double mTreeElerhoEG;

  int    mTreeNPFEle;
  int    mTreePFEleTruth[100];
  int    mTreeNPFEletrign[100];
  int    mTreePFEletrig[100][500];
  int    mTreePFEleSC[100];
  double mTreePFEleCharge[100];
  double mTreePFEleEt[100];
  double mTreePFEleCaloEt[100];
  double mTreePFEleP[100];
  double mTreePFElePt[100];
  double mTreePFElePx[100];
  double mTreePFElePy[100];
  double mTreePFElePz[100];
  double mTreePFEleE[100];
  double mTreePFEleEta[100];
  double mTreePFElePhi[100];
  double mTreePFEleSwissCross[100];
  double mTreePFEleTrackptError[100];
  double mTreePFEleTrackpt[100];
  double mTreePFEleSCEta[100];
  double mTreePFEleHCalOverEm[100];
  double mTreePFEleDr03TkSumPt[100];
  double mTreePFEleDr04HCalTowerSumEt[100];
  double mTreePFEleDr03HCalTowerSumEt[100];
  double mTreePFEleDr04ECalRecHitSumEt[100];
  double mTreePFEleDr03ECalRecHitSumEt[100];
  double mTreePFEleParticleIso[100];
  double mTreePFEleChadIso[100];
  double mTreePFEleNhadIso[100];
  double mTreePFEleGamIso[100];

  int    mTreeNmuo;
  int    mTreeNmuotrign[100];
  int    mTreeMuotrig[100][500];
  int    mTreeMuoTruth[100];
  int    mTreeMuoHitsCm[100];
  int    mTreeMuoHitsTk[100];
  int    mTreeMuoGood[100];
  int    mTreeMuoValidMuonHitsCm[100];
  int    mTreeMuoValidTrackerHitsCm[100];
  int    mTreeMuoValidPixelHitsCm[100];
  int    mTreeMuoChambersMatched[100];
  int    mTreeMuoStationsMatched[100];
  int    mTreeMuoTrackerLayersMeasCm[100];
  int    mTreeMuoTrackerLayersNotMeasCm[100];
  int    mTreeMuoValidMuonHitsTk[100];
  int    mTreeMuoValidTrackerHitsTk[100];
  int    mTreeMuoValidPixelHitsTk[100];
  int    mTreeMuoTrackerLayersMeasTk[100];
  int    mTreeMuoTrackerLayersNotMeasTk[100];
  int    mTreeMuoLostHits[100];
  int    mTreeMuoLostHitsTk[100];
  int    mTreeMuoID[100][24];
  int    mTreeMuoIsPF[100];

  
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
  double mTreeMuodzTk[100];
  double mTreeMuocalocomp[100];
  double mTreeMuocaltowe[100];
  double mTreeMuod0OriginCm[100];
  double mTreeMuovx[100];
  double mTreeMuovy[100];
  double mTreeMuovz[100];
  double mTreeMuoValidFraction[100];
  double mTreeMuoPFiso[100][7];
  
  double mTreeMuoCocktailPt[100];
  double mTreeMuoCocktailPhi[100];
  double mTreeMuoCocktailEta[100];
  
  double mTreeMuoTevRecoPtError[100][7];
  double mTreeMuoTevRecoPt[100][7];
  double mTreeMuoTevRecoEta[100][7];
  double mTreeMuoTevRecoPhi[100][7];
  double mTreeMuoTevRecoChi2[100][7];
  double mTreeMuoTevRecoNdof[100][7];
  
  double mTreeMuoPFCandPx[100];
  double mTreeMuoPFCandPy[100];
  double mTreeMuoPFCandPz[100];
  double mTreeMuoPFCandE[100];
  double mTreeMuoPFCandeta[100];
  double mTreeMuoPFCandphi[100];
  int mTreeMuoPFCandpfid[100];
  double mTreeMuoPFCandpfDeltaR[100];
  
  int    mTreeNPFmuons;
  int    mTreePFMuoHitsCm[100];
  int    mTreePFMuoHitsTk[100];
  int    mTreePFMuoValidMuonHitsCm[100];
  int    mTreePFMuoValidTrackerHitsCm[100];
  int    mTreePFMuoValidPixelHitsCm[100];
  int    mTreePFMuoChambersMatched[100];
  int    mTreePFMuoTrackerLayersMeasCm[100];
  int    mTreePFMuoTrackerLayersNotMeasCm[100];
  //int    mTreePFMuoID[100][24];
  
  double mTreePFMuonP[100];
  double mTreePFMuonPt[100];
  double mTreePFMuonE[100];
  double mTreePFMuonEt[100];
  double mTreePFMuonPx[100];
  double mTreePFMuonPy[100];
  double mTreePFMuonPz[100];
  double mTreePFMuonEta[100];
  double mTreePFMuonPhi[100];
  double mTreePFMuonCharge[100];
  double mTreePFMuonParticleIso[100];
  double mTreePFMuonChadIso[100];
  double mTreePFMuonNhadIso[100];
  double mTreePFMuonGamIso[100];
  double mTreePFMuonTrkIso[100];
  double mTreePFMuonRelTrkIso[100];
  double mTreePFMuonECalIso[100];
  double mTreePFMuonHCalIso[100];
  double mTreePFMuonAllIso[100];
  double mTreePFMuoECalIsoDep[100];
  double mTreePFMuoHCalIsoDep[100];
  double mTreePFMuoTrkIsoDep[100];
  double mTreePFMuoTrkChiNormCm[100];
  double mTreePFMuoTrkChiNormTk[100];
  double mTreePFMuoCharge[100];
  double mTreePFMuod0Cm[100];
  double mTreePFMuod0Tk[100];
  double mTreePFMuosd0Cm[100];
  double mTreePFMuosd0Tk[100];
  double mTreePFMuocalocomp[100];
  double mTreePFMuocaltowe[100];
  double mTreePFMuod0OriginCm[100];
  double mTreePFMuoCmvx[100];
  double mTreePFMuoCmvy[100];
  double mTreePFMuoCmvz[100];
  double mTreePFMuoValidFraction[100];
  
  
  int mTreeNtaus;  
  int    mTreeNtautrign[100];
  int    mTreeTautrig[100][500];
  int mTreeTauDecayMode[100];
  int mTreeTauPFChargedHadrCands[100];
  int mTreeTauPFGammaCands[100];
  double mTreeTauID[100][18];
  double mTreeTauP[100];
  double mTreeTauPt[100];
  double mTreeTauE[100];
  double mTreeTauEt[100];
  double mTreeTauMass[100];
  double mTreeTauMt[100];
  double mTreeTauPx[100];
  double mTreeTauPy[100];
  double mTreeTauPz[100];
  double mTreeTauEta[100];
  double mTreeTauPhi[100];
  int mTreeTauCharge[100];
  double mTreeTauvx[100];
  double mTreeTauvy[100];
  double mTreeTauvz[100];
  double mTreeTauvx2[100];
  double mTreeTauvy2[100];
  double mTreeTauvz2[100];
  double mTreeTauParticleIso[100];
  double mTreeTauChadIso[100];
  double mTreeTauNhadIso[100];
  double mTreeTauGamIso[100];


  int mTreeNTruthMatchTaus;
  int mTreePosTruthMatchTaus[100];
  int mTreeGenTauDecay[100][50];

  double mTreeTauGenJetE[100];
  double mTreeTauGenJetEt[100];
  double mTreeTauGenJetEmE[100];
  double mTreeTauGenJetHadE[100];
  double mTreeTauGenJetInvE[100];
  double mTreeTauGenJetEta[100];
  double mTreeTauGenJetPhi[100];
  double mTreeTauGenJetMass[100];
  double mTreeTauGenJetMt[100];
  double mTreeTauGenJetP[100];
  double mTreeTauGenJetPt[100];
  double mTreeTauGenJetPx[100];
  double mTreeTauGenJetPy[100];
  double mTreeTauGenJetPz[100];

  double mTreeTauIsolationPFChargedHadrCandsPtSum[100]; 
  double mTreeTauIsolationPFGammaCandsEtSum[100]; 
  double mTreeTauEcalStripSumEOverPLead[100]; 
  double mTreeTauEMFraction[100]; 
  double mTreeTauHcal3x3OverPLead[100]; 
  double mTreeTauHcalMaxOverPLead[100]; 
  double mTreeTauHcalTotOverPLead[100]; 
  double mTreeTauLeadPFChargedHadrCandsignedSipt[100]; 
  double mTreeTauPhiphiMoment[100]; 
  double mTreeTauEtaphiMoment[100]; 
  double mTreeTauEtaetaMoment[100]; 
  double mTreeTauElectronPreIDOutput[100]; 
  double mTreeTauPFLeadChargedPT[100];
  int mTreeTauNSignalTracks[100];
  double mTreeTauBremsRecoveryEOverPLead[100];
  
  double mTreesusyScanM0;
  double mTreesusyScanM12;
  double mTreesusyScanA0;
  double mTreesusyScanCrossSection;
  double mTreesusyScanMu;
  double mTreesusyScanRun;
  double mTreesusyScantanbeta;
    
  int    mTreeNvtx;
  int    mTreeVtxntr[100];
  int    mTreeVtxfake[100];
  double mTreeVtxndf[100];
  double mTreeVtxx[100];
  double mTreeVtxy[100];
  double mTreeVtxz[100];
  double mTreeVtxchi[100];

  double mTreebspX;
  double mTreebspY;
  double mTreebspZ;

  double mTreeEventWeight;
  int    mTreeProcID;
  double mTreeQscale;
  double mTreebfield;

};
