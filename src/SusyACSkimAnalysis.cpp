//
// Package:    UserCode/aachen3a/ACSusyAnalysis
// Class:      SusyACSkimAnalysis
// 
// Description: Skeleton analysis for SUSY search with Lepton + Jets + MET
//
// Original Author:  Carsten Magass
//          Created: November 2008
//
#include "aachen3a/ACSusyAnalysis/interface/SusyACSkimAnalysis.h"
 
////////////////////////////////
//
// Constructor
//
SusyACSkimAnalysis::SusyACSkimAnalysis(const edm::ParameterSet& iConfig):
  nrEventTotalRaw_(0),
  nrEventPassedPthatRaw_(0),
  nrEventPassedRaw_(0)
{

  // get the data tags
  calojetTag_           = iConfig.getParameter<edm::InputTag>("calojetTag");
  pfjetTag_              = iConfig.getParameter<edm::InputTag>("pfjetTag");
  metTag_               = iConfig.getParameter<edm::InputTag>("metTag");
  metTagPF_           = iConfig.getParameter<edm::InputTag>("metTagPF");
  metTagTC_           = iConfig.getParameter<edm::InputTag>("metTagTC");
  elecTag_              = iConfig.getParameter<edm::InputTag>("elecTag");
  PFelecTag_           = iConfig.getParameter<edm::InputTag>("pfelecTag");
  muonTag_           = iConfig.getParameter<edm::InputTag>("muonTag");
  PFmuonTag_        = iConfig.getParameter<edm::InputTag>("PFmuonTag");
  tauSrc_                = iConfig.getParameter<edm::InputTag>("tauTag");
  metTagPFnoPU_   = iConfig.getParameter<edm::InputTag >("metTagPFnoPU");
  metTagJESCorAK5CaloJetMuons_  =iConfig.getParameter<edm::InputTag >("metTagJESCorAK5CaloJetMuons");
  metTagcorMetGlobalMuons_      = iConfig.getParameter<edm::InputTag >("metTagcorMetGlobalMuons");
  metTagHO_          = iConfig.getParameter<edm::InputTag >("metTagHO");
  metTagNoHF_       = iConfig.getParameter<edm::InputTag >("metTagNoHF");
  genTag_               = iConfig.getParameter<edm::InputTag>("genTag");
  genJetTag_           = iConfig.getParameter<edm::InputTag>("genJetTag");
  vertexTag_           = iConfig.getParameter<edm::InputTag>("vtxTag");
  ebhitsTag_           = iConfig.getParameter<edm::InputTag>("ebhitsTag");

  is_MC                 = iConfig.getParameter<bool>("is_MC");
  is_SHERPA          = iConfig.getParameter<bool>("is_SHERPA");
  do_fatjets           = iConfig.getParameter<bool>("do_fatjets");
  matchAll_           = iConfig.getParameter<bool>("matchAll");
  susyPar_             = iConfig.getParameter<bool>("susyPar");

  edm::LogVerbatim("SusyACSkimAnalysis") << " Running with flag is_MC      = " << is_MC << endl;
  edm::LogVerbatim("SusyACSkimAnalysis") << " Running with flag is_SHERPA  = " << is_SHERPA << endl;
  edm::LogVerbatim("SusyACSkimAnalysis") << " Running with flag do_fatjets = " << do_fatjets << endl;

  btag_ = iConfig.getParameter<std::string>("btag");

  vers_ = iConfig.getParameter<std::string>("jetselvers");
  qual_ = iConfig.getParameter<std::string>("jetselqual");

  edm::LogVerbatim("SusyACSkimAnalysis") << " JetID version [ " << vers_ << " ] quality [ " << qual_ << " ] " << endl;

  // Implementation according to 
  // PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h
  if      ( qual_ == "MINIMAL" )   qual = JetIDSelectionFunctor::MINIMAL;
  else if ( qual_ == "LOOSE_AOD" ) qual = JetIDSelectionFunctor::LOOSE_AOD;
  else if ( qual_ == "LOOSE" )     qual = JetIDSelectionFunctor::LOOSE;
  else if ( qual_ == "TIGHT" )     qual = JetIDSelectionFunctor::TIGHT;
  else       
    throw cms::Exception("InvalidInput") << "Expect quality to be one of MINIMAL, LOOSE_AOD, LOOSE,TIGHT" << std::endl;

  if      ( vers_ == "CRAFT08" ) vers = JetIDSelectionFunctor::CRAFT08;
  else if ( vers_ == "PURE09" )  vers = JetIDSelectionFunctor::PURE09;
  else if ( vers_ == "DQM09" )   vers = JetIDSelectionFunctor::DQM09;
  else 
    throw cms::Exception("InvalidInput") << "Expect version to be one of CRAFT08, PURE09, DQM09" << std::endl;

  // get the cuts
  muopt_      = iConfig.getParameter<double>("muopt");
  muoeta_     = iConfig.getParameter<double>("muoeta");
  elept_      = iConfig.getParameter<double>("elept");
  eleeta_     = iConfig.getParameter<double>("eleeta");
  pfelept_    = iConfig.getParameter<double>("pfelept");
  pfeleeta_   = iConfig.getParameter<double>("pfeleeta");
  calojetpt_  = iConfig.getParameter<double>("calojetpt");
  calojeteta_ = iConfig.getParameter<double>("calojeteta");
  pfjetpt_    = iConfig.getParameter<double>("pfjetpt");
  pfjeteta_   = iConfig.getParameter<double>("pfjeteta");
  metcalo_    = iConfig.getParameter<double>("metcalo");
  metpf_      = iConfig.getParameter<double>("metpf");
  mettc_      = iConfig.getParameter<double>("mettc");
  htc_        = iConfig.getParameter<double>("htc");
  PFhtc_      = iConfig.getParameter<double>("PFhtc");

  nele_     = iConfig.getParameter<int>("nele");
  npfele_   = iConfig.getParameter<int>("npfele");
  nmuo_     = iConfig.getParameter<int>("nmuo");
  ncalojet_ = iConfig.getParameter<int>("ncalojet");
  npfjet_   = iConfig.getParameter<int>("npfjet");

  pthat_low_  = iConfig.getParameter<double>("pthat_low");
  pthat_high_ = iConfig.getParameter<double>("pthat_high");

  localPi = acos(-1.0);

  // Initialize 
  initPlots();

}

////////////////////////////////
//
// Destructor
//
SusyACSkimAnalysis::~SusyACSkimAnalysis() {}

////////////////////////////////
//
// Called in for each run
//
void SusyACSkimAnalysis::beginMyRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  
  // Initialize HLTConfigProvider
  hltConfigInit_ = false;
  bool changed( true );
  if ( ! hltConfig_.init( iRun, iSetup, processName_, changed ) ) {
    edm::LogError( "SusyACSkimAnalysis" ) << "HLT config extraction error with process name " << processName_;
  } 
  else if ( hltConfig_.size() <= 0 ) {
    edm::LogError( "SusyACSkimAnalysis" ) << "HLT config size error";
  } 
  else hltConfigInit_ = true;

}

////////////////////////////////
//
// Called in for each event
//
bool SusyACSkimAnalysis::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  using namespace reco;
  using namespace pat;
  cmuo_     = 0;
  cele_     = 0;
  ccalojet_ = 0;
  cpfjet_   = 0;

  pre1 = 0;
  pre2 = 0;
  
  for (int i=0; i<20; i++) {
    mTreetrighltname[i] = 0;

    for (int j=0; j<1000; j++) {
      mTreetrigname[j][i] = 0;
      mTreefiltname[j][i] = 0;
    }
  }

  JetIDSelectionFunctor jetId( vers, qual );

  mTreeNtrig = 0;

  // procedure for getting auto process trigger process name
  // and accessing prescales
  Handle< trigger::TriggerEvent > handleTriggerEvent;
  iEvent.getByLabel( "hltTriggerSummaryAOD", handleTriggerEvent );
  const edm::Provenance *meta = handleTriggerEvent.provenance();
  processName_ = meta->processName();

  beginMyRun(iEvent.getRun(), iSetup);

  Handle< pat::TriggerEvent > myTriggerEvent;
  iEvent.getByLabel( "patTriggerEvent", myTriggerEvent );

  if ( !myTriggerEvent.isValid() )
    edm::LogWarning("SusyACSkimAnalysis") << "No pat::TriggerEvent found for InputTag patTriggerEvent";
  else {

    int *tlname = pack(myTriggerEvent->nameHltTable().c_str());
    for (int l=0; l<get_size(tlname); l++) mTreetrighltname[l] = tlname[l];
    
    const TriggerPathRefVector mypath = myTriggerEvent->acceptedPaths();
    for (TriggerPathRefVector::const_iterator it = mypath.begin();
	 it != mypath.end(); ++it) {
      
      TString ttname = (*it)->name();
      std::string tname   = (*it)->name();

      // save all accepted HLT_ trigger results
      if ( myTriggerEvent->path(tname)->wasAccept() &&
	   ttname.Contains("HLT_") &&
	   !ttname.Contains("AlCa") &&
	   !ttname.Contains("L1Tech") &&
	   !ttname.Contains("NZS") &&
	   !ttname.Contains("_step") &&
	   !ttname.Contains("DQM")) {

	pre1 = -1;
	pre2 = -1;
	if ( hltConfigInit_ ) {

	  l1GtUtils_.retrieveL1EventSetup(iSetup);

	  if (hltConfig_.hltL1GTSeeds(tname).size() == 1) {

	    l1error = 0;
	    l1GtUtils_.prescaleFactor(iEvent,
				      hltConfig_.hltL1GTSeeds(tname).at(0).second,
				      l1error);
	    // 	    cout << hltConfig_.hltL1GTSeeds(tname).size() << " <-- " << l1error << endl;

	    if (!l1error) {
	      const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,tname));
	      pre1 = prescales.first;
	      pre2 = prescales.second;
	    }
	  }
	}
	else
	  edm::LogWarning("SusyACSkimAnalysis") << "HLT configuration not properly extracted";
	
	//	cout << "   --> " << ttname << endl;
	
	int *tempname = pack(tname.c_str());
	
	const TriggerFilterRefVector mpf = myTriggerEvent->pathFilters(tname,false);
	for ( TriggerFilterRefVector::const_iterator ll=mpf.begin(); ll!=mpf.end(); ++ll ) {
	  TriggerObjectRefVector torv = myTriggerEvent->filterObjects((*ll)->label());   
	  for ( TriggerObjectRefVector::const_iterator itt = torv.begin(); 
		itt != torv.end(); ++itt ) {
	    const TriggerObjectRef objRef( *itt );
	    
	    //	    cout << "  " << (*ll)->label() << "  " << (*itt)->pt() << "  "  << (*itt)->eta() << endl;
	    int *filtname = pack((*ll)->label().c_str());

	    for (int l=0; l<get_size(tempname); l++) mTreetrigname[mTreeNtrig][l] = tempname[l];
	    for (int l=0; l<get_size(filtname); l++) mTreefiltname[mTreeNtrig][l] = filtname[l];
	    
	    mTreetrigL1pre[mTreeNtrig]  = pre1; 
	    mTreetrigHLTpre[mTreeNtrig] = pre2; // (*it)->prescale();
	    mTreetrigpt[mTreeNtrig]  = objRef->pt();
	    mTreetrigeta[mTreeNtrig] = objRef->eta();
	    mTreetrigphi[mTreeNtrig] = objRef->phi();
	    mTreeNtrig++;
	    if (mTreeNtrig==1000) break;
	  }
	  if (mTreeNtrig==1000) break;
	}
	if (mTreeNtrig==1000) break;
	
	/*
	const TriggerFilterRefVector mpf = myTriggerEvent->pathFilters(tname);
	for ( TriggerFilterRefVector::const_iterator ll=mpf.begin(); ll!=mpf.end(); ++ll )
	  cout << (*ll)->label() << "  " << (*ll)->type() << endl;

	// object ID according to HLTReco/interface/TriggerTypeDefs.h
	// 81 : TriggerL1Mu
	// 93 : TriggerMuon

	const TriggerObjectRefVector triggerMuonsL2L3 = myTriggerEvent->objects( 93 ); 
	for ( TriggerObjectRefVector::const_iterator tt = triggerMuonsL2L3.begin(); 
	      tt != triggerMuonsL2L3.end(); ++tt ) {
	  const TriggerObjectRef objRef( *tt );

	  if ( myTriggerEvent->objectInPath(*tt, tname) ) 
	    cout << "  " << (*tt)->pt() << "  "  << (*tt)->eta() << endl;
	}
	
	const TriggerObjectRefVector triggerMuonsL1 = myTriggerEvent->objects( 81 ); 
	for ( TriggerObjectRefVector::const_iterator tt = triggerMuonsL1.begin(); 
	      tt != triggerMuonsL1.end(); ++tt ) {
	  const TriggerObjectRef objRef( *tt );
	  
	  if ( myTriggerEvent->objectInPath(*tt, tname) ) 
	    cout << "  " << (*tt)->pt() << "  "  << (*tt)->eta() << endl;
	}
	*/
      }
    }
  }
  // Count all events
  nrEventTotalRaw_++;

  // Global event variables
  mTreerun     = iEvent.id().run();
  mTreeevent   = iEvent.id().event();
  mTreelumiblk = iEvent.luminosityBlock();
  mTreebx      = iEvent.bunchCrossing();
  mTreeorbit   = iEvent.orbitNumber();
  mTreeexp     = iEvent.experimentType();
  mTreedata    = iEvent.isRealData();
  mTreestore   = iEvent.eventAuxiliary().storeNumber();

  // Luminosity information
  LuminosityBlock const& lumiBlock = iEvent.getLuminosityBlock();

  // iEvent.luminosityBlock() = lumiBlock.luminosityBlock()

  Handle<LumiSummary> ls;
  lumiBlock.getByLabel("lumiProducer", ls);

  mTreedellumi = mTreereclumi = mTreedellumierr = mTreereclumierr = 0.;
  if ( !ls.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No lumiProducer found ";
  else {
    if (!is_MC) {
      mTreedellumi    = ls->avgInsDelLumi();
      mTreedellumierr = ls->avgInsDelLumiErr();
      mTreereclumi    = ls->avgInsRecLumi();
      mTreereclumierr = ls->avgInsRecLumiErr();
    }
  }
  // Get some MC event information (process ID, weight, PDF)
  mTreeProcID       = -1;
  mTreeEventWeight  = 1.;
  mTreePthat        = -999.;

  //  if (mTreeevent!=9703142 && mTreeevent!=4215340 && mTreeevent!=12799746) return 0;

  if (is_MC) {

    Handle<GenEventInfoProduct> evt_info;
    iEvent.getByType(evt_info);
    
    if ( !evt_info.isValid() )
      edm::LogWarning("SusyACSkimAnalysis") << "No GenEventInfoProduct found ";
    else {
      mTreeProcID      = evt_info->signalProcessID();
      mTreeEventWeight = evt_info->weight();
      mTreePthat       = evt_info->qScale(); 

      const GenEventInfoProduct::PDF *tpdf = evt_info->pdf();
      mTreepdfid1   = (int)tpdf->id.first;
      mTreepdfid2   = (int)tpdf->id.second;
      mTreepdfx1    = tpdf->x.first;
      mTreepdfx2    = tpdf->x.second;
      mTreepdff1    = tpdf->xPDF.first;
      mTreepdff2    = tpdf->xPDF.second;
      mTreepdfscale = tpdf->scalePDF;

    }
    // IMPORTANT for QCD : avoid overlap of samples
    if (mTreePthat>-990. && pthat_low_ > -1. && pthat_high_> -1.) {
      
      if (mTreePthat < pthat_low_ || mTreePthat > pthat_high_)
    
	return 0;

    }
    nrEventPassedPthatRaw_++;

  }
  
  
  nTreePileUp=0;
  if(is_MC){
    //PileUp Information
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
    if(!PupInfo.isValid()) {
    edm::LogWarning("SusyACSkimAnalysis") << " could not find PileupSummeryInfo.\n";
    }
    else{
      nTreePileUp=PupInfo->size();
      
      if(nTreePileUp>20) nTreePileUp=20;
      
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      int i=0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
        mTreePUbx                 = PVI->getBunchCrossing();
        
        mTreePUNumInteractions[i]       = PVI->getPU_NumInteractions();
        mTreePUN[i]                     = PVI->getPU_instLumi().size();
        
        for(int j=0; j<mTreePUN[i];j++){
          mTreePUInstLumi[i][j]        = PVI->getPU_instLumi()[j];
          mTreePUzPosi[i][j]           = PVI->getPU_zpositions()[j];
          mTreePUsumPthi[i][j]         = PVI->getPU_sumpT_highpT()[j];
          mTreePUsumPtlo[i][j]         = PVI->getPU_sumpT_lowpT()[j];
        }
        i++;
        
      }

    }
  }



  
  /*
  edm::LogInfo("SusyACSkimAnalysis") << " Event properties: EvtWeight = " << mTreeEventWeight 
				     << "  ProcID = " << mTreeProcID 
				     << "  Pthat = " <<  mTreePthat << std::endl;
  */

  // HCAL Noise
  Handle<HcalNoiseSummary> noiseHandle;
  iEvent.getByLabel("hcalnoise", noiseHandle);
  
  if ( !noiseHandle.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No HcalNoiseSummary found  ";
  }
  else {
    const HcalNoiseSummary noisesummary = *noiseHandle;
    mTreenoisel = noisesummary.passLooseNoiseFilter();
    mTreenoiset = noisesummary.passTightNoiseFilter();
    mTreenoiseh = noisesummary.passHighLevelNoiseFilter();

  }

  //HCAL Noise advanced!!
  edm::Handle<HcalNoiseSummary> summary_h;
  iEvent.getByType(summary_h);
  if(!summary_h.isValid()) {
    edm::LogWarning("SusyACSkimAnalysis") << " could not find HcalNoiseSummary.\n";
  }
  else{
    const HcalNoiseSummary summary = *summary_h;
    
    //these are to study HCAL noise in 50ns bunches can be removed if not needed
    mTreenoiseHCALeventChargeFraction   = summary.eventChargeFraction();
    mTreenoiseHCALeventEMEnergy         = summary.eventEMEnergy();
    mTreenoiseHCALeventEMFraction       = summary.eventEMFraction();
    mTreenoiseHCALeventHadEnergy        = summary.eventHadEnergy();
    mTreenoiseHCALeventTrackEnergy      = summary.eventTrackEnergy();
    mTreenoiseHCALflatNoiseSumE         = summary.flatNoiseSumE();
    mTreenoiseHCALflatNoiseSumEt        = summary.flatNoiseSumEt();
    mTreenoiseHCALHasBadRBXTS4TS5       = summary.HasBadRBXTS4TS5();
    //mTreenoiseHCALhighLevelNoiseTowers = summary.highLevelNoiseTowers();
    mTreenoiseHCALisolatedNoiseSumE     = summary.isolatedNoiseSumE();
    mTreenoiseHCALisolatedNoiseSumEt    = summary.isolatedNoiseSumEt();
    //mTreenoiseHCALlooseNoiseTowers      = summary.looseNoiseTowers();
    mTreenoiseHCALmax10GeVHitTime       = summary.max10GeVHitTime();
    mTreenoiseHCALmax25GeVHitTime       = summary.max25GeVHitTime();
    mTreenoiseHCALmaxE10TS              = summary.maxE10TS();
    mTreenoiseHCALmaxE2Over10TS         = summary.maxE2Over10TS();
    mTreenoiseHCALmaxE2TS               = summary.maxE2TS();
    mTreenoiseHCALmaxHPDHits            = summary.maxHPDHits();
    mTreenoiseHCALmaxHPDNoOtherHits     = summary.maxHPDNoOtherHits();
    mTreenoiseHCALmaxRBXHits            = summary.maxRBXHits();
    mTreenoiseHCALmaxZeros              = summary.maxZeros();
    mTreenoiseHCALmin10GeVHitTime       = summary.min10GeVHitTime();
    mTreenoiseHCALmin25GeVHitTime       = summary.min25GeVHitTime();
    mTreenoiseHCALminE10TS              = summary.minE10TS();
    mTreenoiseHCALminE2Over10TS         = summary.minE2Over10TS();
    mTreenoiseHCALminE2TS               = summary.minE2TS();
    mTreenoiseHCALminHPDEMF             = summary.minHPDEMF();
    mTreenoiseHCALminRBXEMF             = summary.minRBXEMF();
    mTreenoiseHCALnoiseFilterStatus     = summary.noiseFilterStatus();
    mTreenoiseHCALnoiseType             = summary.noiseType();
    mTreenoiseHCALnum10GeVHits          = summary.num10GeVHits();
    mTreenoiseHCALnum25GeVHits          = summary.num25GeVHits();
    mTreenoiseHCALnumFlatNoiseChannels  = summary.numFlatNoiseChannels();
    mTreenoiseHCALnumIsolatedNoiseChannels = summary.numIsolatedNoiseChannels();
    mTreenoiseHCALnumProblematicRBXs    = summary.numProblematicRBXs();
    mTreenoiseHCALnumSpikeNoiseChannels = summary.numSpikeNoiseChannels();
    mTreenoiseHCALnumTriangleNoiseChannels = summary.numTriangleNoiseChannels();
    mTreenoiseHCALnumTS4TS5NoiseChannels = summary.numTS4TS5NoiseChannels();
    mTreenoiseHCALpassHighLevelNoiseFilter = summary.passHighLevelNoiseFilter();
    mTreenoiseHCALpassLooseNoiseFilter  = summary.passLooseNoiseFilter();
    mTreenoiseHCALpassTightNoiseFilter  = summary.passTightNoiseFilter();
    //mTreenoiseHCALproblematicJets       = summary.problematicJets();
    mTreenoiseHCALrms10GeVHitTime       = summary.rms10GeVHitTime();
    mTreenoiseHCALrms25GeVHitTime       = summary.rms25GeVHitTime();
    mTreenoiseHCALspikeNoiseSumE        = summary.spikeNoiseSumE();
    mTreenoiseHCALspikeNoiseSumEt       = summary.spikeNoiseSumEt();
    //mTreenoiseHCALtightNoiseTowers      = summary.tightNoiseTowers();
    mTreenoiseHCALtriangleNoiseSumE     = summary.triangleNoiseSumE();
    mTreenoiseHCALtriangleNoiseSumEt    = summary.triangleNoiseSumEt();
    mTreenoiseHCALTS4TS5NoiseSumE       = summary.TS4TS5NoiseSumE();
    mTreenoiseHCALTS4TS5NoiseSumEt      = summary.TS4TS5NoiseSumEt();

  }
  
  
  //this does not work for the moment!!
  
  //Hcal Noise cut variable:
  //edm::Handle<bool> hcalnoiseFlag;
  //iEvent.getByLabel("HBHENoiseFilterResult",hcalnoiseFlag);
  //if(!hcalnoiseFlag.isValid()) {
    //edm::LogWarning("SusyACSkimAnalysis") << " could not find HcalNoiseFlag.\n";
    //mTreenoiseHCALhcalnoiseFlag=true;
  //}
  //else{
    //mTreenoiseHCALhcalnoiseFlag=*hcalnoiseFlag;
  //}
  // Vertex
  edm::Handle< reco::VertexCollection > vertexHandle;
  iEvent.getByLabel(vertexTag_, vertexHandle);

  math::XYZPoint vtxPoint(0,0,0);

  if ( !vertexHandle.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::Vertex found for InputTag " << vertexTag_;
  }
  else {

    mTreeNvtx = vertexHandle->size();
    if ( mTreeNvtx > 100 ) mTreeNvtx = 100;

    for (int i=0; i<mTreeNvtx; ++ i ) {

      mTreeVtxfake[i] = (*vertexHandle)[i].isFake();
      mTreeVtxchi[i]  = (*vertexHandle)[i].normalizedChi2();
      mTreeVtxntr[i]  = (*vertexHandle)[i].tracksSize();
      mTreeVtxndf[i]  = (*vertexHandle)[i].ndof();
      mTreeVtxx[i]    = (*vertexHandle)[i].x();
      mTreeVtxy[i]    = (*vertexHandle)[i].y();
      mTreeVtxz[i]    = (*vertexHandle)[i].z();
    }
  }
  if (mTreeNvtx>0) vtxPoint.SetXYZ(mTreeVtxx[0], mTreeVtxy[0], mTreeVtxz[0]);

  // Beam Spot
  edm::Handle< reco::BeamSpot > BeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", BeamSpotHandle);

  math::XYZPoint bsPoint(0,0,0);

  if ( !BeamSpotHandle.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::BeamSpot found for InputTag offlineBeamSpot";
  }
  else {
  
    const reco::BeamSpot offlineBeamSpot = *BeamSpotHandle;

    mTreebspX = BeamSpotHandle->x0();
    mTreebspY = BeamSpotHandle->y0();
    mTreebspZ = BeamSpotHandle->z0();
    bsPoint.SetXYZ(mTreebspX, mTreebspY, mTreebspZ);
  }
    
  // Track information
  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel("generalTracks", tkRef);   

  if ( !tkRef.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::TrackCollection found for InputTag generalTracks";
  }
  else {

    const reco::TrackCollection* tkColl = tkRef.product();

    mTreeNtracks = tkColl->size();

    int numhighpurity = 0;
    reco::TrackBase::TrackQuality  _trackQuality = reco::TrackBase::qualityByName("highPurity");

    reco::TrackCollection::const_iterator itk   = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();

    for ( ; itk!=itk_e; ++itk) {

      if (itk->quality(_trackQuality)) numhighpurity++;
    }
    mTreetrackshqf = (double)numhighpurity/(double)mTreeNtracks;
  }
  // Magnet Fiels
  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel("scalersRawToDigi", dcsHandle);

  bfield = -1.;

  if ( !dcsHandle.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No DcsStatusCollection found for InputTag scalersRawToDigi";
  }
  else {

    if (is_MC) {
      
      ESHandle<MagneticField> magneticField;
      iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
      
      bfield = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
      
    }
    else {
      
      float currentToBFieldScaleFactor = 2.09237036221512717e-04;
      float current = (*dcsHandle)[0].magnetCurrent();
      bfield = current*currentToBFieldScaleFactor;
    }
  }

  mTreebfield = bfield;

  // Truth tree
  mTreeNtruth  = 0;
  mTreeNtruthl = 0;

  std::vector<const reco::Candidate*> truthl;
  std::vector<const reco::Candidate*> truth;

  if (is_MC) {
    
    Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(genTag_, genParticles);    
    
    if ( !genParticles.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenParticleCollection found with InputTag genParticles ";
    else {
  
      std::vector<const reco::Candidate*> part1;
      std::vector<const reco::Candidate*> part2;
      std::vector<int> bvtx1;
      std::vector<int> bvtx2;
      int vtxid = 2;
      int ntruth = 0;
      int evtxid = -999;

      for( size_t i = 0; i < genParticles->size(); ++ i ) {
	const reco::Candidate& p = (*genParticles)[ i ];
       
	//if (p.status() == 1 && (abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)) {
	if ((!matchAll_&&(p.status() == 1 && (abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)))||( matchAll_&& (abs(p.pdgId())< 19&&abs(p.pdgId())> 10))) {
	  mTreetruthlpdgid[mTreeNtruthl] = p.pdgId();
	  mTreetruthlE[mTreeNtruthl]     = p.p4().E();
	  mTreetruthlEt[mTreeNtruthl]    = p.p4().Et();
	  mTreetruthlp[mTreeNtruthl]     = p.p4().P();
	  mTreetruthlpt[mTreeNtruthl]    = p.p4().Pt();
	  mTreetruthlpx[mTreeNtruthl]    = p.p4().px();
	  mTreetruthlpy[mTreeNtruthl]    = p.p4().py();
	  mTreetruthlpz[mTreeNtruthl]    = p.p4().pz();
	  mTreetruthleta[mTreeNtruthl]   = p.p4().eta();
	  mTreetruthlphi[mTreeNtruthl]   = p.p4().phi();

	  truthl.push_back(&p);

	  mTreeNtruthl++;
	  if (mTreeNtruthl == 100) break;
	}
      }

      // (1) Reconstruct 2 -> 2 primary process, if t, W, Z, tau or SUSY particles are involved
      int count_sherpa = 0; 
      for( size_t i = 0; i < genParticles->size(); ++ i ) {
	const reco::Candidate& p = (*genParticles)[ i ];
      
	if (is_SHERPA && p.status() == 3) {

	  count_sherpa++;  
	  //	  cout << i << "  PDGID : " << p.pdgId() << "  " << p.numberOfMothers() << endl;
	}

	if (p.numberOfMothers()==2) {
	  
	  if ( (!is_SHERPA && p.status() == 3) || // && isDecaying(p.pdgId())
	       ( is_SHERPA && count_sherpa == 3) ) {
	    
	    const reco::Candidate &m0 = *(p.mother(0));
	    const reco::Candidate &m1 = *(p.mother(1));
	    
	    mTreetruthpdgid[0]  = m0.pdgId();
	    mTreetruthbvtxid[0] = 0;
	    mTreetruthevtxid[0] = 1;
	    mTreetruthE[0]      = m0.p4().E();
	    mTreetruthEt[0]     = m0.p4().Et();
	    mTreetruthp[0]      = m0.p4().P();
	    mTreetruthpt[0]     = m0.p4().Pt();
	    mTreetruthpx[0]     = m0.p4().px();
	    mTreetruthpy[0]     = m0.p4().py();
	    mTreetruthpz[0]     = m0.p4().pz();
	    mTreetrutheta[0]    = m0.p4().eta();
	    mTreetruthphi[0]    = m0.p4().phi();
	    mTreetruthm[0]      = m0.p4().M();
	    
	    mTreetruthpdgid[1]  = m1.pdgId();
	    mTreetruthbvtxid[1] = 0;
	    mTreetruthevtxid[1] = 1;
	    mTreetruthE[1]      = m1.p4().E();
	    mTreetruthEt[1]     = m1.p4().Et();
	    mTreetruthp[1]      = m1.p4().P();
	    mTreetruthpt[1]     = m1.p4().Pt();
	    mTreetruthpx[1]     = m1.p4().px();
	    mTreetruthpy[1]     = m1.p4().py();
	    mTreetruthpz[1]     = m1.p4().pz();
	    mTreetrutheta[1]    = m1.p4().eta();
	    mTreetruthphi[1]    = m1.p4().phi();
	    mTreetruthm[1]      = m1.p4().M();
	    
	    truth.push_back(p.mother(0));
	    truth.push_back(p.mother(1));
	    
	    for (unsigned int j=0; j<m0.numberOfDaughters(); j++) {
	      part1.push_back(m0.daughter(j));
	      bvtx1.push_back(1);
	    }
	    ntruth = 2;
	    break;
	  }
	}
      }

      // (2) Reconstruct decay chains
      if (part1.size()>0) {
	while (part1.size()>0) {
	
	  for (unsigned int d=0; d<part1.size(); d++) {

	    //	    cout << (*part1[d]).pdgId() << "  " << (*part1[d]).status() << endl;

	    //	    if (!isStable((*part1[d]).pdgId())) {
	    if (isDecaying((*part1[d]).pdgId())) {

	      for (unsigned int j=0; j<(*part1[d]).numberOfDaughters(); j++) {
		//const reco::Candidate &t = *((*part1[d]).daughter(j));

		//		if (abs((*part1[d]).pdgId()) == 15) { // tau
		//		  cout << " found tau with daughter " << t.pdgId() << " and status " << t.status() << endl;
		//		}
		//if (t.status()!=3 && abs((*part1[d]).pdgId()) != 15) continue;

		part2.push_back((*part1[d]).daughter(j));
		bvtx2.push_back(vtxid);
	      }
	      evtxid = vtxid;
	      vtxid++;
	    }
	    else
	      evtxid = -1;
	  
	    mTreetruthpdgid[ntruth]  = (*part1[d]).pdgId();
	    mTreetruthbvtxid[ntruth] = bvtx1[d];
	    mTreetruthevtxid[ntruth] = evtxid;
	    mTreetruthE[ntruth]      = (*part1[d]).p4().E();
	    mTreetruthEt[ntruth]     = (*part1[d]).p4().Et();
	    mTreetruthp[ntruth]      = (*part1[d]).p4().P();
	    mTreetruthpt[ntruth]     = (*part1[d]).p4().Pt();
	    mTreetruthpx[ntruth]     = (*part1[d]).p4().px();
	    mTreetruthpy[ntruth]     = (*part1[d]).p4().py();
	    mTreetruthpz[ntruth]     = (*part1[d]).p4().pz();
	    mTreetrutheta[ntruth]    = (*part1[d]).p4().eta();
	    mTreetruthphi[ntruth]    = (*part1[d]).p4().phi();
	    mTreetruthm[ntruth]      = (*part1[d]).p4().M();

	    truth.push_back(&(*part1[d]));

	    ntruth++;
	    if (ntruth == 100) break;
	  }
	  part1.clear();
	  bvtx1.clear();
	  part1 = part2;
	  bvtx1 = bvtx2;
	  part2.clear();
	  bvtx2.clear();
	  if (ntruth == 100) break;
	}
      }
      mTreeNtruth = ntruth;
    }
  }  
  
  // is the truthl stemming from truth ?
  //  -> match even if bremsstrahlung 

  for (unsigned int t=0; t<truthl.size(); t++) {

    std::vector<const reco::Candidate*> part;
    std::vector<const reco::Candidate*> parth;
    part.push_back(truthl[t]);

    while (part.size()>0) {
      if ((*part[0]).numberOfMothers()==0) break;
      parth.clear();
      for (unsigned int m=0; m<(*part[0]).numberOfMothers(); m++) {
	const reco::Candidate& p2 = *((*part[0]).mother(m));
	
	if (p2.pdgId() == (*truthl[t]).pdgId()) {
	  parth.push_back(&p2);
	  break;
	}
      }
      part.clear();
      part = parth;
    }

    mTreetruthlori[t] = -1;
    for (unsigned int j=0; j<truth.size(); j++) {
      
      if (fabs((*truth[j]).p4().E() - (*part[0]).p4().E()) < 1e-5 &&
	  (*part[0]).pdgId() == (*truth[j]).pdgId()) {
	mTreetruthlori[t] = j;
	break;
      }      
    }
  }

  // Truth jets
  mTreeNtruthjet = 0;

  if (is_MC) {
    
    Handle<reco::GenJetCollection> genjet;
    iEvent.getByLabel(genJetTag_, genjet);
    
    std::vector<reco::GenJet> genjets;
    
    if ( !genjet.isValid() )
      edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenJet found for InputTag " << genJetTag_;
    else {
      
      mTreeNtruthjet = genjet->size();
      
      for (int i=0; i<mTreeNtruthjet; i++) {
        genjets.push_back((*genjet)[i]);
      }
      
      sort(genjets.begin(), genjets.end(), ptcomp_genjet);  
      
      if ( mTreeNtruthjet > 100 ) mTreeNtruthjet = 100;
      
      for (int i=0; i<mTreeNtruthjet; ++ i ) {
	
	mTreetruthJetP[i]    = genjets[i].p();
	mTreetruthJetPt[i]   = genjets[i].pt();
	mTreetruthJetE[i]    = genjets[i].energy();
	mTreetruthJetEt[i]   = genjets[i].et();
	mTreetruthJetPx[i]   = genjets[i].momentum().X();
	mTreetruthJetPy[i]   = genjets[i].momentum().Y();
	mTreetruthJetPz[i]   = genjets[i].momentum().Z();
	mTreetruthJetEta[i]  = genjets[i].eta();
	mTreetruthJetPhi[i]  = genjets[i].phi();
      }
    }
  }
  // Get Hybrid SuperClusters Barrel
  mTreeNSC = 0;

  edm::Handle<reco::SuperClusterCollection> SuperClusterHandle;
  iEvent.getByLabel("correctedHybridSuperClusters", SuperClusterHandle);

  if (!SuperClusterHandle.isValid()) 
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::SuperClusterCollection found for InputTag correctedHybridSuperClusters";
  else {
    const reco::SuperClusterCollection *scCollection = SuperClusterHandle.product();
    
    for(reco::SuperClusterCollection::const_iterator scIt = scCollection->begin();   
	scIt != scCollection->end(); scIt++){

      mTreeSCE[mTreeNSC]   = scIt->energy();
      mTreeSCPhi[mTreeNSC] = scIt->phi();
      mTreeSCEta[mTreeNSC] = scIt->eta();
      mTreeNSC++;
      if (mTreeNSC==200) break;
    }
  }
  
  // Get Island SuperClusters Endcap
  
  edm::Handle<reco::SuperClusterCollection> SuperClusterHandle1;
  iEvent.getByLabel( "correctedMulti5x5SuperClustersWithPreshower", SuperClusterHandle1);

  if (!SuperClusterHandle1.isValid()) 
    edm::LogWarning("SusyACSkimAnalysis") <<  "No reco::SuperClusterCollection found for InputTag correctedMulti5x5SuperClustersWithPreshower";
  else {
    if (mTreeNSC<200) {
      const reco::SuperClusterCollection *scCollection1 = SuperClusterHandle1.product();
      
      for(reco::SuperClusterCollection::const_iterator scIt1 = scCollection1->begin();   
	  scIt1 != scCollection1->end(); scIt1++) {
	
	mTreeSCE[mTreeNSC]   = scIt1->energy();
	mTreeSCPhi[mTreeNSC] = scIt1->phi();
	mTreeSCEta[mTreeNSC] = scIt1->eta();
	mTreeNSC++;
	if (mTreeNSC==200) break;
      }
    }
  }  

  // match SC to truthl
  for(int l=0; l<mTreeNSC; l++) {
    mTreeSCTruth[l] = -1;
    for (int k=0; k<mTreeNtruthl; k++) {
      if (abs(mTreetruthlpdgid[k]) == 11) {
	if (DeltaR(mTreetruthleta[k],mTreeSCEta[l],mTreetruthlphi[k],mTreeSCPhi[l]) < 0.2 ) {
	  mTreeSCTruth[l] = k;
	  break;
	} 
      } 
    } 
  }
  // Electrons
  mTreeNele = 0;

  edm::Handle< std::vector<pat::Electron> > elecHandle;
  iEvent.getByLabel(elecTag_, elecHandle);

  std::vector<pat::Electron> eles;

  if ( !elecHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Electron results found for InputTag " << elecTag_;
  else {
    
    int countele = 0;
    
    mTreeNele= elecHandle->size();

    for (int i=0; i<mTreeNele; i++) {
      eles.push_back((*elecHandle)[i]);
    }

    sort(eles.begin(), eles.end(), ptcomp_ele);  

    if ( mTreeNele > 100 ) mTreeNele = 100;

    // Electron ID :
    //  0  -  eidLoose
    //  1  -  eidTight
    //  2  -  eidRobustLoose
    //  3  -  eidRobustTight
    //  4  -  eidRobustHighEnergy
    
    for (int i=0; i<mTreeNele; i++) {

      if (eles[i].pt() > elept_ && fabs(eles[i].eta()) < eleeta_) 
	cele_++;

      const std::vector<IdPair> &  electronIDs_ = eles[i].electronIDs();

      for (std::vector<IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
	if      (it->first == "eidLoose")             mTreeEleID[countele][0] = (int)it->second;
	else if (it->first == "eidTight")             mTreeEleID[countele][1] = (int)it->second;
	else if (it->first == "eidRobustLoose")       mTreeEleID[countele][2] = (int)it->second;
	else if (it->first == "eidRobustTight")       mTreeEleID[countele][3] = (int)it->second;
	else if (it->first == "eidRobustHighEnergy")  mTreeEleID[countele][4] = (int)it->second;
	else
	  edm::LogWarning("SusyACSkimAnalysis") << "Unknown ElectronID : " << it->first;
      }

      mTreeEleTruth[countele] = -1;
      const reco::GenParticle * gl = eles[i].genLepton();

      if (gl) {
	for (int k=0; k<mTreeNtruthl; k++) {
	  if (mTreetruthlpdgid[k] == gl->pdgId() &&
	      fabs(mTreetruthlE[k]   - gl->energy()) < 1e-5 &&
	      fabs(mTreetruthleta[k] - gl->eta()) < 1e-5 &&
	      fabs(mTreetruthlphi[k] - gl->phi()) < 1e-5) {
	    
	    mTreeEleTruth[countele] = k;
	    break;
	  }
	} 
      }

      mTreeNeletrign[countele] = 0;
      for (int k=0; k<mTreeNtrig; k++) {
	if (mTreeNeletrign[countele]<500 &&
	    DeltaPhi(eles[i].phi(), mTreetrigphi[k])<0.2 &&
	    fabs(eles[i].eta()-mTreetrigeta[k])<0.2) {
	  mTreeEletrig[countele][mTreeNeletrign[countele]] = k;
// 	  	  cout << " ELE " << eles[i].eta() << " TRIG " << mTreetrigeta[k] << "  " << unpack(mTreetrigname[k]) << endl;
	  mTreeNeletrign[countele]++;
	}
      }

      // Match reco electrons to SC
      mTreeEleSC[countele] = -1;
      for (int k=0; k<mTreeNSC; k++) {
	if (fabs(mTreeSCEta[k] - eles[i].superCluster().get()->eta()) < 1e-2 &&
	    DeltaPhi(mTreeSCPhi[k], eles[i].superCluster().get()->phi()) < 1e-2 ) {
	  mTreeEleSC[countele] = k;
	  break;
	}
      } 

      /* DEBUG
      cout << " Electron ID name : ";
      for (std::vector<IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
        cout << " - " << it->first << " [ " << it->second << " ] ";
      }
      cout << endl;
      */

      mTreeEleP[countele]          = eles[i].p();
      mTreeElePt[countele]         = eles[i].pt();
      mTreeEleE[countele]          = eles[i].energy();
      mTreeEleEt[countele]         = eles[i].et();
      mTreeElePx[countele]         = eles[i].momentum().X();
      mTreeElePy[countele]         = eles[i].momentum().Y();
      mTreeElePz[countele]         = eles[i].momentum().Z();
      mTreeEleEta[countele]        = eles[i].eta();
      mTreeElePhi[countele]        = eles[i].phi();
      mTreeEleCharge[countele]     = eles[i].charge();
      mTreeEleCaloEt[countele]     = eles[i].caloEnergy()*sin(eles[i].p4().theta());
      mTreeEleSCEta[countele]      = eles[i].caloPosition().eta();

      mTreeEleHCalOverEm[countele]          = eles[i].hadronicOverEm(); 
      mTreeEleDr03TkSumPt[countele]         = eles[i].dr03TkSumPt(); 
      mTreeEleDr04HCalTowerSumEt[countele]  = eles[i].dr04HcalTowerSumEt();
      mTreeEleDr03HCalTowerSumEt[countele]  = eles[i].dr03HcalTowerSumEt();
      mTreeEleDr04ECalRecHitSumEt[countele] = eles[i].dr04EcalRecHitSumEt(); 
      mTreeEleDr03ECalRecHitSumEt[countele] = eles[i].dr03EcalRecHitSumEt(); 

      // weighted cluster rms along eta and inside 5x5
      mTreeEleSigmaIetaIeta[countele] = eles[i].sigmaIetaIeta(); 

      // the supercluster eta - track eta position at calo extrapolated from innermost track state
      mTreeEleDeltaEtaSuperClusterTrackAtVtx[countele] = eles[i].deltaEtaSuperClusterTrackAtVtx(); 

      // the supercluster phi - track phi position at calo extrapolated from the innermost track state
      mTreeEleDeltaPhiSuperClusterTrackAtVtx[countele] = eles[i].deltaPhiSuperClusterTrackAtVtx(); 

      mTreeElefbrem[countele] = eles[i].fbrem();

      // iso for HEEP cuts
      mTreeEledr03HcalDepth1[countele] = eles[i].dr03HcalDepth1TowerSumEt();
      mTreeEledr03HcalDepth2[countele] = eles[i].dr03HcalDepth2TowerSumEt();
      mTreeElee2x5Max[countele]        = eles[i].e2x5Max();
      mTreeElee5x5[countele]           = eles[i].e5x5();
      mTreeElee1x5[countele]           = eles[i].e1x5();

      if (eles[i].ecalDrivenSeed()==true) 
	mTreeEleisECal[countele] = 1;
      else 
	mTreeEleisECal[countele] = 0;


      if (eles[i].trackerDrivenSeed()==true) 
	mTreeEleisTracker[countele] = 1;
      else 
	mTreeEleisTracker[countele] = 0;

      // Check for conversions
      reco::GsfElectron el = eles[i];
      
      ConversionFinder convFinder;
      ConversionInfo convInfo = convFinder.getConversionInfo(el, tkRef, bfield);
      
      mTreeEleConvdist[countele] = convInfo.dist();
      mTreeEleConvdcot[countele] = convInfo.dcot();
      mTreeEleConvr[countele]    = convInfo.radiusOfConversion();

      if (eles[i].gsfTrack().isNonnull()) {

	mTreeEleHits[countele]  = eles[i].gsfTrack().get()->numberOfValidHits();
  //missing Hits for WPXX
  mTreeElenumberOfHits[countele] = eles[i].gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	mTreeEled0vtx[countele] = (-1.)* eles[i].gsfTrack().get()->dxy(vtxPoint);
	mTreeElesd0[countele]   = eles[i].gsfTrack().get()->d0Error();

	mTreeEled0bs[countele]  = fabs(eles[i].gsfTrack()->dxy(bsPoint));
	// hit pattern used for expected crossed layers after the last trajectory's hit 
	mTreeEleTrkExpHitsInner[countele] = eles[i].gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 
	// has valid hit in Pixel layer 1
	mTreeEleValidHitFirstPxlB[countele] = eles[i].gsfTrack()->hitPattern().hasValidHitInFirstPixelBarrel(); 

	if (eles[i].gsfTrack().get()->ndof() > 0)
	  mTreeEleTrkChiNorm[countele] = eles[i].gsfTrack().get()->chi2()/ eles[i].gsfTrack().get()->ndof();
	else  mTreeEleTrkChiNorm[countele] = 999.;
      }
      else {
	mTreeEleHits[countele]              = -1;
	mTreeEleTrkChiNorm[countele]        = 999.;
 	mTreeEled0vtx[countele]             = 999.;
	mTreeElesd0[countele]               = 999.;
 	mTreeEled0bs[countele]              = 999.;
	mTreeEleTrkExpHitsInner[countele]   = -1;
	mTreeEleValidHitFirstPxlB[countele] = -1;
      }

      // CM Trigger not working
      /*
      const std::vector<pat::TriggerPrimitive> & eletrig = eles[i].triggerMatches();
      cout << "eletrig " << eletrig.size() << endl;
      for(unsigned int j=0; j<eletrig.size(); ++j) {
	cout << "triggerPrimitive["<< j << "] :  filterName = " << eletrig[j].filterName() 
	     << " , triggerObjId = "<< eletrig[j].triggerObjectId() << endl;
      }
      */
      countele++;
    }
    mTreeNele = countele;
  }
  // PF Electrons
  mTreeNPFEle = 0;

  edm::Handle< std::vector<pat::Electron> > PFelecHandle;
  iEvent.getByLabel(PFelecTag_, PFelecHandle);

  std::vector<pat::Electron> PFeles;

  if ( !PFelecHandle.isValid() )
    edm::LogWarning("SusyACSkimAnalysis") << "No PF Electron results found for InputTag " << PFelecTag_;
  else {

    mTreeNPFEle= PFelecHandle->size();

    for (int i=0; i<mTreeNPFEle; i++) {
      PFeles.push_back((*PFelecHandle)[i]);
    }

    sort(PFeles.begin(), PFeles.end(), ptcomp_ele);

    if ( mTreeNPFEle > 100 ) mTreeNPFEle = 100;

    int countPFele = 0;

    for (int i=0; i<mTreeNPFEle; i++) {
     
      if (PFeles[i].pt() > pfelept_ && fabs(PFeles[i].eta()) < pfeleeta_) 
	cpfele_++;

      mTreePFEleP[countPFele]          = PFeles[i].p();
      mTreePFElePt[countPFele]         = PFeles[i].pt();
      mTreePFEleE[countPFele]          = PFeles[i].energy();
      mTreePFEleEt[countPFele]         = PFeles[i].et();
      mTreePFElePx[countPFele]         = PFeles[i].momentum().X();
      mTreePFElePy[countPFele]         = PFeles[i].momentum().Y();
      mTreePFElePz[countPFele]         = PFeles[i].momentum().Z();
      mTreePFEleEta[countPFele]        = PFeles[i].eta();
      mTreePFElePhi[countPFele]        = PFeles[i].phi();
      mTreePFEleCharge[countPFele]     = PFeles[i].charge();

      // truth matching
      mTreePFEleTruth[countPFele] = -1;
      const reco::GenParticle * PFgl = PFeles[i].genLepton();
      
      if (PFgl) {
	for (int k=0; k<mTreeNtruthl; k++) {
	  if (mTreetruthlpdgid[k] == PFgl->pdgId() &&
	      fabs(mTreetruthlE[k]   - PFgl->energy()) < 1e-5 &&
	      fabs(mTreetruthleta[k] - PFgl->eta()) < 1e-5 &&
	      fabs(mTreetruthlphi[k] - PFgl->phi()) < 1e-5) {
	    
	    mTreePFEleTruth[countPFele] = k;
	    break;
	  }
	}
      }

      mTreeNPFEletrign[countPFele] = 0;
      for (int k=0; k<mTreeNtrig; k++) {
	if (mTreeNPFEletrign[countPFele]<500 &&
	    DeltaPhi(PFeles[i].phi(), mTreetrigphi[k])<0.2 &&
	    fabs(PFeles[i].eta()-mTreetrigeta[k])<0.2) {
	  mTreePFEletrig[countPFele][mTreeNPFEletrign[countPFele]] = k;
	  //              cout << " ELE " << eles[i].eta() << " TRIG " << mTreetrigeta[k] << "  " << unpack(mTreetrigname[k]) << endl;
	  mTreeNPFEletrign[countPFele]++;
	}
      }

      mTreePFEleSC[countPFele] = -1;
      for (int k=0; k<mTreeNSC; k++) {
	if (fabs(mTreeSCEta[k] - PFeles[i].superCluster().get()->eta()) < 1e-2 &&
	    DeltaPhi(mTreeSCPhi[k], PFeles[i].superCluster().get()->phi()) < 1e-2 ) {
	  mTreePFEleSC[countPFele] = k;
	  break;
	}
      }
      
      
      countPFele++;
    }
    mTreeNPFEle = countPFele;
  }

  // Muons
  mTreeNmuo = 0;

  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonTag_, muonHandle);

  std::vector<pat::Muon> muons;

  // prepare for cocktail muons, get refit types
  Handle <reco::TrackToTrackMap> tevMap1;
  Handle <reco::TrackToTrackMap> tevMap2;
  Handle <reco::TrackToTrackMap> tevMap3;
  iEvent.getByLabel("tevMuons", "default",  tevMap1);
  iEvent.getByLabel("tevMuons", "firstHit", tevMap2);
  iEvent.getByLabel("tevMuons", "picky",    tevMap3);

  if ( !muonHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Muon results found for InputTag " << muonTag_;
  else {
    
    int countmuo = 0;

    mTreeNmuo= muonHandle->size();
    
    for (int i=0; i<mTreeNmuo; i++) {
      muons.push_back((*muonHandle)[i]);
    }

    sort(muons.begin(), muons.end(), ptcomp_muo);  

    if ( mTreeNmuo > 100 ) mTreeNmuo = 100;
    
    for (int i=0; i<mTreeNmuo; i++) {

      if (muons[i].pt() > muopt_ && fabs(muons[i].eta()) < muoeta_ ) 
	cmuo_++;

      mTreeMuoTruth[countmuo] = -1;
      const reco::GenParticle * gl = muons[i].genLepton();
      
      if (gl) {
	for (int k=0; k<mTreeNtruthl; k++) {
	  if (mTreetruthlpdgid[k] == gl->pdgId() &&
	      fabs(mTreetruthlE[k]   - gl->energy()) < 1e-5 &&
	      fabs(mTreetruthleta[k] - gl->eta()) < 1e-5 &&
	      fabs(mTreetruthlphi[k] - gl->phi()) < 1e-5) {
	    
	    mTreeMuoTruth[countmuo] = k;
	    break;
	  }
	}
      }
      
      mTreeNmuotrign[countmuo] = 0;
      for (int k=0; k<mTreeNtrig; k++) {
	if (mTreeNmuotrign[countmuo]<500 && 
	    DeltaPhi(muons[i].phi(), mTreetrigphi[k])<0.1 &&
	    fabs(muons[i].eta()-mTreetrigeta[k])<0.1) {
	  mTreeMuotrig[countmuo][mTreeNmuotrign[countmuo]] = k;
	  //	  cout << " MUO " << muons[i].eta() << " TRIG " << mTreetrigeta[k] << "  " << unpack(mTreetrigname[k]) << endl;
	  mTreeNmuotrign[countmuo]++;
	}
      }

      mTreeMuoTrkIsoDep[countmuo]  = -1.;
      mTreeMuoECalIsoDep[countmuo] = -1.;
      mTreeMuoHCalIsoDep[countmuo] = -1.;
 
      for (int j=0; j<24; j++) {
	mTreeMuoID[countmuo][j] = muons[i].muonID(ACmuonID[j].Data());
      }

      mTreeMuoP[countmuo]          = muons[i].p();
      mTreeMuoPt[countmuo]         = muons[i].pt();
      mTreeMuoE[countmuo]          = muons[i].energy();
      mTreeMuoEt[countmuo]         = muons[i].et();
      mTreeMuoPx[countmuo]         = muons[i].momentum().X();
      mTreeMuoPy[countmuo]         = muons[i].momentum().Y();
      mTreeMuoPz[countmuo]         = muons[i].momentum().Z();
      mTreeMuoEta[countmuo]        = muons[i].eta();
      mTreeMuoPhi[countmuo]        = muons[i].phi();
      mTreeMuoRelTrkIso[countmuo]  = ( muons[i].trackIso()+muons[i].et() )/muons[i].et();
      mTreeMuoTrkIso[countmuo]     = muons[i].trackIso();
      mTreeMuoCharge[countmuo]     = muons[i].charge();
      mTreeMuoECalIso[countmuo]    = muons[i].ecalIso();
      mTreeMuoHCalIso[countmuo]    = muons[i].hcalIso() ;
      mTreeMuoAllIso[countmuo]     = muons[i].caloIso() ;
      mTreeMuoGood[countmuo]       = muons[i].isGood("GlobalMuonPromptTight");

      if (muons[i].trackIsoDeposit()) 
	mTreeMuoTrkIsoDep[countmuo]  = muons[i].trackIsoDeposit()->candEnergy();
      if (muons[i].ecalIsoDeposit())     
	mTreeMuoECalIsoDep[countmuo] = muons[i].ecalIsoDeposit()->candEnergy();
      if (muons[i].hcalIsoDeposit())     
	mTreeMuoHCalIsoDep[countmuo] = muons[i].hcalIsoDeposit()->candEnergy();

      mTreeMuocalocomp[countmuo]   = muons[i].caloCompatibility();
      mTreeMuocaltowe[countmuo]    = muons[i].calEnergy().tower;

      mTreeMuoChambersMatched[countmuo] = muons[i].numberOfMatches();

      // combined muon
      if ( muons[i].combinedMuon().isNonnull() ) {
  mTreeMuoValidFraction[countmuo]         = muons[i].combinedMuon().get()->validFraction();
	mTreeMuoHitsCm[countmuo]     = muons[i].combinedMuon().get()->numberOfValidHits();
	mTreeMuod0Cm[countmuo]       = (-1.)* muons[i].combinedMuon().get()->dxy(vtxPoint);
	mTreeMuosd0Cm[countmuo]      = muons[i].combinedMuon().get()->d0Error();
	mTreeMuod0bsCm[countmuo]     = (-1.)* muons[i].combinedMuon().get()->dxy(bsPoint); 
	mTreeMuod0OriginCm[countmuo] = (-1.)* muons[i].combinedMuon().get()->dxy();
	mTreeMuodzbsCm[countmuo]     = (-1.)* muons[i].combinedMuon().get()->dz(bsPoint);
	/*
	  /// dxy parameter. (This is the transverse impact parameter w.r.t. to (0,0,0) 
	  ONLY if refPoint is close to (0,0,0): see parametrization definition above for details). 
	  See also function dxy(myBeamSpot) below.
	  
	  double dxy() const { return ( - vx() * py() + vy() * px() ) / pt(); }

	  /// dxy parameter in perigee convention (d0 = - dxy)
	  double d0() const { return - dxy(); }

	  /// dxy parameter with respect to a user-given beamSpot (WARNING: 
	  this quantity can only be interpreted as a minimum transverse distance 
	  if beamSpot, if the beam spot is reasonably close to the refPoint, since 
	  linear approximations are involved). This is a good approximation for Tracker tracks.
	  double dxy(const Point& myBeamSpot) const { 
	  return ( - (vx()-myBeamSpot.x()) * py() + (vy()-myBeamSpot.y()) * px() ) / pt(); 
	  }
	*/
  
	mTreeMuoValidMuonHitsCm[countmuo]        = muons[i].combinedMuon().get()->hitPattern().numberOfValidMuonHits();
	mTreeMuoValidTrackerHitsCm[countmuo]     = muons[i].combinedMuon().get()->hitPattern().numberOfValidTrackerHits();
	mTreeMuoValidPixelHitsCm[countmuo]       = muons[i].combinedMuon().get()->hitPattern().numberOfValidPixelHits();
	mTreeMuoTrackerLayersMeasCm[countmuo]    = muons[i].combinedMuon().get()->hitPattern().trackerLayersWithMeasurement();
	mTreeMuoTrackerLayersNotMeasCm[countmuo] = muons[i].combinedMuon().get()->hitPattern().trackerLayersWithoutMeasurement();

	if (muons[i].combinedMuon().get()->ndof() > 0) {
	  mTreeMuoTrkChiNormCm[countmuo] = muons[i].combinedMuon().get() ->chi2()/ muons[i].combinedMuon().get()->ndof();
	  /*
	  cout << " offline pt " << muons[i].pt() << "  eta " << muons[i].eta() 
	       << "  phi " << muons[i].phi() << "  trackchi " 
	       << muons[i].combinedMuon().get() ->chi2()/ muons[i].combinedMuon().get()->ndof() 
	       << "  trackiso " << muons[i].trackIso() << endl;
	  */
	}
	else  
	  mTreeMuoTrkChiNormCm[countmuo] = 999.;
      }
      else {
	mTreeMuoTrkChiNormCm[countmuo] = 999.;
	mTreeMuoHitsCm[countmuo]       = -1;
	mTreeMuod0Cm[countmuo]         = 999.;
	mTreeMuosd0Cm[countmuo]        = 999.;
	mTreeMuod0OriginCm[countmuo]   = 999.;
	mTreeMuod0bsCm[countmuo]       = 999.;
	mTreeMuoValidMuonHitsCm[countmuo]        = -1;
	mTreeMuoValidTrackerHitsCm[countmuo]     = -1;
	mTreeMuoValidPixelHitsCm[countmuo]       = -1;
	mTreeMuoTrackerLayersMeasCm[countmuo]    = -1;
	mTreeMuoTrackerLayersNotMeasCm[countmuo] = -1;
      }

      // tracker muon
      if ( muons[i].innerTrack().isNonnull() ) {

	mTreeMuoHitsTk[countmuo] = muons[i].innerTrack().get()->numberOfValidHits();
	mTreeMuod0Tk[countmuo]   = (-1.)* muons[i].innerTrack().get()->dxy(vtxPoint);
	mTreeMuosd0Tk[countmuo]  = muons[i].innerTrack().get()->d0Error();

	if (muons[i].innerTrack().get()->ndof() > 0) {
	  mTreeMuoTrkChiNormTk[countmuo] = muons[i].innerTrack().get() ->chi2()/ muons[i].innerTrack().get()->ndof();
	  /*
	  cout << " offline pt " << muons[i].pt() << "  eta " << muons[i].eta() 
	       << "  phi " << muons[i].phi() << "  trackchi " 
	       << muons[i].innerTrack().get() ->chi2()/ muons[i].innerTrack().get()->ndof() 
	       << "  trackiso " << muons[i].trackIso() << endl;
	  */
	}
	else  
	  mTreeMuoTrkChiNormTk[countmuo] = 999.;
      }
      else {
	mTreeMuoTrkChiNormTk[countmuo] = 999.;
	mTreeMuoHitsTk[countmuo]       = -1;
	mTreeMuod0Tk[countmuo]         = 999.;
	mTreeMuosd0Tk[countmuo]        = 999.;
      }

      // cocktail muons
      if (muons[i].combinedMuon().isNonnull() && 
	  (!tevMap1.failedToGet() && !tevMap2.failedToGet() && !tevMap3.failedToGet())) {
	
	reco::TrackRef cocktailTrack = muon::tevOptimized(muons[i], *tevMap1, *tevMap2, *tevMap3);
	
	if ( cocktailTrack.isAvailable() ) {
	    mTreeMuoCocktailPt[countmuo] = cocktailTrack.get()->pt(); 
	    mTreeMuoCocktailEta[countmuo] = cocktailTrack.get()->eta();
	    mTreeMuoCocktailPhi[countmuo] = cocktailTrack.get()->phi(); 
	}
	else {
	  mTreeMuoCocktailPt[countmuo] = -1.; 
	  mTreeMuoCocktailEta[countmuo] = 999.;
	  mTreeMuoCocktailPhi[countmuo] = 999.;	
	}
      }
      else {
	mTreeMuoCocktailPt[countmuo] = -1.;
	mTreeMuoCocktailEta[countmuo] = 999.;
	mTreeMuoCocktailPhi[countmuo] = 999.;
      }

      countmuo++;
    }
    mTreeNmuo = countmuo;
  }
  
  // PF Muons
  mTreeNPFmuons = 0;

  edm::Handle< std::vector<pat::Muon> > PFmuonHandle;
  iEvent.getByLabel(PFmuonTag_, PFmuonHandle);
  std::vector<pat::Muon> PFmuons;
  if ( !PFmuonHandle.isValid() )
    edm::LogWarning("SusyACSkimAnalysis") << "No PF Muon results found for InputTag " << PFmuonTag_; 
  else{
      int countPFmuon = 0;
      mTreeNPFmuons= PFmuonHandle->size();
      
      for (int i=0; i<mTreeNPFmuons; i++) {
        PFmuons.push_back((*PFmuonHandle)[i]);
      }
      sort(PFmuons.begin(), PFmuons.end(), ptcomp_muo);  

      if ( mTreeNPFmuons > 100 ) mTreeNPFmuons = 100;
      
      for (int i=0; i<mTreeNPFmuons; i++) {
        //For the Moment no PFID in Pat
        //for (int j=0; j<24; j++) {
          //mTreePFMuoID[countPFmuon][j] = 0;//PFmuons[i].muonID(ACmuonID[j].Data());
        //}
        
        mTreePFMuonP[countPFmuon]          = PFmuons[i].p();
        mTreePFMuonPt[countPFmuon]         = PFmuons[i].pt();
        mTreePFMuonE[countPFmuon]          = PFmuons[i].energy();
        mTreePFMuonEt[countPFmuon]         = PFmuons[i].et();
        mTreePFMuonPx[countPFmuon]         = PFmuons[i].momentum().X();
        mTreePFMuonPy[countPFmuon]         = PFmuons[i].momentum().Y();
        mTreePFMuonPz[countPFmuon]         = PFmuons[i].momentum().Z();
        mTreePFMuonEta[countPFmuon]        = PFmuons[i].eta();
        mTreePFMuonPhi[countPFmuon]        = PFmuons[i].phi();
        mTreePFMuonTrkIso[countPFmuon]     = PFmuons[i].trackIso();
        mTreePFMuonCharge[countPFmuon]        = PFmuons[i].charge();
        mTreePFMuonECalIso[countPFmuon]       = PFmuons[i].ecalIso();
        mTreePFMuonHCalIso[countPFmuon]       = PFmuons[i].hcalIso();
        mTreePFMuonAllIso[countPFmuon]        = PFmuons[i].caloIso();
        mTreePFMuonParticleIso[countPFmuon]   =   PFmuons[i].userIsolation("pat::PfAllParticleIso");
        mTreePFMuonChadIso[countPFmuon]       =   PFmuons[i].userIsolation("pat::PfChargedHadronIso");
        mTreePFMuonNhadIso[countPFmuon]       =   PFmuons[i].userIsolation("pat::PfNeutralHadronIso");
        mTreePFMuonGamIso[countPFmuon]        =   PFmuons[i].userIsolation("pat::PfGammaIso");
        
        //mTreePFMuonParticleIsoDep[countPFmuon]   =   PFmuons[i].userIsoDeposit();
        //mTreePFMuonChadIsoDep[countPFmuon]       =   PFmuons[i].userIsoDeposit("pat::PfChargedHadronIso");
        //mTreePFMuonNhadIsoDep[countPFmuon]       =   PFmuons[i].userIsoDeposit("pat::PfNeutralHadronIso");
        //mTreePFMuonGamIsoDep[countPFmuon]        =   PFmuons[i].userIsoDeposit("pat::PfGammaIso");
        
      
        
        if ( PFmuons[i].combinedMuon().isNonnull() ) {

          mTreePFMuoHitsCm[countPFmuon]     = PFmuons[i].combinedMuon().get()->numberOfValidHits();
          mTreePFMuoValidFraction[countPFmuon]= PFmuons[i].combinedMuon().get()->validFraction();
          mTreePFMuod0Cm[countPFmuon]       = (-1.)* PFmuons[i].combinedMuon().get()->dxy(vtxPoint);
          mTreePFMuosd0Cm[countPFmuon]      = PFmuons[i].combinedMuon().get()->d0Error();
          mTreePFMuod0bsCm[countPFmuon]     = (-1.)* PFmuons[i].combinedMuon().get()->dxy(bsPoint); 
          mTreePFMuod0OriginCm[countPFmuon] = (-1.)* PFmuons[i].combinedMuon().get()->dxy();
          mTreePFMuodzbsCm[countPFmuon]     = (-1.)* PFmuons[i].combinedMuon().get()->dz(bsPoint);
          /*
            /// dxy parameter. (This is the transverse impact parameter w.r.t. to (0,0,0) 
            ONLY if refPoint is close to (0,0,0): see parametrization definition above for details). 
            See also function dxy(myBeamSpot) below.
            
            double dxy() const { return ( - vx() * py() + vy() * px() ) / pt(); }

            /// dxy parameter in perigee convention (d0 = - dxy)
            double d0() const { return - dxy(); }

            /// dxy parameter with respect to a user-given beamSpot (WARNING: 
            this quantity can only be interpreted as a minimum transverse distance 
            if beamSpot, if the beam spot is reasonably close to the refPoint, since 
            linear approximations are involved). This is a good approximation for Tracker tracks.
            double dxy(const Point& myBeamSpot) const { 
            return ( - (vx()-myBeamSpot.x()) * py() + (vy()-myBeamSpot.y()) * px() ) / pt(); 
            }
          */

          mTreePFMuoValidMuonHitsCm[countPFmuon]        = PFmuons[i].combinedMuon().get()->hitPattern().numberOfValidMuonHits();
          mTreePFMuoValidTrackerHitsCm[countPFmuon]     = PFmuons[i].combinedMuon().get()->hitPattern().numberOfValidTrackerHits();
          mTreePFMuoValidPixelHitsCm[countPFmuon]       = PFmuons[i].combinedMuon().get()->hitPattern().numberOfValidPixelHits();
          mTreePFMuoTrackerLayersMeasCm[countPFmuon]    = PFmuons[i].combinedMuon().get()->hitPattern().trackerLayersWithMeasurement();
          mTreePFMuoTrackerLayersNotMeasCm[countPFmuon] = PFmuons[i].combinedMuon().get()->hitPattern().trackerLayersWithoutMeasurement();

          if (PFmuons[i].combinedMuon().get()->ndof() > 0) {
            mTreePFMuoTrkChiNormCm[countPFmuon] = PFmuons[i].combinedMuon().get() ->chi2()/ PFmuons[i].combinedMuon().get()->ndof();
            /*
            cout << " offline pt " << PFmuons[i].pt() << "  eta " << PFmuons[i].eta() 
                 << "  phi " << PFmuons[i].phi() << "  trackchi " 
                 << PFmuons[i].combinedMuon().get() ->chi2()/ PFmuons[i].combinedMuon().get()->ndof() 
                 << "  trackiso " << PFmuons[i].trackIso() << endl;
            */
          }
          else  
            mTreePFMuoTrkChiNormCm[countPFmuon] = 999.;
        }
        else {
          mTreePFMuoValidFraction[countPFmuon]= -1;
          mTreePFMuoTrkChiNormCm[countPFmuon] = 999.;
          mTreePFMuoHitsCm[countPFmuon]       = -1;
          mTreePFMuod0Cm[countPFmuon]         = 999.;
          mTreePFMuosd0Cm[countPFmuon]        = 999.;
          mTreePFMuod0OriginCm[countPFmuon]   = 999.;
          mTreePFMuod0bsCm[countPFmuon]       = 999.;
          mTreePFMuoValidMuonHitsCm[countPFmuon]        = -1;
          mTreePFMuoValidTrackerHitsCm[countPFmuon]     = -1;
          mTreePFMuoValidPixelHitsCm[countPFmuon]       = -1;
          mTreePFMuoTrackerLayersMeasCm[countPFmuon]    = -1;
          mTreePFMuoTrackerLayersNotMeasCm[countPFmuon] = -1;
        }

              // tracker PFmuon
        if ( PFmuons[i].innerTrack().isNonnull() ) {

        mTreePFMuoHitsTk[countPFmuon] = PFmuons[i].innerTrack().get()->numberOfValidHits();
        mTreePFMuod0Tk[countPFmuon]   = (-1.)* PFmuons[i].innerTrack().get()->dxy(vtxPoint);
        mTreePFMuosd0Tk[countPFmuon]  = PFmuons[i].innerTrack().get()->d0Error();

        if (PFmuons[i].innerTrack().get()->ndof() > 0) {
          mTreePFMuoTrkChiNormTk[countPFmuon] = PFmuons[i].innerTrack().get() ->chi2()/ PFmuons[i].innerTrack().get()->ndof();
          /*
          cout << " offline pt " << PFmuons[i].pt() << "  eta " << PFmuons[i].eta() 
               << "  phi " << PFmuons[i].phi() << "  trackchi " 
               << PFmuons[i].innerTrack().get() ->chi2()/ PFmuons[i].innerTrack().get()->ndof() 
               << "  trackiso " << PFmuons[i].trackIso() << endl;
          */
        }
        else  
          mTreePFMuoTrkChiNormTk[countPFmuon] = 999.;
        }
        else {
          mTreePFMuoTrkChiNormTk[countPFmuon] = 999.;
          mTreePFMuoHitsTk[countPFmuon]       = -1;
          mTreePFMuod0Tk[countPFmuon]         = 999.;
          mTreePFMuosd0Tk[countPFmuon]        = 999.;
        }
        
        countPFmuon++;
        
      }
      mTreeNPFmuons=countPFmuon;
  }
  
  
  
  // Taus
  
  
  edm::Handle<std::vector<pat::Tau> > tauHandle;
  iEvent.getByLabel(tauSrc_, tauHandle);
  

  mTreeNtaus = 0;
  std::vector<pat::Tau> taus;
  
  if ( !tauHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Tau results found for InputTag"<<tauSrc_;
  else {
    
    mTreeNtaus = tauHandle->size();
    int counttau=0;
    
    for( int i=0;i<mTreeNtaus ; i++){
      taus.push_back((*tauHandle)[i]);
    }
    sort(taus.begin(), taus.end(), ptcomp_tau);  
    
    math::XYZPoint tauVtxPoint(0,0,0);
    if(mTreeNtaus >100) mTreeNtaus =100;
    
    
    for( int i=0;i<mTreeNtaus ; i++){
      for(int j=0; j<16; j++){
        mTreeTauID[counttau][j] = taus[i].tauID(ACtauID[j].Data());
      }
      mTreeTauP[counttau]          = taus[i].p();
      mTreeTauPt[counttau]         = taus[i].pt();
      mTreeTauE[counttau]          = taus[i].energy();
      mTreeTauEt[counttau]         = taus[i].et();
      mTreeTauPx[counttau]         = taus[i].momentum().X();
      mTreeTauPy[counttau]         = taus[i].momentum().Y();
      mTreeTauPz[counttau]         = taus[i].momentum().Z();
      mTreeTauEta[counttau]        = taus[i].eta();
      mTreeTauPhi[counttau]        = taus[i].phi();
      mTreeTauDecayMode[counttau]  = taus[i].decayMode();
      mTreeTauvx[counttau]         = taus[i].vx();
      mTreeTauvy[counttau]         = taus[i].vy();
      mTreeTauvz[counttau]         = taus[i].vz();
      tauVtxPoint                  = taus[i].leadPFChargedHadrCand()->vertex();
      mTreeTauvx2[counttau]         = tauVtxPoint.X();
      mTreeTauvy2[counttau]         = tauVtxPoint.Y();
      mTreeTauvz2[counttau]         = tauVtxPoint.Z();
      mTreeTauECalIso[counttau]    = taus[i].ecalIso();
      mTreeTauHCalIso[counttau]    = taus[i].hcalIso() ;
      mTreeTauAllIso[counttau]     = taus[i].caloIso() ;
      mTreeTauTrackIso[counttau]   = taus[i].trackIso() ;
      mTreeTauParticleIso[counttau]= taus[i].userIsolation("pat::PfAllParticleIso");
      mTreeTauChadIso[counttau]    = taus[i].userIsolation("pat::PfChargedHadronIso");
      mTreeTauNhadIso[counttau]    = taus[i].userIsolation("pat::PfNeutralHadronIso");
      mTreeTauGamIso[counttau]     = taus[i].userIsolation("pat::PfGammaIso");      
      
      counttau++;
    }
    mTreeNtaus=counttau;
  }
  
  
  
  
  
  // fatjets
  mTreeNfatjet = 0;

  for (int k=0; k<100; k++) {
    _jeta[k] = -100.;
    _jphi[k] = -100.;
  }

  if (do_fatjets) {

    edm::Handle< std::vector<pat::Jet> > jetHandle2;
    iEvent.getByLabel("patJetsBHS", jetHandle2);

    int mTreeNjet = 0;

    // this is required for JES corrected subjets
    if ( !jetHandle2.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No Jet results found for InputTag allLayer1JetsBHS";
    else {

      mTreeNjet = jetHandle2->size();
      
      if ( mTreeNjet > 100 ) mTreeNjet = 100;
      
      for (int i=0; i<mTreeNjet; i++) {
	
	_jeta[i] = (*jetHandle2)[i].eta();
	_jphi[i] = (*jetHandle2)[i].phi();
	
      }
    }

    // now associate fat jets and subjets
    edm::Handle<vector<reco::BasicJet> > subjets;
    iEvent.getByLabel("BoostedHiggsSubjets", subjets);

    if ( !subjets.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No Jet results found for InputTag BoostedHiggsSubjets";
    else {
      
      const vector<reco::BasicJet>* c_jets = subjets.product();
      
      if (c_jets->size()<100)
	mTreeNfatjet = c_jets->size();
      else
	mTreeNfatjet = 100;
	
      int cjet = 0;

      for(int j=0; j < mTreeNfatjet; j++){
	
	const reco::BasicJet & fat_jet = c_jets->at(j);
	
	mTreefatjetpt[cjet]   = fat_jet.pt();
	mTreefatjetpx[cjet]   = fat_jet.px();
	mTreefatjetpy[cjet]   = fat_jet.py();
	mTreefatjetpz[cjet]   = fat_jet.pz();
	mTreefatjete[cjet ]   = fat_jet.energy();
	mTreefatjeteta[cjet]  = fat_jet.eta();
	mTreefatjetphi[cjet]  = fat_jet.phi();
	mTreefatjetnsub[cjet] = fat_jet.numberOfDaughters();
	
	for (UInt_t k = 0; k < fat_jet.numberOfDaughters(); k++) {

	  mTreefatjetsubeta[cjet][k] = fat_jet.daughter(k)->eta();
	  mTreefatjetsubphi[cjet][k] = fat_jet.daughter(k)->phi();

	  // matching
	  for (int l=0; l<mTreeNjet; l++) {
	    if (fabs(mTreefatjetsubeta[cjet][k]-_jeta[l])<1e-7 && fabs(mTreefatjetsubphi[cjet][k]-_jphi[l])<1e-7) {
 
	      const pat::Jet& jet = (*jetHandle2)[l];

	      mTreefatjetsubpt[cjet][k]   = jet.pt();
	      mTreefatjetsubpx[cjet][k]   = jet.px();
	      mTreefatjetsubpy[cjet][k]   = jet.py();
	      mTreefatjetsubpz[cjet][k]   = jet.pz();
	      mTreefatjetsube[cjet][k]    = jet.energy();
	      mTreefatjetsubfem[cjet][k]  = jet.emEnergyFraction();
	      mTreefatjetsubfhad[cjet][k] = jet.energyFractionHadronic();
	      mTreefatjetsubbtag[cjet][k] = jet.bDiscriminator(btag_);
	      mTreefatjetsubn90[cjet][k]  = jet.jetID().n90Hits;
	      mTreefatjetsubfhpd[cjet][k] = jet.jetID().fHPD;
	      mTreefatjetsubfrbx[cjet][k] = jet.jetID().fRBX;
	      break;
	    }
	  }
	  
	  if (k==10) {
	    mTreefatjetnsub[cjet] = 10;
	    break;
	  }
	}
	
	cjet++;
	if (cjet == 100) break;
      }
    }
  }

  // Calo Jets
  mTreeNCalojet = 0;

  edm::Handle< std::vector<pat::Jet> > calojetHandle;
  iEvent.getByLabel(calojetTag_, calojetHandle);

  std::vector<pat::Jet> calojets;

  if ( !calojetHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Calo Jet results found for InputTag " << calojetTag_;
  else {
    
    int countjets = 0;
    htcutsum_=0;
    
    mTreeNCalojet = calojetHandle->size();

    for (int i=0; i<mTreeNCalojet; i++) {
      calojets.push_back((*calojetHandle)[i]);
      htcutsum_ +=(*calojetHandle)[i].et();
    }

    sort(calojets.begin(), calojets.end(), ptcomp_jet);  

    if ( mTreeNCalojet > 100 ) mTreeNCalojet = 100;
    
    for (int i=0; i<mTreeNCalojet; i++) {
      
      if (jetId(calojets[i])) mTreeCaloJetID[countjets] = 1;
      else                    mTreeCaloJetID[countjets] = 0;

      /*
      // result of all b tagging algorithms
      const vector<pair<string, float> > bvec = jet.getPairDiscri();

      for(unsigned int l=0; l!=bvec.size(); l++){
	cout << " b LABEL " << bvec[l].first << "  " << bvec[l].second << endl;
      }
      */

      // Jet Corrections applied via Python Configuration:
      //  L0 : raw
      //  L1 : off
      //  L2 : rel
      //  L3 : abs  <-- default
      //  L4 : emf
      //  L5 : had  with glu/uds/c/b
      //  L6 : ue   with glu/uds/c/b
      //  L7 : part with glu/uds/c/b
      
      if (calojets[i].pt() > calojetpt_ && fabs(calojets[i].eta()) < calojeteta_) 
	ccalojet_++;
    
      mTreeCaloJetP[countjets]     = calojets[i].p();
      mTreeCaloJetPt[countjets]    = calojets[i].pt();
      mTreeCaloJetE[countjets]     = calojets[i].energy();
      mTreeCaloJetEt[countjets]    = calojets[i].et();
      mTreeCaloJetPx[countjets]    = calojets[i].momentum().X();
      mTreeCaloJetPy[countjets]    = calojets[i].momentum().Y();
      mTreeCaloJetPz[countjets]    = calojets[i].momentum().Z();
      mTreeCaloJetEta[countjets]   = calojets[i].eta();
      mTreeCaloJetPhi[countjets]   = calojets[i].phi();

      mTreeCaloJetPtRaw[countjets] = calojets[i].correctedJet("Uncorrected").pt();

      mTreeCaloJetCharge[countjets]  = calojets[i].jetCharge();
      mTreeCaloJetn90[countjets]     = calojets[i].n90();
      mTreeCaloJetn90hits[countjets] = calojets[i].jetID().n90Hits;
      mTreeCaloJetfhpd[countjets]    = calojets[i].jetID().fHPD;
      mTreeCaloJetfrbx[countjets]    = calojets[i].jetID().fRBX;

      //      if (calojets[i].isCaloJet()) {
      mTreeCaloJetFem[countjets]   = calojets[i].emEnergyFraction();
      mTreeCaloJetFhad[countjets]  = calojets[i].energyFractionHadronic();
      mTreeCaloJetConst[countjets] = calojets[i].getCaloConstituents().size();
      //      }

      // default b tagger
      mTreeCaloJetBtag[countjets]    = calojets[i].bDiscriminator(btag_);

      mTreeCaloJetPart[countjets]    = calojets[i].partonFlavour();
      mTreeCaloJetTruth[countjets]   = -1;
      const reco::GenJet * gl = calojets[i].genJet();
      
      if (gl) {
	for (int k=0; k<mTreeNtruthjet; k++) {
	  if (fabs(mTreetruthJetE[k]   - gl->energy()) < 1e-5 &&
	      fabs(mTreetruthJetEta[k] - gl->eta()) < 1e-5 &&
	      fabs(mTreetruthJetPhi[k] - gl->phi()) < 1e-5) {
	    mTreeCaloJetTruth[countjets] = k;
	    break;
	  }
	}
      }

      countjets++;
    }
    mTreeNCalojet = countjets;
  }

  // Particle Flow Jets
  mTreeNPFjet = 0;

  edm::Handle< std::vector<pat::Jet> > pfjetHandle;
  iEvent.getByLabel(pfjetTag_, pfjetHandle);

  std::vector<pat::Jet> pfjets;

  if ( !pfjetHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No PF Jet results found for InputTag " << pfjetTag_;
  else {
    
    int countjets = 0;

    mTreeNPFjet = pfjetHandle->size();
    
    PFhtcutsum_=0;
    
    for (int i=0; i<mTreeNPFjet; i++) {
      pfjets.push_back((*pfjetHandle)[i]);
      PFhtcutsum_ +=(*pfjetHandle)[i].et();
    }

    sort(pfjets.begin(), pfjets.end(), ptcomp_jet);  
    
      
    
    if ( mTreeNPFjet > 100 ) mTreeNPFjet = 100;
    
    for (int i=0; i<mTreeNPFjet; i++) {
      
      if (pfjets[i].pt() > pfjetpt_ && fabs(pfjets[i].eta()) < pfjeteta_) 
	cpfjet_++;
    
      mTreePFJetP[countjets]     = pfjets[i].p();
      mTreePFJetPt[countjets]    = pfjets[i].pt();
      mTreePFJetE[countjets]     = pfjets[i].energy();
      mTreePFJetEt[countjets]    = pfjets[i].et();
      mTreePFJetPx[countjets]    = pfjets[i].momentum().X();
      mTreePFJetPy[countjets]    = pfjets[i].momentum().Y();
      mTreePFJetPz[countjets]    = pfjets[i].momentum().Z();
      mTreePFJetEta[countjets]   = pfjets[i].eta();
      mTreePFJetPhi[countjets]   = pfjets[i].phi();

      mTreePFJetPtRaw[countjets] = pfjets[i].correctedJet("Uncorrected").pt();

      mTreePFJetCharge[countjets] = pfjets[i].jetCharge();
      mTreePFJetn90[countjets]    = pfjets[i].n90();

      //      if (pfjets[i].isPFJet()) {
      mTreePFJetF[countjets][0]  = pfjets[i].chargedHadronEnergyFraction();
      mTreePFJetN[countjets][0]  = pfjets[i].chargedHadronMultiplicity();  
      mTreePFJetF[countjets][1]  = pfjets[i].neutralHadronEnergyFraction(); 
      mTreePFJetN[countjets][1]  = pfjets[i].neutralHadronMultiplicity(); 
      mTreePFJetF[countjets][2]  = pfjets[i].neutralEmEnergyFraction();
      mTreePFJetN[countjets][2]  = 0.;    
      mTreePFJetF[countjets][3]  = pfjets[i].chargedEmEnergyFraction(); 
      mTreePFJetN[countjets][3]  = 0.;   
      mTreePFJetF[countjets][4]  = pfjets[i].chargedMuEnergyFraction();
      mTreePFJetN[countjets][4]  = pfjets[i].muonMultiplicity();    
      // workaround: missing: pfjets[i].electronEnergyFraction();
      mTreePFJetF[countjets][5]  = pfjets[i].electronEnergy()/(pfjets[i].jecFactor(0)*pfjets[i].energy());
      mTreePFJetN[countjets][5]  = pfjets[i].electronMultiplicity();    
      mTreePFJetF[countjets][6]  = pfjets[i].photonEnergyFraction();
      mTreePFJetN[countjets][6]  = pfjets[i].photonMultiplicity();   
      
      mTreePFJetConst[countjets] = pfjets[i].getPFConstituents().size();
      //      }

      // default b tagger
      mTreePFJetBtag[countjets]    = pfjets[i].bDiscriminator(btag_);

      mTreePFJetPart[countjets]    = pfjets[i].partonFlavour();
      mTreePFJetTruth[countjets]   = -1;
      const reco::GenJet * gl = pfjets[i].genJet();
      
      if (gl) {
	for (int k=0; k<mTreeNtruthjet; k++) {
	  if (fabs(mTreetruthJetE[k]   - gl->energy()) < 1e-5 &&
	      fabs(mTreetruthJetEta[k] - gl->eta()) < 1e-5 &&
	      fabs(mTreetruthJetPhi[k] - gl->phi()) < 1e-5) {
	    mTreePFJetTruth[countjets] = k;
	    break;
	  }
	}
      }

      countjets++;
    }
    mTreeNPFjet = countjets;
  }






  // MET 
  edm::Handle< std::vector<pat::MET> > metHandle;
  iEvent.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTag_;
  if ( metHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
				      << metHandle->size() << " instead of 1";
  if ( metHandle.isValid() && metHandle->size()==1 ) {
    
    // in order to get uncorrected MET use
    // - metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrALL")) 
    // - metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrJES"))
    // - metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrMUON"))
    //  NOT metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrMAXN"))
 
    /*
    const reco::GenMET * gm = metHandle->front().genMET();
    cout << " genMET et " << gm->et() << " x " << gm->momentum().X() 
    	 << " y " << gm->momentum().Y() << endl;
    */

    /*
    cout << "MET " << metHandle->front().et()
	 << " uncorrALL " << metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrALL"))
	 << " uncorrJES " << metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrJES"))
	 << " uncorrMUON " << metHandle->front().uncorrectedPt(pat::uncorrectionTypeMET("uncorrMUON")) << endl;
    */
    int metn=0;
    mTreeMET[metn]         = metHandle->front().et();
    mTreeMEX[metn]         = metHandle->front().momentum().X();
    mTreeMEY[metn]         = metHandle->front().momentum().Y();
    mTreeSumET[metn]       = metHandle->front().sumEt();
    mTreeMETeta[metn]      = metHandle->front().eta();
    mTreeMETphi[metn]      = metHandle->front().phi();
    mTreeSumETSignif[metn] = metHandle->front().mEtSig();
    mTreeMETSignif[metn]   = -1;
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = metHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = metHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = metHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = metHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = metHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = metHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = metHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= metHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = metHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = metHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = metHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = metHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = metHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = metHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = metHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = metHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = metHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = -1; //metHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = -1; //metHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = -1; //metHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = -1; //metHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = -1; //metHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = -1; //metHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = -1; //metHandle->front().Type7EtFraction();
    
    
    
  }

  if (is_MC) {
    
    // genMET 
    edm::Handle< std::vector<reco::GenMET> > genmetHandle;
    iEvent.getByLabel("genMetCalo", genmetHandle);
    if ( !genmetHandle.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenMET results found for InputTag genMetCalo";
    if ( genmetHandle->size()!=1 ) 
      edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					    << genmetHandle->size() << " instead of 1";
    if ( genmetHandle.isValid() && genmetHandle->size()==1 ) {
      
      int metn=1;
      mTreeMET[metn]         = metHandle->front().et();
      mTreeMEX[metn]         = metHandle->front().momentum().X();
      mTreeMEY[metn]         = metHandle->front().momentum().Y();
      mTreeSumET[metn]       = metHandle->front().sumEt();
      mTreeMETeta[metn]      = metHandle->front().eta();
      mTreeMETphi[metn]      = metHandle->front().phi();
      mTreeSumETSignif[metn] = metHandle->front().mEtSig();
      mTreeMETSignif[metn]   = -1;
      
      //calo spesific 
      mTreeMETCaloMETInmHF[metn]     = -1; //metHandle->front().CaloMETInmHF();
      mTreeMETCaloMETInpHF[metn]     = -1; //metHandle->front().CaloMETInpHF();
      mTreeMETCaloMETPhiInmHF[metn]  = -1; //metHandle->front().CaloMETPhiInmHF();
      mTreeMETCaloMETPhiInpHF[metn]  = -1; //metHandle->front().CaloMETPhiInpHF();
      mTreeMETCaloSETInmHF[metn]     = -1; //metHandle->front().CaloSETInmHF();
      mTreeMETCaloSETInpHF[metn]     = -1; //metHandle->front().CaloSETInpHF();
      mTreeMETemEtFraction[metn]      = -1; //metHandle->front().emEtFraction();
      mTreeMETetFractionHadronic[metn]= -1; //metHandle->front().etFractionHadronic();
      mTreeMETmaxEtInEmTowers[metn]   = -1; //metHandle->front().maxEtInEmTowers();
      mTreeMETmaxEtInHadTowers[metn]  = -1; //metHandle->front().maxEtInHadTowers();
      mTreeMETemEtInHF[metn]          = -1; //metHandle->front().emEtInHF();
      mTreeMETemEtInEE[metn]          = -1; //metHandle->front().emEtInEE();
      mTreeMETemEtInEB[metn]          = -1; //metHandle->front().emEtInEB();
      mTreeMEThadEtInHF[metn]         = -1; //metHandle->front().hadEtInHF();
      mTreeMEThadEtInHE[metn]         = -1; //metHandle->front().hadEtInHE();
      mTreeMEThadEtInHO[metn]         = -1; //metHandle->front().hadEtInHO();
      mTreeMEThadEtInHB[metn]         = -1; //metHandle->front().hadEtInHB();
      
      //pf spesific
      mTreeMETChargedEMEtFraction[metn]  = -1; //metHandle->front().ChargedEMEtFraction();
      mTreeMETChargedHadEtFraction[metn] = -1; //metHandle->front().ChargedHadEtFraction();
      mTreeMETMuonEtFraction[metn]       = -1; //metHandle->front().MuonEtFraction();
      mTreeMETNeutralEMFraction[metn]  = -1; //metHandle->front().NeutralEMFraction();
      mTreeMETNeutralHadEtFraction[metn] = -1; //metHandle->front().NeutralHadEtFraction();
      mTreeMETType6EtFraction[metn]      = -1; //metHandle->front().Type6EtFraction();
      mTreeMETType7EtFraction[metn]      = -1; //metHandle->front().Type7EtFraction();
    }
    
    edm::Handle< std::vector<reco::GenMET> > genmet2Handle;
    iEvent.getByLabel("genMetCaloAndNonPrompt", genmet2Handle);
    if ( !genmet2Handle.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenMET results found for InputTag genMetCaloAndNonPrompt";
    if ( genmet2Handle->size()!=1 ) 
      edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					    << genmet2Handle->size() << " instead of 1";
    if ( genmet2Handle.isValid() && genmet2Handle->size()==1 ) {
      
      int metn=2;
      mTreeMET[metn]         = metHandle->front().et();
      mTreeMEX[metn]         = metHandle->front().momentum().X();
      mTreeMEY[metn]         = metHandle->front().momentum().Y();
      mTreeSumET[metn]       = metHandle->front().sumEt();
      mTreeMETeta[metn]      = metHandle->front().eta();
      mTreeMETphi[metn]      = metHandle->front().phi();
      mTreeSumETSignif[metn] = metHandle->front().mEtSig();
      mTreeMETSignif[metn]   = -1;
      
      //calo spesific 
      mTreeMETCaloMETInmHF[metn]     = -1; //metHandle->front().CaloMETInmHF();
      mTreeMETCaloMETInpHF[metn]     = -1; //metHandle->front().CaloMETInpHF();
      mTreeMETCaloMETPhiInmHF[metn]  = -1; //metHandle->front().CaloMETPhiInmHF();
      mTreeMETCaloMETPhiInpHF[metn]  = -1; //metHandle->front().CaloMETPhiInpHF();
      mTreeMETCaloSETInmHF[metn]     = -1; //metHandle->front().CaloSETInmHF();
      mTreeMETCaloSETInpHF[metn]     = -1; //metHandle->front().CaloSETInpHF();
      mTreeMETemEtFraction[metn]      = -1; //metHandle->front().emEtFraction();
      mTreeMETetFractionHadronic[metn]= -1; //metHandle->front().etFractionHadronic();
      mTreeMETmaxEtInEmTowers[metn]   = -1; //metHandle->front().maxEtInEmTowers();
      mTreeMETmaxEtInHadTowers[metn]  = -1; //metHandle->front().maxEtInHadTowers();
      mTreeMETemEtInHF[metn]          = -1; //metHandle->front().emEtInHF();
      mTreeMETemEtInEE[metn]          = -1; //metHandle->front().emEtInEE();
      mTreeMETemEtInEB[metn]          = -1; //metHandle->front().emEtInEB();
      mTreeMEThadEtInHF[metn]         = -1; //metHandle->front().hadEtInHF();
      mTreeMEThadEtInHE[metn]         = -1; //metHandle->front().hadEtInHE();
      mTreeMEThadEtInHO[metn]         = -1; //metHandle->front().hadEtInHO();
      mTreeMEThadEtInHB[metn]         = -1; //metHandle->front().hadEtInHB();
      
      //pf spesific
      mTreeMETChargedEMEtFraction[metn]  = -1; //metHandle->front().ChargedEMEtFraction();
      mTreeMETChargedHadEtFraction[metn] = -1; //metHandle->front().ChargedHadEtFraction();
      mTreeMETMuonEtFraction[metn]       = -1; //metHandle->front().MuonEtFraction();
      mTreeMETNeutralEMFraction[metn]  = -1; //metHandle->front().NeutralEMFraction();
      mTreeMETNeutralHadEtFraction[metn] = -1; //metHandle->front().NeutralHadEtFraction();
      mTreeMETType6EtFraction[metn]      = -1; //metHandle->front().Type6EtFraction();
      mTreeMETType7EtFraction[metn]      = -1; //metHandle->front().Type7EtFraction();
    }
  }
  else {
    for (int k=1; k<=2; k++) {
      mTreeMET[k]         = 0.;
      mTreeMEX[k]         = 0.;
      mTreeMEY[k]         = 0.;
      mTreeSumET[k]       = 0.;
      mTreeMETeta[k]      = 0.;
      mTreeMETphi[k]      = 0.;
      mTreeSumETSignif[k] = 0.;
      mTreeMETSignif[k]   = 0.;
      mTreeMETCaloMETInmHF[k]     = -1; //metHandle->front().CaloMETInmHF();
      mTreeMETCaloMETInpHF[k]     = -1; //metHandle->front().CaloMETInpHF();
      mTreeMETCaloMETPhiInmHF[k]  = -1; //metHandle->front().CaloMETPhiInmHF();
      mTreeMETCaloMETPhiInpHF[k]  = -1; //metHandle->front().CaloMETPhiInpHF();
      mTreeMETCaloSETInmHF[k]     = -1; //metHandle->front().CaloSETInmHF();
      mTreeMETCaloSETInpHF[k]     = -1; //metHandle->front().CaloSETInpHF();
      mTreeMETemEtFraction[k]      = -1; //metHandle->front().emEtFraction();
      mTreeMETetFractionHadronic[k]= -1; //metHandle->front().etFractionHadronic();
      mTreeMETmaxEtInEmTowers[k]   = -1; //metHandle->front().maxEtInEmTowers();
      mTreeMETmaxEtInHadTowers[k]  = -1; //metHandle->front().maxEtInHadTowers();
      mTreeMETemEtInHF[k]          = -1; //metHandle->front().emEtInHF();
      mTreeMETemEtInEE[k]          = -1; //metHandle->front().emEtInEE();
      mTreeMETemEtInEB[k]          = -1; //metHandle->front().emEtInEB();
      mTreeMEThadEtInHF[k]         = -1; //metHandle->front().hadEtInHF();
      mTreeMEThadEtInHE[k]         = -1; //metHandle->front().hadEtInHE();
      mTreeMEThadEtInHO[k]         = -1; //metHandle->front().hadEtInHO();
      mTreeMEThadEtInHB[k]         = -1; //metHandle->front().hadEtInHB();
      
      //pf spesific
      mTreeMETChargedEMEtFraction[k]  = -1; //metHandle->front().ChargedEMEtFraction();
      mTreeMETChargedHadEtFraction[k] = -1; //metHandle->front().ChargedHadEtFraction();
      mTreeMETMuonEtFraction[k]       = -1; //metHandle->front().MuonEtFraction();
      mTreeMETNeutralEMFraction[k]  = -1; //metHandle->front().NeutralEMFraction();
      mTreeMETNeutralHadEtFraction[k] = -1; //metHandle->front().NeutralHadEtFraction();
      mTreeMETType6EtFraction[k]      = -1; //metHandle->front().Type6EtFraction();
      mTreeMETType7EtFraction[k]      = -1; //metHandle->front().Type7EtFraction();
    }
  }


  // ParticleFlow MET
  iEvent.getByLabel(metTagPF_, metHandle);
  if ( !metHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagPF_;
  if ( metHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << metHandle->size() << " instead of 1";
  if ( metHandle.isValid() && metHandle->size()==1 ) {
    
    int metn=3;
    mTreeMET[metn]         = metHandle->front().et();
    mTreeMEX[metn]         = metHandle->front().momentum().X();
    mTreeMEY[metn]         = metHandle->front().momentum().Y();
    mTreeSumET[metn]       = metHandle->front().sumEt();
    mTreeMETeta[metn]      = metHandle->front().eta();
    mTreeMETphi[metn]      = metHandle->front().phi();
    mTreeSumETSignif[metn] = metHandle->front().mEtSig();
    //if(metHandle->front().getSignificanceMatrix());
    //mTreeMETSignif[metn]   = metHandle->front().significance();
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = -1; //metHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = -1; //metHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = -1; //metHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = -1; //metHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = -1; //metHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = -1; //metHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = -1; //metHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= -1; //metHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = -1; //metHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = -1; //metHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = -1; //metHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = -1; //metHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = -1; //metHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = -1; //metHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = -1; //metHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = -1; //metHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = -1; //metHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = metHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = metHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = metHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = metHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = metHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = metHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = metHandle->front().Type7EtFraction();
  }

  // TrackCorrected MET
  iEvent.getByLabel(metTagTC_, metHandle);
  if ( !metHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagTC_;
  if ( metHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << metHandle->size() << " instead of 1";
  if ( metHandle.isValid() && metHandle->size()==1 ) {
    
    int metn=4;
    mTreeMET[metn]         = metHandle->front().et();
    mTreeMEX[metn]         = metHandle->front().momentum().X();
    mTreeMEY[metn]         = metHandle->front().momentum().Y();
    mTreeSumET[metn]       = metHandle->front().sumEt();
    mTreeMETeta[metn]      = metHandle->front().eta();
    mTreeMETphi[metn]      = metHandle->front().phi();
    mTreeSumETSignif[metn] = metHandle->front().mEtSig();
    mTreeMETSignif[metn]   = -1;
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = -1; //metHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = -1; //metHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = -1; //metHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = -1; //metHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = -1; //metHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = -1; //metHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = -1; //metHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= -1; //metHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = -1; //metHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = -1; //metHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = -1; //metHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = -1; //metHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = -1; //metHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = -1; //metHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = -1; //metHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = -1; //metHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = -1; //metHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = -1; //metHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = -1; //metHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = -1; //metHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = -1; //metHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = -1; //metHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = -1; //metHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = -1; //metHandle->front().Type7EtFraction();
  }
  
  
  //PFMETnoPU
  edm::Handle< std::vector<reco::PFMET> > PFmetHandle;
  iEvent.getByLabel(metTagPFnoPU_, PFmetHandle);
  if ( !PFmetHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagPF_;
  if ( PFmetHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << PFmetHandle->size() << " instead of 1";
  if ( PFmetHandle.isValid() && PFmetHandle->size()==1 ) {
    
    int metn=5;
    mTreeMET[metn]         = PFmetHandle->front().et();
    mTreeMEX[metn]         = PFmetHandle->front().momentum().X();
    mTreeMEY[metn]         = PFmetHandle->front().momentum().Y();
    mTreeSumET[metn]       = PFmetHandle->front().sumEt();
    mTreeMETeta[metn]      = PFmetHandle->front().eta();
    mTreeMETphi[metn]      = PFmetHandle->front().phi();
    mTreeSumETSignif[metn] = PFmetHandle->front().mEtSig();
    //mTreeMETSignif[metn]   = PFmetHandle->front().significance();
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = -1; //PFmetHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = -1; //PFmetHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = -1; //PFmetHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = -1; //PFmetHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = -1; //PFmetHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = -1; //PFmetHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = -1; //PFmetHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= -1; //PFmetHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = -1; //PFmetHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = -1; //PFmetHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = -1; //PFmetHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = -1; //PFmetHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = -1; //PFmetHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = -1; //PFmetHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = -1; //PFmetHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = -1; //PFmetHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = -1; //PFmetHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = PFmetHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = PFmetHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = PFmetHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = PFmetHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = PFmetHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = PFmetHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = PFmetHandle->front().Type7EtFraction();
  }
  
  
  // metTagJESCorAK5CaloJetMuons
  edm::Handle< std::vector<reco::CaloMET> > CalometHandle;
  iEvent.getByLabel(metTagJESCorAK5CaloJetMuons_, CalometHandle);
  if ( !CalometHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagJESCorAK5CaloJetMuons_;
  if ( CalometHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << CalometHandle->size() << " instead of 1";
  if ( CalometHandle.isValid() && CalometHandle->size()==1 ) {
    
    int metn=6;
    mTreeMET[metn]         = CalometHandle->front().et();
    mTreeMEX[metn]         = CalometHandle->front().momentum().X();
    mTreeMEY[metn]         = CalometHandle->front().momentum().Y();
    mTreeSumET[metn]       = CalometHandle->front().sumEt();
    mTreeMETeta[metn]      = CalometHandle->front().eta();
    mTreeMETphi[metn]      = CalometHandle->front().phi();
    mTreeSumETSignif[metn] = CalometHandle->front().mEtSig();
    mTreeMETSignif[metn]   = -1;
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = CalometHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = CalometHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = CalometHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = CalometHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = CalometHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = CalometHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = CalometHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= CalometHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = CalometHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = CalometHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = CalometHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = CalometHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = CalometHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = CalometHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = CalometHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = CalometHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = CalometHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = -1; //CalometHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = -1; //CalometHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = -1; //CalometHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = -1; //CalometHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = -1; //CalometHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = -1; //CalometHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = -1; //CalometHandle->front().Type7EtFraction();
  }
  
  
  // metTagcorMetGlobalMuons
  iEvent.getByLabel(metTagcorMetGlobalMuons_, CalometHandle);
  if ( !CalometHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagcorMetGlobalMuons_;
  if ( CalometHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << CalometHandle->size() << " instead of 1";
  if ( CalometHandle.isValid() && CalometHandle->size()==1 ) {
    
    int metn=7;
    mTreeMET[metn]         = CalometHandle->front().et();
    mTreeMEX[metn]         = CalometHandle->front().momentum().X();
    mTreeMEY[metn]         = CalometHandle->front().momentum().Y();
    mTreeSumET[metn]       = CalometHandle->front().sumEt();
    mTreeMETeta[metn]      = CalometHandle->front().eta();
    mTreeMETphi[metn]      = CalometHandle->front().phi();
    mTreeSumETSignif[metn] = CalometHandle->front().mEtSig();
    mTreeMETSignif[metn]   = -1;
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = CalometHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = CalometHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = CalometHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = CalometHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = CalometHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = CalometHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = CalometHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= CalometHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = CalometHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = CalometHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = CalometHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = CalometHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = CalometHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = CalometHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = CalometHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = CalometHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = CalometHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = -1; //CalometHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = -1; //CalometHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = -1; //CalometHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = -1; //CalometHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = -1; //CalometHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = -1; //CalometHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = -1; //CalometHandle->front().Type7EtFraction();
  }
  
  
  // metTagHO
  iEvent.getByLabel(metTagHO_, CalometHandle);
  if ( !CalometHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagHO_;
  if ( CalometHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << CalometHandle->size() << " instead of 1";
  if ( CalometHandle.isValid() && CalometHandle->size()==1 ) {
    
    int metn=8;
    mTreeMET[metn]         = CalometHandle->front().et();
    mTreeMEX[metn]         = CalometHandle->front().momentum().X();
    mTreeMEY[metn]         = CalometHandle->front().momentum().Y();
    mTreeSumET[metn]       = CalometHandle->front().sumEt();
    mTreeMETeta[metn]      = CalometHandle->front().eta();
    mTreeMETphi[metn]      = CalometHandle->front().phi();
    mTreeSumETSignif[metn] = CalometHandle->front().mEtSig();
    mTreeMETSignif[metn]   = -1;
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = CalometHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = CalometHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = CalometHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = CalometHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = CalometHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = CalometHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = CalometHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= CalometHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = CalometHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = CalometHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = CalometHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = CalometHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = CalometHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = CalometHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = CalometHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = CalometHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = CalometHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = -1; //CalometHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = -1; //CalometHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = -1; //CalometHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]    = -1; //CalometHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = -1; //CalometHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = -1; //CalometHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = -1; //CalometHandle->front().Type7EtFraction();
  }
  
  // metTagNoHF
  iEvent.getByLabel(metTagNoHF_, CalometHandle);
  if ( !CalometHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagNoHF_;
  if ( CalometHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << CalometHandle->size() << " instead of 1";
  if ( CalometHandle.isValid() && CalometHandle->size()==1 ) {
    
    int metn=9;
    mTreeMET[metn]         = CalometHandle->front().et();
    mTreeMEX[metn]         = CalometHandle->front().momentum().X();
    mTreeMEY[metn]         = CalometHandle->front().momentum().Y();
    mTreeSumET[metn]       = CalometHandle->front().sumEt();
    mTreeMETeta[metn]      = CalometHandle->front().eta();
    mTreeMETphi[metn]      = CalometHandle->front().phi();
    mTreeSumETSignif[metn] = CalometHandle->front().mEtSig();
    mTreeMETSignif[metn]   = -1;
    
    //calo spesific 
    mTreeMETCaloMETInmHF[metn]     = CalometHandle->front().CaloMETInmHF();
    mTreeMETCaloMETInpHF[metn]     = CalometHandle->front().CaloMETInpHF();
    mTreeMETCaloMETPhiInmHF[metn]  = CalometHandle->front().CaloMETPhiInmHF();
    mTreeMETCaloMETPhiInpHF[metn]  = CalometHandle->front().CaloMETPhiInpHF();
    mTreeMETCaloSETInmHF[metn]     = CalometHandle->front().CaloSETInmHF();
    mTreeMETCaloSETInpHF[metn]     = CalometHandle->front().CaloSETInpHF();
    mTreeMETemEtFraction[metn]      = CalometHandle->front().emEtFraction();
    mTreeMETetFractionHadronic[metn]= CalometHandle->front().etFractionHadronic();
    mTreeMETmaxEtInEmTowers[metn]   = CalometHandle->front().maxEtInEmTowers();
    mTreeMETmaxEtInHadTowers[metn]  = CalometHandle->front().maxEtInHadTowers();
    mTreeMETemEtInHF[metn]          = CalometHandle->front().emEtInHF();
    mTreeMETemEtInEE[metn]          = CalometHandle->front().emEtInEE();
    mTreeMETemEtInEB[metn]          = CalometHandle->front().emEtInEB();
    mTreeMEThadEtInHF[metn]         = CalometHandle->front().hadEtInHF();
    mTreeMEThadEtInHE[metn]         = CalometHandle->front().hadEtInHE();
    mTreeMEThadEtInHO[metn]         = CalometHandle->front().hadEtInHO();
    mTreeMEThadEtInHB[metn]         = CalometHandle->front().hadEtInHB();
    
    //pf spesific
    mTreeMETChargedEMEtFraction[metn]  = -1; //CalometHandle->front().ChargedEMEtFraction();
    mTreeMETChargedHadEtFraction[metn] = -1; //CalometHandle->front().ChargedHadEtFraction();
    mTreeMETMuonEtFraction[metn]       = -1; //CalometHandle->front().MuonEtFraction();
    mTreeMETNeutralEMFraction[metn]  = -1; //CalometHandle->front().NeutralEMFraction();
    mTreeMETNeutralHadEtFraction[metn] = -1; //CalometHandle->front().NeutralHadEtFraction();
    mTreeMETType6EtFraction[metn]      = -1; //CalometHandle->front().Type6EtFraction();
    mTreeMETType7EtFraction[metn]      = -1; //CalometHandle->front().Type7EtFraction();
  }
  
  // This filter
  if (nele_     > 0 && cele_     < nele_)     return 0;
  if (npfele_   > 0 && cpfele_   < npfele_)   return 0;
  if (nmuo_     > 0 && cmuo_     < nmuo_)     return 0;
  if (ncalojet_ > 0 && ccalojet_ < ncalojet_) return 0;
  if (npfjet_   > 0 && cpfjet_   < npfjet_)   return 0;

  if (metcalo_ > 0 && mTreeMET[0] < metcalo_) return 0;
  if (metpf_ > 0   && mTreeMET[3] < metpf_)   return 0;
  if (mettc_ > 0   && mTreeMET[4] < mettc_)   return 0;
  //cut on HT
  if(ncalojet_ > 0 && PFhtc_  > 0   && PFhtcutsum_<PFhtc_ )  return 0;
  if(npfjet_   > 0 && htc_  > 0     && htcutsum_  <htc_   )  return 0;

  nrEventPassedRaw_++;

  // ECal Noise
  edm::Handle<EcalRecHitCollection> pEBRecHits;
  //  iEvent.getByLabel( "ecalRecHit", "EcalRecHitsEB", pEBRecHits );
  //  iEvent.getByLabel( "reducedEcalRecHitsEB", pEBRecHits );
  iEvent.getByLabel(ebhitsTag_, pEBRecHits);


  if ( !pEBRecHits.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No EcalRecHitCollection found for InputTag " << ebhitsTag_;
  }
  else {

    const EcalRecHitCollection *ebRecHits = pEBRecHits.product();
    
    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry* geo = pG.product();

    edm::ESHandle<CaloTopology> pTopology;
    iSetup.get<CaloTopologyRecord>().get(pTopology);

    if ( pTopology.isValid() && pG.isValid() ) {

      const CaloTopology *topology = pTopology.product();
      
      Double_t energy = 0.;
      EcalRecHit maxhit;

      for (EcalRecHitCollection::const_iterator elem = ebRecHits->begin(); 
	   elem != ebRecHits->end(); ++elem) {
	
	if (elem->energy() > energy) {
          energy = elem->energy();
	  maxhit = *elem;
	}        
      }   

      EBDetId EBId = maxhit.id();  
      GlobalPoint pos = geo->getPosition(maxhit.detid());

      const reco::BasicCluster *cluster = 0;  // dummy
      Double_t e3x3 = EcalClusterTools::matrixEnergy( *cluster, ebRecHits, topology, maxhit.id(), -1, 1, -1, 1 );

      if (e3x3>1e-6) {
	
	double ctgTheta = (pos.z() - vtxPoint.Z())/pos.perp();
	double newEta = asinh(ctgTheta);  
	double pf = 1.0/cosh(newEta);
	
	mTreeecalr9   = energy/e3x3;
	mTreeecale    = energy;
	mTreeecalpt   = energy*pf ;
	mTreeecaleta  = newEta;
	mTreeecalphi  = pos.phi();
	mTreeecalpx   = mTreeecalpt*cos(mTreeecalphi);
	mTreeecalpy   = mTreeecalpt*sin(mTreeecalphi);
	mTreeecalpz   = mTreeecale*TMath::TanH(mTreeecaleta);
	mTreeecaltime = maxhit.time();
	mTreeecalchi  = maxhit.chi2();
	mTreeecalflag = maxhit.recoFlag();
	mTreeecalieta = EBId.ieta();
	mTreeecaliphi = EBId.iphi();
      }
      else {

	mTreeecalr9   = -1.;
	mTreeecale    = -1.;
	mTreeecalpt   = -1.;
	mTreeecaleta  = -999.;
	mTreeecalphi  = -999.;
	mTreeecalpx   = -999.;
	mTreeecalpy   = -999.;
	mTreeecalpz   = -999.;
	mTreeecaltime = -999.;
	mTreeecalchi  = -999.;
	mTreeecalflag = -999;
	mTreeecalieta = -999;
	mTreeecaliphi = -999;

      }
      //      cout << iEvent.time().value() << " -> " << energy << "  " << energy/e3x3 << "  " << maxhit.time() << "  " 
      //	   << maxhit.outOfTimeEnergy() << "  " << maxhit.chi2Prob() << endl;
      //      cout << " -> " << energy << "  " << e3x3 << "  " << pos.eta() << "  " << pos.phi() 
      //   << "  " << energy*pf << "  " << newEta  << "  " << mTreeecalpx << "  " 
      //   << mTreeecalpy << "  " << mTreeecalpz << endl;
    }
  }
  
  if(susyPar_){
    // SUSY variables

    Handle<double> susyScanM0;
    iEvent.getByLabel("susyScanM0", susyScanM0);
    mTreesusyScanM0 = *susyScanM0;

    Handle<double> susyScanM12;
    iEvent.getByLabel("susyScanM12", susyScanM12);
    mTreesusyScanM12 = *susyScanM12;

    Handle<double> susyScanA0;
    iEvent.getByLabel("susyScanA0", susyScanA0);
    mTreesusyScanA0 = *susyScanA0;

    Handle<double> susyScanCrossSection;
    iEvent.getByLabel("susyScanCrossSection", susyScanCrossSection);
    mTreesusyScanCrossSection = *susyScanCrossSection;

    Handle<double> susyScanMu;
    iEvent.getByLabel("susyScanCrossSection", susyScanMu);
    mTreesusyScanMu = *susyScanMu;

    Handle<double> susyScanRun;
    iEvent.getByLabel("susyScanRun", susyScanRun);
    mTreesusyScanRun = *susyScanRun;

    Handle<double> susyScantanbeta;
    iEvent.getByLabel("susyScantanbeta", susyScantanbeta);
    mTreesusyScantanbeta = *susyScantanbeta;
  }else{
    mTreesusyScanM0 =-1;
    mTreesusyScanM12 =-1;
    mTreesusyScanA0 = -1;
    mTreesusyScanCrossSection = -1;
    mTreesusyScanMu = -1;
    mTreesusyScanRun = -1;
    mTreesusyScantanbeta = 99999;
  }
  // Fill the tree
  mAllData->Fill();
  return 1;

}
////////////////////////////////
//
// Helper
//
bool SusyACSkimAnalysis::isSUSY(int pdgid) {

  int tid = abs(pdgid);

  if ( (tid>1000000 && tid<1000017) ||
       (tid>1000020 && tid<1000040) ||
       (tid>2000000 && tid<2000016) )
    return true;
  else
    return false;
}
bool SusyACSkimAnalysis::isStable(int pdgid) {

  int tid = abs(pdgid);

  if ( ( tid>0  && tid<6 ) ||     // quarks
       ( tid>10 && tid<15 ) ||    // leptons
       tid == 16 ||               // lepons
       tid == 21 || tid == 22 ||  // bosons
       tid == 1000022)            // neutralino
    return true;
  else
    return false;
}

bool SusyACSkimAnalysis::isDecaying(int pdgid) {

  int tid = abs(pdgid);

  if ( (tid>1000000 && tid<1000017) ||
       (tid>1000020&&  tid<1000040) || //&&  tid!=1000022) || 
       (tid>2000000 && tid<2000016) ||
       tid==25 || tid==24 || tid==23 || tid==6 || tid==15 )
    return true;
  else
    return false;
}
  
////////////////////////////////
//
// Begin of Job
//
void SusyACSkimAnalysis::beginJob() {}

////////////////////////////////
//
// End of Job
//
void SusyACSkimAnalysis::endJob() {

  h_counters->SetBinContent(1, nrEventTotalRaw_);
  h_counters->SetBinContent(2, nrEventPassedPthatRaw_);
  h_counters->SetBinContent(3, nrEventPassedRaw_);

  printSummary();

}

////////////////////////////////
//
// Summary
//
void SusyACSkimAnalysis::printSummary( void ) {
  
  edm::LogInfo("SusyACSkimAnalysis") << "*** Summary of counters: ";
  edm::LogVerbatim("SusyACSkimAnalysis") << "Total number of events = " << nrEventTotalRaw_
					 << " ; selected = " << nrEventPassedRaw_ << endl; 

}

////////////////////////////////
//
// Initialize Plots for ntuple
//
void SusyACSkimAnalysis::initPlots() {

  // Register this ntuple
  edm::Service<TFileService> fs;

  // Now we add some additional ones for the dijet analysis
  mAllData = fs->make<TTree>( "allData", "data after cuts" );

  h_counters = fs->make<TH1F>("h_counters", "Event Counter", 10, 0, 10);

  mAllData->SetAutoSave(10000);

  // Add the branches

  // Event
  mAllData->Branch("global_weight",  &mTreeEventWeight, "global_weight/double");
  mAllData->Branch("global_procID",  &mTreeProcID,      "global_procID/I");
  mAllData->Branch("global_pthat",   &mTreePthat,       "global_pthat/double");
  mAllData->Branch("global_bfield",  &mTreebfield,      "global_bfield/double");
  mAllData->Branch("global_store",   &mTreestore,       "global_store/I");
  mAllData->Branch("global_run",     &mTreerun,         "global_run/I");
  mAllData->Branch("global_event",   &mTreeevent,       "global_event/I");
  mAllData->Branch("global_bx",      &mTreebx,          "global_bx/I");
  mAllData->Branch("global_orbit",   &mTreeorbit,       "global_orbit/I");
  mAllData->Branch("global_exp",     &mTreeexp,         "global_exp/I");
  mAllData->Branch("global_isdata",  &mTreedata,        "global_isdata/I");

  mAllData->Branch("lumi_section", &mTreelumiblk,    "lumi_section/I");
  mAllData->Branch("lumi_del",     &mTreedellumi,    "lumi_del/double");
  mAllData->Branch("lumi_rec",     &mTreereclumi,    "lumi_rec/double");
  mAllData->Branch("lumi_delerr",  &mTreedellumierr, "lumi_delerr/double");
  mAllData->Branch("lumi_recerr",  &mTreereclumierr, "lumi_recerr/double");
  
  
  
  mAllData->Branch("pu_bunchx",  &mTreePUbx, "pu_bunchx/I");
  mAllData->Branch("pu_n",  &nTreePileUp, "pu_n/I");
  mAllData->Branch("pu_vtxn",  &mTreePUN, "pu_vtxn/I");
  mAllData->Branch("pu_num_int", mTreePUNumInteractions, "pu_num_int[pu_n]/I");
  mAllData->Branch("pu_inst_Lumi", mTreePUInstLumi, "pu_inst_Lumi[pu_n][100]/double");
  mAllData->Branch("pu_zPos", mTreePUzPosi , "pu_zPos[pu_n][100]/double");
  mAllData->Branch("pu_sumPthi",mTreePUsumPthi, "pu_sumPthi[pu_n][100]/double");
  mAllData->Branch("pu_sumPtlo",mTreePUsumPtlo, "pu_sumPtlo[pu_n][100]/double");

  mAllData->Branch("noise_pLoose", &mTreenoisel, "noise_pLoose/I");
  mAllData->Branch("noise_pTight", &mTreenoiset, "noise_pTight/I");
  mAllData->Branch("noise_pHigh",  &mTreenoiseh, "noise_pHigh/I");

  mAllData->Branch("noise_ecal_r9",   &mTreeecalr9,   "noise_ecal_r9/double");
  mAllData->Branch("noise_ecal_E",    &mTreeecale,    "noise_ecal_E/double");
  mAllData->Branch("noise_ecal_pt",   &mTreeecalpt,   "noise_ecal_pt/double");
  mAllData->Branch("noise_ecal_px",   &mTreeecalpx,   "noise_ecal_px/double");
  mAllData->Branch("noise_ecal_py",   &mTreeecalpy,   "noise_ecal_py/double");
  mAllData->Branch("noise_ecal_pz" ,  &mTreeecalpz,   "noise_ecal_pz/double");
  mAllData->Branch("noise_ecal_eta",  &mTreeecaleta,  "noise_ecal_eta/double");
  mAllData->Branch("noise_ecal_phi",  &mTreeecalphi,  "noise_ecal_phi/double");
  mAllData->Branch("noise_ecal_time", &mTreeecaltime, "noise_ecal_time/double");
  mAllData->Branch("noise_ecal_chi",  &mTreeecalchi,  "noise_ecal_chi/double");
  mAllData->Branch("noise_ecal_flag", &mTreeecalflag, "noise_ecal_flag/I");
  mAllData->Branch("noise_ecal_ieta", &mTreeecalieta, "noise_ecal_ieta/I");
  mAllData->Branch("noise_ecal_iphi", &mTreeecaliphi, "noise_ecal_iphi/I");
  
  mAllData->Branch("noise_hcal_eventChargeFraction",  &mTreenoiseHCALeventChargeFraction, "noise_hcal_eventChargeFraction/double");
  mAllData->Branch("noise_hcal_eventEMEnergy",        &mTreenoiseHCALeventEMEnergy, "noise_hcal_eventEMEnergy/double");
  mAllData->Branch("noise_hcal_eventEMFraction",      &mTreenoiseHCALeventEMFraction, "noise_hcal_eventEMFraction/double");
  mAllData->Branch("noise_hcal_eventHadEnergy",       &mTreenoiseHCALeventHadEnergy, "noise_hcal_eventHadEnergy/double");
  mAllData->Branch("noise_hcal_eventTrackEnergy",     &mTreenoiseHCALeventTrackEnergy, "noise_hcal_eventTrackEnergy/double");
  mAllData->Branch("noise_hcal_flatNoiseSumE",        &mTreenoiseHCALflatNoiseSumE, "noise_hcal_flatNoiseSumE/double");
  mAllData->Branch("noise_hcal_flatNoiseSumEt",       &mTreenoiseHCALflatNoiseSumEt, "noise_hcal_flatNoiseSumEt/double");
  mAllData->Branch("noise_hcal_HasBadRBXTS4TS5",      &mTreenoiseHCALHasBadRBXTS4TS5, "noise_hcal_HasBadRBXTS4TS5/I");
  mAllData->Branch("noise_hcal_isolatedNoiseSumE",    &mTreenoiseHCALisolatedNoiseSumE, "noise_hcal_isolatedNoiseSumE/double");
  mAllData->Branch("noise_hcal_isolatedNoiseSumEt",   &mTreenoiseHCALisolatedNoiseSumEt, "noise_hcal_isolatedNoiseSumEt/double");
  mAllData->Branch("noise_hcal_max10GeVHitTime",      &mTreenoiseHCALmax10GeVHitTime, "noise_hcal_max10GeVHitTime/double");
  mAllData->Branch("noise_hcal_max25GeVHitTime",      &mTreenoiseHCALmax25GeVHitTime, "noise_hcal_max25GeVHitTime/double");
  mAllData->Branch("noise_hcal_maxE10TS",             &mTreenoiseHCALmaxE10TS, "noise_hcal_maxE10TS/double");
  mAllData->Branch("noise_hcal_maxE2Over10TS",        &mTreenoiseHCALmaxE2Over10TS, "noise_hcal_maxE2Over10TS/double");
  mAllData->Branch("noise_hcal_maxE2TS",              &mTreenoiseHCALmaxE2TS, "noise_hcal_maxE2TS/double");
  mAllData->Branch("noise_hcal_maxHPDHits",           &mTreenoiseHCALmaxHPDHits, "noise_hcal_maxHPDHits/I");
  mAllData->Branch("noise_hcal_maxHPDNoOtherHits",    &mTreenoiseHCALmaxHPDNoOtherHits, "noise_hcal_maxHPDNoOtherHits/I");
  mAllData->Branch("noise_hcal_maxRBXHits",           &mTreenoiseHCALmaxRBXHits, "noise_hcal_maxRBXHits/I");
  mAllData->Branch("noise_hcal_maxZeros",             &mTreenoiseHCALmaxZeros, "noise_hcal_maxZeros/I");
  mAllData->Branch("noise_hcal_min10GeVHitTime",      &mTreenoiseHCALmin10GeVHitTime, "noise_hcal_min10GeVHitTime/double");
  mAllData->Branch("noise_hcal_min25GeVHitTime",      &mTreenoiseHCALmin25GeVHitTime, "noise_hcal_min25GeVHitTime/double");
  mAllData->Branch("noise_hcal_minE10TS",             &mTreenoiseHCALminE10TS, "noise_hcal_minE10TS/double");
  mAllData->Branch("noise_hcal_minE2Over10TS",        &mTreenoiseHCALminE2Over10TS, "noise_hcal_minE2Over10TS/double");
  mAllData->Branch("noise_hcal_minE2TS",              &mTreenoiseHCALminE2TS, "noise_hcal_minE2TS/double");
  mAllData->Branch("noise_hcal_minHPDEMF",            &mTreenoiseHCALminHPDEMF, "noise_hcal_minHPDEMF/double");
  mAllData->Branch("noise_hcal_minRBXEMF",            &mTreenoiseHCALminRBXEMF, "noise_hcal_minRBXEMF/double");
  mAllData->Branch("noise_hcal_noiseFilterStatus",    &mTreenoiseHCALnoiseFilterStatus, "noise_hcal_noiseFilterStatus/I");
  mAllData->Branch("noise_hcal_noiseType",            &mTreenoiseHCALnoiseType, "noise_hcal_noiseType/I");
  mAllData->Branch("noise_hcal_num10GeVHits",         &mTreenoiseHCALnum10GeVHits, "noise_hcal_num10GeVHits/I");
  mAllData->Branch("noise_hcal_num25GeVHits",         &mTreenoiseHCALnum25GeVHits, "noise_hcal_num25GeVHits/I");
  mAllData->Branch("noise_hcal_numFlatNoiseChannels", &mTreenoiseHCALnumFlatNoiseChannels, "noise_hcal_numFlatNoiseChannels/I");
  mAllData->Branch("noise_hcal_numIsolatedNoiseChannels",&mTreenoiseHCALnumIsolatedNoiseChannels, "noise_hcal_numIsolatedNoiseChannels/I");
  mAllData->Branch("noise_hcal_numProblematicRBXs",   &mTreenoiseHCALnumProblematicRBXs, "noise_hcal_numProblematicRBXs/I");
  mAllData->Branch("noise_hcal_numSpikeNoiseChannels",&mTreenoiseHCALnumSpikeNoiseChannels, "noise_hcal_numSpikeNoiseChannels/I");
  mAllData->Branch("noise_hcal_numTriangleNoiseChannels",&mTreenoiseHCALnumTriangleNoiseChannels, "noise_hcal_numTriangleNoiseChannels/I");
  mAllData->Branch("noise_hcal_numTS4TS5NoiseChannels",&mTreenoiseHCALnumTS4TS5NoiseChannels, "noise_hcal_numTS4TS5NoiseChannels/I");
  mAllData->Branch("noise_hcal_passHighLevelNoiseFilter",&mTreenoiseHCALpassHighLevelNoiseFilter, "noise_hcal_passHighLevelNoiseFilter/I");
  mAllData->Branch("noise_hcal_passLooseNoiseFilter", &mTreenoiseHCALpassLooseNoiseFilter, "noise_hcal_passLooseNoiseFilter/I");
  mAllData->Branch("noise_hcal_passTightNoiseFilter", &mTreenoiseHCALpassTightNoiseFilter, "noise_hcal_passTightNoiseFilter/I");
  mAllData->Branch("noise_hcal_rms10GeVHitTime",      &mTreenoiseHCALrms10GeVHitTime, "noise_hcal_rms10GeVHitTime/double");
  mAllData->Branch("noise_hcal_rms25GeVHitTime",      &mTreenoiseHCALrms25GeVHitTime, "noise_hcal_rms25GeVHitTime/double");
  mAllData->Branch("noise_hcal_spikeNoiseSumE",       &mTreenoiseHCALspikeNoiseSumE, "noise_hcal_spikeNoiseSumE/double");
  mAllData->Branch("noise_hcal_spikeNoiseSumEt",      &mTreenoiseHCALspikeNoiseSumEt, "noise_hcal_spikeNoiseSumEt/double");
  mAllData->Branch("noise_hcal_triangleNoiseSumE",    &mTreenoiseHCALtriangleNoiseSumE, "noise_hcal_triangleNoiseSumE/double");
  mAllData->Branch("noise_hcal_triangleNoiseSumEt",   &mTreenoiseHCALtriangleNoiseSumEt, "noise_hcal_triangleNoiseSumEt/double");
  mAllData->Branch("noise_hcal_TS4TS5NoiseSumE",      &mTreenoiseHCALTS4TS5NoiseSumE, "noise_hcal_TS4TS5NoiseSumE/double");
  mAllData->Branch("noise_hcal_TS4TS5NoiseSumEt",     &mTreenoiseHCALTS4TS5NoiseSumEt, "noise_hcal_TS4TS5NoiseSumEt/double");
  
  
  
  //mAllData->Branch("noise_hcal_flag",     &mTreenoiseHCALhcalnoiseFlag, "noise_hcal_flag/bool");
  

  mAllData->Branch("trig_HLTName",    &mTreetrighltname, "trig_HLTName[20]/I");
  mAllData->Branch("trig_n",          &mTreeNtrig,       "trig_n/I");
  mAllData->Branch("trig_L1prescale",  mTreetrigL1pre,   "trig_L1prescale[trig_n]/I");
  mAllData->Branch("trig_HLTprescale", mTreetrigHLTpre,  "trig_HLTprescale[trig_n]/I");
  mAllData->Branch("trig_name",        mTreetrigname,    "trig_name[trig_n][20]/I");
  mAllData->Branch("trig_filter",      mTreefiltname,    "trig_filter[trig_n][20]/I");
  mAllData->Branch("trig_pt",          mTreetrigpt,      "trig_pt[trig_n]/double");
  mAllData->Branch("trig_eta",         mTreetrigeta,     "trig_eta[trig_n]/double");
  mAllData->Branch("trig_phi",         mTreetrigphi,     "trig_phi[trig_n]/double");

  // Truth
  mAllData->Branch("truth_n",      &mTreeNtruth,      "truth_n/I");
  mAllData->Branch("truth_pdgid",   mTreetruthpdgid,  "truth_pdgid[truth_n]/I");
  mAllData->Branch("truth_bvtxid",  mTreetruthbvtxid, "truth_bvtxid[truth_n]/I");
  mAllData->Branch("truth_evtxid",  mTreetruthevtxid, "truth_evtxid[truth_n]/I");
  mAllData->Branch("truth_E",       mTreetruthE,      "truth_E[truth_n]/double");
  mAllData->Branch("truth_Et",      mTreetruthEt,     "truth_Et[truth_n]/double");
  mAllData->Branch("truth_p",       mTreetruthp,      "truth_p[truth_n]/double");
  mAllData->Branch("truth_pt",      mTreetruthpt,     "truth_pt[truth_n]/double");
  mAllData->Branch("truth_px",      mTreetruthpx,     "truth_px[truth_n]/double");
  mAllData->Branch("truth_py",      mTreetruthpy,     "truth_py[truth_n]/double");
  mAllData->Branch("truth_pz",      mTreetruthpz,     "truth_pz[truth_n]/double");
  mAllData->Branch("truth_eta",     mTreetrutheta,    "truth_eta[truth_n]/double");
  mAllData->Branch("truth_phi",     mTreetruthphi,    "truth_phi[truth_n]/double");
  mAllData->Branch("truth_m",       mTreetruthm,      "truth_m[truth_n]/double");

  // Truth stable e/mu in final state
  mAllData->Branch("truthl_n",      &mTreeNtruthl,      "truthl_n/I");
  mAllData->Branch("truthl_ori",     mTreetruthlori,    "truthl_ori[truthl_n]/I");
  mAllData->Branch("truthl_pdgid",   mTreetruthlpdgid,  "truthl_pdgid[truthl_n]/I");
  mAllData->Branch("truthl_E",       mTreetruthlE,      "truthl_E[truthl_n]/double");
  mAllData->Branch("truthl_Et",      mTreetruthlEt,     "truthl_Et[truthl_n]/double");
  mAllData->Branch("truthl_p",       mTreetruthlp,      "truthl_p[truthl_n]/double");
  mAllData->Branch("truthl_pt",      mTreetruthlpt,     "truthl_pt[truthl_n]/double");
  mAllData->Branch("truthl_px",      mTreetruthlpx,     "truthl_px[truthl_n]/double");
  mAllData->Branch("truthl_py",      mTreetruthlpy,     "truthl_py[truthl_n]/double");
  mAllData->Branch("truthl_pz",      mTreetruthlpz,     "truthl_pz[truthl_n]/double");
  mAllData->Branch("truthl_eta",     mTreetruthleta,    "truthl_eta[truthl_n]/double");
  mAllData->Branch("truthl_phi",     mTreetruthlphi,    "truthl_phi[truthl_n]/double");

  // PDF Info
  mAllData->Branch("pdf_id1",   &mTreepdfid1,   "pdf_id1/I");
  mAllData->Branch("pdf_id2",   &mTreepdfid2,   "pdf_id2/I");
  mAllData->Branch("pdf_x1",    &mTreepdfx1,    "pdf_x1/F");
  mAllData->Branch("pdf_x2",    &mTreepdfx2,    "pdf_x2/F");
  mAllData->Branch("pdf_f1",    &mTreepdff1,    "pdf_f1/F");
  mAllData->Branch("pdf_f2",    &mTreepdff2,    "pdf_f2/F");
  mAllData->Branch("pdf_scale", &mTreepdfscale, "pdf_scale/F");

  // Vertex
  mAllData->Branch("vtx_n",   &mTreeNvtx,   "vtx_n/I");
  mAllData->Branch("vtx_ntr",  mTreeVtxntr, "vtx_ntr[vtx_n]/I");
  mAllData->Branch("vtx_fake", mTreeVtxfake,"vtx_fake[vtx_n]/I");
  mAllData->Branch("vtx_ndof", mTreeVtxndf, "vtx_ndof[vtx_n]/double");
  mAllData->Branch("vtx_x",    mTreeVtxx,   "vtx_x[vtx_n]/double");
  mAllData->Branch("vtx_y",    mTreeVtxy,   "vtx_y[vtx_n]/double");
  mAllData->Branch("vtx_z",    mTreeVtxz,   "vtx_z[vtx_n]/double");
  mAllData->Branch("vtx_chi",  mTreeVtxchi, "vtx_chi[vtx_n]/double");

  // Beam Spot
  mAllData->Branch("bs_x", &mTreebspX, "bs_x/double");
  mAllData->Branch("bs_y", &mTreebspY, "bs_y/double");
  mAllData->Branch("bs_z", &mTreebspZ, "bs_z/double");

  // Tracks
  mAllData->Branch("tracks_n",   &mTreeNtracks,   "tracks_n/I");
  mAllData->Branch("tracks_hqf", &mTreetrackshqf, "tracks_hqf/double");

  // MET
  mAllData->Branch("met_et",       &mTreeMET,         "met_et[10]/double");
  mAllData->Branch("met_ex",       &mTreeMEX,         "met_ex[10]/double");
  mAllData->Branch("met_ey",       &mTreeMEY,         "met_ey[10]/double");
  mAllData->Branch("met_eta",      &mTreeMETeta,      "met_eta[10]/double");
  mAllData->Branch("met_phi",      &mTreeMETphi,      "met_phi[10]/double");
  mAllData->Branch("met_sumet",    &mTreeSumET,       "met_sumet[10]/double");
  mAllData->Branch("met_sumetsig", &mTreeSumETSignif, "met_sumetsig[10]/double");
  mAllData->Branch("met_etsignif", &mTreeMETSignif,   "met_etsignif[10]/double");
  mAllData->Branch("met_CaloMETInmHF",      &mTreeMETCaloMETInmHF,   "met_CaloMETInmHF[10]/double");
  mAllData->Branch("met_CaloMETInpHF",      &mTreeMETCaloMETInpHF,   "met_CaloMETInpHF[10]/double");
  mAllData->Branch("met_CaloMETPhiInmHF",   &mTreeMETCaloMETPhiInmHF,   "met_CaloMETPhiInmHF[10]/double");
  mAllData->Branch("met_CaloMETPhiInpHF",   &mTreeMETCaloMETPhiInpHF,   "met_CaloMETPhiInpHF[10]/double");
  mAllData->Branch("met_CaloSETInmHF",      &mTreeMETCaloSETInmHF,   "met_CaloSETInmHF[10]/double");
  mAllData->Branch("met_CaloSETInpHF",      &mTreeMETCaloSETInpHF,   "met_CaloSETInpHF[10]/double");
  mAllData->Branch("met_emEtFraction",      &mTreeMETemEtFraction,   "met_emEtFraction[10]/double");
  mAllData->Branch("met_etFractionHadronic",&mTreeMETetFractionHadronic,   "met_etFractionHadronic[10]/double");
  mAllData->Branch("met_maxEtInEmTowers",   &mTreeMETmaxEtInEmTowers,   "met_maxEtInEmTowers[10]/double");
  mAllData->Branch("met_maxEtInHadTowers",  &mTreeMETmaxEtInHadTowers,   "met_maxEtInHadTowers[10]/double");
  mAllData->Branch("met_emEtInHF",          &mTreeMETemEtInHF,   "met_emEtInHF[10]/double");
  mAllData->Branch("met_emEtInEE",          &mTreeMETemEtInEE,   "met_emEtInEE[10]/double");
  mAllData->Branch("met_emEtInEB",          &mTreeMETemEtInEB,   "met_emEtInEB[10]/double");
  mAllData->Branch("met_hadEtInHF",         &mTreeMEThadEtInHF,   "met_hadEtInHF[10]/double");
  mAllData->Branch("met_hadEtInHE",         &mTreeMEThadEtInHE,   "met_hadEtInHE[10]/double");
  mAllData->Branch("met_hadEtInHO",         &mTreeMEThadEtInHO,   "met_hadEtInHO[10]/double");
  mAllData->Branch("met_hadEtInHB",         &mTreeMEThadEtInHB,   "met_hadEtInHB[10]/double");
  mAllData->Branch("met_ChargedEMEtFraction", &mTreeMETChargedEMEtFraction,   "met_ChargedEMEtFraction[10]/double");
  mAllData->Branch("met_ChargedHadEtFraction",&mTreeMETChargedHadEtFraction,   "met_ChargedHadEtFraction[10]/double");
  mAllData->Branch("met_MuonEtFraction",      &mTreeMETMuonEtFraction,   "met_MuonEtFraction[10]/double");
  mAllData->Branch("met_NeutralEMFraction", &mTreeMETNeutralEMFraction,   "met_NeutralEMFraction[10]/double");
  mAllData->Branch("met_NeutralHadEtFraction",&mTreeMETNeutralHadEtFraction,   "met_NeutralHadEtFraction[10]/double");
  mAllData->Branch("met_Type6EtFraction",     &mTreeMETType6EtFraction,   "met_Type6EtFraction[10]/double");
  mAllData->Branch("met_Type7EtFraction",     &mTreeMETType7EtFraction,   "met_Type7EtFraction[10]/double");

  // Calo Jets
  mAllData->Branch("calojet_n",     &mTreeNCalojet,   "calojet_n/I");  
  mAllData->Branch("calojet_E" ,     mTreeCaloJetE,      "calojet_E[calojet_n]/double");
  mAllData->Branch("calojet_Et",     mTreeCaloJetEt,     "calojet_Et[calojet_n]/double");
  mAllData->Branch("calojet_p",      mTreeCaloJetP,      "calojet_p[calojet_n]/double");
  mAllData->Branch("calojet_pt",     mTreeCaloJetPt,     "calojet_pt[calojet_n]/double");
  mAllData->Branch("calojet_pt_raw", mTreeCaloJetPtRaw,  "calojet_pt_raw[calojet_n]/double");
  mAllData->Branch("calojet_px",     mTreeCaloJetPx,     "calojet_px[calojet_n]/double");
  mAllData->Branch("calojet_py",     mTreeCaloJetPy,     "calojet_py[calojet_n]/double");
  mAllData->Branch("calojet_pz",     mTreeCaloJetPz,     "calojet_pz[calojet_n]/double");
  mAllData->Branch("calojet_eta",    mTreeCaloJetEta,    "calojet_eta[calojet_n]/double");
  mAllData->Branch("calojet_phi",    mTreeCaloJetPhi,    "calojet_phi[calojet_n]/double");
  mAllData->Branch("calojet_fem",    mTreeCaloJetFem,    "calojet_fem[calojet_n]/double");
  mAllData->Branch("calojet_fhad",   mTreeCaloJetFhad,   "calojet_fhad[calojet_n]/double");
  mAllData->Branch("calojet_btag",   mTreeCaloJetBtag,   "calojet_btag[calojet_n]/double");
  mAllData->Branch("calojet_charge", mTreeCaloJetCharge, "calojet_charge[calojet_n]/double");
  mAllData->Branch("calojet_fHPD",   mTreeCaloJetfhpd,   "calojet_fHPD[calojet_n]/double");
  mAllData->Branch("calojet_fRBX",   mTreeCaloJetfrbx,   "calojet_fRBX[calojet_n]/double");
  mAllData->Branch("calojet_n90hits",mTreeCaloJetn90hits,"calojet_n90hits[calojet_n]/I");
  mAllData->Branch("calojet_n90",    mTreeCaloJetn90,    "calojet_n90[calojet_n]/I");
  mAllData->Branch("calojet_flav",   mTreeCaloJetPart,   "calojet_flav[calojet_n]/I");
  mAllData->Branch("calojet_truth",  mTreeCaloJetTruth,  "calojet_truth[calojet_n]/I");
  mAllData->Branch("calojet_const",  mTreeCaloJetConst,  "calojet_const[calojet_n]/I");
  mAllData->Branch("calojet_ID",     mTreeCaloJetID,     "calojet_ID[calojet_n]/I");

  // PF Jets
  mAllData->Branch("pfjet_n",     &mTreeNPFjet,      "pfjet_n/I");  
  mAllData->Branch("pfjet_E" ,     mTreePFJetE,      "pfjet_E[pfjet_n]/double");
  mAllData->Branch("pfjet_Et",     mTreePFJetEt,     "pfjet_Et[pfjet_n]/double");
  mAllData->Branch("pfjet_p",      mTreePFJetP,      "pfjet_p[pfjet_n]/double");
  mAllData->Branch("pfjet_pt",     mTreePFJetPt,     "pfjet_pt[pfjet_n]/double");
  mAllData->Branch("pfjet_pt_raw", mTreePFJetPtRaw,  "pfjet_pt_raw[pfjet_n]/double");
  mAllData->Branch("pfjet_px",     mTreePFJetPx,     "pfjet_px[pfjet_n]/double");
  mAllData->Branch("pfjet_py",     mTreePFJetPy,     "pfjet_py[pfjet_n]/double");
  mAllData->Branch("pfjet_pz",     mTreePFJetPz,     "pfjet_pz[pfjet_n]/double");
  mAllData->Branch("pfjet_eta",    mTreePFJetEta,    "pfjet_eta[pfjet_n]/double");
  mAllData->Branch("pfjet_phi",    mTreePFJetPhi,    "pfjet_phi[pfjet_n]/double");
  mAllData->Branch("pfjet_btag",   mTreePFJetBtag,   "pfjet_btag[pfjet_n]/double");
  mAllData->Branch("pfjet_charge", mTreePFJetCharge, "pfjet_charge[pfjet_n]/double");
  mAllData->Branch("pfjet_n90",    mTreePFJetn90,    "pfjet_n90[pfjet_n]/I");
  mAllData->Branch("pfjet_flav",   mTreePFJetPart,   "pfjet_flav[pfjet_n]/I");
  mAllData->Branch("pfjet_truth",  mTreePFJetTruth,  "pfjet_truth[pfjet_n]/I");
  mAllData->Branch("pfjet_const",  mTreePFJetConst,  "pfjet_const[pfjet_n]/I");
  mAllData->Branch("pfjet_PFN",    mTreePFJetN,      "pfjet_PFN[pfjet_n][7]/I");
  mAllData->Branch("pfjet_PFF",    mTreePFJetF,      "pfjet_PFF[pfjet_n][7]/double");

  // Generator Jets
  mAllData->Branch("truthjet_n",   &mTreeNtruthjet,    "truthjet_n/I");  
  mAllData->Branch("truthjet_E" ,   mTreetruthJetE,    "truthjet_E[truthjet_n]/double");
  mAllData->Branch("truthjet_Et",   mTreetruthJetEt,   "truthjet_Et[truthjet_n]/double");
  mAllData->Branch("truthjet_p",    mTreetruthJetP,    "truthjet_p[truthjet_n]/double");
  mAllData->Branch("truthjet_pt",   mTreetruthJetPt,   "truthjet_pt[truthjet_n]/double");
  mAllData->Branch("truthjet_px",   mTreetruthJetPx,   "truthjet_px[truthjet_n]/double");
  mAllData->Branch("truthjet_py",   mTreetruthJetPy,   "truthjet_py[truthjet_n]/double");
  mAllData->Branch("truthjet_pz",   mTreetruthJetPz,   "truthjet_pz[truthjet_n]/double");
  mAllData->Branch("truthjet_eta",  mTreetruthJetEta,  "truthjet_eta[truthjet_n]/double");
  mAllData->Branch("truthjet_phi",  mTreetruthJetPhi,  "truthjet_phi[truthjet_n]/double");

  // Fat Jets
  mAllData->Branch("fatjet_n",       &mTreeNfatjet,       "fatjet_n/I");
  mAllData->Branch("fatjet_nsub",     mTreefatjetnsub,    "fatjet_nsub[fatjet_n]/I");
  mAllData->Branch("fatjet_pt",       mTreefatjetpt,      "fatjet_pt[fatjet_n]/double");
  mAllData->Branch("fatjet_px",       mTreefatjetpx,      "fatjet_px[fatjet_n]/double");
  mAllData->Branch("fatjet_py",       mTreefatjetpy,      "fatjet_py[fatjet_n]/double");
  mAllData->Branch("fatjet_pz",       mTreefatjetpz,      "fatjet_pz[fatjet_n]/double");
  mAllData->Branch("fatjet_E",        mTreefatjete,       "fatjet_E[fatjet_n]/double");
  mAllData->Branch("fatjet_eta",      mTreefatjeteta,     "fatjet_eta[fatjet_n]/double");
  mAllData->Branch("fatjet_phi",      mTreefatjetphi,     "fatjet_phi[fatjet_n]/double");
  mAllData->Branch("fatjet_sub_pt",   mTreefatjetsubpt,   "fatjet_sub_pt[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_px",   mTreefatjetsubpx,   "fatjet_sub_px[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_py",   mTreefatjetsubpy,   "fatjet_sub_py[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_pz",   mTreefatjetsubpz,   "fatjet_sub_pz[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_E",    mTreefatjetsube,    "fatjet_sub_E[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_eta",  mTreefatjetsubeta,  "fatjet_sub_eta[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_phi",  mTreefatjetsubphi,  "fatjet_sub_phi[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_fem",  mTreefatjetsubfem,  "fatjet_sub_fem[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_fhad", mTreefatjetsubfhad, "fatjet_sub_fhad[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_btag", mTreefatjetsubbtag, "fatjet_sub_btag[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_n90",  mTreefatjetsubn90,  "fatjet_sub_n90[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_fHPD", mTreefatjetsubfhpd, "fatjet_sub_fHPD[fatjet_n][10]/double");
  mAllData->Branch("fatjet_sub_fRBX", mTreefatjetsubfrbx, "fatjet_sub_fRBX[fatjet_n][10]/double");

  // SuperClusters
  mAllData->Branch("SC_n",         &mTreeNSC,          "SC_n/I");
  mAllData->Branch("SC_truth",      mTreeSCTruth,      "SC_truth[SC_n]/I");
  mAllData->Branch("SC_E",          mTreeSCE,          "SC_E[SC_n]/double");
  mAllData->Branch("SC_phi",        mTreeSCPhi,        "SC_phi[SC_n]/double");
  mAllData->Branch("SC_eta",        mTreeSCEta,        "SC_eta[SC_n]/double");
  
  // Electrons
  mAllData->Branch("ele_n",         &mTreeNele,          "ele_n/I");  
  mAllData->Branch("ele_E",          mTreeEleE,          "ele_E[ele_n]/double");
  mAllData->Branch("ele_Et",         mTreeEleEt,         "ele_Et[ele_n]/double");
  mAllData->Branch("ele_p",          mTreeEleP,          "ele_p[ele_n]/double");
  mAllData->Branch("ele_pt",         mTreeElePt,         "ele_pt[ele_n]/double");
  mAllData->Branch("ele_px",         mTreeElePx,         "ele_px[ele_n]/double");
  mAllData->Branch("ele_py",         mTreeElePy,         "ele_py[ele_n]/double");
  mAllData->Branch("ele_pz",         mTreeElePz,         "ele_pz[ele_n]/double");
  mAllData->Branch("ele_eta",        mTreeEleEta,        "ele_eta[ele_n]/double");
  mAllData->Branch("ele_phi",        mTreeElePhi,        "ele_phi[ele_n]/double");
  mAllData->Branch("ele_charge",     mTreeEleCharge,     "ele_charge[ele_n]/double");
  mAllData->Branch("ele_TrkChiNorm", mTreeEleTrkChiNorm, "ele_TrkChiNorm[ele_n]/double");
  mAllData->Branch("ele_d0vtx",      mTreeEled0vtx,      "ele_d0vtx[ele_n]/double");
  mAllData->Branch("ele_d0bs",       mTreeEled0bs,       "ele_d0bs[ele_n]/double");
  mAllData->Branch("ele_sd0",        mTreeElesd0,        "ele_sd0[ele_n]/double");
  mAllData->Branch("ele_hits",       mTreeEleHits,       "ele_hits[ele_n]/I");
  mAllData->Branch("ele_truth",      mTreeEleTruth,      "ele_truth[ele_n]/I");
  mAllData->Branch("ele_isECal",     mTreeEleisECal,     "ele_isECal[ele_n]/I");
  mAllData->Branch("ele_isTracker",  mTreeEleisTracker,  "ele_isTracker[ele_n]/I");
  mAllData->Branch("ele_ID",         mTreeEleID,         "ele_ID[ele_n][5]/I");
  mAllData->Branch("ele_ValidHitFirstPxlB",mTreeEleValidHitFirstPxlB,              "ele_ValidHitFirstPxlB[ele_n]/I");
  mAllData->Branch("ele_TrkExpHitsInner",  mTreeEleTrkExpHitsInner,                "ele_EleTrkExpHitsInner[ele_n]/I");
  mAllData->Branch("ele_HCalOverEm",       mTreeEleHCalOverEm,                     "ele_HCalOverEm[ele_n]/double");
  mAllData->Branch("ele_Dr03TkSumPt",      mTreeEleDr03TkSumPt,                    "ele_Dr03TkSumPt[ele_n]/double");
  mAllData->Branch("ele_Dr04HCalSumEt",    mTreeEleDr04HCalTowerSumEt,             "ele_Dr04HCalSumEt[ele_n]/double");
  mAllData->Branch("ele_Dr03HCalSumEt",    mTreeEleDr03HCalTowerSumEt,             "ele_Dr03HCalSumEt[ele_n]/double");
  mAllData->Branch("ele_Dr04ECalSumEt",    mTreeEleDr04ECalRecHitSumEt,            "ele_Dr04ECalSumEt[ele_n]/double");
  mAllData->Branch("ele_Dr03ECalSumEt",    mTreeEleDr03ECalRecHitSumEt,            "ele_Dr03ECalSumEt[ele_n]/double");
  mAllData->Branch("ele_SigmaIetaIeta",    mTreeEleSigmaIetaIeta,                  "ele_SigmaIetaIeta[ele_n]/double");
  mAllData->Branch("ele_dEtaSCTrackAtVtx", mTreeEleDeltaEtaSuperClusterTrackAtVtx, "ele_dEtaSCTrackAtVtx[ele_n]/double");
  mAllData->Branch("ele_dPhiSCTrackAtVtx", mTreeEleDeltaPhiSuperClusterTrackAtVtx, "ele_dPhiSCTrackAtVtx[ele_n]/double");
  mAllData->Branch("ele_dr03HcalDepth1",   mTreeEledr03HcalDepth1,                 "ele_dr03HcalDepth1[ele_n]/double");
  mAllData->Branch("ele_dr03HcalDepth2",   mTreeEledr03HcalDepth2,                 "ele_dr03HcalDepth2[ele_n]/double");
  mAllData->Branch("ele_e2x5Max",          mTreeElee2x5Max,                        "ele_e2x5Max[ele_n]/double");
  mAllData->Branch("ele_e5x5",             mTreeElee5x5,                           "ele_e5x5[ele_n]/double");
  mAllData->Branch("ele_e1x5",             mTreeElee1x5,                           "ele_e1x5[ele_n]/double");
  mAllData->Branch("ele_caloEt",           mTreeEleCaloEt,                         "ele_caloEt[ele_n]/double");
  mAllData->Branch("ele_SCeta",            mTreeEleSCEta,                          "ele_SCeta[ele_n]/double");
  mAllData->Branch("ele_convdist", mTreeEleConvdist, "ele_convdist[ele_n]/double");
  mAllData->Branch("ele_convdcot", mTreeEleConvdcot, "ele_convdcot[ele_n]/double");
  mAllData->Branch("ele_convr",    mTreeEleConvr,    "ele_convr[ele_n]/double");
  mAllData->Branch("ele_fbrem",    mTreeElefbrem,    "ele_fbrem[ele_n]/double");
  mAllData->Branch("ele_trign",    mTreeNeletrign,   "ele_trign[ele_n]/I");
  mAllData->Branch("ele_trig",     mTreeEletrig,     "ele_trig[ele_n][500]/I");  
  mAllData->Branch("ele_SC" ,      mTreeEleSC,       "ele_SC[ele_n]/I");
  mAllData->Branch("ele_numberOfHits",mTreeElenumberOfHits, "ele_numberOfHits[ele_n]/I");

  // PF Electrons
  mAllData->Branch("pfele_n",      &mTreeNPFEle,      "pfele_n/I");
  mAllData->Branch("pfele_p",       mTreePFEleP,      "pfele_p[pfele_n]/double");
  mAllData->Branch("pfele_E",       mTreePFEleE,      "pfele_E[pfele_n]/double");
  mAllData->Branch("pfele_Et",      mTreePFEleEt,     "pfele_Et[pfele_n]/double");
  mAllData->Branch("pfele_pt",      mTreePFElePt,     "pfele_pt[pfele_n]/double");
  mAllData->Branch("pfele_px",      mTreePFElePx,     "pfele_px[pfele_n]/double");
  mAllData->Branch("pfele_py",      mTreePFElePy,     "pfele_py[pfele_n]/double");
  mAllData->Branch("pfele_pz",      mTreePFElePz,     "pfele_pz[pfele_n]/double");
  mAllData->Branch("pfele_eta",     mTreePFEleEta,    "pfele_eta[pfele_n]/double");
  mAllData->Branch("pfele_phi",     mTreePFElePhi,    "pfele_phi[pfele_n]/double");
  mAllData->Branch("pfele_charge",  mTreePFEleCharge, "pfele_charge[pfele_n]/double");
  mAllData->Branch("pfele_truth",   mTreePFEleTruth,  "pfele_truth[pfele_n]/I");
  mAllData->Branch("pfele_trign",   mTreeNPFEletrign, "pfele_trign[pfele_n]/I");
  mAllData->Branch("pfele_trig",    mTreePFEletrig,   "pfele_trig[pfele_n][500]/I");
  mAllData->Branch("pfele_SC",      mTreePFEleSC,     "pfele_SC[pfele_n]/I");

  // Muons
  mAllData->Branch("muo_n" ,           &mTreeNmuo,            "muo_n/I");  
  mAllData->Branch("muo_E" ,            mTreeMuoE,            "muo_E[muo_n]/double");
  mAllData->Branch("muo_Et",            mTreeMuoEt,           "muo_Et[muo_n]/double");
  mAllData->Branch("muo_p",             mTreeMuoP,            "muo_p[muo_n]/double");
  mAllData->Branch("muo_pt",            mTreeMuoPt,           "muo_pt[muo_n]/double");
  mAllData->Branch("muo_px",            mTreeMuoPx,           "muo_px[muo_n]/double");
  mAllData->Branch("muo_py",            mTreeMuoPy,           "muo_py[muo_n]/double");
  mAllData->Branch("muo_pz",            mTreeMuoPz,           "muo_pz[muo_n]/double");
  mAllData->Branch("muo_eta",           mTreeMuoEta,          "muo_eta[muo_n]/double");
  mAllData->Branch("muo_phi",           mTreeMuoPhi,          "muo_phi[muo_n]/double");
  mAllData->Branch("muo_charge",        mTreeMuoCharge,       "muo_charge[muo_n]/double");
  mAllData->Branch("muo_RelTrkIso",     mTreeMuoRelTrkIso,    "muo_RelTrkIso[muo_n]/double");
  mAllData->Branch("muo_TrkIso",        mTreeMuoTrkIso,       "muo_TrkIso[muo_n]/double");
  mAllData->Branch("muo_ECalIso",       mTreeMuoECalIso,      "muo_ECalIso[muo_n]/double");
  mAllData->Branch("muo_HCalIso",       mTreeMuoHCalIso,      "muo_HCalIso[muo_n]/double");
  mAllData->Branch("muo_TrkIsoDep",     mTreeMuoTrkIsoDep,    "muo_TrkIsoDep[muo_n]/double");
  mAllData->Branch("muo_ECalIsoDep",    mTreeMuoECalIsoDep,   "muo_ECalIsoDep[muo_n]/double");
  mAllData->Branch("muo_HCalIsoDep",    mTreeMuoHCalIsoDep,   "muo_HCalIsoDep[muo_n]/double");
  mAllData->Branch("muo_AllIso",        mTreeMuoAllIso,       "muo_AllIso[muo_n]/double");
  mAllData->Branch("muo_TrkChiNormCm",  mTreeMuoTrkChiNormCm, "muo_TrkChiNormCm[muo_n]/double");
  mAllData->Branch("muo_TrkChiNormTk",  mTreeMuoTrkChiNormTk, "muo_TrkChiNormTk[muo_n]/double");
  mAllData->Branch("muo_d0Cm",          mTreeMuod0Cm,         "muo_d0Cm[muo_n]/double");
  mAllData->Branch("muo_d0Tk",          mTreeMuod0Tk,         "muo_d0Tk[muo_n]/double");
  mAllData->Branch("muo_sd0Cm",         mTreeMuosd0Cm,        "muo_sd0Cm[muo_n]/double");
  mAllData->Branch("muo_sd0Tk",         mTreeMuosd0Tk,        "muo_sd0Tk[muo_n]/double");
  mAllData->Branch("muo_calocomp",      mTreeMuocalocomp,     "muo_calocomp[muo_n]/double");
  mAllData->Branch("muo_calotower_e",   mTreeMuocaltowe,      "muo_calotower_e[muo_n]/double");
  mAllData->Branch("muo_prompttight",   mTreeMuoGood,         "muo_prompttight[muo_n]/I");
  mAllData->Branch("muo_hitsCm",        mTreeMuoHitsCm,       "muo_hitsCm[muo_n]/I");
  mAllData->Branch("muo_hitsTk",        mTreeMuoHitsTk,       "muo_hitsTk[muo_n]/I");
  mAllData->Branch("muo_truth",         mTreeMuoTruth,        "muo_truth[muo_n]/I");
  mAllData->Branch("muo_trign",         mTreeNmuotrign,       "muo_trign[muo_n]/I");
  mAllData->Branch("muo_trig" ,         mTreeMuotrig,         "muo_trig[muo_n][500]/I");  
  mAllData->Branch("muo_ID",            mTreeMuoID,           "muo_ID[muo_n][24]/I");
  mAllData->Branch("muo_ValidMuonHitsCm",        mTreeMuoValidMuonHitsCm,       "muo_ValidMuonHitsCm[muo_n]/I");
  mAllData->Branch("muo_ValidTrackerHitsCm",     mTreeMuoValidTrackerHitsCm,    "muo_ValidTrackerHitsCm[muo_n]/I");
  mAllData->Branch("muo_ValidPixelHitsCm",       mTreeMuoValidPixelHitsCm,      "muo_ValidPixelHitsCm[muo_n]/I");
  mAllData->Branch("muo_ChambersMatched",        mTreeMuoChambersMatched,       "muo_ChambersMatched[muo_n]/I");
  mAllData->Branch("muo_d0bsCm",                 mTreeMuod0bsCm,                "muo_d0bsCm[muo_n]/double");
  mAllData->Branch("muo_d0OriginCm",             mTreeMuod0OriginCm,            "muo_d0OriginCm[muo_n]/double");
  mAllData->Branch("muo_dzbsCm",                 mTreeMuodzbsCm,                "muo_dzbsCm[muo_n]/double");
  mAllData->Branch("muo_TrackerLayersMeasCm",    mTreeMuoTrackerLayersMeasCm,   "muo_TrackerLayersMeasCm[muo_n]/I");
  mAllData->Branch("muo_TrackerLayersNotMeasCm", mTreeMuoTrackerLayersNotMeasCm,"muo_TrackerLayersNotMeasCm[muo_n]/I");  
  mAllData->Branch("muo_Cocktail_pt",  mTreeMuoCocktailPt,  "muo_Cocktail_pt[muo_n]/double");   
  mAllData->Branch("muo_Cocktail_phi", mTreeMuoCocktailPhi, "muo_Cocktail_phi[muo_n]/double");
  mAllData->Branch("muo_Cocktail_eta", mTreeMuoCocktailEta, "muo_Cocktail_eta[muo_n]/double");  
  mAllData->Branch("muo_Valid_fraction",mTreeMuoValidFraction,"muo_Valid_fraction[muo_n]/double");


  mAllData->Branch("PFmuo_n",   &mTreeNPFmuons,    "PFmuo_n/I");  
  
  mAllData->Branch("PFmuo_p",       mTreePFMuonP,  "PFmuo_p[PFmuo_n]/double");
  mAllData->Branch("PFmuo_pt",      mTreePFMuonPt,  "PFmuo_pt[PFmuo_n]/double");
  mAllData->Branch("PFmuo_E",       mTreePFMuonE,  "PFmuo_E[PFmuo_n]/double");
  mAllData->Branch("PFmuo_Et",      mTreePFMuonEt,  "PFmuo_Et[PFmuo_n]/double");
  mAllData->Branch("PFmuo_px",      mTreePFMuonPx,  "PFmuo_px[PFmuo_n]/double");
  mAllData->Branch("PFmuo_py",      mTreePFMuonPy,  "PFmuo_py[PFmuo_n]/double");
  mAllData->Branch("PFmuo_pz",      mTreePFMuonPz,  "PFmuo_pz[PFmuo_n]/double");
  mAllData->Branch("PFmuo_eta",     mTreePFMuonEta,  "PFmuo_eta[PFmuo_n]/double");
  mAllData->Branch("PFmuo_phi",     mTreePFMuonPhi,  "PFmuo_phi[PFmuo_n]/double");
  mAllData->Branch("PFmuo_Charge",  mTreePFMuonCharge,  "PFmuo_Charge[PFmuo_n]/double");
  mAllData->Branch("PFmuo_particleIso",   mTreePFMuonParticleIso,  "PFmuo_particleIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_chadIso",       mTreePFMuonChadIso,      "PFmuo_chadIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_nhadIso",       mTreePFMuonNhadIso,      "PFmuo_nhadIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_gamIso",        mTreePFMuonGamIso,       "PFmuo_gamIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_RelTrkIso",     mTreePFMuonRelTrkIso,    "PFmuo_RelTrkIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_TrkIso",        mTreePFMuonTrkIso,       "PFmuo_TrkIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_ECalIso",       mTreePFMuonECalIso,      "PFmuo_ECalIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_HCalIso",       mTreePFMuonHCalIso,      "PFmuo_HCalIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_TrkIsoDep",     mTreePFMuoTrkIsoDep,    "PFmuo_TrkIsoDep[PFmuo_n]/double");
  mAllData->Branch("PFmuo_ECalIsoDep",    mTreePFMuoECalIsoDep,   "PFmuo_ECalIsoDep[PFmuo_n]/double");
  mAllData->Branch("PFmuo_HCalIsoDep",    mTreePFMuoHCalIsoDep,   "PFmuo_HCalIsoDep[PFmuo_n]/double");
  mAllData->Branch("PFmuo_AllIso",        mTreePFMuonAllIso,       "PFmuo_AllIso[PFmuo_n]/double");
  mAllData->Branch("PFmuo_TrkChiNormCm",  mTreePFMuoTrkChiNormCm, "PFmuo_TrkChiNormCm[PFmuo_n]/double");
  mAllData->Branch("PFmuo_TrkChiNormTk",  mTreePFMuoTrkChiNormTk, "PFmuo_TrkChiNormTk[PFmuo_n]/double");
  mAllData->Branch("PFmuo_d0Cm",          mTreePFMuod0Cm,         "PFmuo_d0Cm[PFmuo_n]/double");
  mAllData->Branch("PFmuo_d0Tk",          mTreePFMuod0Tk,         "PFmuo_d0Tk[PFmuo_n]/double");
  mAllData->Branch("PFmuo_sd0Cm",         mTreePFMuosd0Cm,        "PFmuo_sd0Cm[PFmuo_n]/double");
  mAllData->Branch("PFmuo_sd0Tk",         mTreePFMuosd0Tk,        "PFmuo_sd0Tk[PFmuo_n]/double");
  mAllData->Branch("PFmuo_calocomp",      mTreePFMuocalocomp,     "PFmuo_calocomp[PFmuo_n]/double");
  mAllData->Branch("PFmuo_calotower_e",   mTreePFMuocaltowe,      "PFmuo_calotower_e[PFmuo_n]/double");
  mAllData->Branch("PFmuo_hitsCm",        mTreePFMuoHitsCm,       "PFmuo_hitsCm[PFmuo_n]/I");
  mAllData->Branch("PFmuo_hitsTk",        mTreePFMuoHitsTk,       "PFmuo_hitsTk[PFmuo_n]/I");
  mAllData->Branch("PFmuo_ValidMuonHitsCm",        mTreePFMuoValidMuonHitsCm,       "PFmuo_ValidMuonHitsCm[PFmuo_n]/I");
  mAllData->Branch("PFmuo_ValidTrackerHitsCm",     mTreePFMuoValidTrackerHitsCm,    "PFmuo_ValidTrackerHitsCm[PFmuo_n]/I");
  mAllData->Branch("PFmuo_ValidPixelHitsCm",       mTreePFMuoValidPixelHitsCm,      "PFmuo_ValidPixelHitsCm[PFmuo_n]/I");
  mAllData->Branch("PFmuo_ChambersMatched",        mTreePFMuoChambersMatched,       "PFmuo_ChambersMatched[PFmuo_n]/I");
  mAllData->Branch("PFmuo_d0bsCm",                 mTreePFMuod0bsCm,                "PFmuo_d0bsCm[PFmuo_n]/double");
  mAllData->Branch("PFmuo_d0OriginCm",             mTreePFMuod0OriginCm,            "PFmuo_d0OriginCm[PFmuo_n]/double");
  mAllData->Branch("PFmuo_dzbsCm",                 mTreePFMuodzbsCm,                "PFmuo_dzbsCm[PFmuo_n]/double");
  mAllData->Branch("PFmuo_TrackerLayersMeasCm",    mTreePFMuoTrackerLayersMeasCm,   "PFmuo_TrackerLayersMeasCm[PFmuo_n]/I");
  mAllData->Branch("PFmuo_TrackerLayersNotMeasCm", mTreePFMuoTrackerLayersNotMeasCm,"PFmuo_TrackerLayersNotMeasCm[PFmuo_n]/I");  
  mAllData->Branch("PFmuo_Valid_fraction",         mTreePFMuoValidFraction,         "PFmuo_Valid_fraction[PFmuo_n]/double");
  
  
  //Taus
  mAllData->Branch("tau_n",   &mTreeNtaus,    "tau_n/I");
  
  mAllData->Branch("tau_p",   mTreeTauP,      "tau_p[tau_n]/double");
  mAllData->Branch("tau_pt", mTreeTauPt, "tau_pt[tau_n]/double");
  mAllData->Branch("tau_E", mTreeTauE, "tau_E[tau_n]/double");
  mAllData->Branch("tau_Et", mTreeTauEt, "tau_Et[tau_n]/double"); 
  mAllData->Branch("tau_Px", mTreeTauPx, "tau_Px[tau_n]/double"); 
  mAllData->Branch("tau_Py", mTreeTauPy, "tau_Py[tau_n]/double"); 
  mAllData->Branch("tau_Pz", mTreeTauPz, "tau_Pz[tau_n]/double"); 
  mAllData->Branch("tau_Eta", mTreeTauEta, "tau_Eta[tau_n]/double"); 
  mAllData->Branch("tau_Phi", mTreeTauPhi, "tau_Phi[tau_n]/double"); 
  mAllData->Branch("tau_DecayMode",   mTreeTauDecayMode,    "tau_DecayMode[tau_n]/I");   
  mAllData->Branch("tau_vx", mTreeTauvx, "tau_vx[tau_n]/double"); 
  mAllData->Branch("tau_vy", mTreeTauvy, "tau_vy[tau_n]/double"); 
  mAllData->Branch("tau_vz", mTreeTauvz, "tau_vz[tau_n]/double"); 
  mAllData->Branch("tau_vx2", mTreeTauvx2, "tau_vx2[tau_n]/double"); 
  mAllData->Branch("tau_vy2", mTreeTauvy2, "tau_vy2[tau_n]/double"); 
  mAllData->Branch("tau_vz2", mTreeTauvz2, "tau_vz2[tau_n]/double"); 
  mAllData->Branch("tau_ECalIso", mTreeTauECalIso, "tau_ECalIso[tau_n]/double"); 
  mAllData->Branch("tau_HCalIso", mTreeTauHCalIso, "tau_HCalIso[tau_n]/double"); 
  mAllData->Branch("tau_AllIso", mTreeTauAllIso, "tau_AllIso[tau_n]/double"); 
  mAllData->Branch("tau_TrackIso", mTreeTauTrackIso, "tau_TrackIso[tau_n]/double"); 
  mAllData->Branch("tau_ParticleIso", mTreeTauParticleIso, "tau_ParticleIso[tau_n]/double"); 
  mAllData->Branch("tau_ChadIso", mTreeTauChadIso, "tau_ChadIso[tau_n]/double"); 
  mAllData->Branch("tau_NhadIso", mTreeTauNhadIso, "tau_NhadIso[tau_n]/double"); 
  mAllData->Branch("tau_GamIso", mTreeTauGamIso, "tau_GamIso[tau_n]/double"); 
  mAllData->Branch("tau_id", mTreeTauID, "tau_id[tau_n][16]/double");
  
  // SUSY
  mAllData->Branch("susyScanM0",  &mTreesusyScanM0, "susyScanM0/double");
  mAllData->Branch("susyScanM12",  &mTreesusyScanM12, "susyScanM12/double");
  mAllData->Branch("susyScanA0",  &mTreesusyScanA0, "susyScanA0/double");
  mAllData->Branch("susyScanCrossSection",  &mTreesusyScanCrossSection, "susyScanCrossSection/double");
  mAllData->Branch("susyScanMu",  &mTreesusyScanMu, "susyScanMu/double");
  mAllData->Branch("susyScanRun",  &mTreesusyScanRun, "susyScanRun/double");
  mAllData->Branch("susyScantanbeta",  &mTreesusyScantanbeta, "susyScantanbeta/double");
     

  // DataFormats/MuonReco/interface/MuonSelectors.h
  ACmuonID[0]  = "All";                     // dummy options - always true
  ACmuonID[1]  = "AllGlobalMuons";          // checks isGlobalMuon flag
  ACmuonID[2]  = "AllStandAloneMuons";      // checks isStandAloneMuon flag
  ACmuonID[3]  = "AllTrackerMuons";         // checks isTrackerMuon flag
  ACmuonID[4]  = "TrackerMuonArbitrated";   // resolve ambiguity of sharing segments
  ACmuonID[5]  = "AllArbitrated";           // all muons with the tracker muon arbitrated
  ACmuonID[6]  = "GlobalMuonPromptTight";   // global muons with tighter fit requirements
  ACmuonID[7]  = "TMLastStationLoose";      // penetration depth loose selector
  ACmuonID[8]  = "TMLastStationTight";      // penetration depth tight selector
  ACmuonID[9]  = "TM2DCompatibilityLoose";  // likelihood based loose selector
  ACmuonID[10] = "TM2DCompatibilityTight";  // likelihood based tight selector
  ACmuonID[11] = "TMOneStationLoose";       // require one well matched segment
  ACmuonID[12] = "TMOneStationTight";       // require one well matched segment
  ACmuonID[13] = "TMLastStationOptimizedLowPtLoose"; // combination of TMLastStation and TMOneStation
  ACmuonID[14] = "TMLastStationOptimizedLowPtTight"; // combination of TMLastStation and TMOneStation
  ACmuonID[15] = "GMTkChiCompatibility";    // require tk stub have good chi2 relative to glb track
  ACmuonID[16] = "GMStaChiCompatibility";   // require sta stub have good chi2 compatibility relative to glb track
  ACmuonID[17] = "GMTkKinkTight";           // require a small kink value in the tracker stub
  ACmuonID[18] = "TMLastStationAngLoose";   // TMLastStationLoose with additional angular cuts
  ACmuonID[19] = "TMLastStationAngTight";   // TMLastStationTight with additional angular cuts
  ACmuonID[20] = "TMOneStationAngLoose";    // TMOneStationLoose with additional angular cuts
  ACmuonID[21] = "TMOneStationAngTight";    // TMOneStationTight with additional angular cuts
  ACmuonID[22] = "TMLastStationOptimizedBarrelLowPtLoose"; 
			  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  ACmuonID[23] = "TMLastStationOptimizedBarrelLowPtTight";
			  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
        
        
  ACtauID[0]="againstElectron";
  ACtauID[1]="againstMuon";
  ACtauID[2]="byIsolation";
  ACtauID[3]="byIsolationUsingLeadingPion";
  ACtauID[4]="byTaNC";
  ACtauID[5]="byTaNCfrHalfPercent";
  ACtauID[6]="byTaNCfrOnePercent";
  ACtauID[7]="byTaNCfrQuarterPercent";
  ACtauID[8]="byTaNCfrTenthPercent";
  ACtauID[9]="ecalIsolation";
  ACtauID[10]="ecalIsolationUsingLeadingPion";
  ACtauID[11]="leadingPionPtCut";
  ACtauID[12]="leadingTrackFinding";
  ACtauID[13]="leadingTrackPtCut";
  ACtauID[14]="trackIsolation";
  ACtauID[15]="trackIsolationUsingLeadingPion";


}

////////////////////////////////
//
// Helpers
//
double SusyACSkimAnalysis::DeltaPhi(double a, double b) {
  double temp = fabs(a-b);
  if (temp <= localPi)
    return temp;
  else
    return  2.*localPi - temp;
}
double SusyACSkimAnalysis::DeltaR(double a, double b, double c, double d ) {
  double dr = sqrt( (a-b)*(a-b)+ SusyACSkimAnalysis::DeltaPhi(c,d)*SusyACSkimAnalysis::DeltaPhi(c,d));
  return dr; 
}

// Define this as a plug-in
DEFINE_FWK_MODULE(SusyACSkimAnalysis);
