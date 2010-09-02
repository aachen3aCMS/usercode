//
// Package:    UserCode/aachen3a/ACSusyAnalysis
// Class:      SusyACSkimAnalysis
// 
// Description: Skeleton analysis for SUSY search with Lepton + Jets + MET
//
// Original Author:  Carsten Magass
//         Created:  November 2008
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
  jetTag_    = iConfig.getParameter<edm::InputTag>("jetTag");
  metTag_    = iConfig.getParameter<edm::InputTag>("metTag");
  metTagPF_  = iConfig.getParameter<edm::InputTag>("metTagPF");
  metTagTC_  = iConfig.getParameter<edm::InputTag>("metTagTC");
  elecTag_   = iConfig.getParameter<edm::InputTag>("elecTag");
  muonTag_   = iConfig.getParameter<edm::InputTag>("muonTag");
  genTag_    = iConfig.getParameter<edm::InputTag>("genTag");
  genJetTag_ = iConfig.getParameter<edm::InputTag>("genJetTag");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vtxTag");

  is_MC      = iConfig.getParameter<bool>("is_MC");
  is_SHERPA  = iConfig.getParameter<bool>("is_SHERPA");
  do_fatjets = iConfig.getParameter<bool>("do_fatjets");

  edm::LogVerbatim("SusyACSkimAnalysis") << " Running with flag is_MC      = " << is_MC << endl;
  edm::LogVerbatim("SusyACSkimAnalysis") << " Running with flag is_SHERPA  = " << is_SHERPA << endl;
  edm::LogVerbatim("SusyACSkimAnalysis") << " Running with flag do_fatjets = " << do_fatjets << endl;

  cor_  = iConfig.getParameter<std::string>("correction");
  flav_ = iConfig.getParameter<std::string>("flavour");
  btag_ = iConfig.getParameter<std::string>("btag");

  correction_ = pat::JetCorrFactors::corrStep(cor_, flav_);
 
  edm::LogVerbatim("SusyACSkimAnalysis") << " Using the following Jet correction: " << cor_ << ", flavour: " << flav_ << endl;

  // get the cuts
  muopt_  = iConfig.getParameter<double>("muopt");
  muoeta_ = iConfig.getParameter<double>("muoeta");
  elept_  = iConfig.getParameter<double>("elept");
  eleeta_ = iConfig.getParameter<double>("eleeta");
  jetpt_  = iConfig.getParameter<double>("jetpt");
  jeteta_ = iConfig.getParameter<double>("jeteta");
  met_    = iConfig.getParameter<double>("met");

  nele_ = iConfig.getParameter<int>("nele");
  nmuo_ = iConfig.getParameter<int>("nmuo");
  njet_ = iConfig.getParameter<int>("njet");

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
// Called in for each event
//
bool SusyACSkimAnalysis::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace pat;
  
  mTreeNtrig = 0;

  Handle< pat::TriggerEvent > myTriggerEvent;
  iEvent.getByLabel( "patTriggerEvent", myTriggerEvent );

  string trigname =":";

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

      trigname = trigname + tname + ":";

      // save only muon trigger results
      if ( (ttname.Contains("HLT_Mu") || ttname.Contains("HLT_IsoMu") || 
	    ttname.Contains("HLT_L1Mu") || ttname.Contains("HLT_L2Mu") ||
	    ttname.Contains("HLT_DoubleMu")  || 
	    ttname.Contains("HLT_Ele") || ttname.Contains("HLT_L1DoubleEG") ||
	    ttname.Contains("HLT_Photon") || ttname.Contains("HLT_L1SingleEG") ||
	    ttname.Contains("HLT_DoubleEle") || ttname.Contains("HLT_DoublePhoton")) &&
	   !ttname.Contains("Jet") && !ttname.Contains("Tau") &&
	   myTriggerEvent->path(tname)->wasAccept() ) {

	//	cout << "   --> " << ttname << endl;
	
	int *tempname = pack(tname.c_str());

	const TriggerFilterRefVector mpf = myTriggerEvent->pathFilters(tname);
	for ( TriggerFilterRefVector::const_iterator ll=mpf.begin(); ll!=mpf.end(); ++ll ) {
	  TriggerObjectRefVector torv = myTriggerEvent->filterObjects((*ll)->label());   
	  for ( TriggerObjectRefVector::const_iterator itt = torv.begin(); 
		itt != torv.end(); ++itt ) {
	    const TriggerObjectRef objRef( *itt );
	    
	    //	    cout << "  " << (*ll)->label() << "  " << (*itt)->pt() << "  "  << (*itt)->eta() << endl;
	    int *filtname = pack((*ll)->label().c_str());

	    for (int l=0; l<get_size(tempname); l++) mTreetrigname[mTreeNtrig][l] = tempname[l];
	    for (int l=0; l<get_size(filtname); l++) mTreefiltname[mTreeNtrig][l] = filtname[l];
	    
	    mTreetrigpre[mTreeNtrig] = (*it)->prescale();
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
  strcpy(mTreeHLT, trigname.c_str());

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


  /*
  edm::LogInfo("SusyACSkimAnalysis") << " Event properties: EvtWeight = " << mTreeEventWeight 
				     << "  ProcID = " << mTreeProcID 
				     << "  Pthat = " <<  mTreePthat << std::endl;
  */

  // HCAL Noise
  Handle<HcalNoiseSummary> noiseHandle;
  iEvent.getByLabel("hcalnoise", noiseHandle);
  const HcalNoiseSummary noisesummary = *noiseHandle;

  if ( !noiseHandle.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No HcalNoiseSummary found  ";
  }
  else {

    mTreenoisel = noisesummary.passLooseNoiseFilter();
    mTreenoiset = noisesummary.passTightNoiseFilter();
    mTreenoiseh = noisesummary.passHighLevelNoiseFilter();

  }

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

  if (tkRef.isValid()) {
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
       
	if (p.status() == 1 && (abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)) {
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
		const reco::Candidate &t = *((*part1[d]).daughter(j));

		//		if (abs((*part1[d]).pdgId()) == 15) { // tau
		//		  cout << " found tau with daughter " << t.pdgId() << " and status " << t.status() << endl;
		//		}
		if (t.status()!=3 && abs((*part1[d]).pdgId()) != 15) continue;

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
    edm::LogWarning("SusyACSkimAnalysis") <<  "no correctedHybridSuperClusters" << endl;
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
    edm::LogWarning("SusyACSkimAnalysis") <<  "no correctedMulti5x5SuperClustersWithPreshower" << endl;
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

      if (eles[i].pt() < elept_ || fabs(eles[i].eta()) > eleeta_) continue;

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
	if (DeltaPhi(eles[i].phi(), mTreetrigphi[k])<0.2 &&
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

      // Check for conversions
      reco::GsfElectron el = eles[i];
      
      ConversionFinder convFinder;
      ConversionInfo convInfo = convFinder.getConversionInfo(el, tkRef, bfield);
      
      mTreeEleConvdist[countele] = convInfo.dist();
      mTreeEleConvdcot[countele] = convInfo.dcot();
      mTreeEleConvr[countele]    = convInfo.radiusOfConversion();

      if (eles[i].gsfTrack().isNonnull()) {

	mTreeEleHits[countele]  = eles[i].gsfTrack().get()->numberOfValidHits();
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

  // Muons
  mTreeNmuo = 0;

  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonTag_, muonHandle);

  std::vector<pat::Muon> muons;

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

      if (muons[i].pt() < muopt_ || fabs(muons[i].eta()) > muoeta_ ) continue;

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
	if (DeltaPhi(muons[i].phi(), mTreetrigphi[k])<0.1 &&
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

      // combined muon
      if ( muons[i].combinedMuon().isNonnull() ) {

	mTreeMuoHitsCm[countmuo] = muons[i].combinedMuon().get()->numberOfValidHits();
	mTreeMuod0Cm[countmuo]   = (-1.)* muons[i].combinedMuon().get()->dxy(vtxPoint);
	mTreeMuosd0Cm[countmuo]  = muons[i].combinedMuon().get()->d0Error();

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



      countmuo++;
    }
    mTreeNmuo = countmuo;
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

    mTreeNjet = 0;

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

  // Jets
  mTreeNjet = 0;

  edm::Handle< std::vector<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetTag_, jetHandle);

  std::vector<pat::Jet> jets;

  if ( !jetHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Jet results found for InputTag " << jetTag_;
  else {
    
    int countjets = 0;

    mTreeNjet = jetHandle->size();

    for (int i=0; i<mTreeNjet; i++) {
      jets.push_back((*jetHandle)[i]);
    }

    sort(jets.begin(), jets.end(), ptcomp_jet);  

    if ( mTreeNjet > 100 ) mTreeNjet = 100;
    
    for (int i=0; i<mTreeNjet; i++) {
      
      const pat::Jet& jet = jets[i];
      float correction = 1.;

      /*
      // result of all b tagging algorithms
      const vector<pair<string, float> > bvec = jet.getPairDiscri();

      for(unsigned int l=0; l!=bvec.size(); l++){
	cout << " b LABEL " << bvec[l].first << "  " << bvec[l].second << endl;
      }
      */

      if ( cor_ != jet.corrStep() || flav_ != jet.corrFlavour())  
	correction = jet.corrFactor(cor_, flav_);
      
      // Jet Corrections:
      //  L0 : raw
      //  L1 : off
      //  L2 : rel
      //  L3 : abs  <-- default
      //  L4 : emf
      //  L5 : had  with glu/uds/c/b
      //  L6 : ue   with glu/uds/c/b
      //  L7 : part with glu/uds/c/b
      
      //cout << " Jet Correction [ " << i << " ]  (" << jet.corrStep() << ", " << jet.corrFlavour() 
      //   << " ) -> (" << cor_ << ", " << flav_ << ") : " << correction << endl;
      
      if ((jets[i].pt() * correction) < jetpt_ || fabs(jets[i].eta()) > jeteta_) continue;
    
      mTreeJetP[countjets]      = jets[i].p() * correction;
      mTreeJetPt[countjets]     = jets[i].pt() * correction;
      mTreeJetE[countjets]      = jets[i].energy() * correction;
      mTreeJetEt[countjets]     = jets[i].et() * correction;
      mTreeJetPx[countjets]     = jets[i].momentum().X() * correction;
      mTreeJetPy[countjets]     = jets[i].momentum().Y() * correction;
      mTreeJetPz[countjets]     = jets[i].momentum().Z() * correction;
      mTreeJetEta[countjets]    = jets[i].eta();
      mTreeJetPhi[countjets]    = jets[i].phi();

      if (jets[i].isCaloJet()) {
	mTreeJetFem[countjets]   = jets[i].emEnergyFraction();
	mTreeJetFhad[countjets]  = jets[i].energyFractionHadronic();
	mTreeJetConst[countjets] = jets[i].getCaloConstituents().size();
      }
      if (jets[i].isPFJet()) {
	mTreeJetF[countjets][0]  = jets[i].chargedHadronEnergyFraction();
	mTreeJetN[countjets][0]  = jets[i].chargedMultiplicity();     // pfSpecific().chargedHadronMultiplicity();
	mTreeJetF[countjets][1]  = jets[i].neutralHadronEnergyFraction(); 
	mTreeJetN[countjets][1]  = jets[i].neutralMultiplicity();     // pfSpecific().neutralHadronMultiplicity();
	mTreeJetF[countjets][2]  = jets[i].neutralEmEnergyFraction(); // photonEnergyFraction();
	mTreeJetN[countjets][2]  = 0;                                 // jets[i].pfSpecific().photonMultiplicity();
	mTreeJetF[countjets][3]  = jets[i].chargedEmEnergyFraction(); // electronEnergyFraction();
	mTreeJetN[countjets][3]  = 0;                                 // jets[i].pfSpecific().electronMultiplicity();
	mTreeJetF[countjets][4]  = jets[i].chargedMuEnergyFraction(); // muonEnergyFraction();
	mTreeJetN[countjets][4]  = jets[i].muonMultiplicity();        // pfSpecific().muonMultiplicity();

	mTreeJetConst[countjets] = jets[i].getPFConstituents().size();
      }
      mTreeJetPart[countjets]    = jets[i].partonFlavour();
      // default b tagger
      mTreeJetBtag[countjets]    = jets[i].bDiscriminator(btag_);
      mTreeJetCharge[countjets]  = jets[i].jetCharge();
      mTreeJetn90[countjets]     = jets[i].n90();
      mTreeJetn90hits[countjets] = jets[i].jetID().n90Hits;
      mTreeJetfhpd[countjets]    = jets[i].jetID().fHPD;
      mTreeJetfrbx[countjets]    = jets[i].jetID().fRBX;

      mTreeJetTruth[countjets] = -1;
      const reco::GenJet * gl = jets[i].genJet();
      
      if (gl) {
	for (int k=0; k<mTreeNtruthjet; k++) {
	  if (fabs(mTreetruthJetE[k]   - gl->energy()) < 1e-5 &&
	      fabs(mTreetruthJetEta[k] - gl->eta()) < 1e-5 &&
	      fabs(mTreetruthJetPhi[k] - gl->phi()) < 1e-5) {
	    mTreeJetTruth[countjets] = k;
	    break;
	  }
	}
      }

      countjets++;
    }
    mTreeNjet = countjets;
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

    mTreeMET[0]         = metHandle->front().et();
    mTreeMEX[0]         = metHandle->front().momentum().X();
    mTreeMEY[0]         = metHandle->front().momentum().Y();
    mTreeSumET[0]       = metHandle->front().sumEt();
    mTreeMETeta[0]      = metHandle->front().eta();
    mTreeMETphi[0]      = metHandle->front().phi();
    mTreeSumETSignif[0] = metHandle->front().mEtSig();
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
      
      mTreeMET[1]         = genmetHandle->front().et();
      mTreeMEX[1]         = genmetHandle->front().momentum().X();
      mTreeMEY[1]         = genmetHandle->front().momentum().Y();
      mTreeSumET[1]       = genmetHandle->front().sumEt();
      mTreeMETeta[1]      = genmetHandle->front().eta();
      mTreeMETphi[1]      = genmetHandle->front().phi();
      mTreeSumETSignif[1] = genmetHandle->front().mEtSig();
      
    }
    
    edm::Handle< std::vector<reco::GenMET> > genmet2Handle;
    iEvent.getByLabel("genMetCaloAndNonPrompt", genmet2Handle);
    if ( !genmet2Handle.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenMET results found for InputTag genMetCaloAndNonPrompt";
    if ( genmet2Handle->size()!=1 ) 
      edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					    << genmet2Handle->size() << " instead of 1";
    if ( genmet2Handle.isValid() && genmet2Handle->size()==1 ) {
      
      mTreeMET[2]         = genmet2Handle->front().et();
      mTreeMEX[2]         = genmet2Handle->front().momentum().X();
      mTreeMEY[2]         = genmet2Handle->front().momentum().Y();
      mTreeSumET[2]       = genmet2Handle->front().sumEt();
      mTreeMETeta[2]      = genmet2Handle->front().eta();
      mTreeMETphi[2]      = genmet2Handle->front().phi();
      mTreeSumETSignif[2] = genmet2Handle->front().mEtSig();
      
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
    
    mTreeMET[3]         = metHandle->front().et();
    mTreeMEX[3]         = metHandle->front().momentum().X();
    mTreeMEY[3]         = metHandle->front().momentum().Y();
    mTreeSumET[3]       = metHandle->front().sumEt();
    mTreeMETeta[3]      = metHandle->front().eta();
    mTreeMETphi[3]      = metHandle->front().phi();
    mTreeSumETSignif[3] = metHandle->front().mEtSig();
  }

  // TrackCorrected MET
  iEvent.getByLabel(metTagTC_, metHandle);
  if ( !metHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No Met results found for InputTag " << metTagTC_;
  if ( metHandle->size()!=1 ) 
    edm::LogWarning("SusyACSkimAnalysis") << "MET collection size is "
					  << metHandle->size() << " instead of 1";
  if ( metHandle.isValid() && metHandle->size()==1 ) {
    
    mTreeMET[4]         = metHandle->front().et();
    mTreeMEX[4]         = metHandle->front().momentum().X();
    mTreeMEY[4]         = metHandle->front().momentum().Y();
    mTreeSumET[4]       = metHandle->front().sumEt();
    mTreeMETeta[4]      = metHandle->front().eta();
    mTreeMETphi[4]      = metHandle->front().phi();
    mTreeSumETSignif[4] = metHandle->front().mEtSig();
  }
  

  // This filter
  if (nele_>0 && mTreeNele<nele_) return 0;
  if (nmuo_>0 && mTreeNmuo<nmuo_) return 0;
  if (njet_>0 && mTreeNjet<njet_) return 0;
  if (mTreeMET[0] < met_) return 0;

  nrEventPassedRaw_++;

  // ECal Noise
  edm::Handle<EcalRecHitCollection> pEBRecHits;
  iEvent.getByLabel( "ecalRecHit", "EcalRecHitsEB", pEBRecHits );
  
  if ( pEBRecHits.isValid() ) {
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

      //      cout << iEvent.time().value() << " -> " << energy << "  " << energy/e3x3 << "  " << maxhit.time() << "  " 
      //	   << maxhit.outOfTimeEnergy() << "  " << maxhit.chi2Prob() << endl;
      //      cout << " -> " << energy << "  " << e3x3 << "  " << pos.eta() << "  " << pos.phi() 
      //   << "  " << energy*pf << "  " << newEta  << "  " << mTreeecalpx << "  " 
      //   << mTreeecalpy << "  " << mTreeecalpz << endl;
    }
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
       (tid>1000020 && tid<1000040 && tid!=1000022) ||
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
  mAllData->Branch("global_store",   &mTreestore,       "global_store/I");
  mAllData->Branch("global_run",     &mTreerun,         "global_run/I");
  mAllData->Branch("global_event",   &mTreeevent,       "global_event/I");
  mAllData->Branch("global_bx",      &mTreebx,          "global_bx/I");
  mAllData->Branch("global_orbit",   &mTreeorbit,       "global_orbit/I");
  mAllData->Branch("global_exp",     &mTreeexp,         "global_exp/I");
  mAllData->Branch("global_isdata",  &mTreedata,        "global_isdata/I");
  mAllData->Branch("global_HLT",     &mTreeHLT,         "global_HLT/C");

  mAllData->Branch("lumi_section", &mTreelumiblk,    "lumi_section/I");
  mAllData->Branch("lumi_del",     &mTreedellumi,    "lumi_del/double");
  mAllData->Branch("lumi_rec",     &mTreereclumi,    "lumi_rec/double");
  mAllData->Branch("lumi_delerr",  &mTreedellumierr, "lumi_delerr/double");
  mAllData->Branch("lumi_recerr",  &mTreereclumierr, "lumi_recerr/double");

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

  mAllData->Branch("trig_HLTName",   &mTreetrighltname, "trig_HLTName[50]/I");
  mAllData->Branch("trig_n",         &mTreeNtrig,       "trig_n/I");
  mAllData->Branch("trig_prescale",   mTreetrigpre,     "trig_prescale[trig_n]/I");
  mAllData->Branch("trig_name",       mTreetrigname,    "trig_name[trig_n][100]/I");
  mAllData->Branch("trig_filter",     mTreefiltname,    "trig_filter[trig_n][100]/I");
  mAllData->Branch("trig_pt",         mTreetrigpt,      "trig_pt[trig_n]/double");
  mAllData->Branch("trig_eta",        mTreetrigeta,     "trig_eta[trig_n]/double");
  mAllData->Branch("trig_phi",        mTreetrigphi,     "trig_phi[trig_n]/double");

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
  mAllData->Branch("met_et",       &mTreeMET,         "met_et[5]/double");
  mAllData->Branch("met_ex",       &mTreeMEX,         "met_ex[5]/double");
  mAllData->Branch("met_ey",       &mTreeMEY,         "met_ey[5]/double");
  mAllData->Branch("met_eta",      &mTreeMETeta,      "met_eta[5]/double");
  mAllData->Branch("met_phi",      &mTreeMETphi,      "met_phi[5]/double");
  mAllData->Branch("met_sumet",    &mTreeSumET,       "met_sumet[5]/double");
  mAllData->Branch("met_sumetsig", &mTreeSumETSignif, "met_sumetsig[5]/double");

  // Jets
  mAllData->Branch("jet_n",     &mTreeNjet,      "jet_n/I");  
  mAllData->Branch("jet_E" ,     mTreeJetE,      "jet_E[jet_n]/double");
  mAllData->Branch("jet_Et",     mTreeJetEt,     "jet_Et[jet_n]/double");
  mAllData->Branch("jet_p",      mTreeJetP,      "jet_p[jet_n]/double");
  mAllData->Branch("jet_pt",     mTreeJetPt,     "jet_pt[jet_n]/double");
  mAllData->Branch("jet_px",     mTreeJetPx,     "jet_px[jet_n]/double");
  mAllData->Branch("jet_py",     mTreeJetPy,     "jet_py[jet_n]/double");
  mAllData->Branch("jet_pz",     mTreeJetPz,     "jet_pz[jet_n]/double");
  mAllData->Branch("jet_eta",    mTreeJetEta,    "jet_eta[jet_n]/double");
  mAllData->Branch("jet_phi",    mTreeJetPhi,    "jet_phi[jet_n]/double");
  mAllData->Branch("jet_fem",    mTreeJetFem,    "jet_fem[jet_n]/double");
  mAllData->Branch("jet_fhad",   mTreeJetFhad,   "jet_fhad[jet_n]/double");
  mAllData->Branch("jet_btag",   mTreeJetBtag,   "jet_btag[jet_n]/double");
  mAllData->Branch("jet_charge", mTreeJetCharge, "jet_charge[jet_n]/double");
  mAllData->Branch("jet_fHPD",   mTreeJetfhpd,   "jet_fHPD[jet_n]/double");
  mAllData->Branch("jet_fRBX",   mTreeJetfrbx,   "jet_fRBX[jet_n]/double");
  mAllData->Branch("jet_n90hits",mTreeJetn90hits,"jet_n90hits[jet_n]/I");
  mAllData->Branch("jet_n90",    mTreeJetn90,    "jet_n90[jet_n]/I");
  mAllData->Branch("jet_flav",   mTreeJetPart,   "jet_flav[jet_n]/I");
  mAllData->Branch("jet_truth",  mTreeJetTruth,  "jet_truth[jet_n]/I");
  mAllData->Branch("jet_const",  mTreeJetConst,  "jet_const[jet_n]/I");
  mAllData->Branch("jet_PFN",    mTreeJetN,      "jet_PFN[jet_n][5]/I");
  mAllData->Branch("jet_PFF",    mTreeJetF,      "jet_PFF[jet_n][5]/double");

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
  mAllData->Branch("ele_convdist", mTreeEleConvdist, "ele_convdist[ele_n]/double");
  mAllData->Branch("ele_convdcot", mTreeEleConvdcot, "ele_convdcot[ele_n]/double");
  mAllData->Branch("ele_convr",    mTreeEleConvr,    "ele_convr[ele_n]/double");
  mAllData->Branch("ele_fbrem",    mTreeElefbrem,    "ele_fbrem[ele_n]/double");
  mAllData->Branch("ele_trign",    mTreeNeletrign,   "ele_trign[ele_n]/I");
  mAllData->Branch("ele_trig" ,    mTreeEletrig,     "ele_trig[ele_n][100]/I");  
  mAllData->Branch("ele_SC"   ,    mTreeEleSC,       "ele_SC[ele_n]/I");

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
  mAllData->Branch("muo_trig" ,         mTreeMuotrig,         "muo_trig[muo_n][100]/I");  
  mAllData->Branch("muo_ID",            mTreeMuoID,           "muo_ID[muo_n][24]/I");

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
