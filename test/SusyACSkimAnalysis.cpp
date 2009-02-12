//
// Package:    SusyACSkimAnalysis
// Class:      SusyACSkimAnalysis
// 
// Description: Skeleton analysis for SUSY search with Lepton + Jets + MET
//
// Original Author:  Carsten Magass
//         Created:  November 2008
//
#include "SusyAnalysis/ACSusyAnalysis/test/SusyACSkimAnalysis.h"

////////////////////////////////
//
// Constructor
//
SusyACSkimAnalysis::SusyACSkimAnalysis(const edm::ParameterSet& iConfig):
  sequence_( iConfig.getParameter<edm::ParameterSet>("selections") ),
  plotSelection_( iConfig.getParameter<std::vector<std::string> >("plotSelection") ),
  nrEventTotalRaw_(0),
  nrEventPassedRaw_(0),
  nrEventPassedACRaw_(0),
  nrEventTotalWeighted_(0.0)
{
  // Translate plotSelection strings to indices
  plotSelectionIndices_.reserve(plotSelection_.size());
  for ( size_t i=0; i<plotSelection_.size(); ++i )
    plotSelectionIndices_.push_back(sequence_.selectorIndex(plotSelection_[i]));
    
  // List all selectors and selection variables
  edm::LogVerbatim("SusyACSkimAnalysis") << "Selectors are:" << std::endl;
  for ( std::vector<const SusyEventSelector*>::const_iterator it = sequence_.selectors().begin();
        it != sequence_.selectors().end(); ++it ) {

    edm::LogVerbatim("SusyACSkimAnalysis") << " * " << (*it)->name()
				       << " selects on following " 
				       << (*it)->numberOfVariables() << " variable(s):";
    for ( unsigned int i=0; i<(*it)->numberOfVariables(); ++i )
      edm::LogVerbatim("SusyACSkimAnalysis") << "    - " << (*it)->variableNames()[i];
    edm::LogVerbatim("SusyACSkimAnalysis") << std::endl;
  }

  mSelectorResults = new unsigned int[sequence_.size()];

  // get the data tags
  jetTag_    = iConfig.getParameter<edm::InputTag>("jetTag");
  metTag_    = iConfig.getParameter<edm::InputTag>("metTag");
  elecTag_   = iConfig.getParameter<edm::InputTag>("elecTag");
  muonTag_   = iConfig.getParameter<edm::InputTag>("muonTag");
  genJetTag_ = iConfig.getParameter<edm::InputTag>("genJetTag");
  trigTag_   = iConfig.getParameter<edm::InputTag>("trigTag");
  vertexTag_ = iConfig.getParameter<edm::InputTag>("vtxTag");

  gen  = iConfig.getParameter<std::string>("generator");

  string cor  = iConfig.getParameter<std::string>("correction");
  string flav = iConfig.getParameter<std::string>("flavour");

  correction_ = pat::JetCorrFactors::corrStep(cor, flav);
 
  edm::LogVerbatim("SusyACSkimAnalysis") << " Using the following Jet correction: " << cor << " , flavour: " << flav << endl;

  // get the cuts
  muopt_  = iConfig.getParameter<double>("muopt");
  muoeta_ = iConfig.getParameter<double>("muoeta");
  elept_  = iConfig.getParameter<double>("elept");
  eleeta_ = iConfig.getParameter<double>("eleeta");
  jetpt_  = iConfig.getParameter<double>("jetpt");
  jeteta_ = iConfig.getParameter<double>("jeteta");
  jetfem_ = iConfig.getParameter<double>("jetfem");
  met_    = iConfig.getParameter<double>("met");

  nele_ = iConfig.getParameter<int>("nele");
  nmuo_ = iConfig.getParameter<int>("nmuo");
  njet_ = iConfig.getParameter<int>("njet");

  // Initialise counters
  nrEventSelected_.resize( sequence_.size(), 0.0 );
  nrEventSelectedRaw_.resize( sequence_.size(), 0.0 );
  nrEventAllButOne_.resize( sequence_.size(), 0.0 );
  nrEventAllButOneRaw_.resize( sequence_.size(), 0.0 );
  nrEventCumulative_.resize( sequence_.size(), 0.0 );
  nrEventCumulativeRaw_.resize( sequence_.size(), 0.0 );

  localPi = acos(-1.0);

  // Initialise plots [should improve in the future]
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

  // Get some event information (process ID, weight)

  mTreeProcID       = 0;
  mTreeEventWeight  = 1.;
  mTreePthat        = -999.;

  Handle<int> myProcess;
  iEvent.getByLabel("genEventProcID",myProcess);
  
  if (myProcess.isValid())
    mTreeProcID = (*myProcess);

  Handle<double> genEventScale;
  iEvent.getByLabel("genEventScale", genEventScale);

  if (genEventScale.isValid())
    mTreePthat = (*genEventScale);

  Handle<double> genEventWeight;
  iEvent.getByLabel("genEventWeight", genEventWeight);

  if (genEventWeight.isValid())
    mTreeEventWeight = (*genEventWeight);

  /*
  edm::LogInfo("SusyACSkimAnalysis") << " Event properties: EvtWeight = " << mTreeEventWeight 
				     << "  ProcID = " << mTreeProcID 
				     << "  Pthat = " <<  mTreePthat << std::endl;
  */

  // for ntuple
  weight_       = mTreeEventWeight;
  processId_    = mTreeProcID;

  // Retrieve the decision of each selector module
  SelectorDecisions decisions = sequence_.decisions(iEvent);

  // Count all events
  nrEventTotalRaw_++;
  nrEventTotalWeighted_ += weight_;
  
  // Fill plots with all variables
  bool dec(true);
  for ( size_t i=0; i<sequence_.size(); ++i ) {
    dec = dec && decisions.decision(i);
    /*
    edm::LogVerbatim("SusyACSkimAnalysis") 
      << " " << sequence_.selectorName(i)
      << " " << decisions.decision(i)
      << " " << decisions.complementaryDecision(i)
      << " " << decisions.cumulativeDecision(i)
      << " " << dec << std::endl;*/
    
    // Add the decision to the tree
    mSelectorResults[i] = (decisions.decision(i)?1:0);
    
    // Update counters
    if ( decisions.decision(i) ) {
      nrEventSelected_[i] += weight_;
      nrEventSelectedRaw_[i] ++;
    }
    if ( decisions.complementaryDecision(i) ) {
      nrEventAllButOne_[i] += weight_;
      nrEventAllButOneRaw_[i] ++;
    }
    if ( decisions.cumulativeDecision(i) ) {
      nrEventCumulative_[i] += weight_;
      nrEventCumulativeRaw_[i] ++;
    }
  }

  // Fill some plots (only if some selections passed, as configured)
  dec = true;
  for ( size_t i=0; i<plotSelectionIndices_.size(); ++i )
    dec = dec&&decisions.decision(plotSelectionIndices_[i]);
  //  fillPlots( iEvent, decisions );
  if ( dec ) fillPlots( iEvent, decisions );

  // Print summary so far (every 10 till 100, every 100 till 1000, etc.)
  for ( unsigned int i=10; i<nrEventTotalRaw_; i*=10 )
    if ( nrEventTotalRaw_<=10*i && (nrEventTotalRaw_%i)==0 )
      printSummary();

  if ( !sequence_.decisions(iEvent).globalDecision() ) {
    // Just fill the selector data
    mGlobalDecision = 0;
    //    mSelectorData->Fill();
    return mGlobalDecision;
  }

  // Fill all event information
  mGlobalDecision = 1;
  //  mSelectorData->Fill();

  nrEventPassedRaw_++;

  mTreerun     = iEvent.id().run();
  mTreeevent   = iEvent.id().event();
  mTreelumiblk = iEvent.luminosityBlock();
  mTreebx      = iEvent.bunchCrossing();
  mTreeorbit   = iEvent.orbitNumber();
  mTreeexp     = iEvent.experimentType();
  mTreedata    = iEvent.isRealData();
  mTreestore   = 0;

  // Global Event information
  // does not work at all !
  //  Handle<EventAuxiliary> evt;
  //  iEvent.getByLabel("EventAuxiliary", evt);

  //  cout << " fooled you " << (*evt).storeNumber() << endl;
    /* 
  if ( !evt.isValid() )
    edm::LogWarning("SusyACSkimAnalysis") << "No edm::EventAuxiliary is found with tag EventAuxiliary ";
  else {
    mTreestore   = (*evt).storeNumber();
    cout << " fooled you " << (*evt).storeNumber() << endl;
    }*/
  
  // Trigger information
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(trigTag_,triggerResults);
  
  string tempname =":";

  if ( !triggerResults.isValid() )
    edm::LogWarning("SusyACSkimAnalysis") << "No trigger::TriggerResults found for InputTag " << trigTag_;
  else {
    
    edm::TriggerNames triggerNames( *triggerResults );

    /*    
    cout << " found triggers : " << (*triggerResults).size() << endl;
    if((*triggerResults).accept()) 
      cout << "accepted" << endl;
    */

    // loop over all paths, get trigger decision
    for (unsigned i = 0; i != (*triggerResults).size(); ++i) {
      std::string name = triggerNames.triggerName(i);
      //      cout << " -> " << i << " " << name << " ";
      if ( (*triggerResults).accept(i) ) {
	//	cout << "1" << endl;
	tempname = tempname + name + ":";
      }
      //      else
      //	cout << endl;
    }
  }
  strcpy(_mTreeHLT, tempname.c_str());

  // PDF information 
  Handle<reco::PdfInfo> pdfi;
  iEvent.getByLabel("genEventPdfInfo", pdfi);

  if ( !pdfi.isValid() )
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::PdfInfo found for Tag genEventPdfInfo ";
  else {
    mTreepdfid1   = (int)(*pdfi).id1;
    mTreepdfid2   = (int)(*pdfi).id2;
    mTreepdfx1    = (*pdfi).x1;
    mTreepdfx2    = (*pdfi).x2;
    mTreepdff1    = (*pdfi).pdf1;
    mTreepdff2    = (*pdfi).pdf2;
    mTreepdfscale = (*pdfi).scalePDF;
  }

  // Vertex
  edm::Handle< reco::VertexCollection > vertexHandle;
  //  Handle< reco::Vertex > vertexHandle;
  iEvent.getByLabel(vertexTag_, vertexHandle);

  if ( !vertexHandle.isValid() ) {
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::Vertex found for InputTag " << vertexTag_;
  }
  else {

    mTreeNvtx = vertexHandle->size();
    if ( mTreeNvtx > 50 ) mTreeNvtx = 50;

    for (int i=0; i<mTreeNvtx; ++ i ) {

      mTreeVtxchi[i] = (*vertexHandle)[i].normalizedChi2();
      mTreeVtxntr[i] = (*vertexHandle)[i].tracksSize();
      mTreeVtxx[i]   = (*vertexHandle)[i].x();
      mTreeVtxy[i]   = (*vertexHandle)[i].y();
      mTreeVtxz[i]   = (*vertexHandle)[i].z();
    }
  }


  // Truth tree
  mTreeNtruth  = 0;
  mTreeNtruthl = 0;

  if (gen == "source") {
    
    Handle<HepMCProduct> mc;
    iEvent.getByLabel("source", mc);
    
    if ( !mc.isValid() ) 
      edm::LogWarning("SusyACSkimAnalysis") << "No HepMCProduct found with InputTag source ";
    else {
      
      HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(mc->GetEvent()));
      
      //  myGenEvent->print();
      //      cout << " ProcessID  : " << myGenEvent->signal_process_id() << endl;
      //      cout << " Eventscale : " << myGenEvent->event_scale() << endl; 
    
      bool primary = true;
      int pvid = 0;
    
      int ntruth = 0;
    
      for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	    p != myGenEvent->particles_end(); ++p ) {
      
	int id = (*p)->pdg_id();
	int st = (*p)->status();
      
	// stable e/mu in final state
	if (st == 1 && (abs(id) == 11 || abs(id) == 13)) {
	  mTreetruthlpdgid[mTreeNtruthl] = (*p)->pdg_id();
	  mTreetruthlE[mTreeNtruthl]     = (*p)->momentum().e();
	  mTreetruthlp[mTreeNtruthl]     = (*p)->momentum().mag();
	  mTreetruthlpt[mTreeNtruthl]    = (*p)->momentum().perp();
	  mTreetruthlpx[mTreeNtruthl]    = (*p)->momentum().px();
	  mTreetruthlpy[mTreeNtruthl]    = (*p)->momentum().py();
	  mTreetruthlpz[mTreeNtruthl]    = (*p)->momentum().pz();
	  mTreetruthleta[mTreeNtruthl]   = (*p)->momentum().eta();
	  mTreetruthlphi[mTreeNtruthl]   = (*p)->momentum().phi();
	  mTreeNtruthl++;
	}

	HepMC::GenVertex* bVertex = (*p)->production_vertex();
	HepMC::GenVertex* eVertex = (*p)->end_vertex();
      
	int bvtxid = -999;
	int evtxid = -999;
      
	if (bVertex) bvtxid = bVertex->barcode();
	if (eVertex) evtxid = eVertex->barcode();
      
	// (1) primary interaction vertex
	if ( primary && st == 3 && isDecaying(id) ) {
	
	  //	  cout << " Found primary interaction vertex ! " << endl;
	
	  if (bVertex) {
	  
	    pvid = bvtxid;
	  
	    //	    cout << " + Now analysing incoming particles at vtx " << bvtxid << endl;
	  
	    for (HepMC::GenVertex::particles_in_const_iterator iter = bVertex->particles_in_const_begin();
		 iter != bVertex->particles_in_const_end(); iter++) {
	    
	      if ((*iter)->status() != 3) continue;
	    
	      mTreetruthpdgid[ntruth] = (*iter)->pdg_id();
	      mTreetruthbvtxid[ntruth] = 0;
	      mTreetruthevtxid[ntruth] = bvtxid;

	      mTreetruthE[ntruth]   = (*iter)->momentum().e();
	      mTreetruthp[ntruth]   = (*iter)->momentum().mag();
	      mTreetruthpt[ntruth]  = (*iter)->momentum().perp();
	      mTreetruthpx[ntruth]  = (*iter)->momentum().px();
	      mTreetruthpy[ntruth]  = (*iter)->momentum().py();
	      mTreetruthpz[ntruth]  = (*iter)->momentum().pz();
	      mTreetrutheta[ntruth] = (*iter)->momentum().eta();
	      mTreetruthphi[ntruth] = (*iter)->momentum().phi();
	      mTreetruthm[ntruth]  = (*iter)->momentum().m();

	      //	      cout << ntruth << "   --> in " << " [" << (*iter)->barcode() << "] " << mTreetruthpdgid[ntruth] << endl;
	    
	      ntruth++;
	    
	    }
	  
	    //	    cout << " + Now analysing outgoing particles at vtx " << bvtxid << endl;
	  
	    for (HepMC::GenVertex::particles_out_const_iterator iter = bVertex->particles_out_const_begin();
		 iter != bVertex->particles_out_const_end(); iter++) {
	    
	      if ((*iter)->status() != 3) continue;
	    
	      HepMC::GenVertex* eVertex2 = (*iter)->end_vertex();
	      int ev = -999;
	      if (eVertex) ev = eVertex2->barcode();
	    
	      if (isStable((*iter)->pdg_id())) ev = -999;
	    
	      mTreetruthpdgid[ntruth] = (*iter)->pdg_id();
	      mTreetruthbvtxid[ntruth] = bvtxid;
	      mTreetruthevtxid[ntruth] = ev;
	    
	      mTreetruthE[ntruth]   = (*iter)->momentum().e();
	      mTreetruthp[ntruth]   = (*iter)->momentum().mag();
	      mTreetruthpt[ntruth]  = (*iter)->momentum().perp();
	      mTreetruthpx[ntruth]  = (*iter)->momentum().px();
	      mTreetruthpy[ntruth]  = (*iter)->momentum().py();
	      mTreetruthpz[ntruth]  = (*iter)->momentum().pz();
	      mTreetrutheta[ntruth] = (*iter)->momentum().eta();
	      mTreetruthphi[ntruth] = (*iter)->momentum().phi();
	      mTreetruthm[ntruth]   = (*iter)->momentum().m();

	      //	      cout << ntruth << "   --> out " << " [" << (*iter)->barcode() << "] " 
	      //		   << mTreetruthpdgid[ntruth] << " --> " << mTreetruthevtxid[ntruth] << endl;
	    
	      ntruth++;
	    }
	  }
	  //	  cout << endl;
	  primary = false;
	}
      
	// (2) Now decay cascade of primary particles
	if ( !primary && st == 3 && isDecaying(id) ) {
	
	  //	  cout << " + Now analysing outgoing particles at vtx " << evtxid << endl;
	
	  for (HepMC::GenVertex::particles_out_const_iterator iter = eVertex->particles_out_const_begin();
	       iter != eVertex->particles_out_const_end(); iter++) {
	  
	    if ((*iter)->status() != 3) continue;
	    HepMC::GenVertex* eVertex2 = (*iter)->end_vertex();
	    int ev = -999;
	    if (eVertex) ev = eVertex2->barcode();
	  
	    if (isStable((*iter)->pdg_id())) ev = -999;
	  
	    mTreetruthpdgid[ntruth] = (*iter)->pdg_id();
	    mTreetruthbvtxid[ntruth] = evtxid;
	    mTreetruthevtxid[ntruth] = ev;
	  
	    mTreetruthE[ntruth]   = (*iter)->momentum().e();
	    mTreetruthp[ntruth]   = (*iter)->momentum().mag();
	    mTreetruthpt[ntruth]  = (*iter)->momentum().perp();
	    mTreetruthpx[ntruth]  = (*iter)->momentum().px();
	    mTreetruthpy[ntruth]  = (*iter)->momentum().py();
	    mTreetruthpz[ntruth]  = (*iter)->momentum().pz();
	    mTreetrutheta[ntruth] = (*iter)->momentum().eta();
	    mTreetruthphi[ntruth] = (*iter)->momentum().phi();
	    mTreetruthm[ntruth]   = (*iter)->momentum().m();
	  
	    //	    cout << ntruth << "   --> out " << " [" << (*iter)->barcode() << "] " 
	    //		 << mTreetruthpdgid[ntruth] << "  --> " << mTreetruthevtxid[ntruth] << endl;
	  
	    ntruth++;
	  }
	}
      }
      mTreeNtruth = ntruth;
    }
  }

  std::vector<const reco::Candidate*> truthl;
  std::vector<const reco::Candidate*> truth;

  if (gen == "genParticles") {
    
    Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);    
    
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
	  if (mTreeNtruthl == 50) break;
	}
      }

      // (1) Reconstruct 2 -> 2 primary process, if t, W, Z, tau or SUSY particles are involved
      for( size_t i = 0; i < genParticles->size(); ++ i ) {
	const reco::Candidate& p = (*genParticles)[ i ];
      
	if (isDecaying(p.pdgId()) && p.status() == 3) {
	
	  if (p.numberOfMothers()==2) {
	  
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

    if ( mTreeNtruthjet > 50 ) mTreeNtruthjet = 50;

    for (int i=0; i<mTreeNtruthjet; ++ i ) {

      mTreetruthJetP[i]   = genjets[i].p();
      mTreetruthJetPt[i]  = genjets[i].pt();
      mTreetruthJetE[i]   = genjets[i].energy();
      mTreetruthJetEt[i]  = genjets[i].et();
      mTreetruthJetPx[i]  = genjets[i].momentum().X();
      mTreetruthJetPy[i]  = genjets[i].momentum().Y();
      mTreetruthJetPz[i]  = genjets[i].momentum().Z();
      mTreetruthJetEta[i] = genjets[i].eta();
      mTreetruthJetPhi[i] = genjets[i].phi();
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

    if ( mTreeNele > 50 ) mTreeNele = 50;

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

      /* DEBUG
      cout << " Electron ID name : ";
      for (std::vector<IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
        cout << " - " << it->first << " [ " << it->second << " ] ";
      }
      cout << endl;
      */

      mTreeEleP[countele]       = eles[i].p();
      mTreeElePt[countele]      = eles[i].pt();
      mTreeEleE[countele]       = eles[i].energy();
      mTreeEleEt[countele]      = eles[i].et();
      mTreeElePx[countele]      = eles[i].momentum().X();
      mTreeElePy[countele]      = eles[i].momentum().Y();
      mTreeElePz[countele]      = eles[i].momentum().Z();
      mTreeEleEta[countele]     = eles[i].eta();
      mTreeElePhi[countele]     = eles[i].phi();
      mTreeEleTrkIso[countele]  = ( eles[i].trackIso()+eles[i].et() )/eles[i].et();
      mTreeEleCharge[countele]  = eles[i].charge();
      
      mTreeEleECalIso[countele] = eles[i].ecalIso();
      mTreeEleHCalIso[countele] = eles[i].hcalIso();
      mTreeEleAllIso[countele]  = eles[i].caloIso();
      
      if (eles[i].gsfTrack().isNonnull()) {
	if (eles[i].gsfTrack().get()->ndof() > 0)
	  mTreeEleTrkChiNorm[countele] = eles[i].gsfTrack().get()->chi2()/ eles[i].gsfTrack().get()->ndof();
	else  mTreeEleTrkChiNorm[countele] = 999.;
      }
      else  mTreeEleTrkChiNorm[countele] = 999.;
      
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

    if ( mTreeNmuo > 50 ) mTreeNmuo = 50;
    
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

      mTreeMuoP[countmuo]       = muons[i].p();
      mTreeMuoPt[countmuo]      = muons[i].pt();
      mTreeMuoE[countmuo]       = muons[i].energy();
      mTreeMuoEt[countmuo]      = muons[i].et();
      mTreeMuoPx[countmuo]      = muons[i].momentum().X();
      mTreeMuoPy[countmuo]      = muons[i].momentum().Y();
      mTreeMuoPz[countmuo]      = muons[i].momentum().Z();
      mTreeMuoEta[countmuo]     = muons[i].eta();
      mTreeMuoPhi[countmuo]     = muons[i].phi();
      mTreeMuoTrkIso[countmuo]  = ( muons[i].trackIso()+muons[i].et() )/muons[i].et();
      mTreeMuoCharge[countmuo]  = muons[i].charge();
      mTreeMuoECalIso[countmuo] = muons[i].ecalIso();
      mTreeMuoHCalIso[countmuo] = muons[i].hcalIso() ;
      mTreeMuoAllIso[countmuo]  = muons[i].caloIso() ;

      // globalTrack() returns combinedMuon()
      if ( muons[i].globalTrack().isNonnull() ) {
	if (muons[i].globalTrack().get()->ndof() > 0) 
	  mTreeMuoTrkChiNorm[countmuo] = muons[i].globalTrack().get() ->chi2()/ muons[i].globalTrack().get()->ndof();
	else  
	  mTreeMuoTrkChiNorm[countmuo] = 999.;
      }
      else  
	mTreeMuoTrkChiNorm[countmuo] = 999.;

      if ( muons[i].combinedMuon().isNonnull() ) 
	mTreeMuoHits[countmuo] = muons[i].combinedMuon().get()->numberOfValidHits();
      else
	mTreeMuoHits[countmuo] = -1;

      countmuo++;
    }
    mTreeNmuo = countmuo;
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

    if ( mTreeNjet >50 ) mTreeNjet = 50;
    
    for (int i=0; i<mTreeNjet; i++) {
      
      const pat::Jet& jet = jets[i];
      float correction = 1.;
      if ( correction_ != jet.jetCorrStep() )
	correction = jet.jetCorrFactors().correction(correction_);

      // Jet Corrections:
      //  L0 : RAW
      //  L1 : OFF
      //  L2 : REL
      //  L3 : ABS  <-- default
      //  L4 : EMF
      //  L5 : HAD  with GLU/UDS/C/B
      //  L6 : UE   with GLU/UDS/C/B
      //  L7 : PART with GLU/UDS/C/B

      //      cout << " Jet Correction [ " << i << " ] " << correction << endl;

      if ((jets[i].pt() * correction) < jetpt_ || fabs(jets[i].eta()) > jeteta_ ||
	  jets[i].emEnergyFraction() > jetfem_) continue;
     
      mTreeJetP[countjets]   = jets[i].p() * correction;
      mTreeJetPt[countjets]  = jets[i].pt() * correction;
      mTreeJetE[countjets]   = jets[i].energy() * correction;
      mTreeJetEt[countjets]  = jets[i].et() * correction;
      mTreeJetPx[countjets]  = jets[i].momentum().X() * correction;
      mTreeJetPy[countjets]  = jets[i].momentum().Y() * correction;
      mTreeJetPz[countjets]  = jets[i].momentum().Z() * correction;
      mTreeJetEta[countjets] = jets[i].eta();
      mTreeJetPhi[countjets] = jets[i].phi();
      mTreeJetFem[countjets] = jets[i].emEnergyFraction();

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

  // genMET 
  edm::Handle< std::vector<reco::GenMET> > genmetHandle;
  iEvent.getByLabel("genMet", genmetHandle);
  if ( !genmetHandle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenMET results found for InputTag genMet";
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
  iEvent.getByLabel("genMetNoNuBSM", genmet2Handle);
  if ( !genmet2Handle.isValid() ) 
    edm::LogWarning("SusyACSkimAnalysis") << "No reco::GenMET results found for InputTag genMetNoNuBSM";
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

  // This filter
  if (nele_>0 && mTreeNele<nele_) return 0;
  if (nmuo_>0 && mTreeNmuo<nmuo_) return 0;
  if (njet_>0 && mTreeNjet<njet_) return 0;
  if (mTreeMET[0] < met_) return 0;

  nrEventPassedACRaw_++;

  // Fill the tree
  mAllData->Fill();

  return mGlobalDecision;
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
       tid==24 || tid==23 || tid==6 || tid==15 )
    return true;
  else
    return false;
}
  
////////////////////////////////
//
// Begin of Job
//
void SusyACSkimAnalysis::beginJob(const edm::EventSetup&) {}

////////////////////////////////
//
// End of Job
//

void SusyACSkimAnalysis::endJob() {

  h_counters->SetBinContent(1, nrEventTotalRaw_);
  h_counters->SetBinContent(2, nrEventPassedRaw_);
  h_counters->SetBinContent(3, nrEventPassedACRaw_);

  printSummary();

}

////////////////////////////////
//
// Summary
//
void SusyACSkimAnalysis::printSummary( void ) {
  
  edm::LogInfo("SusyACSkimAnalysis") << "*** Summary of counters: ";
  edm::LogVerbatim("SusyACSkimAnalysis") << "Total number of events = " << nrEventTotalRaw_
					 << " ; selected = " << nrEventPassedRaw_  
					 << " ; selectedAC = " << nrEventPassedACRaw_ << endl;

  std::ostringstream summary;
  summary << std::setw(21) << std::left << "Name"
          << std::setw(35) << "Selected"
          << std::setw(35) << "AllButOne"
          << std::setw(35) << "Cumulative" << std::endl;
         
  for ( size_t i=0; i<sequence_.size(); ++i ) {
    summary << std::setw(15) << std::left << sequence_.selectorName(i) << std::right;
    summary << std::setw(10) << std::setprecision(2) << std::fixed
            << nrEventSelected_[i] 
	    << " (" << std::setw(10) << std::setprecision(0) << nrEventSelectedRaw_[i] << ")" 
            << " [" << std::setw(6) << std::setprecision(2) << std::fixed
            << (nrEventSelected_[i]/nrEventTotalWeighted_)*100. << "%]  ";
    summary << std::setw(10) << nrEventAllButOne_[i] 
	    << " (" << std::setw(10) << std::setprecision(0) << nrEventAllButOneRaw_[i] << ")" 
            << " [" << std::setw(6) << std::setprecision(2) << std::fixed
            << (nrEventAllButOne_[i]/nrEventTotalWeighted_)*100. << "%]  ";
    summary << std::setw(10) << nrEventCumulative_[i] 
	    << " (" << std::setw(10) << std::setprecision(0) << nrEventCumulativeRaw_[i] << ")"
            << " [" << std::setw(6) << std::setprecision(2) << std::fixed
            << (nrEventCumulative_[i]/nrEventTotalWeighted_)*100. << "%]  ";
    summary << std::endl; 
  }
  edm::LogVerbatim("SusyACSkimAnalysis") << summary.str();

}

////////////////////////////////
//
// Initialize Plots for ntuple
//
void SusyACSkimAnalysis::initPlots() {

  std::ostringstream variables; // Container for all variables

  // 1. Event variables
  //  variables << "evtWeight:procID";

  // 2. Decision from all selectors
  for ( std::vector<const SusyEventSelector*>::const_iterator it = sequence_.selectors().begin();
        it != sequence_.selectors().end(); ++it ) {
    std::string var( (*it)->name() );
    var += "_result";
    // Push to list of variables
    variables << var << ":";
  }
  variables << "all_result"; // Also store global decision
    
  // 3. All variables from sequence
  for ( std::vector<const SusyEventSelector*>::const_iterator it = sequence_.selectors().begin();
        it != sequence_.selectors().end(); ++it ) {
    for ( unsigned int i=0; i<(*it)->numberOfVariables(); ++i ) {
      std::string var( (*it)->name() ); // prefix variable with selector name
      var += "." + (*it)->variableNames()[i];
      // Push to list of variables
      variables << ":" << var;
    }
  }

  // Register this ntuple
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TNtuple>( "ntuple","SusyACSkimAnalysis variables",
                               variables.str().c_str() );

  // Now we add some additional ones for the dijet analysis
  mAllData = fs->make<TTree>( "allData", "data after cuts" );
  //  mSelectorData = fs->make<TTree>( "selectorData" , "Bit results for selectors");

  h_counters = fs->make<TH1F>("h_counters", "Event Counter", 10, 0, 10);

  mAllData->SetAutoSave(10000);
  //  mSelectorData->SetAutoSave(10000);

  /*
  std::vector<std::string> names = sequence_.selectorNames();
  for ( size_t i = 0 ; i < sequence_.size() ; ++i ) {
    std::string tempName = names[i] + "/i";
    mSelectorData->Branch(names[i].c_str(),&mSelectorResults[i],tempName.c_str());
  }
  mSelectorData->Branch("globalDecision",&mGlobalDecision,"globalDecision/i");
  */

  // Add the branches

  // Event
  mAllData->Branch("global_weight",  &mTreeEventWeight, "global_weight/double");
  mAllData->Branch("global_procID",  &mTreeProcID,      "global_procID/I");
  mAllData->Branch("global_pthat",   &mTreePthat,       "global_pthat/double");
  mAllData->Branch("global_store",   &mTreestore,       "global_store/I");
  mAllData->Branch("global_run",     &mTreerun,         "global_run/I");
  mAllData->Branch("global_event",   &mTreeevent,       "global_event/I");
  mAllData->Branch("global_bx",      &mTreebx,          "global_bx/I");
  mAllData->Branch("global_lumiblk", &mTreelumiblk,     "global_lumiblk/I");
  mAllData->Branch("global_orbit",   &mTreeorbit,       "global_orbit/I");
  mAllData->Branch("global_exp",     &mTreeexp,         "global_exp/I");
  mAllData->Branch("global_isdata",  &mTreedata,        "global_isdata/I");
  mAllData->Branch("global_HLT",     &_mTreeHLT,        "global_HLT/C");

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
  mAllData->Branch("vtx_n",  &mTreeNvtx,   "vtx_n/int");
  mAllData->Branch("vtx_ntr", mTreeVtxntr, "vtx_ntr[vtx_n]/int");
  mAllData->Branch("vtx_x",   mTreeVtxx,   "vtx_x[vtx_n]/double");
  mAllData->Branch("vtx_y",   mTreeVtxy,   "vtx_y[vtx_n]/double");
  mAllData->Branch("vtx_z",   mTreeVtxz,   "vtx_z[vtx_n]/double");
  mAllData->Branch("vtx_chi", mTreeVtxchi, "vtx_chi[vtx_n]/double");

  // MET
  mAllData->Branch("met",      &mTreeMET,         "met[3]/double");
  mAllData->Branch("mex",      &mTreeMEY,         "mex[3]/double");
  mAllData->Branch("mey",      &mTreeMEX,         "mey[3]/double");
  mAllData->Branch("meteta",   &mTreeMETeta,      "meteta[3]/double");
  mAllData->Branch("metphi",   &mTreeMETphi,      "metphi[3]/double");
  mAllData->Branch("sumet",    &mTreeSumET,       "sumet[3]/double");
  mAllData->Branch("sumetsig", &mTreeSumETSignif, "sumetsig[3]/double");

  // Jets
  mAllData->Branch("jet_n",    &mTreeNjet,     "jet_n/int");  
  mAllData->Branch("jet_E" ,    mTreeJetE,     "jet_E[jet_n]/double");
  mAllData->Branch("jet_Et",    mTreeJetEt,    "jet_Et[jet_n]/double");
  mAllData->Branch("jet_p",     mTreeJetP,     "jet_p[jet_n]/double");
  mAllData->Branch("jet_pt",    mTreeJetPt,    "jet_pt[jet_n]/double");
  mAllData->Branch("jet_px",    mTreeJetPx,    "jet_px[jet_n]/double");
  mAllData->Branch("jet_py",    mTreeJetPy,    "jet_py[jet_n]/double");
  mAllData->Branch("jet_pz",    mTreeJetPz,    "jet_pz[jet_n]/double");
  mAllData->Branch("jet_eta",   mTreeJetEta,   "jet_eta[jet_n]/double");
  mAllData->Branch("jet_phi",   mTreeJetPhi,   "jet_phi[jet_n]/double");
  mAllData->Branch("jet_fem",   mTreeJetFem,   "jet_fem[jet_n]/double");
  mAllData->Branch("jet_truth", mTreeJetTruth, "jet_truth[jet_n]/I");

  // Generator Jets
  mAllData->Branch("truthjet_n",   &mTreeNtruthjet,   "truthjet_n/int");  
  mAllData->Branch("truthjet_E" ,   mTreetruthJetE,   "truthjet_E[truthjet_n]/double");
  mAllData->Branch("truthjet_Et",   mTreetruthJetEt,  "truthjet_Et[truthjet_n]/double");
  mAllData->Branch("truthjet_p",    mTreetruthJetP,   "truthjet_p[truthjet_n]/double");
  mAllData->Branch("truthjet_pt",   mTreetruthJetPt,  "truthjet_pt[truthjet_n]/double");
  mAllData->Branch("truthjet_px",   mTreetruthJetPx,  "truthjet_px[truthjet_n]/double");
  mAllData->Branch("truthjet_py",   mTreetruthJetPy,  "truthjet_py[truthjet_n]/double");
  mAllData->Branch("truthjet_pz",   mTreetruthJetPz,  "truthjet_pz[truthjet_n]/double");
  mAllData->Branch("truthjet_eta",  mTreetruthJetEta, "truthjet_eta[truthjet_n]/double");
  mAllData->Branch("truthjet_phi",  mTreetruthJetPhi, "truthjet_phi[truthjet_n]/double");

  // Electrons
  mAllData->Branch("ele_n",         &mTreeNele,          "ele_n/int");  
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
  mAllData->Branch("ele_TrkIso",     mTreeEleTrkIso,     "ele_TrkIso[ele_n]/double");
  mAllData->Branch("ele_ECalIso",    mTreeEleECalIso,    "ele_ECalIso[ele_n]/double");
  mAllData->Branch("ele_HCalIso",    mTreeEleHCalIso,    "ele_HCalIso[ele_n]/double");
  mAllData->Branch("ele_AllIso",     mTreeEleAllIso,     "ele_AllIso[ele_n]/double");
  mAllData->Branch("ele_TrkChiNorm", mTreeEleTrkChiNorm, "ele_TrkChiNorm[ele_n]/double");
  mAllData->Branch("ele_truth",      mTreeEleTruth,      "ele_truth[ele_n]/I");
  mAllData->Branch("ele_ID",         mTreeEleID,         "ele_ID[ele_n][5]/int");

  // Muons
  mAllData->Branch("muo_n" ,        &mTreeNmuo,          "muo_n/int");  
  mAllData->Branch("muo_E" ,         mTreeMuoE,          "muo_E[muo_n]/double");
  mAllData->Branch("muo_Et",         mTreeMuoEt,         "muo_Et[muo_n]/double");
  mAllData->Branch("muo_p",          mTreeMuoP,          "muo_p[muo_n]/double");
  mAllData->Branch("muo_pt",         mTreeMuoPt,         "muo_pt[muo_n]/double");
  mAllData->Branch("muo_px",         mTreeMuoPx,         "muo_px[muo_n]/double");
  mAllData->Branch("muo_py",         mTreeMuoPy,         "muo_py[muo_n]/double");
  mAllData->Branch("muo_pz",         mTreeMuoPz,         "muo_pz[muo_n]/double");
  mAllData->Branch("muo_eta",        mTreeMuoEta,        "muo_eta[muo_n]/double");
  mAllData->Branch("muo_phi",        mTreeMuoPhi,        "muo_phi[muo_n]/double");
  mAllData->Branch("muo_charge",     mTreeMuoCharge,     "muo_charge[muo_n]/double");
  mAllData->Branch("muo_TrkIso",     mTreeMuoTrkIso,     "muo_TrkIso[muo_n]/double");
  mAllData->Branch("muo_ECalIso",    mTreeMuoECalIso,    "muo_ECalIso[muo_n]/double");
  mAllData->Branch("muo_HCalIso",    mTreeMuoHCalIso,    "muo_HCalIso[muo_n]/double");
  mAllData->Branch("muo_AllIso",     mTreeMuoAllIso,     "muo_AllIso[muo_n]/double");
  mAllData->Branch("muo_TrkChiNorm", mTreeMuoTrkChiNorm, "muo_TrkChiNorm[muo_n]/double");
  mAllData->Branch("muo_hits",       mTreeMuoHits,       "muo_hits[muo_n]/I");
  mAllData->Branch("muo_truth",      mTreeMuoTruth,      "muo_truth[muo_n]/I");

  //  h_counters = new TH1F("h_counters", "Counters", 10, 0, 1);

  edm::LogInfo("SusyACSkimAnalysis") << "Ntuple variables " << variables.str();

}

////////////////////////////////
//
// Fill Plots for ntuple
//
void SusyACSkimAnalysis::fillPlots( const edm::Event& iEvent, const SelectorDecisions& decisions ) {
  
  // Container array
  float* x = new float[ntuple_->GetNbranches()];
  int ivar = 0; 

  // 1. Event variables
  //  x[ivar++] = weight_;
  //  x[ivar++] = processId_;
  
  // 2. Decision from all selectors
  for ( size_t i=0; i<sequence_.size(); ++i ) x[ivar++] = decisions.decision(i);
  x[ivar++] = decisions.globalDecision();

  // 3. All variables from sequence
  std::vector<double> values = sequence_.values();
  for ( size_t i=0; i<values.size(); ++i ) x[ivar++] = values[i];

  if ( ntuple_->Fill( x ) < 0 ) { // Fill returns number of bytes committed, -1 on error
    edm::LogWarning("SusyACSkimAnalysis") << "@SUB=fillPlots()" << "Problem filling ntuple";
  }

  delete [] x; // Important! otherwise we'll leak...

}


// Define this as a plug-in
DEFINE_FWK_MODULE(SusyACSkimAnalysis);
