from PhysicsTools.PatAlgos.tools.coreTools import *

process = cms.Process("ANA")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SusyACSkimAnalysis'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('GR_R_42_V13::All')
#process.GlobalTag.globaltag = cms.string('START42_V12::All')
process.GlobalTag.globaltag = cms.string('')
process.load("Configuration.StandardSequences.MagneticField_cff")

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("RecoVertex.Configuration.RecoVertex_cff")
from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")

#-- To get JEC in 4_2 return rho corrections:----------------------------------------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True

#--To modify noPU needed for METnoPU -------------------------------------------
process.load('CommonTools.ParticleFlow.pfNoPileUp_cff')

# Output (PAT file only created if
# process.outpath = cms.EndPath(process.out)
# is called at the end
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out2 = cms.OutputModule(
    "PoolOutputModule",
    fileName       = cms.untracked.string('dummy.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData   = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring('keep *')
    )


##dummy out to be modifyed by pf2pat
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName       = cms.untracked.string('dummy2.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData   = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring('keep *')
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out.root')
                                   )

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(7.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
)

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.2)
                                     )



# add iso deposits
from PhysicsTools.PatAlgos.tools.muonTools import addMuonUserIsolation
addMuonUserIsolation(process)


process.patMuons.embedCombinedMuon = False;
process.patMuons.embedStandAloneMuon = False;

### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([
    'file:/home/home1/institut_3a/padeken/public/QCD_Pt-20_MuEnrichedPt-10_TuneZ2_7TeV-pythia6.root']
    #'file:/home/home1/institut_3a/padeken/public/datenRecoMay10ReReco.root']
    #'file:/home/home1/institut_3a/jschulte/public/MetSigTest.root']
    #'/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/148/002/8A37418A-DFD9-DF11-91E3-0030487CD77E.root']
    #'file:/net/data_cms/institut_3a/gueth/punch_through.root']
    #'file:/home/home1/institut_3a/magass/SUSY/CMSSW_3_5_6/src/aachen3a/ACSusyAnalysis/test/CRAB/DATA_1.root',
    #'file:/home/home1/institut_3a/magass/SUSY/CMSSW_3_5_6/src/aachen3a/ACSusyAnalysis/test/CRAB/DATA_2.root']
    ),
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)

# remove MC matching
#removeMCMatching(process, ['All'])

# for PFnoPU
process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpPFlow.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet = False

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
    rParam = cms.double(0.6),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    voronoiRfact = cms.double(0.9)
    )
process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow", "rho")
process.patJetCorrFactorsPFlow.levels = cms.vstring(
     'L2Relative', 'L3Absolute'  # MC
    #'L1FastJet','L2Relative','L3Absolute' #DATA
    )

# Add anti-kt 5 jets
# this is only for calo met is there an other way?
addJetCollection(process,cms.InputTag('ak5CaloJets'), 'AK5', 'Calo',
                 doJTA                  = True,
                 doBTagging        = True,
                 jetCorrLabel     = ('AK5Calo', cms.vstring(['L2Relative', 'L3Absolute'])),  # MC
                 #jetCorrLabel       = ('AK5Calo', cms.vstring(['L1FastJet','L2Relative','L3Absolute'])),  # DATA
                 doType1MET      = True,
                 doL1Cleaning    = False,
                 doL1Counters    = False,
                 genJetCollection =cms.InputTag("ak5GenJets"),
                 doJetID              = True,
                 jetIdLabel          = "ak5"
                 )

# dummy - IMPORTANT
process.patJetCorrFactors.levels = cms.vstring(
     'L2Relative', 'L3Absolute'  # MC
    #'L1FastJet','L2Relative','L3Absolute' #DATA
    )
    
# Add the PV selector and KT6 producer to the sequence
getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsPFlow )

process.patseq = cms.Sequence(
    process.goodOfflinePrimaryVertices*
    getattr(process,"patPF2PATSequence"+postfix)
    )

# add TrackCorrected  met
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')

# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *

# Add met with NoPU to the Collections
from aachen3a.ACSusyAnalysis.pfMET_cfi import *
process.pfMetPFnoPU         = pfMET.clone()
process.pfMetPFnoPU.alias   = 'pfMetNoPileUp'
process.pfMetPFnoPU.src     = 'pfNoPileUp'




process.patJets.addTagInfos = cms.bool(False)  # AOD only

################################
###                          ###
###  Analysis configuration  ###
###                          ###
################################


### Definition of all tags here
elecTag         = cms.InputTag("patElectrons")
pfelecTag      = cms.InputTag("patElectronsPFlow")
calojetTag     = cms.InputTag("patJetsAK5Calo")
pfjetTag        = cms.InputTag("patJetsPFlow")
muonTag      = cms.InputTag("patMuons")
PFmuonTag  = cms.InputTag("selectedPatMuonsPFlow")
tauTag         = cms.InputTag("patTausPFlow")
metTag        = cms.InputTag("patMETsAK5Calo")
metTagTC    = cms.InputTag("patMETsTC")
metTagPF    = cms.InputTag("patMETsPFlow")
metTagPFnoPU=cms.InputTag("pfMetPFnoPU")
metTagJESCorAK5CaloJetMuons =cms.InputTag("metJESCorAK5CaloJetMuons")
metTagcorMetGlobalMuons     = cms.InputTag("corMetGlobalMuons")
metTagHO    = cms.InputTag("metHO")
metTagNoHF= cms.InputTag("metNoHF")
genTag        = cms.InputTag("genParticles")
genJetTag    = cms.InputTag("ak5GenJets")
vtxTag         = cms.InputTag("offlinePrimaryVertices")
#ebhitsTag  = cms.InputTag("ecalRecHit", "EcalRecHitsEB");  # RECO
ebhitsTag     = cms.InputTag("reducedEcalRecHitsEB");   # AOD

### Cuts and switches ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",

    is_MC      = cms.bool(True),  # set to 'False' for real Data !
    is_SHERPA  = cms.bool(False),  # set to 'True' if running on SHERPA
    do_fatjets = cms.bool(False),  # set to 'True' for fat jets
    matchAll   = cms.bool(False),  # if True all truth leptons are matched else only ele and mu
    susyPar    = cms.bool(False),

    # IMPORTANT for QCD -> configured via ./prepare.sh
    pthat_low  = cms.double(-1.),
    pthat_high = cms.double(-1.),

    calojetTag = calojetTag,
    pfjetTag   = pfjetTag,
    elecTag    = elecTag,
    pfelecTag  = pfelecTag,
    muonTag    = muonTag,
    PFmuonTag  = PFmuonTag,
    tauTag     = tauTag,
    metTag     = metTag,
    metTagTC   = metTagTC,
    metTagPF   = metTagPF,
    metTagPFnoPU   =metTagPFnoPU,
    metTagJESCorAK5CaloJetMuons= metTagJESCorAK5CaloJetMuons,
    metTagcorMetGlobalMuons =metTagcorMetGlobalMuons,
    metTagHO = metTagHO,
    metTagNoHF = metTagNoHF,
    genTag     = genTag,
    genJetTag  = genJetTag,
    vtxTag     = vtxTag,
    ebhitsTag  = ebhitsTag,

    muopt      = cms.double(0.),
    muoeta     = cms.double(25.),
    elept      = cms.double(0.),
    eleeta     = cms.double(25.),
    pfelept    = cms.double(0.),
    pfeleeta   = cms.double(25.),
    calojetpt  = cms.double(0.),
    calojeteta = cms.double(25.),
    pfjetpt    = cms.double(0.),
    pfjeteta   = cms.double(25.),
    metcalo    = cms.double(0.),
    metpf      = cms.double(0.),
    mettc      = cms.double(0.),
    nele       = cms.int32(0),
    npfele     = cms.int32(0),
    nmuo       = cms.int32(0),
    ncalojet   = cms.int32(0),
    npfjet     = cms.int32(0),
    htc        = cms.double(0), 
    PFhtc        = cms.double(0), 

    btag       = cms.string('trackCountingHighEffBJetTags'),

    # Calo Jet ID
    jetselvers = cms.string("PURE09"),
    jetselqual = cms.string("LOOSE")

)

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )
process.patTrigger.processName = "*"
process.patTriggerEvent.processName = "*"

### Define the paths
process.p = cms.Path(
    process.kt6PFJets * 
    process.ak5PFJets *
    process.goodOfflinePrimaryVertices*
    process.patDefaultSequence*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.offlinePrimaryVertices *
    process.pfMetPFnoPU*
    process.ACSkimAnalysis
    )

#process.outpath = cms.EndPath(process.out2)
