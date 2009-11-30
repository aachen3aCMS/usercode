import FWCore.ParameterSet.Config as cms

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
process.GlobalTag.globaltag = cms.string('STARTUP31X_V1::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#-- JES -----------------------------------------------------------------------
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_cff")

#-- Remove all Monte Carlo matching -------------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, 'All')


### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/data/BeamCommissioning09/MinimumBias/RECO/rereco_FIRSTCOLL_v1/0083/EA91CC4D-28D9-DE11-833F-00261894396E.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# this is just to satisfy the parser - NO PAT file is created !
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName       = cms.untracked.string('dummy.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData   = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring('drop *')
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out.root')
                                   )

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,
                           cms.InputTag('ak5CaloJets'),
                           doJTA            = True,           # Run Jet-Track association & JetCharge
                           doBTagging       = True,           # Run b-tagging
                           jetCorrLabel     = ('AK5','Calo'), # jet correction name; set to None for no JEC
                           doType1MET       = True,           # recompute Type1MET using these jets
                           genJetCollection = cms.InputTag("ak5GenJets"))

# Add latest HcalNoiseSummary
process.load("RecoMET.METProducers.hcalnoiseinfoproducer_cfi")
process.hcalnoise.refillRefVectors = True
process.hcalnoise.hcalNoiseRBXCollName = "hcalnoise" # This has changed in 33X

### Definition of all tags here
elecTag   = cms.InputTag("cleanLayer1Electrons")
jetTag    = cms.InputTag("cleanLayer1Jets")
muonTag   = cms.InputTag("cleanLayer1Muons")
metTag    = cms.InputTag("layer1METs")
genTag    = cms.InputTag("genParticles")
genJetTag = cms.InputTag("antikt5GenJets")
vtxTag    = cms.InputTag("offlinePrimaryVertices")

### Analysis configuration ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",
    
    is_MC     = cms.bool(False),   # set to 'False' for real Data !
    is_SHERPA = cms.bool(False),   # set to 'True' if running on SHERPA
    do_fatjets = cms.bool(False),  # set to 'True' for fat jets
                                   # if 'True', include process.BoostedHiggsSubjets (see example)

    # IMPORTANT for QCD -> configured via ./prepare.sh 
    pthat_low  = cms.double(-1.),
    pthat_high = cms.double(-1.),

    jetTag    = jetTag,
    elecTag   = elecTag,
    muonTag   = muonTag,
    metTag    = metTag,
    genTag    = genTag,
    genJetTag = genJetTag,
    vtxTag    = vtxTag,

    muopt  = cms.double(0.),
    muoeta = cms.double(2.5),
    elept  = cms.double(0.),
    eleeta = cms.double(2.5),
    jetpt  = cms.double(0.),
    jeteta = cms.double(2.5),
    jetfem = cms.double(10.),
    met    = cms.double(0.),
    nele   = cms.int32(0),
    nmuo   = cms.int32(0),
    njet   = cms.int32(0),
    
    correction = cms.string('abs'),
    flavour    = cms.string('glu'),

    btag       = cms.string('trackCountingHighEffBJetTags'),
    
)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerSequence.remove( process.patTriggerMatcher )
process.patTriggerEvent.patTriggerMatches  = ()


### Define the paths
process.p = cms.Path(
    process.hcalnoise*
    process.patDefaultSequence*
    process.patTrigger*
    process.patTriggerEvent*
    process.ACSkimAnalysis
    )
