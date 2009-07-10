import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SusyACSkimAnalysis'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet(process,newName='Winter09',oldName='Summer08Redigi') # change from old to Winter08

### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out.root')
                                   )

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,
                           cms.InputTag('sisCone5CaloJets'),
                           doJTA=True,            # Run Jet-Track association & JetCharge
                           doBTagging=True,       # Run b-tagging
                           jetCorrLabel=('SC5','Calo'), # example jet correction name; set to None for no JEC
                           doType1MET=True,       # recompute Type1MET using these jets
                           genJetCollection=cms.InputTag("sisCone5GenJets"))

### Definition of all tags here
elecTag   = cms.InputTag("selectedLayer1Electrons")
jetTag    = cms.InputTag("selectedLayer1Jets")
muonTag   = cms.InputTag("selectedLayer1Muons")
metTag    = cms.InputTag("layer1METs")                         
genTag    = cms.InputTag("genParticles")
genJetTag = cms.InputTag("sisCone5GenJets")
trigTag   = cms.InputTag("TriggerResults::HLT")
vtxTag    = cms.InputTag("offlinePrimaryVertices")

### Analysis configuration ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",
    
    is_MC     = cms.bool(True),  # set to 'False' for real Data !
    is_SHERPA = cms.bool(False),  # set to 'True' if running on SHERPA
    
    jetTag    = jetTag,
    elecTag   = elecTag,
    muonTag   = muonTag,
    metTag    = metTag,
    genTag    = genTag,
    genJetTag = genJetTag,
    trigTag   = trigTag,
    vtxTag    = vtxTag,

    muopt  = cms.double(10.),
    muoeta = cms.double(2.5),
    elept  = cms.double(10.),
    eleeta = cms.double(2.5),
    jetpt  = cms.double(20.),
    jeteta = cms.double(2.5),
    jetfem = cms.double(0.9),
    met    = cms.double(20.),
    nele   = cms.int32(0),
    nmuo   = cms.int32(1),
    njet   = cms.int32(3),
    
    correction = cms.string('abs'),
    flavour    = cms.string('glu'),
    
)

## Necessary fixes to run 2.2.X on 2.1.X data
#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
#run22XonSummer08AODSIM(process)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOffTriggerMatchingOld( process )

##switchOnTrigger( process )
process.patTriggerSequence.remove( process.patTriggerMatcher )
process.patTriggerEvent.patTriggerMatches  = ()

### Define the paths
process.p = cms.Path(
    process.patDefaultSequence*
    process.patTrigger*
    process.patTriggerEvent*
    process.ACSkimAnalysis
    )
