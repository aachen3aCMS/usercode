import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SelectorSequence',
                                         'EventSelectorAND',
                                         'HLTEventSelector',
                                         'JetEventSelector',
                                         'MetEventSelector',
                                         'SusyACSkimAnalysis',
                                         'PATLayer0Summary'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Include PAT Layer 0 & 1 if not running on pattified data
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")

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

### Definition of all tags here
elecTag   = cms.InputTag("selectedLayer1Electrons")
jetTag    = cms.InputTag("selectedLayer1Jets")
muonTag   = cms.InputTag("selectedLayer1Muons")
metTag    = cms.InputTag("selectedLayer1METs")                         
genJetTag = cms.InputTag("iterativeCone5GenJets")
trigTag   = cms.InputTag("TriggerResults::HLT")
vtxTag    = cms.InputTag("offlinePrimaryVertices")

### Analysis configuration ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",
    
    jetTag    = jetTag,
    elecTag   = elecTag,
    muonTag   = muonTag,
    metTag    = metTag,
    genJetTag = genJetTag,
    trigTag   = trigTag,
    vtxTag    = vtxTag,

    muopt  = cms.double(20.),
    muoeta = cms.double(2.5),
    elept  = cms.double(20.),
    eleeta = cms.double(2.5),
    jetpt  = cms.double(20.),
    jeteta = cms.double(2.5),
    jetfem = cms.double(0.9),
    met    = cms.double(20.),
    nele   = cms.int32(0),
    nmuo   = cms.int32(1),
    njet   = cms.int32(3),
    
    generator  = cms.string('genParticles'), # genParticles or source
    correction = cms.string('ABS'),
    flavour    = cms.string(''),
    
)

# Output file
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('test-SKIM.root'),
                               # save only events passing the full path
                               dropMetaDataForDroppedData = cms.untracked.bool(True), # Magic setting to reduce output size
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_selectedLayer1*_*_*',
                                                                      'keep *_TriggerResults_*_HLT',
                                                                      'keep *GenJet*_iterativeCone5GenJets_*_*',
                                                                      'keep *_genEvent*_*_*',
                                                                      'keep *_genMet*_*_*',
                                                                      'keep *_offlinePrimaryVertices_*_*',
                                                                      'keep *_genParticles_*_*'),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               )


## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

### Define the paths
process.p = cms.Path(
    process.patLayer0*process.patLayer1*
    process.ACSkimAnalysis
)
##process.outpath = cms.EndPath(process.out)
