import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SusyACSkimAnalysis'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/home/home1/institut_3a/magass/SUSY/CMSSW_3_1_4/src/SUSYPAT.root'
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
    outputCommands = cms.untracked.vstring('drop *' )
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out.root')
                                   )


### Definition of all tags here
elecTag   = cms.InputTag("cleanLayer1Electrons")
jetTag    = cms.InputTag("cleanLayer1JetsAK5")
muonTag   = cms.InputTag("cleanLayer1Muons")
metTag    = cms.InputTag("layer1METsAK5")
genTag    = cms.InputTag("genParticles")
genJetTag = cms.InputTag("antikt5GenJets")
vtxTag    = cms.InputTag("offlinePrimaryVertices")

### Analysis configuration ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",
    
    is_MC     = cms.bool(True),    # set to 'False' for real Data !
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

    btag       = cms.string('trackCountingHighEffBJetTags'),
    
)

### Define the paths
process.p = cms.Path(
    process.ACSkimAnalysis
    )
