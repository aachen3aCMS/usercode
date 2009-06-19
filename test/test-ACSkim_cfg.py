import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SusyACSkimAnalysis'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V7::All')
#process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("PhysicsTools.PatAlgos.patSequences_cff")

#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet(process,newName='Winter09',oldName='Summer08Redigi') # change from old to Winter08

### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root'
#        '/pnfs/physik.rwth-aachen.de/cms/store/user/ata/test/ata/ADD_150_450_1jet_Mf2p5_2n_GRW_10TeV/ADD_150_450_1jet_Mf2p5_2n_GRW_10TeV_RECO/bfb9b1ca5534929b93d8c7d7bdccf0e7/RECO_2.root'

#    'file:/home/home1/institut_3a/magass/wprime_mu.root'
#    'file:/opt/user/magass/RECO_19.root'
#    '/pnfs/physik.rwth-aachen.de/cms/store/mc/Winter09/Wjets-madgraph/GEN-SIM-DIGI-RECO/IDEAL_V11_FastSim_v1/0043/B2B774A3-A7D3-DD11-A098-0011114FBAD4.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out.root')
                                   )

### Definition of all tags here
elecTag   = cms.InputTag("selectedLayer1Electrons")
jetTag    = cms.InputTag("selectedLayer1Jets")
muonTag   = cms.InputTag("selectedLayer1Muons")
metTag    = cms.InputTag("layer1METs")
genTag    = cms.InputTag("genParticles")
genJetTag = cms.InputTag("iterativeCone5GenJets")
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

    muopt  = cms.double(0.),
    muoeta = cms.double(2.5),
    elept  = cms.double(10.),
    eleeta = cms.double(2.5),
    jetpt  = cms.double(10.),
    jeteta = cms.double(2.5),
    jetfem = cms.double(0.9),
    met    = cms.double(10.),
    nele   = cms.int32(0),
    nmuo   = cms.int32(1),
    njet   = cms.int32(0),

    correction = cms.string('abs'),   # abs
    flavour    = cms.string('glu'),   # glu

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
                                                                      'keep *_genParticles_*_*',
                                                                      'keep *_patTrigger_*_*',
                                                                      'keep *_patTriggerEvent_*_*'
                                                                      ),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               )


## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)

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

##process.outpath = cms.EndPath(process.out)
