import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SusyACSkimAnalysis'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-- Geometry ------------------------------------------------------------------
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
#    'file:/home/home1/institut_3a/magass/SUSY/CMSSW_3_3_4/src/aachen3a/ACSusyAnalysis/test/DATA.root'
#    'file:/user/magass/FIRST_COLLISIONS_23-11-2009.root'
    '/store/data/BeamCommissioning09/MinimumBias/RECO/rereco_FIRSTCOLL_v1/0083/EA91CC4D-28D9-DE11-833F-00261894396E.root'
#       [ '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/F62B040F-6CD8-DE11-9007-001D09F24664.root',
#       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/F4387297-74D8-DE11-996C-001D09F24F1F.root',
#       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/EE7B1AC4-6CD8-DE11-97BB-0030487A1FEC.root',
#       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/9264A28F-87D8-DE11-83A1-001D09F24763.root',
#       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/7C7495C2-71D8-DE11-ACF2-001D09F248F8.root',
#       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/4CAB3B6C-6BD8-DE11-845C-000423D9890C.root',
#       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/314/02E8544C-70D8-DE11-85CF-001617C3B66C.root' ] 
    ),
      duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
      firstRun = cms.untracked.uint32(122314)
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

    is_MC      = cms.bool(False),  # set to 'False' for real Data !
    is_SHERPA  = cms.bool(False),  # set to 'True' if running on SHERPA
    do_fatjets = cms.bool(False),  # set to 'True' for fat jets
                                   # if 'True', include process.CACaloSubjets below

    # IMPORTANT for QCD ! ! !
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
    elept  = cms.double(10.),
    eleeta = cms.double(2.5),
    jetpt  = cms.double(10.),
    jeteta = cms.double(2.5),
    jetfem = cms.double(0.9),
    met    = cms.double(0.),
    nele   = cms.int32(0),
    nmuo   = cms.int32(0),
    njet   = cms.int32(0),

    correction = cms.string('abs'),   # abs
    flavour    = cms.string('glu'),   # glu

    btag       = cms.string('trackCountingHighEffBJetTags'),

)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerSequence.remove( process.patTriggerMatcher )
process.patTriggerEvent.patTriggerMatches  = ()

# Be Careful !!!
# process.patTrigger.processName = "HLT8E29"
# process.patTriggerEvent.processName = "HLT8E29"


#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#run33xOn31xMC(process,
#                  jetSrc = cms.InputTag("antikt5CaloJets"),
#                  jetIdTag = "antikt5" )

### Define the paths
process.p = cms.Path(
    process.hcalnoise*
    process.patDefaultSequence*
    process.patTrigger*
    process.patTriggerEvent*
    process.ACSkimAnalysis
    )

