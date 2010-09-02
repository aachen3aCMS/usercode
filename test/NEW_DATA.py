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
#process.GlobalTag.globaltag = cms.string('GR10_P_V4::All')
process.GlobalTag.globaltag = cms.string('')
process.load("Configuration.StandardSequences.MagneticField_cff")

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#-- Remove all Monte Carlo matching -------------------------------------------
removeMCMatching(process, ['All'])

# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Spring10")

# add iso deposits
from PhysicsTools.PatAlgos.tools.muonTools import addMuonUserIsolation
addMuonUserIsolation(process)

from PhysicsTools.PatAlgos.tools.electronTools import addElectronUserIsolation
addElectronUserIsolation(process)

### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([
    'file:/net/data_cms/institut_3a/gueth/punch_through.root']
#    'file:/home/home1/institut_3a/magass/SUSY/CMSSW_3_5_6/src/aachen3a/ACSusyAnalysis/test/CRAB/DATA_1.root',
#    'file:/home/home1/institut_3a/magass/SUSY/CMSSW_3_5_6/src/aachen3a/ACSusyAnalysis/test/CRAB/DATA_2.root']
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Output (PAT file only created if
# process.outpath = cms.EndPath(process.out)
# is called at the end
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName       = cms.untracked.string('dummy.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData   = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring('keep *')
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('out.root')
                                   )

# add ParticleFlow met
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# add TrackCorrected  met
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')

# run b-tagging sequences
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
run36xOn35xInput( process )

# Boosted Higgs
#addJetCollection(process,
#                 cms.InputTag('BoostedHiggsSubjets'), 'BHS', '',
#                 doJTA            = False,
#                 doBTagging       = True,
#                 jetCorrLabel     = ('AK5','Calo'), # DANGEROUS !!!!!!!!!!!!!!!!
#                 doType1MET       = False,
#                 doJetID          = False,
#                 genJetCollection = cms.InputTag("antikt5GenJets"))

# Add ParticleFlow jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','PF'),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )

# Add anti-kt 5 jets
addJetCollection(process,cms.InputTag('ak5CaloJets'),
                 'AK5', '',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5','Calo'),
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )


################################
###                          ###
###  Analysis configuration  ###
###                          ###
################################


### Definition of all tags here
elecTag   = cms.InputTag("patElectrons")
jetTag    = cms.InputTag("patJetsAK5")   # patJetsAK5PF or patJetsAK5
muonTag   = cms.InputTag("patMuons")
metTag    = cms.InputTag("patMETsAK5")
metTagTC  = cms.InputTag("patMETsTC")
metTagPF  = cms.InputTag("patMETsPF")
genTag    = cms.InputTag("genParticles")
genJetTag = cms.InputTag("ak5GenJets")
vtxTag    = cms.InputTag("offlinePrimaryVertices")

### Cuts and switches ###
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
    metTagTC  = metTagTC,
    metTagPF  = metTagPF,
    genTag    = genTag,
    genJetTag = genJetTag,
    vtxTag    = vtxTag,

    muopt  = cms.double(0.),
    muoeta = cms.double(25.),
    elept  = cms.double(0.),
    eleeta = cms.double(25.),
    jetpt  = cms.double(0.),
    jeteta = cms.double(25.),
    met    = cms.double(0.),
    nele   = cms.int32(0),
    nmuo   = cms.int32(0),
    njet   = cms.int32(0),

    correction = cms.string('abs'),
    flavour    = cms.string('glu'),

    btag       = cms.string('trackCountingHighEffBJetTags'),

)

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

# Boosted Higgs Configuration
from RecoJets.JetProducers.CaloJetParameters_cfi import CaloJetParameters
from RecoJets.JetProducers.AnomalousCellParameters_cfi import AnomalousCellParameters
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.BoostedHiggsParameters_cfi import BoostedHiggsParameters

# these are ignored, but required by VirtualJetProducer:
virtualjet_parameters = cms.PSet(jetAlgorithm=cms.string("SISCone"), rParam=cms.double(0.00001))


process.BoostedHiggsSubjets = cms.EDProducer("BoostedHiggsProducer",
    BoostedHiggsParameters,
    virtualjet_parameters,
    #this is required also for GenJets of PFJets:
    AnomalousCellParameters,
    CaloJetParameters
)

process.BoostedHiggsSubjets.jetSize       = cms.double(1.2)
process.BoostedHiggsSubjets.massThreshold = cms.double(0.667)
process.BoostedHiggsSubjets.rtyCut        = cms.double(0.3)
process.BoostedHiggsSubjets.ptMin         = cms.double(100.0)  # not JES corrected !!!
process.BoostedHiggsSubjets.nSubjets      = cms.int32(3)       # b bbar + radiation

### Define the paths
process.p = cms.Path(
#    process.BoostedHiggsSubjets*
    process.patDefaultSequence*
    process.ACSkimAnalysis
    )
# process.outpath = cms.EndPath(process.out)
