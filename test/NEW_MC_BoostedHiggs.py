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
#process.GlobalTag.globaltag = cms.string('START38_V14::All')
process.GlobalTag.globaltag = cms.string('')
process.load("Configuration.StandardSequences.MagneticField_cff")

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

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

# add iso deposits
from PhysicsTools.PatAlgos.tools.muonTools import addMuonUserIsolation
addMuonUserIsolation(process)

from PhysicsTools.PatAlgos.tools.electronTools import addElectronUserIsolation
addElectronUserIsolation(process)

process.patMuons.embedCombinedMuon = False;
process.patMuons.embedStandAloneMuon = False;

### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'/store/mc/Spring10/LM0/GEN-SIM-RECO/START3X_V26_S09-v1//0025/A6A36FEC-0348-DF11-BA10-E41F13181AB4.root'
    '/store/mc/Spring10/MinBias_7TeV-pythia8/GEN-SIM-RECO/START3X_V26B-v1/0002/2C0B016D-255E-DF11-B47D-0018FE284C82.root'
#    '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0006/36FF3D94-083C-DF11-BF54-0026189437FA.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix) 

# add ParticleFlow met
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')

# add TrackCorrected  met
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')

# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *

# Boosted Higgs
#addJetCollection(process,
#                 cms.InputTag('BoostedHiggsSubjets'), 'BHS', '',
#                 doJTA            = False,
#                 doBTagging       = True,
#                 jetCorrLabel     = ('AK5','Calo'), # DANGEROUS !!!!!!!!!!!!!!!!
#                 doType1MET       = False,
#                 doJetID          = False,
#                 genJetCollection = cms.InputTag("antikt5GenJets"))

# Add ParticleFlow jets - uncleaned
addJetCollection(process,cms.InputTag('ak5PFJets'), 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute'])),   # MC
                 # jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),  # DATA
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5pf"
                 )

# Add anti-kt 5 jets
addJetCollection(process,cms.InputTag('ak5CaloJets'), 'AK5', 'Calo',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5Calo', cms.vstring(['L2Relative', 'L3Absolute'])),  # MC
                 # jetCorrLabel = ('AK5Calo', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),  # DATA
                 doType1MET   = True,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )

# dummy - IMPORTANT
process.patJetCorrFactors.levels = cms.vstring(
    'L2Relative', 'L3Absolute'  # MC
    # 'L2Relative', 'L3Absolute', 'L2L3Residual'  #DATA
    )

# process.patJets.addTagInfos = cms.bool(False)  # AOD only

################################
###                          ###
###  Analysis configuration  ###
###                          ###
################################


### Definition of all tags here
elecTag    = cms.InputTag("patElectrons")
pfelecTag  = cms.InputTag("patElectronsPFlow")
calojetTag = cms.InputTag("patJetsAK5Calo")
pfjetTag   = cms.InputTag("patJetsPFlow")
muonTag    = cms.InputTag("patMuons")
metTag     = cms.InputTag("patMETsAK5Calo")
metTagTC   = cms.InputTag("patMETsTC")
metTagPF   = cms.InputTag("patMETsPFlow")
genTag     = cms.InputTag("genParticles")
genJetTag  = cms.InputTag("ak5GenJets")
vtxTag     = cms.InputTag("offlinePrimaryVertices")
ebhitsTag  = cms.InputTag("ecalRecHit", "EcalRecHitsEB");  # RECO
#ebhitsTag = cms.InputTag("reducedEcalRecHitsEB");   # AOD

### Cuts and switches ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",

    is_MC      = cms.bool(True),   # set to 'False' for real Data !
    is_SHERPA  = cms.bool(False),  # set to 'True' if running on SHERPA
    do_fatjets = cms.bool(False),  # set to 'True' for fat jets
                                   # if 'True', include process.BoostedHiggsSubjets (see example)

    # IMPORTANT for QCD -> configured via ./prepare.sh
    pthat_low  = cms.double(-1.),
    pthat_high = cms.double(-1.),

    calojetTag = calojetTag,
    pfjetTag   = pfjetTag,
    elecTag    = elecTag,
    pfelecTag  = pfelecTag,
    muonTag    = muonTag,
    metTag     = metTag,
    metTagTC   = metTagTC,
    metTagPF   = metTagPF,
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
    getattr(process,"patPF2PATSequence"+postfix)*
    process.ACSkimAnalysis
    )
# process.outpath = cms.EndPath(process.out)
