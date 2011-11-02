from PhysicsTools.PatAlgos.tools.coreTools import *

process = cms.Process("ANA")

isData=@ISDATA@
#isData=False

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
process.GlobalTag.globaltag = cms.string('@GLOBALTAG@')
#process.GlobalTag.globaltag = cms.string('START42_V15B::All')
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
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
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

process.load("aachen3a.ACSusyAnalysis.pfIsoForLeptons_cff")
process.patMuons.embedCombinedMuon = False;
process.patMuons.embedStandAloneMuon = False;

### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/3298E7AE-B19C-E011-86B7-90E6BA442F24.root']
    ),
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
if isData:
    usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix, jetCorrections=('AK5PFchs',['L1FastJet','L2Relative','L3Absolute','L2L3Residual']))
    removeMCMatching(process, ['All'])
else:
     usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix, jetCorrections=('AK5PFchs',['L1FastJet', 'L2Relative', 'L3Absolute']))

adaptPFTaus(process,"hpsPFTau",postfix=postfix)


# for PFnoPU
process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpPFlow.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet = False


from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
    rParam = cms.double(0.6),
    #src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )

process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow", "rho")
if isData:
    process.patJetCorrFactorsPFlow.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual')
else:
    process.patJetCorrFactorsPFlow.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute') # MC


if isData:
    process.patJetCorrFactors.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual') # DATA
else:
    process.patJetCorrFactors.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute') # MC

process.patJetCorrFactors.useRho = True

# Add anti-kt 5 jets
# this is only for calo met is there an other way?
if isData:
    addJetCollection(process,cms.InputTag('ak5CaloJets'), 'AK5', 'Calo',
                 doJTA                  = True,
                 doBTagging        = True,
                 #jetCorrLabel     = ('AK5Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),  # MC
                 jetCorrLabel       = ('AK5Calo', cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])),  # DATA
                 doType1MET      = True,
                 doL1Cleaning    = False,
                 doL1Counters    = False,
                 genJetCollection =cms.InputTag("ak5GenJets"),
                 doJetID              = True,
                 jetIdLabel          = "ak5"
                 )
else:
    addJetCollection(process,cms.InputTag('ak5CaloJets'), 'AK5', 'Calo',
                 doJTA                  = True,
                 doBTagging        = True,
                 jetCorrLabel     = ('AK5Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),  # MC
                 #jetCorrLabel       = ('AK5Calo', cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])),  # DATA
                 doType1MET      = True,
                 doL1Cleaning    = False,
                 doL1Counters    = False,
                 genJetCollection =cms.InputTag("ak5GenJets"),
                 doJetID              = True,
                 jetIdLabel          = "ak5"
                 )



# Add the PV selector and KT6 producer to the sequence
getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsPFlow )


# add TrackCorrected  met
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')

# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *

# Add met with NoPU to the Collections
from aachen3a.ACSusyAnalysis.pfMET_cfi import *
process.pfMetPFnoPU         = pfMET.clone()
process.pfMetPFnoPU.alias   = 'pfMetNoPileUp'
process.pfMetPFnoPU.src     = 'pfNoPileUpPFlow'
process.pfMetPFnoPU.jets = cms.InputTag("ak5PFJets")



from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JetMETCorrections.Type1MET.pfMETCorrections_cff import *

process.pfType1CorrectedMet =   pfType1CorrectedMet.clone()
process.pfCandsNotInJet =       pfCandsNotInJet.clone()
process.pfJetMETcorr =          pfJetMETcorr.clone()
process.pfCandMETcorr =         pfCandMETcorr.clone()

if isData:
    process.pfJetMETcorr.jetCorrLabel = cms.string('ak5PFL1FastL2L3Residual')
else:
    process.pfJetMETcorr.jetCorrLabel = cms.string('ak5PFL1FastL2L3')
    
process.pfType1CorrectedMet.src = cms.InputTag("pfMetPFnoPU")

process.pfType1CorrectedPFMet = pfType1CorrectedMet.clone()
process.pfType1CorrectedPFMet.src     = cms.InputTag("pfMet")


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
metTagJPFnoPUType1 =cms.InputTag("pfType1CorrectedMet")
metTagcorMetGlobalMuons     = cms.InputTag("pfType1CorrectedPFMet")
#these two are ignored for now:
metTagHO    = cms.InputTag("metHO")
metTagNoHF= cms.InputTag("metNoHF")
genTag        = cms.InputTag("genParticles")
genJetTag    = cms.InputTag("ak5GenJets")
vtxTag         = cms.InputTag("offlinePrimaryVertices")
reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")
#ebhitsTag  = cms.InputTag("ecalRecHit", "EcalRecHitsEB");  # RECO
ebhitsTag     = cms.InputTag("reducedEcalRecHitsEB");   # AOD

### Cuts and switches ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",

    is_MC      = cms.bool(not isData),  # set to 'False' for real Data !
    is_SHERPA  = cms.bool(False),  # set to 'True' if running on SHERPA
    do_fatjets = cms.bool(False),  # set to 'True' for fat jets
    matchAll   = cms.bool(False),  # if True all truth leptons are matched else only ele and mu
    susyPar    = cms.bool(False),
    doCaloJet  = cms.bool(False),

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
    metTagJPFnoPUType1= metTagJPFnoPUType1,
    metTagcorMetGlobalMuons =metTagcorMetGlobalMuons,
    metTagHO = metTagHO,
    metTagNoHF = metTagNoHF,
    genTag     = genTag,
    genJetTag  = genJetTag,
    vtxTag     = vtxTag,
    ebhitsTag  = ebhitsTag,
    reducedBarrelRecHitCollection = reducedBarrelRecHitCollection,
    reducedEndcapRecHitCollection = reducedEndcapRecHitCollection,

    muoptfirst = cms.double(0.),
    muoptother = cms.double(0.),
    muoeta     = cms.double(25.),
    elept      = cms.double(0.),
    eleeta     = cms.double(25.),
    pfelept    = cms.double(0.),
    pfeleeta   = cms.double(25.),
    calojetpt  = cms.double(0.),
    calojeteta = cms.double(25.),
    pfjetpt    = cms.double(0.),
    pfjeteta   = cms.double(25.),
    taupt      = cms.double(0.),
    taueta     = cms.double(25.),
    metcalo    = cms.double(0.),
    metpf      = cms.double(0.),
    mettc      = cms.double(0.),
    nele       = cms.int32(0),
    npfele     = cms.int32(0),
    nmuo       = cms.int32(0),
    ncalojet   = cms.int32(0),
    npfjet     = cms.int32(0),
    ntau       = cms.int32(0),
    htc        = cms.double(0), 
    PFhtc      = cms.double(0),
    triggerContains=cms.string('None'),
    muoMinv    = cms.double(0), 
    muoDMinv   = cms.double(0), 

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
    process.PFTau*
    process.patDefaultSequence*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.offlinePrimaryVertices *
    process.pfMetPFnoPU*
    process.pfIsolationAllSequence*
    process.pfCandsNotInJet*
    process.pfJetMETcorr*
    process.pfCandMETcorr*
    process.pfType1CorrectedMet*
    process.pfType1CorrectedPFMet*
    process.ACSkimAnalysis
    )

#process.outpath = cms.EndPath(process.out2)
