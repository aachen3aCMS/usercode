from aachen3a.ACSusyAnalysis.METfilter_cff import *

process = cms.Process("ANA")

isData=@ISDATA@
qscalehigh=@QSCALE_LOW@
qscalelow=@QSCALE_HIGH@
tauSwitch=@DOTAU@
pfeleSwitch=@DOPFELE@
Pythia8Switch=@PYTHIA8@
SherpaSwitch=@SHERPA@
MatchAllSwitch=@MATCHALL@
SusyParSwith=@SUSYPAR@
IsPythiaShowered=@ISPYTHIASHOWERED@
#~ qscalehigh=-1.
#~ qscalelow=-1.
#~ isData=False
#~ tauSwitch=False
#~ IsPythiaShowered=True
# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['SusyACSkimAnalysis'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('@GLOBALTAG@')
#~ process.GlobalTag.globaltag = cms.string('START44_V5::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoVertex.Configuration.RecoVertex_cff")

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi")

# Pythia GEN filter (used to correct for wrong 4-momentum-imbalance
# see https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1489.html
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")

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


### Input / output ###

# Input file
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring([
#    'file:///user/mweber/AODSIM/DoubleMu+Run2012A-13Jul2012-v1.root',
    'file:///user/mweber/AODSIM/Summer12_DR53X+DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph.root'
    ]),
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#PAT Stuff


process.patMuons.embedCombinedMuon = False;
process.patMuons.embedStandAloneMuon = False;



if tauSwitch:
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
if isData:
    _jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])
else:
    _jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
# Removed the pv Collection from the PF2PAT statement because it produced a segfault in some events
#~ usePF2PAT(process,runPF2PAT=True,jetAlgo='AK5', runOnMC=not isData, postfix=postfix,jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))
usePF2PAT(process,runPF2PAT=True,jetAlgo='AK5', runOnMC=not isData, postfix=postfix,jetCorrections=_jetCorrections)

usePFIso( process )


process.patElectrons.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
        )
process.patElectrons.isolationValuesNoPFId = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
        )
          
if isData:
    removeMCMatching(process, ['All'])
    
process.pfPileUpPFlow.checkClosestZVertex = False


process.load("PhysicsTools/PatAlgos/patSequences_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
#
#from aachen3a.ACSusyAnalysis.tauDiscriminator_cff import updateHPSPFTausRelPt
process.load("aachen3a/ACSusyAnalysis/tauDiscriminator_cff")

hpsTauIDSources+=[
    #("decayModeFinding", "DiscriminationByDecayModeFinding"),
    ("ByVLooseChargedIsolation","DiscriminationByVLooseChargedIsolation"),
    ("ByLooseChargedIsolation","DiscriminationByLooseChargedIsolation"),
    ("ByMediumChargedIsolation","DiscriminationByMediumChargedIsolation"),
    ("ByTightChargedIsolation","DiscriminationByTightChargedIsolation"),
    ("ByLooseCombinedIsolationDBRelSumPtCorr","DiscriminationByLooseCombinedIsolationDBRelSumPtCorr"),
    ("ByMediumCombinedIsolationDBRelSumPtCorr","DiscriminationByMediumCombinedIsolationDBRelSumPtCorr"),
    ("ByTightCombinedIsolationDBRelSumPtCorr","DiscriminationByTightCombinedIsolationDBRelSumPtCorr"),
    ("ByRawCombinedIsolationDBRelSumPtCorr","DiscriminationByRawCombinedIsolationDBRelSumPtCorr")
]
switchToPFTauHPS(process)


# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process)



# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *


#~ from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
#~ process.eleIsoSequence = setupPFElectronIso(process, 'patElectrons')
#~ process.phoIsoSequence = setupPFPhotonIso(process, 'patPhotons')


process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfType1CorrectedMet


switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA = True,
    doBTagging = False,
    jetCorrLabel = _jetCorrections,
    doType1MET = False,
    doJetID = True,
    jetIdLabel = "ak5"
)

# apply type I/type I + II PFMEt corrections to pat::MET object 
# and estimate systematic uncertainties on MET
from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
runMEtUncertainties(process,doApplyType0corr=False)

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfType1CorrectedMet
 
process.pfMETType0 = pfType1CorrectedMet.clone()
process.pfMETType0.applyType1Corrections = cms.bool(True)
process.pfMETType0.applyType0Corrections = cms.bool(True)

if tauSwitch:
    tausequence = cms.Sequence( process.PFTau)
else:
    tausequence = cms.Sequence()
if IsPythiaShowered:
    filtersequence = cms.Sequence( process.totalKinematicsFilter )
else:
    filtersequence = cms.Sequence()

# photon
# the isolation has to be reprocessed, cf. https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation#The_25th_May_update
# which require the EGammaAnalysisTools photonIsoProducer from:
# cvs co -r V00-00-21 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
# cvs up -r 1.13 EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h
# cvs up -r 1.20 EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc 
process.load('EGamma.EGammaAnalysisTools.photonIsoProducer_cfi')
process.phoPFIso.verbose = False
IsoValPhotonPF = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
                               cms.InputTag('phoPFIso:phIsoForGsfEle'),
                               cms.InputTag('phoPFIso:nhIsoForGsfEle'))
#it also requires to add phoPFIso in the path

#grrrr this is so anoying
#for isolation correction in HEEP ID
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)


################################
###                          ###
###  Analysis configuration  ###
###                          ###
################################


### Definition of all tags here
elecTag         = cms.InputTag("patElectrons")
pfelecTag      = cms.InputTag("patElectronsPFlow")
gsfelecTag         = cms.InputTag("gsfElectrons")
photonTag         = cms.InputTag("patPhotons")
pfjetTag        = cms.InputTag("patJetsPFlow")
muonTag      = cms.InputTag("patMuons")
PFmuonTag  = cms.InputTag("selectedPatMuonsPFlow")
tauTag         = cms.InputTag("patTaus")
metRAWTag= cms.InputTag("pfMet")
metType1Tag=cms.InputTag("patType1CorrectedPFMet")
metType0Tag=cms.InputTag("pfMETType0")
genTag        = cms.InputTag("genParticles")
genJetTag    = cms.InputTag("ak5GenJets")
vtxTag         = cms.InputTag("offlinePrimaryVertices")
reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE")
#ebhitsTag  = cms.InputTag("ecalRecHit", "EcalRecHitsEB");  # RECO
ebhitsTag     = cms.InputTag("reducedEcalRecHitsEB")   # AOD
HLTInputTag = cms.InputTag('TriggerResults','','HLT')
TriggerSummaryTag = cms.InputTag('hltTriggerSummaryAOD',"","HLT")



### Cuts and switches ###
process.ACSkimAnalysis = cms.EDFilter(
    "SusyACSkimAnalysis",

    is_MC      = cms.bool(not isData),  # set to 'False' for real Data !
    is_PYTHIA8 = cms.bool(Pythia8Switch),  # set to 'True' if running on PYTHIA8
    is_SHERPA  = cms.bool(SherpaSwitch),  # set to 'True' if running on SHERPA
    matchAll   = cms.bool(MatchAllSwitch),  # if True all truth leptons are matched else only ele and mu
    susyPar    = cms.bool(SusyParSwith),
    doTaus     = cms.bool(tauSwitch),
    doPFele    = cms.bool(pfeleSwitch),
    doCommonSkim = cms.bool(@COMMONSKIM@),

    # This is used to access the results of all filters that ran.
    #
    filters = cms.PSet(
        AllFilters = cms.PSet(
            process = cms.string( 'PAT' ),
            results = cms.string( 'TriggerResults' ),
            paths = cms.vstring()
            )
        ),

    pfjetTag   = pfjetTag,
    elecTag    = elecTag,
    gsfelecTag    = gsfelecTag,
    photonTag  = photonTag,
    pfelecTag  = pfelecTag,
    muonTag    = muonTag,
    PFmuonTag  = PFmuonTag,
    tauTag     = tauTag,
    metRAWTag= metRAWTag,
    metType1Tag=metType1Tag,
    metType0Tag=metType0Tag,
    genTag     = genTag,
    genJetTag  = genJetTag,
    vtxTag     = vtxTag,
    ebhitsTag  = ebhitsTag,
    reducedBarrelRecHitCollection = reducedBarrelRecHitCollection,
    reducedEndcapRecHitCollection = reducedEndcapRecHitCollection,
    HLTInputTag = HLTInputTag,
    TriggerSummaryTag = TriggerSummaryTag,

    IsoValPhotonPF = IsoValPhotonPF,


    qscale_low  = cms.double(qscalelow),
    qscale_high = cms.double(qscalehigh),
    muoptfirst = cms.double(@MUOPTFIRST@),
    muoptother = cms.double(@MUOPTOTHER@),
    muoeta     = cms.double(@MUOETA@),
    elept      = cms.double(@ELEPT@),
    eleeta     = cms.double(@ELEETA@),
    phopt      = cms.double(@PHOPT@),
    phoeta     = cms.double(@PHOETA@),
    pfelept    = cms.double(@PFELEPT@),
    pfeleeta   = cms.double(@PFELEETA@),
    pfjetpt    = cms.double(@PFJETPT@),
    pfjeteta   = cms.double(@PFJETETA@),
    taupt      = cms.double(@TAUPT@),
    taueta     = cms.double(@TAUETA@),
    met0       = cms.double(@MET0@),
    met1       = cms.double(@MET1@),
    met2       = cms.double(@MET2@),
    nele       = cms.int32(@NELE@),
    npho       = cms.int32(@NPHO@),
    npfele     = cms.int32(@NPFELE@),
    nmuo       = cms.int32(@NMUO@),
    npfjet     = cms.int32(@NPFJET@),
    ntau       = cms.int32(@NTAU@),
    htc        = cms.double(@HTC@),
    PFhtc      = cms.double(@PFHTC@),
    triggerContains=cms.string(@TRIGGERCONTAINS@),
    muoMinv    = cms.double(@MUOMINV@),
    muoDMinv   = cms.double(@MUODMINV@),

  #  btag       = cms.string('trackCountingHighEffBJetTags'),

    # Calo Jet ID
    jetselvers = cms.string("PURE09"),
    jetselqual = cms.string("LOOSE")



)

# switch off tau top level projection
getattr(process,"pfNoTau"+postfix).enable = False

### Define the paths
process.p = cms.Path(
    filtersequence*
    process.eIdSequence*
    process.goodOfflinePrimaryVertices*
    process.kt6PFJetsForIsolation*
    #tausequence*
    process.recoTauClassicHPSSequence*
    process.updateHPSPFTausRelPt*
    process.patDefaultSequence*
    getattr(process,"patPF2PATSequence"+postfix)*
    #process.pfParticleSelectionSequence*
    #process.eleIsoSequence*
    #process.phoIsoSequence*
    process.phoPFIso*
    process.producePFMETCorrections*
    process.pfMETType0
    #process.pfMETType1
    )



# The skimmer is in the endpath because then the results of all preceding paths
# are available. This is used to access the outcome of filters that ran.
#
process.ACSkimAnalysis.filterlist = cms.vstring()
addScrapingFilter( process )
addCSCHaloFilter( process )
addHCALLaserFilter( process )
addECALDeadCellFilterTP( process )
addECALDeadCellFilterBE( process )
addMuonFailureFilter( process )
addBadSuperCrystalFilter( process )
addTrackingFailureFilter( process )
addHBHENoiseFilter( process )
addecalLaserCorrFilter( process )
addecalLaserCorrFilter( process )
# This only works for datasets reprocessed in 2013 datasets (e.g. 22Jan2013 or Parked datasets)
addHCALLaserFilterFromTriggerResult ( process )


process.ACSkimAnalysis.filters.AllFilters.paths = process.ACSkimAnalysis.filterlist
process.ACSkimAnalysis.filters.AllFilters.process = process.name_()
#process.outpath = cms.EndPath(process.out2)
process.e = cms.EndPath(process.ACSkimAnalysis )
