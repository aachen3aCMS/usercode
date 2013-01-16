from PhysicsTools.PatAlgos.tools.coreTools import *


def addScrapingFilter( process ):
    process.scrapingFilter = cms.EDFilter( 'FilterOutScraping',
                                           applyfilter = cms.untracked.bool( True ),
                                           debugOn = cms.untracked.bool( False ),
                                           numtrack = cms.untracked.uint32( 10 ),
                                           thresh = cms.untracked.double( 0.25 )
                                           )
    process.p_scrapingFilter = cms.Path( process.scrapingFilter )
    process.ACSkimAnalysis.filterlist.append( 'p_scrapingFilter' )
    

def addCSCHaloFilter ( process ):
    process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
    process.p_CSCHaloFilter = cms.Path(process.CSCTightHaloFilter )
    process.ACSkimAnalysis.filterlist.append( 'p_CSCHaloFilter' )

def addHCALLaserFilter ( process ):
    process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cfi")
    process.p_HCALLaserFilter = cms.Path(process.hcallasereventfilter2012)
    process.ACSkimAnalysis.filterlist.append( 'p_HCALLaserFilter' )    

def addECALDeadCellFilterTP ( process ):
    process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
    ## For AOD and RECO recommendation to use recovered rechits
    process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
    process.p_ECALDeadCellFilterTP = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_ECALDeadCellFilterTP' )

def addECALDeadCellFilterBE ( process ):
    process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
    process.p_ECALDeadCellFilterBE = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_ECALDeadCellFilterBE' )

def addTrackingFailureFilter ( process ):
    process.goodVertices = cms.EDFilter(
        "VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )
    process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
    process.p_TrackingFailureFilter = cms.Path(process.goodVertices*process.trackingFailureFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_TrackingFailureFilter' )

def addMuonFailureFilter ( process ):
    process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')
    process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')
    process.p_MuonFailureFilter = cms.Path(process.greedyMuonPFCandidateFilter*process.inconsistentMuonPFCandidateFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_MuonFailureFilter' )


def addBadSuperCrystalFilter( process ):
    process.load('RecoMET.METFilters.eeBadScFilter_cfi')
    process.p_BadSuperCrystalFilter = cms.Path(process.eeBadScFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_BadSuperCrystalFilter' )

def addHBHENoiseFilter( process ):
    #see https://twiki.cern.ch/twiki/bin/view/CMS/HBHEAnomalousSignals2011
    #process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
    ## The iso-based HBHE noise filter ___________________________________________||
    process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
    process.p_HBHENoiseFilter = cms.Path(process.HBHENoiseFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_HBHENoiseFilter' )


def addecalLaserCorrFilter( process):
    ## The ECAL laser correction filter
    process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
    process.p_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
    process.ACSkimAnalysis.filterlist.append( 'p_ecalLaserCorrFilter' )    
    
    
def addtrkPOGFilters( process ):
    # The Tracking POG filters -- Strip Tracker coherent noise filters & log error filters
    ###The filter sets true if the event is bad so be warned
    process.load('RecoMET.METFilters.trackingPOGFilters_cfi')
    process.p_trkPOGFilters = cms.Path( process.trkPOGFilters )
    process.ACSkimAnalysis.filterlist.append( 'p_trkPOGFilters' )
