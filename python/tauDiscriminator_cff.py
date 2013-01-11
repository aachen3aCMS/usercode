import FWCore.ParameterSet.Config as cms
#import copy
from RecoTauTag.Configuration.HPSPFTaus_cff import hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr

hpsPFTauDiscriminationByLooseCombinedIsolationDBRelSumPtCorr = hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr.clone(
    applySumPtCut = False,
    applyRelativeSumPtCut = True,
    maximumSumPtCut = 0.1
    )

hpsPFTauDiscriminationByMediumCombinedIsolationDBRelSumPtCorr = hpsPFTauDiscriminationByLooseCombinedIsolationDBRelSumPtCorr.clone(
    maximumSumPtCut = 0.05
    )

hpsPFTauDiscriminationByTightCombinedIsolationDBRelSumPtCorr = hpsPFTauDiscriminationByLooseCombinedIsolationDBRelSumPtCorr.clone(
    maximumSumPtCut = 0.01
    )
hpsPFTauDiscriminationByRawCombinedIsolationDBRelSumPtCorr =hpsPFTauDiscriminationByLooseCombinedIsolationDBRelSumPtCorr.clone(
    applyRelativeSumPtCut=False,
    storeRelativeRawSumPt=cms.bool(True)
    )
    
updateHPSPFTausRelPt = cms.Sequence(
    hpsPFTauDiscriminationByLooseCombinedIsolationDBRelSumPtCorr*
    hpsPFTauDiscriminationByMediumCombinedIsolationDBRelSumPtCorr*
    hpsPFTauDiscriminationByTightCombinedIsolationDBRelSumPtCorr*
    hpsPFTauDiscriminationByRawCombinedIsolationDBRelSumPtCorr
)
