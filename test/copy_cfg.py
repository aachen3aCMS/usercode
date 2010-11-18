import FWCore.ParameterSet.Config as cms

process = cms.Process("COPY")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/122/318/E0E8C7BF-7DD8-DE11-93F4-001617DC1F70.root")
)


process.display = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('ttbar_AODSIM.root'),
    outputCommands = cms.untracked.vstring('keep *')
    )

process.outpath = cms.EndPath(process.display)
