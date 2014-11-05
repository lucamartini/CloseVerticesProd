import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/p/pchen/public/Run2BphHLTDev/offline/CMSSW_7_2_0_pre6/src/RECO/outputA_RAW2DIGI_L1Reco_RECO.root'
    )
)

process.DimuonPVs = cms.EDProducer('CloseVerticesProd',
	PrimaryVertexCollection = cms.InputTag("pixelVertices"),
	DimuonVertexCollection = cms.InputTag("inclusiveSecondaryVertices"),
	MinCos = cms.untracked.double(-2.),

)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.DimuonPVs)

process.e = cms.EndPath(process.out)
