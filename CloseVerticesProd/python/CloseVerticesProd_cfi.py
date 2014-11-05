import FWCore.ParameterSet.Config as cms

demo = cms.EDProducer('CloseVerticesProd',
	PrimaryVertexCollection = cms.InputTag("pixelVertices"),
	DimuonVertexCollection = cms.InputTag("inclusiveSecondaryVertices"),
	CosMin = cms.untracked.double(0.99),
	MaxPrimaryVerticesPerDimuon = cms.untracked.uint32(3)
)
