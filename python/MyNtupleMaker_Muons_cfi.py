import FWCore.ParameterSet.Config as cms

Muons = cms.EDProducer("MyNtupleMaker_Muons",
    InputTag = cms.InputTag('selectedPatMuons'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1), # turned off if negative
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)
