import FWCore.ParameterSet.Config as cms

Vertices = cms.EDProducer("MyNtupleMaker_Vertices",
    InputTag = cms.InputTag('offlinePrimaryVertices'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1) # turned off if negative
)
