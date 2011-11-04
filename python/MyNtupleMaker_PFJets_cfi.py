import FWCore.ParameterSet.Config as cms

AK5PFJets = cms.EDProducer("MyNtupleMaker_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK5PF'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1), # turned off if negative
    JECUncertainty = cms.string('AK5PF'),
    ReadJECUncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)

AK7PFJets = cms.EDProducer("MyNtupleMaker_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK7PF'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1), # turned off if negative
    JECUncertainty = cms.string('AK7PF'),
    ReadJECUncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('goodOfflinePrimaryVertices')
)
