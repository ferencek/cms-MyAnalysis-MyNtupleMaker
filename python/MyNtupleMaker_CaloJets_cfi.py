import FWCore.ParameterSet.Config as cms

AK5CaloJets = cms.EDProducer("MyNtupleMaker_CaloJets",
    InputTag = cms.InputTag('selectedPatJets'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1), # turned off if negative
    JECUncertainty = cms.string('AK5Calo'),
    ReadJECUncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)

AK7CaloJets = cms.EDProducer("MyNtupleMaker_CaloJets",
    InputTag = cms.InputTag('selectedPatJetsAK7Calo'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1), # turned off if negative
    JECUncertainty = cms.string('AK7Calo'),
    ReadJECUncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)
