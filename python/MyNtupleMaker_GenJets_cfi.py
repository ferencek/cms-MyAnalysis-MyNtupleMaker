import FWCore.ParameterSet.Config as cms

AK5GenJets = cms.EDProducer("MyNtupleMaker_GenJets",
    InputTag = cms.InputTag('ak5GenJets'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1) # turned off if negative
)

AK7GenJets = cms.EDProducer("MyNtupleMaker_GenJets",
    InputTag = cms.InputTag('ak7GenJets'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(-1) # turned off if negative
)
