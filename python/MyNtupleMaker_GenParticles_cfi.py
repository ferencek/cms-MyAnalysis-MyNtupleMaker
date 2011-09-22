import FWCore.ParameterSet.Config as cms

GenParticles = cms.EDProducer("MyNtupleMaker_GenParticles",
    InputTag = cms.InputTag('genParticles'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(30) # turned off if negative
)
