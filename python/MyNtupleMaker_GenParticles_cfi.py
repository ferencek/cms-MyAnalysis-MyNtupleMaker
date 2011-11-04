import FWCore.ParameterSet.Config as cms

GenParticles = cms.EDProducer("MyNtupleMaker_GenParticles",
    InputTag = cms.InputTag('genParticles'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    MaxSize = cms.int32(20), # turned off if negative
    PdgIDsOfInterest = cms.vint32(13) # will be stored regardless of the MaxSize parameter
)
