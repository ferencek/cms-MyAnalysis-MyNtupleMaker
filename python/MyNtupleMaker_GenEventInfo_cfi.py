import FWCore.ParameterSet.Config as cms

GenEventInfo = cms.EDProducer("MyNtupleMaker_GenEventInfo",
    GenEventInfoInputTag = cms.InputTag('generator')
)
