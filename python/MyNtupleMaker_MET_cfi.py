import FWCore.ParameterSet.Config as cms

CaloMET = cms.EDProducer("MyNtupleMaker_MET",
    InputTag = cms.InputTag('patMETs'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(True),
    StoreMETSignificance = cms.bool(False)
)

PFMET = cms.EDProducer("MyNtupleMaker_MET",
    InputTag = cms.InputTag('patMETsPF'),
    Prefix = cms.string(''),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)
