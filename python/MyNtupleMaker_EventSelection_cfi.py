import FWCore.ParameterSet.Config as cms

EventSelection = cms.EDProducer("MyNtupleMaker_EventSelection",
    VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
    VertexMinimumNdof = cms.double(4),
    VertexMaxAbsZ = cms.double(24.),
    VertexMaxd0 = cms.double(2.),
    TracksInputTag = cms.InputTag('generalTracks'),
    NumberOfTracks = cms.uint32(10),
    HPTracksThreshold = cms.double(0.25),
    HcalNoiseInputTag = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
    BeamHaloInputTag = cms.InputTag('BeamHaloSummary'),
    TrackingFailureJets           = cms.InputTag ('patJetsAK5PF'),
    TrackingFailureDzTrkVtzMax     = cms.double(1.0),
    TrackingFailureDxyTrkVtxMax    = cms.double(0.2),
    TrackingFailureMinSumPtOverHT = cms.double(0.10),
    EcalMaskedCellDRFilterInputTag = cms.InputTag('simpleDRfilter','deadCellStatus'),
    CaloBoundaryDRFilterInputTag = cms.InputTag('simpleDRfilter','boundaryStatus')
)
