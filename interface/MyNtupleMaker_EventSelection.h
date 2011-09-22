#ifndef MYNTUPLEMAKER_EVENTSELECTION_H
#define MYNTUPLEMAKER_EVENTSELECTION_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class MyNtupleMaker_EventSelection : public edm::EDProducer {
 public:
  explicit MyNtupleMaker_EventSelection(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup &);
  const edm::InputTag   vtxInputTag;
  const double          vtxMinNdof, vtxMaxAbsZ, vtxMaxd0;
  const edm::InputTag   trkInputTag;
  const unsigned int    nTracks;
  const double          hpTracksThreshold;
  const edm::InputTag   hcalNoiseInputTag;
  const edm::InputTag   beamHaloInputTag;
  const edm::InputTag   trackingFilterJetInputTag;
  const double          trackingFilterDzTrkVtxMax;
  const double          trackingFilterDxyTrkVtxMax;
  const double          trackingFilterMinSumPtOverHT;
  const edm::InputTag   ecalMaskedCellDRFilterInputTag, caloBoundaryDRFilterInputTag;
};

#endif
