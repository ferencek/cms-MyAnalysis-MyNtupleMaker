#ifndef MYNTUPLEMAKER_PFJETS_H
#define MYNTUPLEMAKER_PFJETS_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class MyNtupleMaker_PFJets : public edm::EDProducer {
 public:
  explicit MyNtupleMaker_PFJets(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup &);
  const edm::InputTag     inputTag;
  const std::string       prefix, suffix;
  const int               maxSize;
  const std::string       jecUncTag; 
  const bool              readJECUncertainty;
  const edm::InputTag     vtxInputTag;
};

#endif
