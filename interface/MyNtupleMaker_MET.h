#ifndef MYNTUPLEMAKER_MET_H
#define MYNTUPLEMAKER_MET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class MyNtupleMaker_MET : public edm::EDProducer {
 public:
  explicit MyNtupleMaker_MET(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup &);
  const edm::InputTag   inputTag;
  const std::string     prefix, suffix;
  const bool            store_uncorrected_MET;
  const bool            store_MET_significance;
};

#endif
