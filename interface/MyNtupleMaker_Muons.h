#ifndef MYNTUPLEMAKER_MUONS_H
#define MYNTUPLEMAKER_MUONS_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class MyNtupleMaker_Muons : public edm::EDProducer {
 public:
  explicit MyNtupleMaker_Muons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup &);
  const edm::InputTag   inputTag;
  const std::string     prefix, suffix;
  const int             maxSize;
  const edm::InputTag   vtxInputTag;
};

#endif
