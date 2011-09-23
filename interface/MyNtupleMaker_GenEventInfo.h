#ifndef MYNTUPLEMAKER_GENEVENTINFO_H
#define MYNTUPLEMAKER_GENEVENTINFO_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class MyNtupleMaker_GenEventInfo : public edm::EDProducer {
 public:
  explicit MyNtupleMaker_GenEventInfo(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup &);
  const edm::InputTag   genEvtInfoInputTag, puInfoInputTag;
};

#endif
