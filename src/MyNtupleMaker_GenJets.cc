#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_GenJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

MyNtupleMaker_GenJets::MyNtupleMaker_GenJets(const edm::ParameterSet& iConfig) :

  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize (iConfig.getParameter<int> ("MaxSize"))
  
{
    produces <std::vector<double> > ( prefix + "Eta" + suffix );
    produces <std::vector<double> > ( prefix + "Phi" + suffix );
    produces <std::vector<double> > ( prefix + "P" + suffix );
    produces <std::vector<double> > ( prefix + "Pt" + suffix );
    produces <std::vector<double> > ( prefix + "Energy" + suffix );
    produces <std::vector<double> > ( prefix + "EmEnergy" + suffix );
    produces <std::vector<double> > ( prefix + "HadEnergy" + suffix );
    produces <std::vector<double> > ( prefix + "InvisibleEnergy" + suffix );
}

void
MyNtupleMaker_GenJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  p  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  emEnergy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  hadEnergy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  invisibleEnergy  ( new std::vector<double>()  );
    
    //-----------------------------------------------------------------

    if( !iEvent.isRealData() )
    {
        edm::Handle<reco::GenJetCollection> genJets;
        iEvent.getByLabel(inputTag, genJets);

        if( genJets.isValid() )
        {
            edm::LogInfo("MyNtupleMaker_GenJets") << "Total # GenJets: " << genJets->size();

            for( reco::GenJetCollection::const_iterator it = genJets->begin(); it != genJets->end(); ++it )
            {
                // exit from loop when you reach the required number of jets
                if(maxSize > 0 && eta->size() >= (unsigned int) maxSize)
                    break;

                // fill in all the vectors
                eta->push_back( it->eta() );
                phi->push_back( it->phi() );
                p->push_back( it->p() );
                pt->push_back( it->pt() );
                energy->push_back( it->energy() );
                emEnergy->push_back( it->emEnergy() );
                hadEnergy->push_back( it->hadEnergy() );
                invisibleEnergy->push_back( it->invisibleEnergy() );
            }
        }
        else
        {
            edm::LogError("MyNtupleMaker_GenJets") << "Error! Can't get the product " << inputTag;
        }
    }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( eta, prefix + "Eta" + suffix );
    iEvent.put( phi, prefix + "Phi" + suffix );
    iEvent.put( p, prefix + "P" + suffix );
    iEvent.put( pt, prefix + "Pt" + suffix );
    iEvent.put( energy, prefix + "Energy" + suffix );
    iEvent.put( emEnergy, prefix + "EmEnergy" + suffix );
    iEvent.put( invisibleEnergy, prefix + "InvisibleEnergy" + suffix );
}
