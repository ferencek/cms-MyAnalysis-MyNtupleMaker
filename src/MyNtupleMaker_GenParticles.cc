#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_GenParticles.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

MyNtupleMaker_GenParticles::MyNtupleMaker_GenParticles(const edm::ParameterSet& iConfig) :

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
    produces <std::vector<int> >    ( prefix + "PdgId" + suffix );
    produces <std::vector<int> >    ( prefix + "Status" + suffix );
    produces <std::vector<double> > ( prefix + "VX" + suffix );
    produces <std::vector<double> > ( prefix + "VY" + suffix );
    produces <std::vector<double> > ( prefix + "VZ" + suffix );
    produces <std::vector<int> >    ( prefix + "NumDaughters" + suffix );
    produces <std::vector<int> >    ( prefix + "MotherIndex" + suffix );
}

void
MyNtupleMaker_GenParticles::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  p  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     pdgId ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     status  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  vx  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  vy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  vz  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     numDaughters  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     motherIndex  ( new std::vector<int>()  );
    
    //-----------------------------------------------------------------

    if( !iEvent.isRealData() ) {
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByLabel(inputTag, genParticles);

        if( genParticles.isValid() )
        {
            edm::LogInfo("MyNtupleMaker_GenParticles") << "Total # GenParticles: " << genParticles->size();

            for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it )
            {
                // exit from loop when you reach the required number of GenParticles
                if( maxSize > 0 && eta->size() >= (unsigned int) maxSize )
                    break;

                // fill in all the vectors
                eta->push_back( it->eta() );
                phi->push_back( it->phi() );
                p->push_back( it->p() );
                pt->push_back( it->pt() );
                energy->push_back( it->energy() );
                pdgId->push_back( it->pdgId() );
                status->push_back( it->status() );
                vx->push_back( it->vx() );
                vy->push_back( it->vy() );
                vz->push_back( it->vz() );
                numDaughters->push_back( it->numberOfDaughters() );

                int idx = -1;
                for( reco::GenParticleCollection::const_iterator mIt = genParticles->begin(); mIt != genParticles->end(); ++mIt )
                {
                    if( it->mother() == &(*mIt) )
                    {
                        idx = std::distance(genParticles->begin(),mIt);
                        break;
                    }
                }
                motherIndex->push_back( idx );
            }
        }
        else
        {
            edm::LogError("MyNtupleMaker_GenParticles") << "Error! Can't get the product " << inputTag;
        }
    }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( eta, prefix + "Eta" + suffix );
    iEvent.put( phi, prefix + "Phi" + suffix );
    iEvent.put( p, prefix + "P" + suffix );
    iEvent.put( pt, prefix + "Pt" + suffix );
    iEvent.put( energy, prefix + "Energy" + suffix );
    iEvent.put( pdgId, prefix + "PdgId" + suffix );
    iEvent.put( status, prefix + "Status" + suffix );
    iEvent.put( vx, prefix + "VX" + suffix );
    iEvent.put( vy, prefix + "VY" + suffix );
    iEvent.put( vz, prefix + "VZ" + suffix );
    iEvent.put( numDaughters, prefix + "NumDaughters" + suffix );
    iEvent.put( motherIndex, prefix + "MotherIndex" + suffix );
}
