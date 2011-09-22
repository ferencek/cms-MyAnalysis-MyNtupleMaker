#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_GenEventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

MyNtupleMaker_GenEventInfo::MyNtupleMaker_GenEventInfo(const edm::ParameterSet& iConfig) :

  genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag"))
  
{
    produces <unsigned int> ( "ProcessID" );
    produces <double>       ( "PtHat" );
    produces <double>       ( "Weight" );
}

void
MyNtupleMaker_GenEventInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<unsigned int >         processID   ( new unsigned int() );
    std::auto_ptr<double >               ptHat ( new double() );
    std::auto_ptr<double >               weight ( new double() );
    
    //-----------------------------------------------------------------

    *processID.get() = 0;
    *ptHat.get() = 0.;
    *weight.get() = 0.;

    if( !iEvent.isRealData() ) {

        edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
        iEvent.getByLabel(genEvtInfoInputTag, genEvtInfoProduct);

        if( genEvtInfoProduct.isValid() )
        {
            edm::LogInfo("MyNtupleMaker_GenEventInfo") << "Successfully obtained " << genEvtInfoInputTag;

            *processID.get() = genEvtInfoProduct->signalProcessID();
            *ptHat.get() = ( genEvtInfoProduct->hasBinningValues() ? genEvtInfoProduct->binningValues()[0] : -999. );
            *weight.get() = genEvtInfoProduct->weight();
        }
        else
        {
            edm::LogError("MyNtupleMaker_GenEventInfo") << "Error! Can't get the product " << genEvtInfoInputTag;
        }
    }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( processID, "ProcessID" );
    iEvent.put( ptHat, "PtHat" );
    iEvent.put( weight, "Weight" );
}
