#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_GenEventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

MyNtupleMaker_GenEventInfo::MyNtupleMaker_GenEventInfo(const edm::ParameterSet& iConfig) :

  genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag")),
  puInfoInputTag(iConfig.getParameter<edm::InputTag>("PileupSummaryInfoInputTag"))
  
{
    produces <unsigned int> ( "ProcessID" );
    produces <double>       ( "PtHat" );
    produces <double>       ( "Weight" );
    produces <std::vector<unsigned int> > ( "PileUpNumberOfInteractions" );
    produces <std::vector<int> > ( "PileUpBunchCrossing" ) ;
}

void
MyNtupleMaker_GenEventInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<unsigned int>         processID   ( new unsigned int() );
    std::auto_ptr<double>               ptHat ( new double() );
    std::auto_ptr<double>               weight ( new double() );
    std::auto_ptr<std::vector<unsigned int> > pileupNumberOfInteractions  ( new std::vector<unsigned int>() );
    std::auto_ptr<std::vector<int> >   pileupBX ( new std::vector<int>() );

    *processID.get() = 0;
    *ptHat.get() = 0.;
    *weight.get() = 0.;
    
    //-----------------------------------------------------------------

    if( !iEvent.isRealData() ) {
        // GenEventInfo
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

        // PileupSummaryInfo
       edm::Handle<std::vector<PileupSummaryInfo> >  puInfo;
       iEvent.getByLabel(puInfoInputTag, puInfo);

       if(puInfo.isValid())
       {
           edm::LogInfo("MyNtupleMaker_GenEventInfo") << "Successfully obtained " << puInfoInputTag;
        
           for( std::vector<PileupSummaryInfo>::const_iterator it = puInfo->begin(); it != puInfo->end(); ++it )
           {
               pileupNumberOfInteractions->push_back( it->getPU_NumInteractions() );
               pileupBX->push_back( it -> getBunchCrossing());
           }
       }
       else
       {
           edm::LogError("MyNtupleMaker_GenEventInfo") << "Error! Can't get the product " << puInfoInputTag;
       }
    }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( processID, "ProcessID" );
    iEvent.put( ptHat, "PtHat" );
    iEvent.put( weight, "Weight" );
    iEvent.put( pileupNumberOfInteractions, "PileUpNumberOfInteractions" );
    iEvent.put( pileupBX, "PileUpBunchCrossing" );
}
