#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_Vertices.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

MyNtupleMaker_Vertices::MyNtupleMaker_Vertices(const edm::ParameterSet& iConfig) :

  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize (iConfig.getParameter<int> ("MaxSize"))
  
{
    produces <std::vector<bool> >   ( prefix + "IsFake" + suffix );
    produces <std::vector<double> > ( prefix + "X" + suffix );
    produces <std::vector<double> > ( prefix + "Y" + suffix );
    produces <std::vector<double> > ( prefix + "Z" + suffix );
    produces <std::vector<double> > ( prefix + "XErr" + suffix );
    produces <std::vector<double> > ( prefix + "YErr" + suffix );
    produces <std::vector<double> > ( prefix + "ZErr" + suffix );
    produces <std::vector<double> > ( prefix + "Chi2" + suffix );
    produces <std::vector<double> > ( prefix + "NDF" + suffix );
    produces <std::vector<int> >    ( prefix + "NTracks" + suffix );
    produces <std::vector<int> >    ( prefix + "NTracksW05" + suffix );
}

void
MyNtupleMaker_Vertices::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<bool> >    isFake  ( new std::vector<bool>()  );
    std::auto_ptr<std::vector<double> >  x  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  y  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  z  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  xErr  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  yErr  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  zErr  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  chi2  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  ndf  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     nTracks  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     nTracksW05  ( new std::vector<int>()  );
    
    //-----------------------------------------------------------------

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(inputTag,primaryVertices);

    if(primaryVertices.isValid())
    {
        edm::LogInfo("MyNtupleMaker_Vertices") << "Total # Primary Vertices: " << primaryVertices->size();

        for( reco::VertexCollection::const_iterator it = primaryVertices->begin(); it != primaryVertices->end(); ++it )
        {
            // exit from loop when you reach the required number of vertices
            if( maxSize > 0 && isFake->size() >= (unsigned int) maxSize )
                break;
            
            isFake->push_back( it->isFake() );
            x->push_back( it->x() );
            y->push_back( it->y() );
            z->push_back( it->z() );
            xErr->push_back( it->xError() );
            yErr->push_back( it->yError() );
            zErr->push_back( it->zError() );
            chi2->push_back( it->chi2() );
            ndf->push_back( it->ndof() );
            nTracks->push_back( int(it->tracksSize()) );
            nTracksW05->push_back( it->nTracks(0.5) ); // number of tracks in the vertex with weight above 0.5
        }
    }
    else
    {
        edm::LogError("MyNtupleMaker_Vertices") << "Error! Can't get the product " << inputTag;
    }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( isFake, prefix + "IsFake" + suffix );
    iEvent.put( x, prefix + "X" + suffix );
    iEvent.put( y, prefix + "Y" + suffix );
    iEvent.put( z, prefix + "Z" + suffix );
    iEvent.put( xErr, prefix + "XErr" + suffix );
    iEvent.put( yErr, prefix + "YErr" + suffix );
    iEvent.put( zErr, prefix + "ZErr" + suffix );
    iEvent.put( chi2, prefix + "Chi2" + suffix );
    iEvent.put( ndf, prefix + "NDF" + suffix );
    iEvent.put( nTracks, prefix + "NTracks" + suffix );
    iEvent.put( nTracksW05, prefix + "NTracksW05" + suffix );
}
