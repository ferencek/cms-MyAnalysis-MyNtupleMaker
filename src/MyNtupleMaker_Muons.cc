#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_Muons.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

MyNtupleMaker_Muons::MyNtupleMaker_Muons(const edm::ParameterSet& iConfig) :

  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize (iConfig.getParameter<int> ("MaxSize")),
  vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag"))
  
{
    produces <std::vector<bool> >   ( prefix + "IsGlobal" + suffix );
    produces <std::vector<bool> >   ( prefix + "IsTracker" + suffix );
    produces <std::vector<double> > ( prefix + "Eta" + suffix );
    produces <std::vector<double> > ( prefix + "EtaError" + suffix );
    produces <std::vector<double> > ( prefix + "Phi" + suffix );
    produces <std::vector<double> > ( prefix + "PhiError" + suffix );
    produces <std::vector<double> > ( prefix + "Pt" + suffix );
    produces <std::vector<double> > ( prefix + "PtError" + suffix );
    produces <std::vector<double> > ( prefix + "P" + suffix );
    produces <std::vector<double> > ( prefix + "Energy" + suffix );
    produces <std::vector<int> >    ( prefix + "Charge" + suffix );
    produces <std::vector<double> > ( prefix + "QoverPError" + suffix );
    produces <std::vector<double> > ( prefix + "TrkEta" + suffix );
    produces <std::vector<double> > ( prefix + "TrkPhi" + suffix );
    produces <std::vector<double> > ( prefix + "TrkPt" + suffix );
    produces <std::vector<double> > ( prefix + "TrkEtaError" + suffix );
    produces <std::vector<double> > ( prefix + "TrkPhiError" + suffix );
    produces <std::vector<double> > ( prefix + "TrkPtError" + suffix );
    produces <std::vector<int> >    ( prefix + "NHitsTracker" + suffix );
    produces <std::vector<double> > ( prefix + "TrkValidFraction" + suffix );
    produces <std::vector<int> >    ( prefix + "NHitsPixel" + suffix );
    produces <std::vector<int> >    ( prefix + "NHitsMuon" + suffix );
    produces <std::vector<int> >    ( prefix + "NMatchedChambers" + suffix );
    produces <std::vector<int> >    ( prefix + "NMatchedStations" + suffix );
    produces <std::vector<double> > ( prefix + "TrkChi2" + suffix );
    produces <std::vector<double> > ( prefix + "TrkNdof" + suffix );
    produces <std::vector<double> > ( prefix + "Chi2" + suffix );
    produces <std::vector<double> > ( prefix + "Ndof" + suffix );
    produces <std::vector<double> > ( prefix + "TrkIso" + suffix );
    produces <std::vector<double> > ( prefix + "EcalIso" + suffix );
    produces <std::vector<double> > ( prefix + "HcalIso" + suffix );
    produces <std::vector<double> > ( prefix + "HOIso" + suffix );
    produces <std::vector<double> > ( prefix + "PVXYDistance" + suffix );
    produces <std::vector<double> > ( prefix + "PVXYDistanceError" + suffix );
    produces <std::vector<double> > ( prefix + "PV3DDistance" + suffix );
    produces <std::vector<double> > ( prefix + "PV3DDistanceError" + suffix );
    produces <std::vector<double> > ( prefix + "BeamSpotXYDistance" + suffix );
    produces <std::vector<double> > ( prefix + "BeamSpotXYDistanceError" + suffix );
    produces <std::vector<double> > ( prefix + "BeamSpot3DDistance" + suffix );
    produces <std::vector<double> > ( prefix + "BeamSpot3DDistanceError" + suffix );
    produces <std::vector<double> > ( prefix + "ClosestPV3DDistance" + suffix );
    produces <std::vector<double> > ( prefix + "ClosestPVXYDistance" + suffix );
    produces <std::vector<double> > ( prefix + "ClosestPVZDistance" + suffix );
    produces <std::vector<int> >    ( prefix + "ClosestPV3DIndex" + suffix);
    produces <std::vector<int> >    ( prefix + "ClosestPVXYIndex" + suffix);
    produces <std::vector<int> >    ( prefix + "ClosestPVZIndex" + suffix);
}

void
MyNtupleMaker_Muons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<bool> >    isGlobal   ( new std::vector<bool>()  );
    std::auto_ptr<std::vector<bool> >    isTracker   ( new std::vector<bool>()  );
    std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  etaError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phiError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  ptError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  p  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     charge  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  qoverpError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkEta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkPhi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkPt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkEtaError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkPhiError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkPtError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     nHitsTracker ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  trkValidFraction ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     nHitsPixel ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     nHitsMuon ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     nMatchedChambers ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     nMatchedStations ( new std::vector<int>()  );    
    std::auto_ptr<std::vector<double> >  trkChi2   ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkNdof ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  chi2   ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  ndof   ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trkIso   ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  ecalIso  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  hcalIso  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  hoIso    ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pvXYDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pvXYDistanceError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pv3DDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pv3DDistanceError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  beamSpotXYDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  beamSpotXYDistanceError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  beamSpot3DDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  beamSpot3DDistanceError  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >   closestPV3DDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >   closestPVXYDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >   closestPVZDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >      closestPV3DIndex            ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >      closestPVXYIndex           ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >      closestPVZIndex            ( new std::vector<int>()  );
    
    //-----------------------------------------------------------------

    edm::Handle<std::vector<pat::Muon> > muons;
    iEvent.getByLabel(inputTag, muons);

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(vtxInputTag,primaryVertices);

    if(muons.isValid())
    {
        edm::LogInfo("MyNtupleMaker_Muons") << "Total # Muons: " << muons->size();

        for( std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end(); ++it )
        {
            // exit from loop when you reach the required number of muons
            if( maxSize > 0 && eta->size() >= (unsigned int) maxSize )
                break;

            bool global = it->isGlobalMuon();
            bool tracker = it->isTrackerMuon();

            isGlobal->push_back( global );
            isTracker->push_back( tracker );
            eta->push_back( it->eta() );
            etaError->push_back( ( global ? it->globalTrack()->etaError() : -999. ) );
            phi->push_back( it->phi() );
            phiError->push_back( ( global ? it->globalTrack()->phiError() : -999. ) );
            pt->push_back( it->pt() );
            ptError->push_back( ( global ? it->globalTrack()->ptError() : -999. ) );
            p->push_back( it->p() );
            energy->push_back( it->energy() );
            charge->push_back( it->charge() );
            qoverpError -> push_back ( ( global ? it->globalTrack()->qoverpError() : -999. ) );
            trkEta->push_back( ( tracker ? it->track()->eta() : -999. ) );
            trkEtaError->push_back( ( tracker ? it->track()->etaError() : -999. ) );
            trkPhi->push_back( ( tracker ? it->track()->phi() : -999. ) );
            trkPhiError->push_back( ( tracker ? it->track()->phiError() : -999. ) );
            trkPt->push_back( ( tracker ? it->track()->pt() : -999. ) );
            trkPtError->push_back( ( tracker ? it->track()->ptError() : -999. ) );
            nHitsTracker->push_back( ( tracker ? it->track()->hitPattern().numberOfValidTrackerHits() : -1 ) );
            trkValidFraction->push_back ( ( tracker ? it->track()->validFraction() : -999. ) );
            nHitsPixel->push_back( ( tracker ? it->track()->hitPattern().numberOfValidPixelHits() : -1 ) );
            nHitsMuon->push_back( ( global ? it->globalTrack()->hitPattern().numberOfValidMuonHits() : -1 ) );
            nMatchedChambers->push_back( it->numberOfMatches() );
            nMatchedStations->push_back( it->numberOfMatchedStations() );
            trkChi2->push_back( ( tracker ? it->track()->chi2() : -999. ) );
            trkNdof->push_back( ( tracker ? it->track()->ndof() : -999. ) );
            chi2->push_back( ( global ? it->globalTrack()->chi2() : -999. ) );
            ndof->push_back( ( global ? it->globalTrack()->ndof() : -999. ) );
            trkIso->push_back( it->isolationR03().sumPt );
            ecalIso->push_back( it->isolationR03().emEt );
            hcalIso->push_back( it->isolationR03().hadEt );
            hoIso->push_back( it->isolationR03().hoEt );
            pvXYDistance->push_back( it->dB(pat::Muon::PV2D) );
            pvXYDistanceError->push_back( it->edB(pat::Muon::PV2D) );
            pv3DDistance->push_back( it->dB(pat::Muon::PV3D) );
            pv3DDistanceError->push_back( it->edB(pat::Muon::PV3D) );;
            beamSpotXYDistance->push_back( it->dB(pat::Muon::BS2D) );
            beamSpotXYDistanceError->push_back( it->edB(pat::Muon::BS2D) );
            beamSpot3DDistance->push_back( it->dB(pat::Muon::BS3D) );
            beamSpot3DDistanceError->push_back( it->edB(pat::Muon::BS3D) );

            // vertex association
            double minPVDist3D = 99999.;
            double minPVDistXY = 99999.;
            double minPVDistZ  = 99999.;
            
            int closestPVIndex3D = -1;
            int closestPVIndexXY = -1;
            int closestPVIndexZ = -1;

            if( global || tracker )
            {
                if( primaryVertices.isValid() )
                {
                    edm::LogInfo("MyNtupleMaker_Muons") << "Total # Primary Vertices: " << primaryVertices->size();

                    for( reco::VertexCollection::const_iterator v_it = primaryVertices->begin(); v_it != primaryVertices->end(); ++v_it )
                    {
                        double distXY = - it->track()->dxy( v_it->position() );
                        double distZ = it->track()->dz( v_it->position() );
                        double dist3D = sqrt( pow(distXY,2) + pow(distZ,2) );

                        // find the closest primary vertex
                        if( dist3D < minPVDist3D )
                        {
                            minPVDist3D = dist3D;
                            closestPVIndex3D = int(std::distance(primaryVertices->begin(),v_it));
                        }

                        if( fabs(distXY) < minPVDistXY )
                        {
                            minPVDistXY = fabs(distXY);
                            closestPVIndexXY = int(std::distance(primaryVertices->begin(),v_it));
                        }

                        if( fabs(distZ) < minPVDistZ )
                        {
                            minPVDistZ = fabs(distZ);
                            closestPVIndexZ = int(std::distance(primaryVertices->begin(),v_it));
                        }
                    }
                }
                else
                {
                    edm::LogError("MyNtupleMaker_Muons") << "Error! Can't get the product " << vtxInputTag;
                }
            }

            closestPV3DDistance->push_back( minPVDist3D );
            closestPVXYDistance->push_back( minPVDistXY );
            closestPVZDistance->push_back( minPVDistZ );
            closestPV3DIndex->push_back( closestPVIndex3D );
            closestPVXYIndex->push_back( closestPVIndexXY );
            closestPVZIndex->push_back( closestPVIndexZ );
        }
    }
    else
    {
        edm::LogError("MyNtupleMaker_Muons") << "Error! Can't get the product " << inputTag;
    }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( isGlobal, prefix + "IsGlobal" + suffix );
    iEvent.put( isTracker, prefix + "IsTracker" + suffix );
    iEvent.put( eta, prefix + "Eta" + suffix );
    iEvent.put( etaError, prefix + "EtaError" + suffix );
    iEvent.put( phi, prefix + "Phi" + suffix );
    iEvent.put( phiError, prefix + "PhiError" + suffix );
    iEvent.put( pt, prefix + "Pt" + suffix );
    iEvent.put( ptError, prefix + "PtError" + suffix );
    iEvent.put( p, prefix + "P" + suffix );
    iEvent.put( energy, prefix + "Energy" + suffix );
    iEvent.put( charge, prefix + "Charge" + suffix );
    iEvent.put( qoverpError, prefix + "QoverPError" + suffix );
    iEvent.put( trkEta, prefix + "TrkEta" + suffix );
    iEvent.put( trkPhi, prefix + "TrkPhi" + suffix );
    iEvent.put( trkPt, prefix + "TrkPt" + suffix );
    iEvent.put( trkEtaError, prefix + "TrkEtaError" + suffix );
    iEvent.put( trkPhiError, prefix + "TrkPhiError" + suffix );
    iEvent.put( trkPtError, prefix + "TrkPtError" + suffix );
    iEvent.put( nHitsTracker, prefix + "NHitsTracker" + suffix );
    iEvent.put( trkValidFraction, prefix + "TrkValidFraction" + suffix );
    iEvent.put( nHitsPixel, prefix + "NHitsPixel" + suffix );
    iEvent.put( nHitsMuon, prefix + "NHitsMuon" + suffix );
    iEvent.put( nMatchedChambers, prefix + "NMatchedChambers" + suffix );    
    iEvent.put( nMatchedStations, prefix + "NMatchedStations" + suffix );
    iEvent.put( trkChi2, prefix + "TrkChi2" + suffix );
    iEvent.put( trkNdof, prefix + "TrkNdof" + suffix );
    iEvent.put( chi2, prefix + "Chi2" + suffix );
    iEvent.put( ndof, prefix + "Ndof" + suffix );
    iEvent.put( trkIso, prefix + "TrkIso" + suffix );
    iEvent.put( ecalIso, prefix + "EcalIso" + suffix );
    iEvent.put( hcalIso, prefix + "HcalIso" + suffix );
    iEvent.put( hoIso, prefix + "HOIso" + suffix );
    iEvent.put( pvXYDistance, prefix + "PVXYDistance" + suffix );
    iEvent.put( pvXYDistanceError, prefix + "PVXYDistanceError" + suffix );
    iEvent.put( pv3DDistance, prefix + "PV3DDistance" + suffix );
    iEvent.put( pv3DDistanceError, prefix + "PV3DDistanceError" + suffix );
    iEvent.put( beamSpotXYDistance, prefix + "BeamSpotXYDistance" + suffix );
    iEvent.put( beamSpotXYDistanceError, prefix + "BeamSpotXYDistanceError" + suffix );
    iEvent.put( beamSpot3DDistance, prefix + "BeamSpot3DDistance" + suffix );
    iEvent.put( beamSpot3DDistanceError, prefix + "BeamSpot3DDistanceError" + suffix );
    iEvent.put( closestPV3DDistance,prefix + "ClosestPV3DDistance" + suffix );
    iEvent.put( closestPVXYDistance,prefix + "ClosestPVXYDistance" + suffix );
    iEvent.put( closestPVZDistance,prefix + "ClosestPVZDistance" + suffix );
    iEvent.put( closestPV3DIndex,prefix + "ClosestPV3DIndex" + suffix);
    iEvent.put( closestPVXYIndex,prefix + "ClosestPVXYIndex" + suffix);
    iEvent.put( closestPVZIndex,prefix + "ClosestPVZIndex" + suffix);
}
