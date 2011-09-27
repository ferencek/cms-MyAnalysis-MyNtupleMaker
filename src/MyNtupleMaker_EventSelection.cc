#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_EventSelection.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

MyNtupleMaker_EventSelection::MyNtupleMaker_EventSelection(const edm::ParameterSet& iConfig) :

    vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag")),
    vtxMinNdof(iConfig.getParameter<double>("VertexMinimumNdof")),
    vtxMaxAbsZ(iConfig.getParameter<double>("VertexMaxAbsZ")),
    vtxMaxd0(iConfig.getParameter<double>("VertexMaxd0")),
    trkInputTag(iConfig.getParameter<edm::InputTag>("TracksInputTag")),
    nTracks(iConfig.getParameter<unsigned int>("NumberOfTracks")),
    hpTracksThreshold(iConfig.getParameter<double>("HPTracksThreshold")),
    hcalNoiseInputTag(iConfig.getParameter<edm::InputTag>("HcalNoiseInputTag")),
    beamHaloInputTag(iConfig.getParameter<edm::InputTag>("BeamHaloInputTag")),
    trackingFilterJetInputTag   (iConfig.getParameter<edm::InputTag>("TrackingFailureJets")),
    trackingFilterDzTrkVtxMax    (iConfig.getParameter<double>       ("TrackingFailureDzTrkVtzMax")),
    trackingFilterDxyTrkVtxMax   (iConfig.getParameter<double>       ("TrackingFailureDxyTrkVtxMax")) ,
    trackingFilterMinSumPtOverHT(iConfig.getParameter<double>       ("TrackingFailureMinSumPtOverHT")),
    ecalMaskedCellDRFilterInputTag(iConfig.getParameter<edm::InputTag>("EcalMaskedCellDRFilterInputTag")),
    caloBoundaryDRFilterInputTag(iConfig.getParameter<edm::InputTag>("CaloBoundaryDRFilterInputTag"))
  
{
    produces <bool> ("PassPrimaryVertex");
    produces <bool> ("PassBeamScraping");
    produces <bool> ("PassHBHENoiseFilter");
    produces <bool> ("PassBeamHaloFilterLoose");
    produces <bool> ("PassBeamHaloFilterTight");
    produces <bool> ("PassTrackingFailure");
    produces <bool> ("PassEcalMaskedCellDRFilter");
    produces <bool> ("PassCaloBoundaryDRFilter");
}

void
MyNtupleMaker_EventSelection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<bool> passPrimaryVertex( new bool() );
    std::auto_ptr<bool> passBeamScraping( new bool() );
    std::auto_ptr<bool> passHBHENoiseFilter( new bool() );
    std::auto_ptr<bool> passBeamHaloFilterLoose( new bool() );
    std::auto_ptr<bool> passBeamHaloFilterTight( new bool() );
    std::auto_ptr<bool> passTrackingFailure ( new bool() ) ;
    std::auto_ptr<bool> passEcalMaskedCellDRFilter ( new bool() ) ;
    std::auto_ptr<bool> passCaloBoundaryDRFilter ( new bool() ) ;

    *passPrimaryVertex.get() = false;
    *passBeamScraping.get() = false;
    *passHBHENoiseFilter.get() = true;
    *passBeamHaloFilterLoose.get() = true;
    *passBeamHaloFilterTight.get() = true;
    *passTrackingFailure.get() = true;
    *passEcalMaskedCellDRFilter.get() = true;
    *passCaloBoundaryDRFilter.get() = true;
    
    //-----------------------------------------------------------------
    // good primary vertex filter
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(vtxInputTag,primaryVertices);

    if(primaryVertices.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Total # Primary Vertices: " << primaryVertices->size();

        for( reco::VertexCollection::const_iterator it = primaryVertices->begin(); it != primaryVertices->end(); ++it ) {
          if( !(it->isFake()) && it->ndof() > vtxMinNdof &&
              ( vtxMaxAbsZ<0 || fabs(it->z()) <= vtxMaxAbsZ ) &&
              ( vtxMaxd0<0 || fabs(it->position().rho()) <= vtxMaxd0 )
            ) *passPrimaryVertex.get() = true;
        }
    }
    else
    {
        edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << vtxInputTag;
    }

    // beam scraping filter
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel(trkInputTag,tracks);

    if(tracks.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Total # Tracks: " << tracks->size();

        int nHighPurity = 0;
        double fraction = 1.;
        reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::qualityByName("highPurity");

        if( tracks->size() > nTracks )
        {
            for( reco::TrackCollection::const_iterator it=tracks->begin(); it!=tracks->end(); ++it )
            {
                if( it->quality(trackQuality) ) nHighPurity++;
            }
            fraction = (double)nHighPurity/(double)tracks->size();
            if( fraction > hpTracksThreshold ) *passBeamScraping.get() = true;
        }
    }
    else
    {
        edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << trkInputTag;
    }

    // HBHE noise filter
    edm::Handle<bool> hbheFilterResult;
    iEvent.getByLabel(hcalNoiseInputTag, hbheFilterResult);

    if(hbheFilterResult.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Successfully obtained " << hcalNoiseInputTag;

        *passHBHENoiseFilter.get() = *hbheFilterResult;
    }
    else
    {
        edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << hcalNoiseInputTag;
    }

    // beam halo filter
    edm::Handle<reco::BeamHaloSummary> beamHaloSummary;
    iEvent.getByLabel(beamHaloInputTag, beamHaloSummary);

    if(beamHaloSummary.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Successfully obtained " << beamHaloInputTag;

        *passBeamHaloFilterLoose.get() = beamHaloSummary->CSCLooseHaloId();
        *passBeamHaloFilterTight.get() = beamHaloSummary->CSCTightHaloId();
    }
    else
    {
        edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << beamHaloInputTag;
    }

    // tracking failure filter
    edm::Handle<edm::View<reco::Jet> > jets;
    iEvent.getByLabel(trackingFilterJetInputTag, jets);

    if(jets.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Successfully obtained " << trackingFilterJetInputTag;
     
        double ht = 0;
        for (edm::View<reco::Jet>::const_iterator j = jets->begin(); j != jets->end(); ++j)
        {
            ht += j->pt();
        }

        double sumpt = 0;
        if (primaryVertices->size() > 0)
        {
            for (std::vector<reco::Track>::const_iterator tr = tracks->begin(); tr != tracks->end(); ++tr)
            {
                if (fabs(tr->dz(primaryVertices->at(0).position()))  > trackingFilterDzTrkVtxMax    ) continue;
                if (fabs(tr->dxy(primaryVertices->at(0).position())) > trackingFilterDxyTrkVtxMax   ) continue;
                sumpt += tr->pt();
            }
        }

        *passTrackingFailure.get() = ( (sumpt/ht) > trackingFilterMinSumPtOverHT );
    }
    else
    {
        edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << trackingFilterJetInputTag;
    }

    // ECAL masked cell and calo boundary energy filter ( https://twiki.cern.ch/twiki/bin/view/CMS/SusyEcalMaskedCellSummary )
    edm::Handle<int> ecalMaskedCellDRFilterResult;
    iEvent.getByLabel(ecalMaskedCellDRFilterInputTag, ecalMaskedCellDRFilterResult);

    edm::Handle<int> caloBoundaryDRFilterResult;
    iEvent.getByLabel(caloBoundaryDRFilterInputTag, caloBoundaryDRFilterResult);

    if(ecalMaskedCellDRFilterResult.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Successfully obtained " << ecalMaskedCellDRFilterInputTag;
     
        *passEcalMaskedCellDRFilter.get() = !(*ecalMaskedCellDRFilterResult);
    }
    // LogError call has to be commented out due to a bug/feature in the simpleDRfilter
//     else
//     {
//         edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << ecalMaskedCellDRFilterInputTag;
//     }

    if(caloBoundaryDRFilterResult.isValid())
    {
        edm::LogInfo("MyNtupleMaker_EventSelection") << "Successfully obtained " << caloBoundaryDRFilterInputTag;
     
        *passCaloBoundaryDRFilter.get() = !(*caloBoundaryDRFilterResult);
    }
    // LogError call has to be commented out due to a bug/feature in the simpleDRfilter
//     else
//     {
//         edm::LogError("MyNtupleMaker_EventSelection") << "Error! Can't get the product " << caloBoundaryDRFilterInputTag;
//     }

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put(passPrimaryVertex,"PassPrimaryVertex");
    iEvent.put(passBeamScraping,"PassBeamScraping");
    iEvent.put(passHBHENoiseFilter,"PassHBHENoiseFilter");
    iEvent.put(passBeamHaloFilterLoose,"PassBeamHaloFilterLoose");
    iEvent.put(passBeamHaloFilterTight,"PassBeamHaloFilterTight");
    iEvent.put(passTrackingFailure, "PassTrackingFailure");
    iEvent.put(passEcalMaskedCellDRFilter, "PassEcalMaskedCellDRFilter");
    iEvent.put(passCaloBoundaryDRFilter, "PassCaloBoundaryDRFilter");
}
