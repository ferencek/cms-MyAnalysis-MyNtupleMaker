#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_CaloJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

MyNtupleMaker_CaloJets::MyNtupleMaker_CaloJets(const edm::ParameterSet& iConfig) :

  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize (iConfig.getParameter<int> ("MaxSize")),
  jecUncTag(iConfig.getParameter<std::string>("JECUncertainty")),
  readJECUncertainty (iConfig.getParameter<bool>("ReadJECUncertainty")),
  vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag"))

{
    produces <std::vector<double> > ( prefix + "Eta" + suffix );
    produces <std::vector<double> > ( prefix + "Phi" + suffix );
    produces <std::vector<double> > ( prefix + "Pt" + suffix );
    produces <std::vector<double> > ( prefix + "PtRaw" + suffix );
    produces <std::vector<double> > ( prefix + "Energy" + suffix );
    produces <std::vector<double> > ( prefix + "EnergyRaw" + suffix );
    produces <std::vector<double> > ( prefix + "JECUnc" + suffix );
    produces <std::vector<double> > ( prefix + "L2L3ResJEC" + suffix );
    produces <std::vector<double> > ( prefix + "L3AbsJEC" + suffix );
    produces <std::vector<double> > ( prefix + "L2RelJEC" + suffix );
    produces <std::vector<double> > ( prefix + "L1OffsetJEC" + suffix );
    produces <std::vector<int> >    ( prefix + "PartonFlavor" + suffix );
    produces <std::vector<double> > ( prefix + "CombinedSecondaryVertexBJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "CombinedSecondaryVertexMVABJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "JetBProbabilityBTag" + suffix );
    produces <std::vector<double> > ( prefix + "JetProbabilityBTag" + suffix );
    produces <std::vector<double> > ( prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
    produces <std::vector<double> > ( prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
    produces <std::vector<double> > ( prefix + "SoftElectronByPtBJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "SoftElectronByIP3dBJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "SoftMuonBJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "SoftMuonByPtBJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "SoftMuonByIP3dBJetTag" + suffix );
    produces <std::vector<double> > ( prefix + "TrackCountingHighEffBTag" + suffix );
    produces <std::vector<double> > ( prefix + "TrackCountingHighPurBTag" + suffix );    
    produces <std::vector<int> >    ( prefix + "PassLooseID" + suffix);
    produces <std::vector<int> >    ( prefix + "PassTightID" + suffix);
    produces <std::vector<double> > ( prefix + "ClosestPVWeighted3DDistance" + suffix );
    produces <std::vector<double> > ( prefix + "ClosestPVWeightedXYDistance" + suffix );
    produces <std::vector<double> > ( prefix + "ClosestPVWeightedZDistance" + suffix );
    produces <std::vector<int> >    ( prefix + "ClosestPV3DIndex" + suffix);
    produces <std::vector<int> >    ( prefix + "ClosestPVXYIndex" + suffix);
    produces <std::vector<int> >    ( prefix + "ClosestPVZIndex" + suffix);
    produces <std::vector<double> > ( prefix + "BestPVTrackAssociationFactor" + suffix );
    produces <std::vector<int> >    ( prefix + "BestPVTrackAssociationIndex" + suffix);
}

JetIDSelectionFunctor jetIDLoose( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE );
JetIDSelectionFunctor jetIDTight( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT );

pat::strbitset ret = jetIDLoose.getBitTemplate();

void
MyNtupleMaker_CaloJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  pt_raw  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  energy_raw ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  jecUnc_vec ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  l2l3resJEC_vec ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  l3absJEC_vec ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  l2relJEC_vec ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  l1offsetJEC_vec ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     partonFlavor  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  combinedSecondaryVertexBJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  combinedSecondaryVertexMVABJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  jetBProbabilityBTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  jetProbabilityBTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  simpleSecondaryVertexHighEffBTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  simpleSecondaryVertexHighPurBTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  softElectronByPtBJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  softElectronByIP3dBJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  softMuonBJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  softMuonByPtBJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  softMuonByIP3dBJetTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trackCountingHighEffBTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  trackCountingHighPurBTag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >      passLooseID  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >      passTightID  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >   closestPVWeighted3DDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >   closestPVWeightedXYDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >   closestPVWeightedZDistance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >      closestPV3DIndex            ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >      closestPVXYIndex           ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >      closestPVZIndex            ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >   bestPVTrackAssociationFactor  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >      bestPVTrackAssociationIndex   ( new std::vector<int>()  );
    
    //-----------------------------------------------------------------

    // JEC Uncertainties
    JetCorrectionUncertainty *jecUnc = 0;
    if(readJECUncertainty) {
        // (see https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1075/1.html
        // and https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2367/1.html)
        // handle the jet corrector parameters collection
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
        // get the jet corrector parameters collection from the global tag
        iSetup.get<JetCorrectionsRecord>().get(jecUncTag,JetCorParColl);
        // get the uncertainty parameters from the collection
        JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
        // instantiate the jec uncertainty object
        jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }

    edm::Handle<std::vector<pat::Jet> > jets;
    iEvent.getByLabel(inputTag, jets);
 
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(vtxInputTag,primaryVertices);

    if(jets.isValid())
    {
        edm::LogInfo("MyNtupleMaker_CaloJets") << "Total # CaloJets: " << jets->size();

        for( std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it )
        {
            // exit from loop when you reach the required number of jets
            if( maxSize > 0 && eta->size() >= (unsigned int) maxSize )
                break;

            ret.set(false);
            int passjetLoose = 0;
            if( jetIDLoose( *it, ret ) ) passjetLoose = 1;

            ret.set(false);
            int passjetTight = 0;
            if( jetIDTight( *it, ret ) ) passjetTight = 1;

            if(readJECUncertainty)
            {
                jecUnc->setJetEta( it->eta() );
                jecUnc->setJetPt( it->pt() ); // the uncertainty is a function of the corrected pt
            }

            // status of JEC
            //std::cout << "Calo: currentJECLevel(): " << it->currentJECLevel() << std::endl;
            //std::cout << "Calo: currentJECSet(): " << it->currentJECSet() << std::endl;
            //-------------------

            // vertex association
            int bestPVIndex3Ddist = -1;
            int bestPVIndexXYdist = -1;
            int bestPVIndexZdist = -1;

            int bestPVIndexSharedTracks = -1;

            double minPVDist3D = 999999.;
            double minPVDistXY = 999999.;
            double minPVDistZ  = 999999.;
            double maxTrackAssocRatio = -999999.;


            // loop over primary vertices and jets and perform associations
            if(primaryVertices.isValid())
            {
                edm::LogInfo("MyNtupleMaker_CaloJets") << "Total # Primary Vertices: " << primaryVertices->size();

                // main vertex loop
                for( reco::VertexCollection::const_iterator v_it = primaryVertices->begin() ; v_it != primaryVertices->end() ; ++v_it )
                {

                    double sumweights = 0.0;
                    double dist3Dweighted = 0.0;
                    double distXYweighted = 0.0;
                    double distZweighted = 0.0;
                    double assocsumpttracks = 0.0;
                    double trackassociationratio = 0.000001;

                    // loop over tracks associated to a jet, calculate pT-weighted 3D distance to vertex and pT-weighted shared track ratio
                    const reco::TrackRefVector &jtracks = it->associatedTracks();
                    for(reco::TrackRefVector::const_iterator jtIt = jtracks.begin(); jtIt != jtracks.end(); ++jtIt)
                    {
                        if( jtIt->isNull() ) continue;
                        const reco::Track *jtrack = jtIt->get();
                        double trackptweight = jtrack->pt();
                        sumweights += trackptweight;

                        // weighted distance calculation
                        double distXY = fabs( jtrack->dxy(v_it->position()) );
                        double distZ = fabs( jtrack->dz(v_it->position()) );
                        dist3Dweighted += trackptweight*(sqrt(pow(distXY,2) + pow(distZ,2)));
                        distXYweighted += trackptweight*distXY;
                        distZweighted += trackptweight*distZ;

                        // loop over vertex tracks, find pT-weighted shared tracks
                        for(reco::Vertex::trackRef_iterator vtIt = v_it->tracks_begin(); vtIt != v_it->tracks_end(); ++vtIt)
                        {
                            if( vtIt->isNull() ) continue;
                            const reco::Track *vtrack = vtIt->get();
                            if( vtrack != jtrack ) continue;
                            assocsumpttracks += jtrack->pt();
                            break;
                        }
                    }

                    // divide distances by sum of weights
                    dist3Dweighted = dist3Dweighted / sumweights;
                    distXYweighted = distXYweighted / sumweights;
                    distZweighted  = distZweighted  / sumweights;

                    // calculate track association ratio
                    trackassociationratio = assocsumpttracks/sumweights;
                    
                    // find vertex with minimum weighted distance
                    if( dist3Dweighted < minPVDist3D )
                    {
                        minPVDist3D = dist3Dweighted;
                        bestPVIndex3Ddist = int(std::distance(primaryVertices->begin(),v_it));
                    }

                    if( distXYweighted < minPVDistXY )
                    {
                        minPVDistXY = distXYweighted;
                        bestPVIndexXYdist = int(std::distance(primaryVertices->begin(),v_it));
                    }

                    if( distZweighted < minPVDistZ )
                    {
                        minPVDistZ = distZweighted;
                        bestPVIndexZdist = int(std::distance(primaryVertices->begin(),v_it));
                    }

                    // find vertex with maximum weighted distance
                    if( trackassociationratio > maxTrackAssocRatio )
                    {
                        maxTrackAssocRatio = trackassociationratio ;
                        bestPVIndexSharedTracks = int(std::distance(primaryVertices->begin(),v_it));
                    }

                    //std::cout<<dist3Dweighted<<"  "<<distXYweighted<<"  "<<distZweighted<<"  "<<trackassociationratio<<"  "<<int(std::distance(primaryVertices->begin(),v_it))<<std::endl;
                }
                    //std::cout<<"---------------------"<<std::endl;
            }
            else
            {
                edm::LogError("MyNtupleMaker_CaloJets") << "Error! Can't get the product " << vtxInputTag;
            }

            //std::cout<<bestPVIndex3Ddist<<"  "<<bestPVIndexSharedTracks<<"  "<<minPVDist3D<<"  "<<minPVDistXY<<"  "<<minPVDistZ<<"  "<<maxTrackAssocRatio<<std::endl;
            //std::cout<<"------------------------------------------"<<std::endl;

            eta->push_back( it->eta() );
            phi->push_back( it->phi() );
            pt->push_back( it->pt() );
            pt_raw->push_back( it->correctedJet("Uncorrected").pt() );
            energy->push_back( it->energy() );
            energy_raw->push_back( it->correctedJet("Uncorrected").energy() );
            l2l3resJEC_vec->push_back( it->pt()/it->correctedJet("L3Absolute").pt() );
            l3absJEC_vec->push_back( it->correctedJet("L3Absolute").pt()/it->correctedJet("L2Relative").pt() );
            l2relJEC_vec->push_back( it->correctedJet("L2Relative").pt()/it->correctedJet("L1Offset").pt() );
            l1offsetJEC_vec->push_back( it->correctedJet("L1Offset").pt()/it->correctedJet("Uncorrected").pt() );
            if(readJECUncertainty)
              jecUnc_vec->push_back( jecUnc->getUncertainty(true) );
            else
              jecUnc_vec->push_back( -999. );
            partonFlavor->push_back( it->partonFlavour() );
            combinedSecondaryVertexBJetTag->push_back( it->bDiscriminator("combinedSecondaryVertexBJetTags") );
            combinedSecondaryVertexMVABJetTag->push_back( it->bDiscriminator("combinedSecondaryVertexMVABJetTags") );
            jetBProbabilityBTag->push_back( it->bDiscriminator("jetBProbabilityBJetTags") );
            jetProbabilityBTag->push_back( it->bDiscriminator("jetProbabilityBJetTags") );
            simpleSecondaryVertexHighEffBTag->push_back( it->bDiscriminator("simpleSecondaryVertexHighEffBJetTags") );
            simpleSecondaryVertexHighPurBTag->push_back( it->bDiscriminator("simpleSecondaryVertexHighPurBJetTags") );
            softElectronByPtBJetTag->push_back( it->bDiscriminator("softElectronByPtBJetTags") );
            softElectronByIP3dBJetTag->push_back( it->bDiscriminator("softElectronByIP3dBJetTags") );
            softMuonBJetTag->push_back( it->bDiscriminator("softMuonBJetTags") );
            softMuonByPtBJetTag->push_back( it->bDiscriminator("softMuonByPtBJetTags") );
            softMuonByIP3dBJetTag->push_back( it->bDiscriminator("softMuonByIP3dBJetTags") );
            trackCountingHighEffBTag->push_back( it->bDiscriminator("trackCountingHighEffBJetTags") );
            trackCountingHighPurBTag->push_back( it->bDiscriminator("trackCountingHighPurBJetTags") );
            passLooseID->push_back( passjetLoose );
            passTightID->push_back( passjetTight );
            closestPVWeighted3DDistance ->push_back( minPVDist3D);
            closestPVWeightedXYDistance ->push_back( minPVDistXY );
            closestPVWeightedZDistance ->push_back( minPVDistZ);
            closestPV3DIndex ->push_back( bestPVIndex3Ddist);
            closestPVXYIndex ->push_back( bestPVIndexXYdist);
            closestPVZIndex ->push_back( bestPVIndexZdist);
            bestPVTrackAssociationFactor ->push_back( maxTrackAssocRatio );
            bestPVTrackAssociationIndex ->push_back( bestPVIndexSharedTracks );
        }
    }
    else
    {
        edm::LogError("MyNtupleMaker_CaloJets") << "Error! Can't get the product " << inputTag;
    }

    // necessary to avoid introducing a memory leak
    delete jecUnc;

    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( eta, prefix + "Eta" + suffix );
    iEvent.put( phi, prefix + "Phi" + suffix );
    iEvent.put( pt, prefix + "Pt" + suffix );
    iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
    iEvent.put( energy, prefix + "Energy" + suffix );
    iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
    iEvent.put( jecUnc_vec, prefix + "JECUnc" + suffix );
    iEvent.put( l2l3resJEC_vec, prefix + "L2L3ResJEC" + suffix );
    iEvent.put( l3absJEC_vec, prefix + "L3AbsJEC" + suffix );
    iEvent.put( l2relJEC_vec, prefix + "L2RelJEC" + suffix );
    iEvent.put( l1offsetJEC_vec, prefix + "L1OffsetJEC" + suffix );
    iEvent.put( partonFlavor, prefix + "PartonFlavor" + suffix );
    iEvent.put( combinedSecondaryVertexBJetTag, prefix + "CombinedSecondaryVertexBJetTag" + suffix );
    iEvent.put( combinedSecondaryVertexMVABJetTag, prefix + "CombinedSecondaryVertexMVABJetTag" + suffix );    
    iEvent.put( jetBProbabilityBTag, prefix + "JetBProbabilityBTag" + suffix );
    iEvent.put( jetProbabilityBTag, prefix + "JetProbabilityBTag" + suffix );
    iEvent.put( simpleSecondaryVertexHighEffBTag, prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
    iEvent.put( simpleSecondaryVertexHighPurBTag, prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
    iEvent.put( softElectronByPtBJetTag, prefix + "SoftElectronByPtBJetTag" + suffix );
    iEvent.put( softElectronByIP3dBJetTag, prefix + "SoftElectronByIP3dBJetTag" + suffix );
    iEvent.put( softMuonBJetTag, prefix + "SoftMuonBJetTag" + suffix );
    iEvent.put( softMuonByPtBJetTag, prefix + "SoftMuonByPtBJetTag" + suffix );
    iEvent.put( softMuonByIP3dBJetTag, prefix + "SoftMuonByIP3dBJetTag" + suffix );
    iEvent.put( trackCountingHighEffBTag, prefix + "TrackCountingHighEffBTag" + suffix );
    iEvent.put( trackCountingHighPurBTag, prefix + "TrackCountingHighPurBTag" + suffix );
    iEvent.put( passLooseID, prefix + "PassLooseID" + suffix);
    iEvent.put( passTightID, prefix + "PassTightID" + suffix);
    iEvent.put( closestPVWeighted3DDistance,prefix + "ClosestPVWeighted3DDistance" + suffix );
    iEvent.put( closestPVWeightedXYDistance,prefix + "ClosestPVWeightedXYDistance" + suffix );
    iEvent.put( closestPVWeightedZDistance,prefix + "ClosestPVWeightedZDistance" + suffix );
    iEvent.put( closestPV3DIndex,prefix + "ClosestPV3DIndex" + suffix);
    iEvent.put( closestPVXYIndex,prefix + "ClosestPVXYIndex" + suffix);
    iEvent.put( closestPVZIndex,prefix + "ClosestPVZIndex" + suffix);
    iEvent.put( bestPVTrackAssociationFactor,prefix + "BestPVTrackAssociationFactor" + suffix );
    iEvent.put( bestPVTrackAssociationIndex,prefix + "BestPVTrackAssociationIndex" + suffix);
}
