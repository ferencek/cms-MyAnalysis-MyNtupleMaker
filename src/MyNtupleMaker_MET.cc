#include "MyAnalysis/MyNtupleMaker/interface/MyNtupleMaker_MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/MET.h"

MyNtupleMaker_MET::MyNtupleMaker_MET(const edm::ParameterSet& iConfig) :

  inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  store_uncorrected_MET (iConfig.getParameter<bool>  ("StoreUncorrectedMET")),
  store_MET_significance (iConfig.getParameter<bool>  ("StoreMETSignificance"))
  
{
    produces <std::vector<double> > ( prefix + "Mag" + suffix );
    produces <std::vector<double> > ( prefix + "Phi" + suffix );
    produces <std::vector<double> > ( prefix + "SumET" + suffix );
    if ( store_uncorrected_MET )
    {
        produces <std::vector<double> > ( prefix + "MagUncorr" + suffix );
        produces <std::vector<double> > ( prefix + "PhiUncorr" + suffix );
        produces <std::vector<double> > ( prefix + "SumETUncorr" + suffix );
    }
    if ( store_MET_significance )
    {
        produces <std::vector<double> > ( prefix + "Significance" + suffix );
        produces <std::vector<double> > ( prefix + "SigMatrixDXX" + suffix );
        produces <std::vector<double> > ( prefix + "SigMatrixDXY" + suffix );
        produces <std::vector<double> > ( prefix + "SigMatrixDYX" + suffix );
        produces <std::vector<double> > ( prefix + "SigMatrixDYY" + suffix );
    }
}

void
MyNtupleMaker_MET::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::auto_ptr<std::vector<double> >  mag  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  sumet  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  maguncorr  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  phiuncorr  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  sumetuncorr  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  significance  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  sigmatrixdxx  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  sigmatrixdxy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  sigmatrixdyx  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  sigmatrixdyy  ( new std::vector<double>()  );
    
    //-----------------------------------------------------------------

    edm::Handle<std::vector<pat::MET> > mets;
    iEvent.getByLabel(inputTag, mets);

    if(mets.isValid())
    {
        edm::LogInfo("MyNtupleMaker_MET") << "Successfully obtained " << inputTag;

        for( std::vector<pat::MET>::const_iterator it = mets->begin(); it != mets->end(); ++it )
        {

            // fill in all the vectors
            mag->push_back( it->pt() );
            phi->push_back( it->phi() );
            sumet->push_back( it->sumEt() );

            if ( store_uncorrected_MET )
            {
                maguncorr->push_back( it->uncorrectedPt(pat::MET::uncorrALL) );
                phiuncorr->push_back( it->uncorrectedPhi(pat::MET::uncorrALL) );
                sumetuncorr->push_back( it->sumEt() - it->corSumEt(pat::MET::uncorrALL) );
            }

            if ( store_MET_significance )
            {
                // See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETSignificance#Known_Issues
                double sigmaX2= it->getSignificanceMatrix()(0,0);
                double sigmaY2= it->getSignificanceMatrix()(1,1);
                double sig = 0;
                if(sigmaX2<1.e10 && sigmaY2<1.e10) sig = it->significance();

                significance->push_back( sig );
                sigmatrixdxx->push_back( it->getSignificanceMatrix()(0,0) );
                sigmatrixdxy->push_back( it->getSignificanceMatrix()(0,1) );
                sigmatrixdyx->push_back( it->getSignificanceMatrix()(1,0) );
                sigmatrixdyy->push_back( it->getSignificanceMatrix()(1,1) );
                // See DataFormats/METReco/src/MET.cc
            }

        }
    }
    else
    {
      edm::LogError("MyNtupleMaker_MET") << "Error! Can't get the product " << inputTag;
    }
    
    //-----------------------------------------------------------------
    // put products in the event
    iEvent.put( mag, prefix + "Mag" + suffix );
    iEvent.put( phi, prefix + "Phi" + suffix );
    iEvent.put( sumet, prefix + "SumET" + suffix );
    if ( store_uncorrected_MET )
    {
      iEvent.put( maguncorr, prefix + "MagUncorr" + suffix );
      iEvent.put( phiuncorr, prefix + "PhiUncorr" + suffix );
      iEvent.put( sumetuncorr, prefix + "SumETUncorr" + suffix );
    }
    if ( store_MET_significance )
    {
      iEvent.put( significance, prefix + "Significance" + suffix );
      iEvent.put( sigmatrixdxx, prefix + "SigMatrixDXX" + suffix );
      iEvent.put( sigmatrixdxy, prefix + "SigMatrixDXY" + suffix );
      iEvent.put( sigmatrixdyx, prefix + "SigMatrixDYX" + suffix );
      iEvent.put( sigmatrixdyy, prefix + "SigMatrixDYY" + suffix );
    }
}
