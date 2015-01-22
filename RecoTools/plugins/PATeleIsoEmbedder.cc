/**
* @file PATeleIsoEmbedder.cc
* @author T.M.Perry
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
/**
* @class PATeleIsoEmbedder
* @brief embeds trigger/id scalefactors as userFloats
*/
class PATeleIsoEmbedder : public edm::EDProducer
{

    public:
        PATeleIsoEmbedder(const edm::ParameterSet& pset);
        virtual ~PATeleIsoEmbedder () {}
        void produce( edm::Event& iEvent, const edm::EventSetup& es );

    private:
        typedef reco::Candidate::LorentzVector FourVec;
        typedef reco::IsoDeposit IsoDeposit;
        edm::InputTag src_;
};

PATeleIsoEmbedder::PATeleIsoEmbedder( const edm::ParameterSet& pset)
{
    src_ = pset.getParameter<edm::InputTag>("src");
    produces<pat::ElectronCollection>();
}

void PATeleIsoEmbedder::produce( edm::Event& iEvent, const edm::EventSetup& es )
{
    std::auto_ptr<pat::ElectronCollection> out( new pat::ElectronCollection );
    edm::Handle<pat::ElectronCollection> electrons;
    if( iEvent.getByLabel(src_,electrons) )
    //std::cout<<"RecoTools/plugins/PATeleIsoEmbedder.cc"<<std::endl;
      for(unsigned int i=0;i<electrons->size();++i) {
        pat::Electron electron = electrons->at(i);

        float et;
        float cHI04, nHI04, gI04, pcHI04;
        float cHI03, nHI03, gI03, pcHI03;
        float Iso04, Iso03;
        et  = electron.pt();
        cHI03  = electron.isoDeposit(pat::IsolationKeys(4))->depositWithin(0.3);
        cHI04  = electron.isoDeposit(pat::IsolationKeys(4))->depositWithin(0.4);
        nHI03  = electron.isoDeposit(pat::IsolationKeys(5))->depositWithin(0.3);
        nHI04  = electron.isoDeposit(pat::IsolationKeys(5))->depositWithin(0.4);
        gI03   = electron.isoDeposit(pat::IsolationKeys(6))->depositWithin(0.3);
        gI04   = electron.isoDeposit(pat::IsolationKeys(6))->depositWithin(0.4);
        pcHI03 = electron.isoDeposit(pat::IsolationKeys(12))->depositWithin(0.3);
        pcHI04 = electron.isoDeposit(pat::IsolationKeys(12))->depositWithin(0.4);

        Iso03 = ( cHI04 + std::max( (nHI03 + gI03 - 0.5*pcHI03) ,0.0) ) / et;
        Iso04 = ( cHI04 + std::max( (nHI04 + gI04 - 0.5*pcHI04) ,0.0) ) / et;

        //  std::cout<<"cHI:  "<<cHI04<< " "<<cHI03<<std::endl;
        //  std::cout<<"nHI:  "<<nHI04<< " "<<nHI03<<std::endl;
        //  std::cout<<"gI:   "<<gI04<<  " "<<gI03<<std::endl;
        //  std::cout<<"pcHI: "<<pcHI04<<" "<<pcHI03<<std::endl;
        //  std::cout<<"Iso03: "<<Iso03<<" Iso04: "<<Iso04<<std::endl;

        electron.addUserFloat("Iso03",Iso03);
        electron.addUserFloat("Iso04",Iso04);
        electron.addUserFloat("cHI03",cHI03);
        electron.addUserFloat("cHI04",cHI04);
        electron.addUserFloat("nHI03",nHI03);
        electron.addUserFloat("nHI04",nHI04);
        electron.addUserFloat("gI03",gI03);
        electron.addUserFloat("gI04",gI04);
        electron.addUserFloat("pcHI03",pcHI03);
        electron.addUserFloat("pcHI04",pcHI04);
        out->push_back(electron);
    }

    iEvent.put( out );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATeleIsoEmbedder);
