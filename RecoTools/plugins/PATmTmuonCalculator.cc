// WARNING  never finished writing this
// doesn't check for existance of muon/met and crashes

// system include files
#include <memory>
#include <algorithm>
#include <TLorentzVector.h>
#include <math.h> //fabs

#include "Math/GenVector/VectorUtil.h"
#include "boost/filesystem.hpp"

// user include files
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaPhi.h"

class PATmTmuonCalculator : public edm::EDProducer {
    public:
        PATmTmuonCalculator(const edm::ParameterSet& pset);
        virtual ~PATmTmuonCalculator(){}
        void produce(edm::Event& iEvent, const edm::EventSetup& es);
    private:
	edm::InputTag srcMET_;
        edm::InputTag srcMuon_;
};

PATmTmuonCalculator::PATmTmuonCalculator(const edm::ParameterSet& pset)
{
     srcMET_  = pset.getParameter<edm::InputTag>("srcMET");
     srcMuon_ = pset.getParameter<edm::InputTag>("srcMuon");
     produces<pat::MuonCollection>();
}

void PATmTmuonCalculator::produce(edm::Event& iEvent, const edm::EventSetup& es) {


    using namespace edm;
    using namespace reco;
    double pt_muon,pt_met,phi_muon,phi_met,dphi,mt = -1.;

    // get the collections
    std::auto_ptr<pat::MuonCollection> out(new pat::MuonCollection);

    edm::Handle<pat::MuonCollection> muons;

    edm::Handle<pat::METCollection> mets ;
    if(iEvent.getByLabel( srcMuon_, muons )){
     if(iEvent.getByLabel( srcMET_, mets)){
      if(muons->size()>0){
       if(mets->size()>0){

        pat::MET met=mets->at(0);
        pt_met = met.pt();

        std::string met_names[6] = {
         "eesUp",
         "eesDn",
         "uesUp",
         "uesDn",
         "jesUp",
         "jesDn"
        };
        
        for ( size_t i = 0; i < muons->size(); ++i )
        {
         pat::Muon muon = muons->at(i);
         pt_muon = muon.pt();
         phi_muon = muon.phi();

         dphi = deltaPhi(muon.phi(),met.phi());
         mt = sqrt(2*pt_muon*pt_met*(1. - cos(dphi)));
         muon.addUserFloat("mt",mt);

         for ( size_t j = 0; j < 6; ++j){
          pt_met = met.userFloat(met_names[j]+"_pt");
          phi_met = met.userFloat(met_names[j]+"_phi");

          dphi = deltaPhi(phi_muon,phi_met);
          mt = sqrt(2*pt_muon*pt_met*(1. - cos(dphi)));
          muon.addUserFloat("mt_"+met_names[j],mt);
         }

         out->push_back(muon);
        }
       }
      }
     }
    }
   iEvent.put(out);
   //std::cout<<"muon"<<std::endl<<" "<<pt_muon<<" "<<phi_muon<<std::endl;
   //std::cout<<"met"<<std::endl<<" "<<pt_met<<" "<<phi_met<<std::endl;
   //std::cout<<"dphi = "<<dphi<<" mt = "<<mt<<std::endl;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATmTmuonCalculator);
