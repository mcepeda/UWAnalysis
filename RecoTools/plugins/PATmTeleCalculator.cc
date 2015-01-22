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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/deltaPhi.h"

class PATmTeleCalculator : public edm::EDProducer {
    public:
        PATmTeleCalculator(const edm::ParameterSet& pset);
        virtual ~PATmTeleCalculator(){}
        void produce(edm::Event& iEvent, const edm::EventSetup& es);
    private:
	edm::InputTag srcMET_;
        edm::InputTag srcEle_;
};

PATmTeleCalculator::PATmTeleCalculator(const edm::ParameterSet& pset)
{
     srcMET_  = pset.getParameter<edm::InputTag>("srcMET");
     srcEle_ = pset.getParameter<edm::InputTag>("srcEle");
     produces<pat::ElectronCollection>();
}

void PATmTeleCalculator::produce(edm::Event& iEvent, const edm::EventSetup& es) {


    using namespace edm;
    using namespace reco;
    double pt_ele,pt_met,phi_ele,phi_met,dphi,mt = -1.;

    // get the collections
    std::auto_ptr<pat::ElectronCollection> out(new pat::ElectronCollection);

    edm::Handle<pat::ElectronCollection> electrons;

    edm::Handle<pat::METCollection> mets ;
    if(iEvent.getByLabel( srcEle_, electrons )){
     if(iEvent.getByLabel( srcMET_, mets)){
      if(electrons->size()>0){
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

        for ( size_t i = 0; i < electrons->size(); ++i )
        {
         pat::Electron electron = electrons->at(i);
         pt_ele = electron.pt();
         phi_ele = electron.phi();

         dphi = deltaPhi(electron.phi(),met.phi());
         mt = sqrt(2*pt_ele*pt_met*(1. - cos(dphi)));
         //std::cout<<"Electron: pt "<<pt_ele<<" phi "<<phi_ele<<std::endl;
         //std::cout<<"Raw"<<
         //  " pT: "<<pt_met<<" phi: "<<phi_met<<
         //  " mT: "<<mt<<std::endl;
         electron.addUserFloat("mt",mt);

         for ( size_t j = 0; j < 6; ++j){
          pt_met = met.userFloat(met_names[j]+"_pt");
          phi_met = met.userFloat(met_names[j]+"_phi");

          dphi = deltaPhi(phi_ele,phi_met);
          mt = sqrt(2*pt_ele*pt_met*(1. - cos(dphi)));
          //std::cout<<met_names[j]<<
          //  " pT: "<<pt_met<<" phi: "<<phi_met<<
          //  " mT: "<<mt<<std::endl;
          electron.addUserFloat("mt_"+met_names[j],mt);
         }

         out->push_back(electron);
        }
       }
      }
     }
    }
   iEvent.put(out);
   //std::cout<<"  ele"<<std::endl<<" "<<pt_ele<<" "<<phi_ele<<std::endl;
   //std::cout<<"  met"<<std::endl<<" "<<pt_met<<" "<<phi_met<<std::endl;
   //std::cout<<"  dphi = "<<dphi<<" mt = "<<mt<<std::endl;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATmTeleCalculator);
