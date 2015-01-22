#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// system include files
#include <memory>
#include <algorithm>
#include <TLorentzVector.h>

// user include files

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

class PATMETCorrector : public edm::EDProducer {
    public:

        PATMETCorrector(const edm::ParameterSet& pset);
        virtual ~PATMETCorrector(){}
        void produce(edm::Event& iEvent, const edm::EventSetup& es);
    private:
	edm::InputTag srcMET_;
        bool isMC_;
};

PATMETCorrector::PATMETCorrector(
        const edm::ParameterSet& pset):
     srcMET_(pset.getParameter<edm::InputTag>("srcMET")),
     isMC_(pset.getParameter<bool>( "isMC" ))
     {
      produces<pat::METCollection>();
     }

void PATMETCorrector::produce(edm::Event& iEvent, const edm::EventSetup& es) {

    // get the collections and define the new output MET
    std::auto_ptr<pat::METCollection> outMET(new pat::METCollection);

    edm::Handle<pat::METCollection> inMET;

    edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsRaw;
    edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metType1;
    edm::Handle<pat::METCollection> handle_patPfMetT1xy;
    edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metT0pcT1Txy;

    iEvent.getByLabel("systematicsMET", "metsRaw",  handleLeafCand_metsRaw);
    iEvent.getByLabel("systematicsMET", "metType1", handleLeafCand_metType1);
    iEvent.getByLabel("patPfMetT1xy", handle_patPfMetT1xy);
    iEvent.getByLabel("systematicsMET", "metT0pcT1Txy",  handleLeafCand_metT0pcT1Txy);

    math::PtEtaPhiMLorentzVector rawP4(
     handleLeafCand_metsRaw.product()->front().pt(),
     handleLeafCand_metsRaw.product()->front().eta(),
     handleLeafCand_metsRaw.product()->front().phi(),
     handleLeafCand_metsRaw.product()->front().mass()
    );

    math::PtEtaPhiMLorentzVector type1P4(
     handleLeafCand_metType1.product()->front().pt(),
     handleLeafCand_metType1.product()->front().eta(),
     handleLeafCand_metType1.product()->front().phi(),
     handleLeafCand_metType1.product()->front().mass()
    );

    math::PtEtaPhiMLorentzVector T1TxyP4(
     handle_patPfMetT1xy->at(0).pt(),
     handle_patPfMetT1xy->at(0).eta(),
     handle_patPfMetT1xy->at(0).phi(),
     handle_patPfMetT1xy->at(0).mass()
    );

    math::PtEtaPhiMLorentzVector T0pcT1TxyP4(
     handleLeafCand_metT0pcT1Txy.product()->front().pt(),
     handleLeafCand_metT0pcT1Txy.product()->front().eta(),
     handleLeafCand_metT0pcT1Txy.product()->front().phi(),
     handleLeafCand_metT0pcT1Txy.product()->front().mass()
    );

    if(iEvent.getByLabel(srcMET_,inMET)){
     pat::MET newmet=inMET->at(0);

     //newmet.setP4(rawP4);
     newmet.setP4(type1P4);
     //newmet.setP4(T1TxyP4);
     //newmet.setP4(T0pcT1TxyP4);
     
     if(isMC_){
      edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsFullJESUp;
      edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsFullJESDown;
      edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsFullUESUp;
      edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsFullUESDown;
      edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsEESUp;
      edm::Handle<std::vector<reco::LeafCandidate>> handleLeafCand_metsEESDown;

      iEvent.getByLabel("systematicsMET", "metsFullJESUp",   handleLeafCand_metsFullJESUp);
      iEvent.getByLabel("systematicsMET", "metsFullJESDown", handleLeafCand_metsFullJESDown);
      iEvent.getByLabel("systematicsMET", "metsFullUESUp",   handleLeafCand_metsFullUESUp);
      iEvent.getByLabel("systematicsMET", "metsFullUESDown", handleLeafCand_metsFullUESDown);
      iEvent.getByLabel("systematicsMET", "metsEESUp",   handleLeafCand_metsEESUp);
      iEvent.getByLabel("systematicsMET", "metsEESDown", handleLeafCand_metsEESDown);

      math::PtEtaPhiMLorentzVector eesUpP4(
       handleLeafCand_metsEESUp.product()->front().pt(),
       handleLeafCand_metsEESUp.product()->front().eta(),
       handleLeafCand_metsEESUp.product()->front().phi(),
       handleLeafCand_metsEESUp.product()->front().mass()
      );

      math::PtEtaPhiMLorentzVector eesDnP4(
       handleLeafCand_metsEESDown.product()->front().pt(),
       handleLeafCand_metsEESDown.product()->front().eta(),
       handleLeafCand_metsEESDown.product()->front().phi(),
       handleLeafCand_metsEESDown.product()->front().mass()
      );

      math::PtEtaPhiMLorentzVector jesUpP4(
       handleLeafCand_metsFullJESUp.product()->front().pt(),
       handleLeafCand_metsFullJESUp.product()->front().eta(),
       handleLeafCand_metsFullJESUp.product()->front().phi(),
       handleLeafCand_metsFullJESUp.product()->front().mass()
      );

      math::PtEtaPhiMLorentzVector jesDnP4(
       handleLeafCand_metsFullJESDown.product()->front().pt(),
       handleLeafCand_metsFullJESDown.product()->front().eta(),
       handleLeafCand_metsFullJESDown.product()->front().phi(),
       handleLeafCand_metsFullJESDown.product()->front().mass()
      );

      math::PtEtaPhiMLorentzVector uesUpP4(
       handleLeafCand_metsFullUESUp.product()->front().pt(),
       handleLeafCand_metsFullUESUp.product()->front().eta(),
       handleLeafCand_metsFullUESUp.product()->front().phi(),
       handleLeafCand_metsFullUESUp.product()->front().mass()
      );

      math::PtEtaPhiMLorentzVector uesDnP4(
       handleLeafCand_metsFullUESDown.product()->front().pt(),
       handleLeafCand_metsFullUESDown.product()->front().eta(),
       handleLeafCand_metsFullUESDown.product()->front().phi(),
       handleLeafCand_metsFullUESDown.product()->front().mass()
      );
      newmet.addUserFloat("eesUp_pt",eesUpP4.pt());
      newmet.addUserFloat("eesUp_phi",eesUpP4.phi());
      newmet.addUserFloat("uesUp_pt",uesUpP4.pt());
      newmet.addUserFloat("uesUp_phi",uesUpP4.phi());
      newmet.addUserFloat("jesUp_pt",jesUpP4.pt());
      newmet.addUserFloat("jesUp_phi",jesUpP4.phi());
      newmet.addUserFloat("eesDn_pt",eesDnP4.pt());
      newmet.addUserFloat("eesDn_phi",eesDnP4.phi());
      newmet.addUserFloat("uesDn_pt",uesDnP4.pt());
      newmet.addUserFloat("uesDn_phi",uesDnP4.phi());
      newmet.addUserFloat("jesDn_pt",jesDnP4.pt());
      newmet.addUserFloat("jesDn_phi",jesDnP4.phi());
     }
     outMET->push_back(newmet);
    }
    iEvent.put(outMET);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMETCorrector);
