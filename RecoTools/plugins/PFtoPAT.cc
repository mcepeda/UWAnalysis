// system include files
#include <memory>
#include <algorithm>
#include <TLorentzVector.h>

#include "boost/filesystem.hpp"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "Math/GenVector/VectorUtil.h"

class PFtoPAT : public edm::EDProducer{
 public:
  PFtoPAT(const edm::ParameterSet& pset);
  virtual ~PFtoPAT(){}
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);

 private:
   edm::InputTag src_;

};

PFtoPAT::PFtoPAT(const edm::ParameterSet& pset)
{
   src_ = pset.getParameter<edm::InputTag>("src");
   produces<pat::JetCollection>();
}

void PFtoPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
 {
  
  edm::Handle<edm::View<reco::Jet> > jets;
  //edm::Handle<edm::View<reco::PFJet> > jets;
  iEvent.getByLabel(src_, jets);
  std::auto_ptr<pat::JetCollection> patJets( new pat::JetCollection );

  for (edm::View<reco::Jet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++){
  //for (edm::View<reco::PFJet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++){

   unsigned int idx = itJet - jets->begin();
   edm::RefToBase<reco::Jet> jetRef = jets->refAt(idx);
   //edm::RefToBase<reco::PFJet> jetRef = jets->refAt(idx);
   //pat::Jet ajet = dynamic_cast<pat::Jet *>(jetRef.get());
   pat::Jet ajet(jetRef);
   std::cout<<"PFtoPAT"<<
                " Jet pT: "<<ajet.pt()<<
                " eta: "<<ajet.eta()<<
                " phi: "<<ajet.phi()<<std::endl;
   patJets->push_back(ajet);
  }

  iEvent.put(patJets);
 }

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFtoPAT);
