//
// Adds a variety of userFloats to jets
// Authors: M.Cepeda, I.Ojalvo, T.M.Perry, M.Bachtis
//
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Math/GenVector/VectorUtil.h"
//
// class decleration
#include <boost/foreach.hpp>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"

class PATJetReader : public edm::EDProducer {
   public:
  
  explicit PATJetReader(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src"))
  {
//    produces<pat::JetCollection>();
  }
  
  ~PATJetReader() {}
   private:

  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    using namespace edm;
    using namespace reco;
    //std::auto_ptr<pat::JetCollection> jets(new pat::JetCollection);
    Handle<pat::JetCollection > cands;

    if(iEvent.getByLabel(src_,cands)) 
      for(unsigned int  l=0;l!=cands->size();++l){

	pat::Jet jet = cands->at(l);
        std::cout<<"PATJetReader"<<
                     " Jet pT: "<<jet.pt()<<
                     " eta: "<<jet.eta()<<
                     " phi: "<<jet.phi()<<std::endl;

	//jets->push_back(jet);
      }
    //iEvent.put(jets);
  } 
  
      // ----------member data ---------------------------
  edm::InputTag src_;
};

