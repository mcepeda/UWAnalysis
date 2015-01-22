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

class PATJetOverloader : public edm::EDProducer {
   public:
  
  explicit PATJetOverloader(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src")),
    leptons_(iConfig.getParameter<edm::InputTag>("leptons")),
    vertex_(iConfig.getParameter<edm::InputTag>("vertices"))
  {
    produces<pat::JetCollection>();
  }
  
  ~PATJetOverloader() {}
   private:

  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    using namespace edm;
    using namespace reco;
    std::auto_ptr<pat::JetCollection> jets(new pat::JetCollection);
    Handle<pat::JetCollection > cands;
    Handle<pat::MuonCollection > leptons; ///adding muons! :D
    Handle<reco::VertexCollection> vertices;
    bool verticesExist = iEvent.getByLabel(vertex_,vertices);

    if(verticesExist)
      verticesExist*=(vertices->size()>0)&&(vertices->at(0).isValid());


    if(iEvent.getByLabel(src_,cands)) 
      for(unsigned int  l=0;l!=cands->size();++l){

	pat::Jet jet = cands->at(l);


        //std::cout<<"PATJetOverloader"<<
        //             " Jet pT: "<<jet.pt()<<
        //             " eta: "<<jet.eta()<<
        //             " phi: "<<jet.phi()<<std::endl;

        // Jet Energy Corrections
        math::PtEtaPhiMLorentzVector p4_levelOne;
        math::PtEtaPhiMLorentzVector p4_levelTwo;
        math::PtEtaPhiMLorentzVector p4_levelThree;
        math::PtEtaPhiMLorentzVector p4_levelTwoThree;
        p4_levelOne      = jet.correctedP4("L1FastJet");
        p4_levelTwo      = jet.correctedP4("L2Relative");
        p4_levelThree    = jet.correctedP4("L3Absolute");
        //p4_levelTwoThree = jet.correctedP4("L2L3Residual");

        float pt_L1  = p4_levelOne.pt();
        float eta_L1 = p4_levelOne.eta();
        float phi_L1 = p4_levelOne.phi();
        float pt_L2  = p4_levelTwo.pt();
        float eta_L2 = p4_levelTwo.eta();
        float phi_L2 = p4_levelTwo.phi();
        float pt_L3  = p4_levelThree.pt();
        float eta_L3 = p4_levelThree.eta();
        float phi_L3 = p4_levelThree.phi();
        //float pt_L23  = p4_levelTwoThree.pt();
        //float eta_L23 = p4_levelTwoThree.eta();
        //float phi_L23 = p4_levelTwoThree.phi();

        //std::cout<<"L1 Jet Pt: "<<pt_L1<<" Eta: "<<eta_L1<<" Phi: "<<phi_L1<<std::endl;
        //std::cout<<"L2 Jet Pt: "<<pt_L2<<" Eta: "<<eta_L2<<" Phi: "<<phi_L2<<std::endl;
        //std::cout<<"L3 Jet Pt: "<<pt_L3<<" Eta: "<<eta_L3<<" Phi: "<<phi_L3<<std::endl;
        //std::cout<<"L23 Jet Pt: "<<pt_L23<<" Eta: "<<eta_L23<<" Phi: "<<phi_L23<<std::endl;

        jet.addUserFloat("pt_L1",pt_L1);
        jet.addUserFloat("eta_L1",eta_L1);
        jet.addUserFloat("phi_L1",phi_L1);
        jet.addUserFloat("pt_L2",pt_L2);
        jet.addUserFloat("eta_L2",eta_L2);
        jet.addUserFloat("phi_L2",phi_L2);
        jet.addUserFloat("pt_L3",pt_L3);
        jet.addUserFloat("eta_L3",eta_L3);
        jet.addUserFloat("phi_L3",phi_L3);
        //jet.addUserFloat("pt_L23",pt_L23);
        //jet.addUserFloat("eta_L23",eta_L23);
        //jet.addUserFloat("phi_L23",phi_L23);

	float pt=0.0;
	float sumPt=0.0;
	float sumPt2=0.0;
	float maxpt=0.0;
	float nTrack = 0.0;
	reco::PFCandidatePtr highptpfcand;
	reco::PFCandidatePtr pfcand;
	reco::PFCandidatePtr pfcandMuon1;
	reco::PFCandidatePtr pfcandMuon2;
	bool charged = false;
	bool muon1found;
	bool muon2found;
	float muonpt1 = 0;
	float muonpt2 = 0;
	float maxPFdist = 0;
	float Mu1ZMasstemp=0;
	float Mu1ZMass =0;
	float Mu2ZMasstemp=0;
	float Mu2ZMass =0;
	//int pdgid = 0;
	float dxy =0;
	float dz = 0;

	std::vector<reco::PFCandidatePtr> chargedParticles;
	std::vector<reco::PFCandidatePtr> muons;
	std::vector<reco::PFCandidatePtr> elecs;
	int t = 0; ///number of charged particles
	for(unsigned int j=0;j<jet.getPFConstituents().size();++j) {
	  pt = jet.getPFConstituents().at(j)->pt();
	  //pdgid = jet.getPFConstituents().at(j)->pdgId();

	  if(ROOT::Math::VectorUtil::DeltaR(jet.p4(),jet.getPFConstituents().at(j)->p4())>maxPFdist){
	    maxPFdist = ROOT::Math::VectorUtil::DeltaR(jet.p4(),jet.getPFConstituents().at(j)->p4());
	  }

	  pfcand = jet.getPFConstituents().at(j);
	  sumPt+=pt;
	  sumPt2+=pt*pt;
	  if(pt>maxpt){
	    maxpt = pt;
	    if(jet.getPFConstituents().at(j)->charge()!=0){
	      highptpfcand = jet.getPFConstituents().at(j);
	      charged = true;
	    }
	  }
	if(abs(pfcand->pdgId())==211&&pfcand->pt()>1){
	  t = t+1;
	  chargedParticles.push_back(pfcand);
	  
	}
	
	if(pfcand->trackRef().isNonnull())
	  nTrack +=1;
	  
	  if(abs(pfcand->pdgId())==13){
	    if(pfcand->pt()>muonpt1){
	      if(pfcandMuon1.isNonnull())
		pfcandMuon2 = pfcandMuon1;

	      muon1found = true;
	      pfcandMuon1 = jet.getPFConstituents().at(j);
	      muonpt1 = pfcandMuon1->pt();

	      if(pfcandMuon2.isNonnull())
		muonpt2 = pfcandMuon2->pt();
	      muon2found = true;
	    }
	    else if(pfcand->pt()>muonpt2&&pfcand->pt()<muonpt1){
	      pfcandMuon2 = jet.getPFConstituents().at(j);
	      muonpt2 = pfcandMuon2->pt();
	    }
	  }
	}
	if(iEvent.getByLabel(leptons_,leptons)){
	  for(unsigned int k =0; k!=leptons->size();k++){

	    if(pfcandMuon1.isNonnull()){
	      Mu1ZMasstemp=(pfcandMuon1->p4()+leptons->at(k).p4()).M();
	      if(Mu1ZMasstemp>Mu1ZMass)
		Mu1ZMass = Mu1ZMasstemp;
	    }
	    
	    if(pfcandMuon2.isNonnull()){
	      Mu1ZMasstemp=(pfcandMuon2->p4()+leptons->at(k).p4()).M();
	      if(Mu2ZMasstemp>Mu2ZMass)
		Mu2ZMass = Mu2ZMasstemp;
	    }
	  }
	}

	if(pfcandMuon1.isNonnull()){
	  TLorentzVector muVec(pfcandMuon1->p4().px(),pfcandMuon1->p4().py(),pfcandMuon1->p4().pz(),pfcandMuon1->energy());
	  TLorentzVector jet2(jet.px(),jet.py(),jet.pz());
	  //double jetpx=-jet.px();
	  //double jetpy=-jet.py();
	  //double jetpz=-jet.pz();
	  muVec.Boost(-jet2.BoostVector());
	  jet.addUserFloat("Mu1Boost",muVec.Pt());
	}

	if(pfcandMuon2.isNonnull()){
	  TLorentzVector muVec(pfcandMuon2->p4().px(),pfcandMuon2->p4().py(),pfcandMuon2->p4().pz(),pfcandMuon2->energy());
	  muVec.Boost(-jet.px(),-jet.py(),-jet.pz());
	  jet.addUserFloat("Mu2Boost",muVec.P());
	}

	  if(verticesExist&&maxpt>0&&charged) {
	    if(highptpfcand->trackRef().isNonnull()){
	      dxy = highptpfcand->trackRef()->dxy(vertices->at(0).position());
	      dz = highptpfcand->trackRef()->dz(vertices->at(0).position());
	      jet.addUserFloat("dxy_track",dxy);
	      jet.addUserFloat("dz_track",dz);
	    }
	    dxy = 0;
	    dz = 0;
	    if(pfcandMuon1.isNonnull())
	      if(pfcandMuon1->trackRef().isNonnull()&&verticesExist&&muon1found)
		if(pfcandMuon1->charge()!=0){
		  dxy = pfcandMuon1->trackRef()->dxy(vertices->at(0).position());

		  jet.addUserFloat("dxyMu1",dxy);

		  dz = pfcandMuon1->trackRef()->dz(vertices->at(0).position());
		  jet.addUserFloat("dzMu1",dz);
		}

	    dxy = 0;
	    dz = 0;
	    if(pfcandMuon2.isNonnull())
	      if(pfcandMuon2->trackRef().isNonnull()&&verticesExist&&muon2found&&pfcandMuon2->pt()>0)
		if(pfcandMuon2->charge()!=0){
		  dxy = pfcandMuon2->trackRef()->dxy(vertices->at(0).position());
		  jet.addUserFloat("dxyMu2",dxy);
		  dz = pfcandMuon2->trackRef()->dz(vertices->at(0).position());
		  jet.addUserFloat("dzMu2",dz);
	     }
	  }
        double mass_PV = -777; 
        double pt_PV = -777; 
        double eta_PV = -777; 
        double phi_PV = -777; 
        mass_PV = vertices->at(0).p4().M();
        pt_PV = vertices->at(0).p4().Pt();
        eta_PV = vertices->at(0).p4().Eta();
        phi_PV = vertices->at(0).p4().Phi();
        //std::cout<<"PV "<<mass_PV<<" "<<pt_PV<<" "<<eta_PV<<" "<<phi_PV<<std::endl;
	jet.addUserFloat("muon1ZMass",Mu1ZMass);
	jet.addUserFloat("muon2ZMass",Mu2ZMass);
	jet.addUserFloat("muon1pt",muonpt1);
	jet.addUserFloat("muon2pt",muonpt2);
	jet.addUserFloat("nTracks",nTrack);
	jet.addUserFloat("DR_PV_jet",maxPFdist);
	jet.addUserFloat("ptRMS_PV",sqrt(sumPt2/(sumPt*sumPt)));
        jet.addUserFloat("mass_PV",mass_PV); 
        jet.addUserFloat("pt_PV",pt_PV); 
        jet.addUserFloat("eta_PV",eta_PV); 
        jet.addUserFloat("phi_PV",phi_PV); 
	jets->push_back(jet);
      }
    iEvent.put(jets);
  } 
  
      // ----------member data ---------------------------
  edm::InputTag src_;
  edm::InputTag leptons_;
  edm::InputTag vertex_;
};

