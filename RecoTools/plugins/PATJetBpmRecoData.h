// -*- C++ -*-
//
// Package:    PATMuonTrackVetoSelector
// Class:      PATMuonTrackVetoSelector
// 
/**\class PATMuonTrackVetoSelector PATMuonTrackVetoSelector.cc UWAnalysis/PATMuonTrackVetoSelector/src/PATMuonTrackVetoSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Nov 16
// $Id: PATJetBRreco.h
//
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

#include <boost/foreach.hpp>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"



class PATJetBpmRecoData : public edm::EDProducer {
   public:
  
  

  explicit PATJetBpmRecoData(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src")),
    leptons_(iConfig.getParameter<edm::InputTag>("leptons")),
    vertex_(iConfig.getParameter<edm::InputTag>("vertices"))
  {
    produces<pat::JetCollection>();
  }
  
  ~PATJetBpmRecoData() {}
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

    //printf("Gets Here1\n");
    if(verticesExist)
      verticesExist*=(vertices->size()>0)&&(vertices->at(0).isValid());

    //printf("Gets Here2\n");

    if(iEvent.getByLabel(src_,cands)) 
      for(unsigned int  l=0;l!=cands->size();++l){
	pat::Jet jet = cands->at(l);
	reco::PFCandidatePtr pfcand;
	std::vector<reco::PFCandidatePtr> chargedParticles;
	std::vector<reco::PFCandidatePtr> chargedParticlesNoD0;
	std::vector<reco::PFCandidatePtr> muons;
	std::vector<reco::PFCandidatePtr> elecs;
	std::vector<reco::PFCandidatePtr> secVertexBpm;
	std::vector<reco::PFCandidatePtr> secVertexBpmNoD0;
    //printf("Gets Here3\n");

	////adding D meson Code Here--> Take Vector of particles and loop
	const reco::SecondaryVertexTagInfo& secVertInfo = *jet.tagInfoSecondaryVertex("secondaryVertex");
        //printf("Gets Here3.1\n");
        //const reco::SecondaryVertexTagInfo* secInfo = jet.tagInfoSecondaryVertex("secondaryVertex");
	//printf("Gets Here3.1\n");
	//printf("nVert: ");
        //std::cout<<secVertInfo.nVertices()<<std::endl;
	if(secVertInfo.nVertices()>0){
        //printf("Gets Here3.2\n");
        //if (secInfo && secInfo->vertexTracks().size()>0) {
//std::cout<<"OOOOOOOOOOOOOOOOOOOOOOO"<<std::endl;
	  //float pt = secVertInfo.taggingVariables().size();//
	  //printf("track IP 2d above Charm: %f\n",pt);
	  //const reco::Vertex&sv= secInfo->secondaryVertex(0);
	  const reco::Vertex&sv= secVertInfo.secondaryVertex(0);
        //printf("Gets Here3.3\n");
	  //float pt = secVertInfo.taggingVariables().get(TaggingVariableName(41));//
	  //printf("number sec vertices %i\n",secVertInfo.secondaryVertex(0));
	  if ( sv.tracksSize()>1){
    //printf("Gets Here4\n");
	    for(reco::Vertex::trackRef_iterator track = sv.tracks_begin(); track!= sv.tracks_end(); ++track){
	      for(unsigned int pfj=0;pfj<jet.getPFConstituents().size();++pfj) {
		pfcand = jet.getPFConstituents().at(pfj);
		if(pfcand.isNonnull()){
    //printf("Gets Here5\n");
		  if(abs(pfcand->pdgId())==211&&pfcand->pt()>0.25&&pfcand->charge()!=0){ ///create charged particle collection
		    if((*track).key() == pfcand->trackRef().key()){
		      secVertexBpm.push_back(pfcand);
                    }
                  }
		}
	      }
            }
	  }
	}
	if(jet.pt()>20&&secVertexBpm.size()>1){
	  math::PtEtaPhiMLorentzVector total(0,0,0,0);
	  for(unsigned int k = 0; k<secVertexBpm.size()-1;k++){
	    math::PtEtaPhiMLorentzVector pion(secVertexBpm.at(k)->pt(),
					      secVertexBpm.at(k)->eta(),
					      secVertexBpm.at(k)->phi(),
					      0.140);

	    total +=pion;
    //printf("Gets Here6\n");
	  }
	  jet.addUserFloat("SV_M",total.M());
	  jet.addUserFloat("SV_pt",total.pt());
	  jet.addUserFloat("SV_eta",total.eta());
	  jet.addUserFloat("SV_phi",total.phi());
    //printf("Gets Here7\n");
	}


//	int finalpion;
//	int finalkaon;
//	double massD0 = 0;
//	math::PtEtaPhiMLorentzVector D0;
//	if(jet.pt()>20&&secVertexBpm.size()>1)
//	  for(unsigned int k = 0; k<secVertexBpm.size()-1;k++)
//	    for(unsigned int j = k+1; j<secVertexBpm.size();j++)
//	      if(secVertexBpm.at(k)->charge()+secVertexBpm.at(j)->charge()==0){
//		//printf("chargedparticles size %i\n",(int)secVertexBpm.size());
//		math::PtEtaPhiMLorentzVector kaon(secVertexBpm.at(k)->pt(),
//						  secVertexBpm.at(k)->eta(),
//						  secVertexBpm.at(k)->phi(),
//						  0.494);
//		math::PtEtaPhiMLorentzVector pion(secVertexBpm.at(j)->pt(),
//						  secVertexBpm.at(j)->eta(),
//						  secVertexBpm.at(j)->phi(),
//						  0.140);
//		
//		if(secVertexBpm.at(k)->charge()>0){
//		  pion =math::PtEtaPhiMLorentzVector(secVertexBpm.at(k)->pt(),
//						     secVertexBpm.at(k)->eta(),
//						     secVertexBpm.at(k)->phi(),
//						     0.140);
//		  kaon =math::PtEtaPhiMLorentzVector(secVertexBpm.at(j)->pt(),
//						     secVertexBpm.at(j)->eta(),
//						     secVertexBpm.at(j)->phi(),
//						     0.494);
//		}
//		D0 = pion+kaon;
//		double mass = (pion+kaon).M();
//
//		if(fabs(mass-1.864)<fabs(massD0-1.864)){ //find the candidate closest to D0 mass
//		  massD0 = mass;
//		  finalpion = k;
//		  finalkaon = j;
//		  //printf("k= %i th particle pt %f, j= %i th particle pt %f \n",k,secVertexBpm.at(k)->pt(),j,secVertexBpm.at(j)->pt());
//		}
//	      }
//
//	jet.addUserFloat("sec_massD0",massD0);
//	////try d0 pi pi pi decay mode
//	int charge = 0;
//	double massBpm2 = 0;
//	math::PtEtaPhiMLorentzVector Bpm2;
//	if(jet.pt()>20&&secVertexBpm.size()>4&&massD0>0){
//	  secVertexBpmNoD0 = secVertexBpm;
//	  secVertexBpmNoD0.erase(secVertexBpmNoD0.begin()+finalpion);
//	  
//	  if(finalpion>finalkaon)
//	    secVertexBpmNoD0.erase(secVertexBpmNoD0.begin()+finalkaon);
//	  else
//	    secVertexBpmNoD0.erase(secVertexBpmNoD0.begin()+(finalkaon-1));
//	  
//	  for(unsigned int k = 0; k<secVertexBpmNoD0.size()-1;k++)
//	    for(unsigned int j = k+1; j<secVertexBpmNoD0.size();j++)
//	      for(unsigned int n = j+1; n<secVertexBpmNoD0.size();n++)
//		if(abs(secVertexBpmNoD0.at(k)->charge()
//		       +secVertexBpmNoD0.at(j)->charge()
//		       +secVertexBpmNoD0.at(n)->charge())==1){
//		  //printf("chargedparticles size %i\n",(int)secVertexBpm.size());
//		  math::PtEtaPhiMLorentzVector pion1(secVertexBpmNoD0.at(k)->pt(),
//						    secVertexBpmNoD0.at(k)->eta(),
//						    secVertexBpmNoD0.at(k)->phi(),
//						    0.140);
//		  math::PtEtaPhiMLorentzVector pion2(secVertexBpmNoD0.at(j)->pt(),
//						    secVertexBpmNoD0.at(j)->eta(),
//						    secVertexBpmNoD0.at(j)->phi(),
//						    0.140);
//		  math::PtEtaPhiMLorentzVector pion3(secVertexBpmNoD0.at(j)->pt(),
//						    secVertexBpmNoD0.at(j)->eta(),
//						    secVertexBpmNoD0.at(j)->phi(),
//						    0.140);		
//
//		  Bpm2 = pion1+pion2+pion3;
//		  double mass = (pion1+pion2+pion3+D0).M();
//		  charge = secVertexBpmNoD0.at(k)->charge()
//		          +secVertexBpmNoD0.at(j)->charge()
//		          +secVertexBpmNoD0.at(n)->charge();
//
//		  if(fabs(mass-5.279)<fabs(massD0-5.279)){ //find the candidate closest to D0 mass
//		    massBpm2 = mass;
////		    printf("BPM mass %f\n",massBpm2);
//		  }
//		}
//	}
//	jet.addUserFloat("sec_massBpm_charge",charge);
//	jet.addUserFloat("sec_massBpm",massBpm2);
	/////////////

    //printf("Gets Here8\n");
	jets->push_back(jet);
    //printf("Gets Here9\n");

      }
    //printf("Gets Here10\n");
    iEvent.put(jets);
    //printf("Gets Here11\n");
  }
  
  // ----------member data ---------------------------
  edm::InputTag src_;
  edm::InputTag leptons_;
  edm::InputTag vertex_;
    
};

