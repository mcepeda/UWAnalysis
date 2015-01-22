// Fast implementation of Beta and Beta*
// To be revised

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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Math/GenVector/VectorUtil.h"
#include "TLorentzVector.h"

class PATSimplePUID : public edm::EDProducer {
   public:
  
  explicit PATSimplePUID(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src")),
    vertex_(iConfig.getParameter<edm::InputTag>("vertices"))
  {
    produces<pat::JetCollection>();
  }
  
  ~PATSimplePUID() {}
   private:

  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    using namespace edm;
    using namespace reco;
    std::auto_ptr<pat::JetCollection> jets(new pat::JetCollection);
    Handle<pat::JetCollection > cands;
    Handle<reco::VertexCollection> vertices;

    if(!iEvent.getByLabel(vertex_,vertices)) {edm::LogError("")<<"NO PVs!"; return;}
    if(!iEvent.getByLabel(src_,cands)) {edm::LogError("")<<"NO Jets!"; return;}

    for(unsigned int  l=0;l!=cands->size();++l){

	pat::Jet jet = cands->at(l);
       // std::cout<<
       //             " Jet pT: "<<jet.pt()<<
       //             " eta: "<<jet.eta()<<
       //             " phi: "<<jet.phi()<<std::endl;
       //if(jet.hasTagInfo("trackIP"))printf("has Track IP TagInfo\n");

        double nTrack=0;
	double sumPt=0;
	double sumPt2=0;
        double sumTkPt=0;
        double sumTkPt2=0;
        double sumDR=0;
        double sumDR2=0;
	double closestDxy=1e6;
        double closestDz=1e6;
	double betaClassic=0;
        double betaStarClassic=0;
        double beta=0;
        double betaStar=0;


	for(unsigned int j=0;j<jet.getPFConstituents().size();++j) {
          reco::PFCandidatePtr pfcand = jet.getPFConstituents().at(j);
	  
	  double pdgid = pfcand->pdgId();
          bool  hasTrack=(pfcand->trackRef().isNonnull()&&pfcand->trackRef().isAvailable());

	  double dR=ROOT::Math::VectorUtil::DeltaR(jet.p4(),pfcand->p4());
	  double ptdR=pfcand->pt()*dR;

	  sumPt+=pfcand->pt(); 
	  sumPt2+=pfcand->pt()*pfcand->pt();

	  sumDR+=ptdR;
	  sumDR2+=ptdR*ptdR; 

	  if (!hasTrack) continue;
          double pt = pfcand->trackRef()->pt();

	  sumTkPt+=pt;
	  sumTkPt2+=pt*pt;
          nTrack++;

	  double dxy = pfcand->trackRef()->dxy(vertices->at(0).position());
	  if(fabs(dxy)<fabs(closestDxy)) closestDxy=dxy;
	
          double dz = pfcand->trackRef()->dz(vertices->at(0).position());
          if(fabs(dz)<fabs(closestDz)) closestDz=dz;

	  double dz0=dz; 

	  bool inVtx0 = find( vertices->at(0).tracks_begin(), vertices->at(0).tracks_end(), reco::TrackBaseRef(pfcand->trackRef()) ) != vertices->at(0).tracks_end();
 	  bool inAnyOther = false;

	  // 'classic' beta definition based on track-vertex association
	  // alternative beta definition based on track-vertex distance of closest approach
	  for (unsigned int nVtx=1;nVtx<vertices->size();++nVtx) {
		reco::Vertex iv=vertices->at(nVtx);
		if( iv.isFake() || iv.ndof() < 4 ) { continue; }
		// 'classic' beta definition: check if the track is associated with any vertex other than the primary one
		if(!inAnyOther){
			inAnyOther = find( iv.tracks_begin(), iv.tracks_end(), reco::TrackBaseRef(pfcand->trackRef()) ) != iv.tracks_end();
		}
		// alternative beta: find closest vertex to the track
		double dz_v=pfcand->trackRef()->dz(iv.position());
		if( fabs(dz_v)<fabs(dz)) dz=dz_v;
	  }
	  // classic beta/betaStar
	  if( inVtx0 && ! inAnyOther ) {
		betaClassic += pt;
	  } else if( ! inVtx0 && inAnyOther ) {
		betaStarClassic += pt;
	  }
	  // alternative beta/betaStar
	  if( dz0 < 0.2 ) {
		beta += pt;
	  } else if( dz < 0.2 ) {
		betaStar += pt;
	  }

       }

       double ptRMS= sqrt(sumPt2/(sumPt*sumPt));	
       double ptTKRMS= -1; 
		
       double dRMean = sumDR / jet.pt();
       double dR2Mean = sumDR2 / sumPt2;

       if(sumTkPt!=0){
       ptTKRMS=sqrt(sumTkPt2/(sumTkPt*sumTkPt)); 	
       beta = beta / sumTkPt;
       betaStar = betaStar / sumTkPt;
       betaClassic = betaClassic / sumTkPt;
       betaStarClassic = betaStarClassic / sumTkPt;
       }
       else {
	beta=-1; betaStar=-1; betaClassic=-1; betaStarClassic=-1;
       }		

       double betaStarClassicModified = betaStarClassic/log(vertices->size()-0.64);

	jet.addUserFloat("idPtRMS",ptRMS);
        jet.addUserFloat("idPtTKRMS",ptTKRMS);
        jet.addUserFloat("id_dRMean",dRMean);
        jet.addUserFloat("id_dR2Mean",dR2Mean);
        jet.addUserFloat("idBeta",beta);
        jet.addUserFloat("idBetaClassic",betaClassic);
        jet.addUserFloat("idBetaStar",betaStar);
        jet.addUserFloat("idBetaStarClassic",betaStarClassic);
        jet.addUserFloat("idBetaStarClassicMod",betaStarClassicModified);
        jet.addUserFloat("idClosestDxy",closestDxy);
        jet.addUserFloat("idClosestDz",closestDz);
        jet.addUserFloat("id_nTrack",nTrack);


	jets->push_back(jet);
      }


    iEvent.put(jets);
  } 
  
      // ----------member data ---------------------------
  edm::InputTag src_;
  edm::InputTag vertex_;
};


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(PATSimplePUID);
