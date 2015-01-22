#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "UWAnalysis/DataFormats/interface/BCandFinalState.h"
#include <TLorentzVector.h>

struct Bhadron {
  TLorentzVector        kinematics;
  reco::SecondaryVertex vertex;
  GlobalVector          flightdir;
  int                   isSelected;
};

struct compareBCandbyPt {  
  bool operator()( const Bhadron BCand1, const Bhadron BCand2 ) const {
    return BCand1.kinematics.Pt() > BCand2.kinematics.Pt();
  }
};

class BCProducer : public edm::EDProducer {
   public:
      explicit BCProducer(const edm::ParameterSet&);
      ~BCProducer();

   private:
      virtual void beginJob();
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      GlobalVector flightDirection(reco::Vertex &, reco::Vertex &);
      void addTracks(std::set<reco::TrackRef> &, reco::Vertex &, double);	
      edm::InputTag vertexsrc_,offlinePrimaryVertices_,PFJet_,PFGenJet_;	
      BCandFinalState BCFS;
      compareBCandbyPt sortB;
	bool isMC_;
};

BCProducer::BCProducer(const edm::ParameterSet& iConfig):
  vertexsrc_(iConfig.getParameter<edm::InputTag>("src")),
  offlinePrimaryVertices_(iConfig.getParameter<edm::InputTag>("primaryVertices")),
  PFJet_(iConfig.getParameter<edm::InputTag>("jetSource")),
  PFGenJet_(iConfig.getParameter<edm::InputTag>("jetGenSource")),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
// produces<std::vector<reco::LeafCandidate> >();
produces<BCandFinalStateCollection>("BCandFinalState");
}


BCProducer::~BCProducer()
{}

void BCProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	BCFS.nBC(-777);
	BCFS.BC1PT(-777);
        BCFS.BC2PT(-777);
        BCFS.BC1ETA(-777);
        BCFS.BC2ETA(-777);
        BCFS.BC1PHI(-777);
        BCFS.BC2PHI(-777);
        BCFS.BC1nTracks(-777);
        BCFS.BC2nTracks(-777);
        BCFS.BCEnergySum(-777);
        BCFS.BCGenEnergySum(-777);
        BCFS.nJetAssociatedToBC1(-777);
        BCFS.nJetAssociatedToBC2(-777);
        BCFS.nJetAssociatedToBC(-777);
        BCFS.Jet1BC1PT(-777);
        BCFS.Jet2BC1PT(-777);
        BCFS.Jet3BC1PT(-777);
        BCFS.Jet4BC1PT(-777);
        BCFS.Jet1BC1ETA(-777);
        BCFS.Jet2BC1ETA(-777);
        BCFS.Jet3BC1ETA(-777);
        BCFS.Jet4BC1ETA(-777);
        BCFS.Jet1BC1PHI(-777);
        BCFS.Jet2BC1PHI(-777);
        BCFS.Jet3BC1PHI(-777);
        BCFS.Jet4BC1PHI(-777);
        BCFS.Jet1BC2PT(-777);
        BCFS.Jet2BC2PT(-777);
        BCFS.Jet3BC2PT(-777) ;
        BCFS.Jet4BC2PT(-777);
        BCFS.Jet1BC2ETA(-777);
        BCFS.Jet2BC2ETA(-777);
        BCFS.Jet3BC2ETA(-777);
        BCFS.Jet4BC2ETA(-777);
        BCFS.Jet1BC2PHI(-777);
        BCFS.Jet2BC2PHI(-777);
        BCFS.Jet3BC2PHI(-777);
        BCFS.Jet4BC2PHI(-777);
        BCFS.BCDeltaR(-777);
        BCFS.BDeltaPHI(-777);
	BCFS.BC1MASS(-777);
	BCFS.BC2MASS(-777);
	//std::cout<<"In BCProducer"<<std::endl;
	// std::vector<reco::LeafCandidate> *blc = new std::vector<reco::LeafCandidate>;
	
	std::auto_ptr<BCandFinalStateCollection> BCandFS(new BCandFinalStateCollection);
	//Getting the PV collection
	edm::Handle<reco::VertexCollection> primaryVertexCollection;
	iEvent.getByLabel(offlinePrimaryVertices_, primaryVertexCollection);
	const reco::VertexCollection pvc = *(primaryVertexCollection.product());
	
	//Getting the PFJetCollection
	edm::Handle<edm::View<reco::Jet> > inputJets;
	  iEvent.getByLabel(PFJet_, inputJets);

        edm::Handle<edm::View<reco::GenJet> > inputGenJets;
	
	//Getting the Bcand collection (reco::Vertex)
	edm::Handle<std::vector<reco::Vertex> > SVC;
	iEvent.getByLabel(vertexsrc_, SVC);
	const std::vector<reco::Vertex> svc = *(SVC.product());
	//Define pv from the first PV (ordered by sum of track pt)
	reco::Vertex pv = pvc[0];
	//std::cout<<"PV x="<<pv.x()<<" y="<<pv.y()<<" z="<<pv.z()<<std::endl;
	std::vector<std::set<reco::TrackRef> > bCandtrackSets;
	float charge=0;	
	//std::cout<<"Running on SV, size="<<svc.size()<<std::endl;
	int nSV=0;
	int nJets=0,nGJets=0;
	std::vector<Bhadron> selectedBcands;
	TLorentzVector TLV_jet[100],TLV_Gjet[100];
	for(unsigned int j=0; j<100; j++){
	         TLV_jet[j].SetPxPyPzE(0,0,0,0);
       	}
	for(edm::View<reco::Jet>::const_iterator jet = inputJets->begin(); jet!=inputJets->end(); ++jet){
		if(jet->pt()>10){
                	TLV_jet[nJets].SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),jet->mass());
			nJets++;
		}
        }
	if(isMC_){	
        	for(edm::View<reco::GenJet>::const_iterator Gjet = inputGenJets->begin(); Gjet!=inputGenJets->end(); ++Gjet){
        	        if(Gjet->pt()>2){
        	                TLV_Gjet[nGJets].SetPtEtaPhiM(Gjet->pt(),Gjet->eta(),Gjet->phi(),Gjet->mass());
        	                nGJets++;
        	        }
        	}
	}	
	int index[nJets],Gindex[nGJets];
	for(int j=0; j<nJets; j++){
		index[j]=0;
		Gindex[j]=0;
	}
	float sumEnergy=0,sumGenEnergy=0;
        int nJetAssociated=0,nJetTotalAssociated=0;	

	//Selecting the good BCandidate
	for(unsigned int i=0; i<svc.size(); i++){
		//std::cout<<"Testing the quality of the SV "<<i<<std::endl;		
		reco::Vertex svertex = svc[i];
		GlobalVector flightaxis = flightDirection(pv,svertex);
		reco::SecondaryVertex sv(pv,svertex,flightaxis,true);
		Bhadron oneBcand;
		oneBcand.isSelected = 0;
		if(svc[i].p4().M()>1.4 && sv.dist3d().significance()>5.0 && fabs(flightaxis.eta())<2.1 && svc[i].p4().Pt()>8.0) { 
	  		oneBcand.isSelected = 1; 
		}
		//std::cout<<"    Mass="<<svc[i].p4().M()
		//	<<" sv.dist3d().significance()="<<sv.dist3d().significance()
		//	<<" flightaxis.eta()="<<flightaxis.eta()
		//	<<" svc[i].p4().Pt()="<<svc[i].p4().Pt()
		//	<<" ===> is Selected? "<<oneBcand.isSelected<<std::endl;
		oneBcand.vertex = sv;
		reco::Vertex::Point p(svc[i].x(),svc[i].y(),svc[i].z());
		math::XYZTLorentzVectorD bCandmom;
		bCandmom = svc[i].p4(0.13957,0.5);
                reco::LeafCandidate lcc(charge, bCandmom, p,0);
		TLorentzVector BLorentzV;  BLorentzV.SetPtEtaPhiM(lcc.p4().Pt(),lcc.p4().Eta(),lcc.p4().Phi(),lcc.p4().M());
		oneBcand.kinematics = BLorentzV;
		oneBcand.flightdir  = flightaxis;
		//std::cout<<"	Computing the kinematics:"<<std::endl;
                //std::cout<<"         lc.Eta= "<<lcc.eta()<<"      lcc.Phi="<<lcc.phi()<<std::endl;
                //std::cout<<"    (pv-SV).Eta= "<<flightaxis.eta()<<" (pv-sv).Phi="<<flightaxis.phi()<<std::endl;

		if( oneBcand.isSelected==1 ) selectedBcands.push_back( oneBcand );
	}
	//std::cout<<"After Selection, "<<selectedBcands.size()<<" remains over the "<<svc.size()<<" ones present"<<std::endl;
	sort( selectedBcands.begin(), selectedBcands.end(), sortB );

	BCFS.nBC(selectedBcands.size());
	for(unsigned int i=0; i<selectedBcands.size(); i++){
		float etabc=selectedBcands[i].flightdir.eta();
                float phibc=selectedBcands[i].flightdir.phi();

		std::set<reco::TrackRef> bCandtracks;
		addTracks(bCandtracks, (reco::Vertex&) svc[i], 0.5);
		TLorentzVector SumFourVec;	
		for(std::set<reco::TrackRef>::const_iterator it=bCandtracks.begin(); it!=bCandtracks.end(); it++){
			charge += (*it)->charge();
			TLorentzVector FourVec;
			FourVec.SetPtEtaPhiM((*it)->pt(),(*it)->eta(),(*it)->phi(),0.135);
			SumFourVec+=FourVec;
		}
		nJetAssociated=0;	
		if(nSV==0){
			BCFS.BC1nTracks(selectedBcands[i].vertex.nTracks());
			BCFS.BC1MASS(selectedBcands[i].kinematics.M());
                	BCFS.BC1PT(selectedBcands[i].kinematics.Pt());
                 	BCFS.BC1ETA(selectedBcands[i].flightdir.eta());
                 	BCFS.BC1PHI(selectedBcands[i].flightdir.phi());
			//std::cout<<"	---> First BC, pt="<<selectedBcands[i].kinematics.Pt()
			//	 <<"	eta="<<selectedBcands[i].flightdir.eta()
			//	 <<"	phi="<<selectedBcands[i].flightdir.phi()
			//	 <<"	nTrack="<<bCandtracks.size()
			//	 <<"	3D sig="<<selectedBcands[i].vertex.dist3d().significance()
			//	 <<"	mass="<<selectedBcands[i].kinematics.M()<<std::endl;	
			//std::cout<<"		Running jets to count the total energy associated to the first b-hadron"<<std::endl;
			for(int j=0; j<nJets; j++){
//				std::cout<<"      --- > jet number "<<j<<std::endl;
				float etaj=TLV_jet[j].Eta();
				float phij=TLV_jet[j].Phi();
                               	float dEtaBJ=fabs(etaj-etabc);
				float dPhiBJ=fabs(phij-phibc)?fabs(phij-phibc):6.2831-fabs(phij-phibc);
				float dRBJ=sqrt(pow(dEtaBJ,2)+pow(dPhiBJ,2));
//				std::cout<<"         --- > not used yet pt="<<TLV_jet[j].Pt()<<" eta="<<etaj<<" phi="<<phij<<std::endl;
	                        if(dRBJ<0.5){
					if(index[j]==0){
						//std::cout<<"			jet "<<j<<"    --- > dr<0.5, adding the energy to sumEnergy: "<<sumEnergy
						//	<<"---->"<<sumEnergy+TLV_jet[j].E()
                                                //        <<" , nJetAssociated="<<nJetAssociated
						//	<<" --->"<<nJetAssociated+1<<std::endl;	
						sumEnergy+=TLV_jet[j].E();
						index[j]=1;
						if(nJetAssociated==0){
								BCFS.Jet1BC1PT(TLV_jet[j].Pt());BCFS.Jet1BC1ETA(TLV_jet[j].Eta());BCFS.Jet1BC1ETA(TLV_jet[j].Phi());
						}
                                                if(nJetAssociated==1){
							BCFS.Jet2BC1PT(TLV_jet[j].Pt());BCFS.Jet2BC1ETA(TLV_jet[j].Eta());BCFS.Jet2BC1ETA(TLV_jet[j].Phi());
						}
                                                if(nJetAssociated==2){
								BCFS.Jet3BC1PT(TLV_jet[j].Pt());BCFS.Jet3BC1ETA(TLV_jet[j].Eta());BCFS.Jet3BC1ETA(TLV_jet[j].Phi());
						}
                                                if(nJetAssociated==3){
							BCFS.Jet4BC1PT(TLV_jet[j].Pt());BCFS.Jet4BC1ETA(TLV_jet[j].Eta());BCFS.Jet4BC1ETA(TLV_jet[j].Phi());
						}
						nJetTotalAssociated++;
					}
					//else{std::cout<<"                    jet "<<j<<" already considered"<<std::endl;}	
					nJetAssociated++;
					//std::cout<<"          =============> nJetAssociated="<<nJetAssociated<<" nJetTotalAssociated="<<nJetTotalAssociated<<std::endl;

				}
                        }
			if(isMC_){
				for(int j=0; j<nGJets; j++){
  //                    	          std::cout<<"      --- >Gen  jet number "<<j<<std::endl;
                        	        if(Gindex[j]==0){
                        	                float etaj=TLV_Gjet[j].Eta();
                        	                float phij=TLV_Gjet[j].Phi();
                        	                float dEtaBJ=fabs(etaj-etabc);
                        	                float dPhiBJ=fabs(phij-phibc)?fabs(phij-phibc):6.2831-fabs(phij-phibc);
                        	                float dRBJ=sqrt(pow(dEtaBJ,2)+pow(dPhiBJ,2));
    //                  	                  std::cout<<"         --- > GenJet not used yet pt="<<TLV_Gjet[j].Pt()<<" eta="<<TLV_Gjet[j].Eta()<<" phi="<<TLV_Gjet[j].Phi()<<std::endl;
                        	                 if(dRBJ<0.5){
                        	                        //std::cout<<"			genjet "<<j<<"     --- > dr<0.5, adding the energy to sumGenEnergy: "<<sumGenEnergy
							//	<<"---->"<<sumGenEnergy+TLV_Gjet[j].E()<<std::endl;
                        	                        sumGenEnergy+=TLV_Gjet[j].E();
                        	                        Gindex[j]=1;
                        	                }
                        	        }
                        	}
			}
			BCFS.nJetAssociatedToBC1(nJetAssociated);
		}
		nJetAssociated=0;
		if(nSV==1){
			BCFS.BC2nTracks(selectedBcands[i].vertex.nTracks());
                        BCFS.BC2MASS(selectedBcands[i].kinematics.M());
                        BCFS.BC2PT(selectedBcands[i].kinematics.Pt());
                        BCFS.BC2ETA(selectedBcands[i].flightdir.eta());
                        BCFS.BC2PHI(selectedBcands[i].flightdir.phi());
                        //std::cout<<"	---> Second BC, pt="<<selectedBcands[i].kinematics.Pt()
                        //        <<" eta="<<selectedBcands[i].flightdir.eta()
                        //        <<" phi="<<selectedBcands[i].flightdir.phi()
                        //        <<" nTrack="<<bCandtracks.size()
                        //        <<" 3D sig="<<selectedBcands[i].vertex.dist3d().significance()
                        //        <<" mass="<<selectedBcands[i].kinematics.M()<<std::endl;    				
	
			float dphi=fabs(BCFS.BC1PHI()-BCFS.BC2PHI())<3.14159?fabs(BCFS.BC1PHI()-BCFS.BC2PHI()):6.28318-fabs(BCFS.BC1PHI()-BCFS.BC2PHI());
			BCFS.BDeltaPHI(dphi);
                	BCFS.BCDeltaR(sqrt(pow(BCFS.BC1ETA()-BCFS.BC2ETA(),2)+pow(BCFS.BDeltaPHI(),2)));
			//std::cout<<"			DeltaR="<<sqrt(pow(BCFS.BC1ETA()-BCFS.BC2ETA(),2)+pow(BCFS.BDeltaPHI(),2))<<std::endl; 
			//std::cout<<"			Running jets to count the total energy associated to the second b-hadron"<<std::endl;
			for(int j=0; j<nJets; j++){
//				std::cout<<"      --- > jet number "<<j<<std::endl,
                               float etaj=TLV_jet[j].Eta();
	                       float phij=TLV_jet[j].Phi();
                               float dEtaBJ=fabs(etaj-etabc);
                               float dPhiBJ=fabs(phij-phibc)?fabs(phij-phibc):6.2831-fabs(phij-phibc);
                               float dRBJ=sqrt(pow(dEtaBJ,2)+pow(dPhiBJ,2));
//					std::cout<<"         --- > not used yet pt="<<TLV_jet[j].Pt()<<" eta="<<TLV_jet[j].Eta()<<" phi="<<TLV_jet[j].Phi()<<std::endl;
                	  	if(dRBJ<0.5){
					if(index[j]==0){
						//std::cout<<"			jet "<<j<<"     --- > dr<0.5, adding the energy to sumEnergy: "<<sumEnergy
						//	<<"---->"<<sumEnergy+TLV_jet[j].E()
                                                //        <<" , nJetAssociated="<<nJetAssociated
						//	<<" --->"<<nJetAssociated+1<<std::endl;
						sumEnergy+=TLV_jet[j].E();
						index[j]=1;
						if(nJetAssociated==0){
                                                                BCFS.Jet1BC2PT(TLV_jet[j].Pt());BCFS.Jet1BC2ETA(TLV_jet[j].Eta());BCFS.Jet1BC2ETA(TLV_jet[j].Phi());
                                                }
                                                if(nJetAssociated==1){
                                                        BCFS.Jet2BC2PT(TLV_jet[j].Pt());BCFS.Jet2BC2ETA(TLV_jet[j].Eta());BCFS.Jet2BC2ETA(TLV_jet[j].Phi());
                                                }
                                                if(nJetAssociated==2){
                                                                BCFS.Jet3BC2PT(TLV_jet[j].Pt());BCFS.Jet3BC2ETA(TLV_jet[j].Eta());BCFS.Jet3BC2ETA(TLV_jet[j].Phi());
                                                }
                                                if(nJetAssociated==3){
                                                        BCFS.Jet4BC2PT(TLV_jet[j].Pt());BCFS.Jet4BC2ETA(TLV_jet[j].Eta());BCFS.Jet4BC2ETA(TLV_jet[j].Phi());
                                                }
						nJetTotalAssociated++;
					}
					//else{std::cout<<"                    jet  "<<j<<" already considered"<<std::endl;}
					nJetAssociated++;
					//std::cout<<"          =============> nJetAssociated="<<nJetAssociated<<" nJetTotalAssociated="<<nJetTotalAssociated<<std::endl;
				}
			}
			if(isMC_){
				for(int j=0; j<nGJets; j++){
                //      	          std::cout<<"      --- > Gen jet number "<<j<<std::endl;
                        	        if(Gindex[j]==0){
                        	                float etaj=TLV_Gjet[j].Eta();
                        	                float phij=TLV_Gjet[j].Phi();
                        	                float dEtaBJ=fabs(etaj-etabc);
                        	                float dPhiBJ=fabs(phij-phibc)?fabs(phij-phibc):6.2831-fabs(phij-phibc);
                        	                float dRBJ=sqrt(pow(dEtaBJ,2)+pow(dPhiBJ,2));
                  //    	                  std::cout<<"         --- > GenJet not used yet pt="<<TLV_Gjet[j].Pt()<<" eta="<<TLV_Gjet[j].Eta()<<" phi="<<TLV_Gjet[j].Phi()<<std::endl;
                        	                 if(dRBJ<0.5){
        	                                        sumGenEnergy+=TLV_Gjet[j].E();
        	                                        Gindex[j]=1;
        	                                }
        	                        }
        	                }
			}
			BCFS.nJetAssociatedToBC2(nJetAssociated);
                }
		BCFS.nJetAssociatedToBC(nJetTotalAssociated);
		BCFS.BCEnergySum(sumEnergy);
                BCFS.BCGenEnergySum(sumGenEnergy);
		nSV++;	
	}
        BCandFS->push_back(BCFS);
	//std::cout<<"--->"<<sqrt(pow(BCFS.BC1ETA()-BCFS.BC2ETA(),2)+pow(BCFS.BDeltaPHI(),2))<<std::endl;	
      	std::auto_ptr<BCandFinalStateCollection > bfsColl(BCandFS);
	iEvent.put(bfsColl,"BCandFinalState");	
	
}	

void BCProducer::addTracks(std::set<reco::TrackRef> &trackrefs, reco::Vertex &v, double weight){

  for(reco::Vertex::trackRef_iterator ti = v.tracks_begin(); ti!=v.tracks_end(); ti++){
    	if(v.trackWeight(*ti)>weight){
      		reco::TrackRef t = ti->castTo<reco::TrackRef>();
      		trackrefs.insert(t);
    	}

  }
}

GlobalVector BCProducer::flightDirection(reco::Vertex &pv, reco::Vertex &sv){
  GlobalVector res(sv.position().X() - pv.position().X(),sv.position().Y() - pv.position().Y(),sv.position().Z() - pv.position().Z());
  return res;
}

// ------------ method called once each job just before starting event loop ------------
void
BCProducer::beginJob()
{
 

}

// ------------ method called once each job just after ending the event loop ------------
void
BCProducer::endJob() {

}


DEFINE_FWK_MODULE(BCProducer);
