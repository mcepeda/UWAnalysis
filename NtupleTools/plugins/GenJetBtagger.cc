// -*- C++ -*-
//
// Package: PATMuonTrackVetoSelector
// Class: PATMuonTrackVetoSelector
//
/**\class PATMuonTrackVetoSelector PATMuonTrackVetoSelector.cc UWAnalysis/PATMuonTrackVetoSelector/src/PATMuonTrackVetoSelector.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author: Michail Bachtis
// Created: Sun Jan 31 15:04:57 CST 2010
// $Id: GENJetBtagger.h,v 1.1 2011/06/18 00:41:29 bachtis Exp $
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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "Math/GenVector/VectorUtil.h"
//
// class decleration
#include <boost/foreach.hpp>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/getRef.h"
//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "ZSV/BAnalysis/interface/SimBHadron.h"


class GENJetBtagger : public edm::EDProducer {
 public:

 explicit GENJetBtagger(const edm::ParameterSet& iConfig):
 src_(iConfig.getParameter<edm::InputTag>("src")),
 srcSIMB_(iConfig.getParameter<edm::InputTag>("srcSIMB"))
 {
 produces<pat::JetCollection>();
 }

 ~GENJetBtagger() {}
 private:

 const reco::Candidate* getBAncestors ( const reco::Candidate* gp) {
 if (gp->numberOfMothers()!=1 || gp->status()==3) return NULL;
 int pid = gp->pdgId();
 if ((abs(pid)/100)%10 == 5 || (abs(pid)/1000)%10 == 5) {return gp;}
 const reco::Candidate* p = gp->mother();
 const reco::Candidate* mom = getBAncestors(p);
 if (mom != NULL) { return mom;}
 return NULL;
 }



 virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
 {
 using namespace edm;
 using namespace reco;
 using namespace std;
 std::auto_ptr<pat::JetCollection> jets(new pat::JetCollection);
 Handle<reco::GenJetCollection > cands;
 edm::Handle <SimBHadronCollection> BHadrons;


 Handle<GenParticleCollection> GPCollection;
 iEvent.getByLabel("genParticles", GPCollection);
 const GenParticleCollection & pGenPart = *(GPCollection.product());

 if(!iEvent.getByLabel(srcSIMB_, BHadrons)) {std::cout<<"You forgot the bhadrons!!!"<<std::endl;}



 if(iEvent.getByLabel(src_,cands))
 for(unsigned int i=0;i!=cands->size();++i){
 pat::Jet jet = cands->at(i);
 reco::GenJet genjet=cands->at(i);

 float dRBHadron=1e6;
 std::vector<double> BHadrons_PT;
 std::vector<double> BHadrons_ETA;
 std::vector<double> BHadrons_PHI;
 double matchedBHadron=0;
 double secondJetBHadronPt=-10;
 double leadJetBHadronPt=-10;
 double secondJetBHadronPhi=-10;
 double leadJetBHadronPhi=-10;
 double secondJetBHadronEta=-10;
 double leadJetBHadronEta=-10;
 double secondJetBHadronPdgId=-10;
 double leadJetBHadronPdgId=-10;
 double DRBHadrons=-10;


 for( size_t i = 0; i < BHadrons->size(); ++ i ) {
 SimBHadron genpart = BHadrons->at(i);

 double dR=deltaR(jet.eta(), jet.phi(), genpart.eta(), genpart.phi());
 if(dR<0.5) {
 matchedBHadron++;
 BHadrons_PT.push_back(genpart.pt());
 BHadrons_ETA.push_back(genpart.eta());
 BHadrons_PHI.push_back(genpart.phi());
 if(genpart.pt()>leadJetBHadronPt) {secondJetBHadronPt=leadJetBHadronPt; leadJetBHadronPt=genpart.pt();
 secondJetBHadronPhi=leadJetBHadronPhi; leadJetBHadronPhi=genpart.phi();
 secondJetBHadronEta=leadJetBHadronEta; leadJetBHadronEta=genpart.eta();
 secondJetBHadronPdgId=leadJetBHadronPdgId; leadJetBHadronPdgId=genpart.pdgId();
 }
 else if(genpart.pt()>secondJetBHadronPt) {
 secondJetBHadronEta=genpart.eta(); secondJetBHadronPhi=genpart.phi(); secondJetBHadronPt=genpart.pt();
 secondJetBHadronPdgId=genpart.pdgId();
 }
 }
 if(dR<dRBHadron) dRBHadron=dR;
 }

 if(leadJetBHadronPt>0 && secondJetBHadronPt>0) DRBHadrons=deltaR(leadJetBHadronEta, leadJetBHadronPhi,secondJetBHadronEta,secondJetBHadronPhi);

 jet.addUserFloat("matchedHadron",matchedBHadron);
 jet.addUserFloat("minDRBHadron",dRBHadron);
 jet.addUserFloat("leadJetBHadronPt",leadJetBHadronPt);
 jet.addUserFloat("secondJetBHadronPt",secondJetBHadronPt);
 jet.addUserFloat("leadJetBHadronEta",leadJetBHadronEta);
 jet.addUserFloat("secondJetBHadronEta",secondJetBHadronEta);
 jet.addUserFloat("leadJetBHadronPhi",leadJetBHadronPhi);
 jet.addUserFloat("secondJetBHadronPhi",secondJetBHadronPhi);
 jet.addUserFloat("leadJetBHadronPdgId",leadJetBHadronPdgId);
 jet.addUserFloat("secondJetBHadronPdgId",secondJetBHadronPdgId);
 jet.addUserFloat("DR2BHadrons",DRBHadrons);


 int originalBAncestor=-10;
 int originalBAncestorAndrea=-10;

 if(matchedBHadron>0) {
 //if(jet.pt()>25) std::cout<<jet.pt()<<" "<<jet.eta()<<" "<<jet.phi()<<" "<<jet.userFloat("matchedHadron")<<jet.userFloat("leadJetBHadronPt")<<endl;
 bool foundB=false;
 bool foundB2=false;

 for( size_t i = 0; i < genjet.getGenConstituents().size() &&!foundB; i++){
 const GenParticle* genpart=genjet.getGenConstituent(i);
 if(abs(genpart->pdgId()/1000)==5 || abs(genpart->pdgId()/100)==5) { foundB=true; originalBAncestor=genpart->pdgId();}

 for( size_t j = 0; j < pGenPart.size(); j++ ) {
 const GenParticle& genpart2=pGenPart.at(j); // For whatever reason my constituents did not store the mother info
 // so I go back to finding that particle in the genParticles list. Plus, pointer fun.
 // This is bad code...
 if(genpart2.pdgId()!=genpart->pdgId()) continue;
 if(genpart2.pt()!=genpart->pt()) continue;
 if(genpart2.eta()!=genpart->eta()) continue;
 if(genpart2.phi()!=genpart->phi()) continue;
 if(genpart2.status()!=genpart->status()) continue;

 if(!foundB2){
 // A la Andrea
 const Candidate* oldie=getBAncestors((reco::Candidate*)genpart2.mother());
 if(oldie!=NULL){
 int bmom=oldie->pdgId();
 if(abs(bmom/100)%10==5 || abs(bmom/1000)%10==5) {originalBAncestorAndrea=bmom;}
 //if(jet.pt()>25) std::cout<<"\t"<<oldie->pt()<<" "<<oldie->eta()<<" "<<oldie->phi()<<" "<<bmom<<std::endl;
 }
 }

 if(!foundB&&genpart2.numberOfMothers()==1){
 reco::Candidate *mother=(reco::Candidate*)genpart2.mother();

 // My loop
 int momPDGID=mother->pdgId();

 if( abs(momPDGID/100)%10==5 || abs(momPDGID/1000)%10==5) {foundB=true; originalBAncestor=momPDGID; }


 while(!foundB && mother->status()!=3 && mother->numberOfMothers()==1){
 reco::Candidate* motherB=(reco::Candidate*)mother->mother();
 mother=motherB;
 momPDGID=mother->pdgId();
 if( abs(momPDGID/100)%10==5 || abs(momPDGID/1000)%10==5) {foundB=true; originalBAncestor=momPDGID; }
 }

 }

 }

 }

 }

 jet.addUserFloat("originalBAncestor", originalBAncestor);
 jet.addUserFloat("originalBAncestorAndrea", originalBAncestorAndrea);



 jets->push_back(jet);

 }

 iEvent.put(jets);

 }

 // ----------member data ---------------------------
 edm::InputTag src_;
 edm::InputTag srcSIMB_;

};

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(GENJetBtagger);
