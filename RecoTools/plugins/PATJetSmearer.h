// system include files
#include <memory>
#include <algorithm>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

class PATJetSmearer : public edm::EDProducer{
 public:

  explicit PATJetSmearer(const edm::ParameterSet& iConfig):
   src_(iConfig.getParameter<edm::InputTag>("src")),
   srcGenJet_(iConfig.getParameter<edm::InputTag>("genJet"))
  {
   produces<pat::JetCollection>();
  }
   ~PATJetSmearer() {}


 private:
   double genRecoMatch(pat::Jet recoJet, edm::Handle<reco::GenJetCollection> genJets, double& gen_pt, double& gen_eta, double& gen_phi){
    const reco::GenJet* retVal = 0;
//    double gen_pt = -10;
//    double gen_eta = -99;
//    double gen_phi = -99;
    //double gen_kin[3] = {0,0,0};
    double dRbestMatch = 0.5;
    for ( reco::GenJetCollection::const_iterator genJet = genJets->begin();
      genJet != genJets->end(); ++genJet ) {
        double dR = deltaR(recoJet.p4(), genJet->p4());
        if ( dR < dRbestMatch ) {
        retVal = &(*genJet) ;
        dRbestMatch = dR;
      }
    }
 if(retVal!=0){
  gen_pt = (retVal->pt());
  gen_eta = (retVal->eta());
  gen_phi = (retVal->phi());
 }
 //double [3] gen_kin = {};
 //gen_kin[0] = gen_pt;
 //gen_kin[1] = gen_eta;
 //gen_kin[2] = gen_phi;
 return gen_pt;
}

//    math::PtEtaPhiMLorentzVector genNoNu_p4;
//    if(retVal!=0){
//     genNoNu_p4->setPtEtaPhiM(retVal->pt(), retVal->eta(), retVal->phi(), retVal->mass());
//    }
//    return genNoNu_p4;
//   }

// double genRecoMatch(pat::Jet recoJet, edm::Handle<reco::GenJetCollection> genJets){
//const reco::GenJet* retVal = 0;
//double PT=-20;	
//double dRbestMatch = 0.5;
//for ( reco::GenJetCollection::const_iterator genJet = genJets->begin();
//genJet != genJets->end(); ++genJet ) {
//double dR = deltaR(recoJet.p4(), genJet->p4());
//if ( dR < dRbestMatch ) {
//retVal = &(*genJet) ;
//dRbestMatch = dR;
//}
//}
//if(retVal!=0) PT=retVal->pt();	
//return PT;
//}


  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
   using namespace edm;
   using namespace reco;

   std::auto_ptr<pat::JetCollection> smeared(new pat::JetCollection);
   edm::Handle<pat::JetCollection> recoJets;
   edm::Handle<reco::GenJetCollection> genJets;
   iEvent.getByLabel(srcGenJet_, genJets);


   double etaRange[8] = {0.0, 0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
   double scale[7]   = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
   double scaleUp[7] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
   double scaleDn[7] = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
   double sf = 1.;
   double sfUp = 1.;
   double sfDown = 1.;

   if(iEvent.getByLabel(src_,recoJets)){
    for(unsigned int i=0;i!=recoJets->size();++i){
     pat::Jet jet = recoJets->at(i);

     double reco_pt = jet.pt();

     double reco_eta = fabs( jet.eta() );
     // matching using pointer back (may not be good for PU jets)
     double gen_pt_Nu = -10;
     double gen_eta_Nu = -99;
     double gen_phi_Nu = -99;
     bool matchedGenJet = jet.genJet();
     if(matchedGenJet) {
      gen_pt_Nu = jet.genJet()->pt();
      gen_eta_Nu = jet.genJet()->eta();
      gen_phi_Nu = jet.genJet()->phi();
     }

     //alternative matching by looping over genJets collections and going for dR
     double gen_pt_NoNu = -10;
     double gen_eta_NoNu = -99;
     double gen_phi_NoNu = -99;
     genRecoMatch(jet,genJets,gen_pt_NoNu,gen_eta_NoNu,gen_phi_NoNu);
     //std::cout<<"jet.genJet "<<gen_pt_Nu<<" "<<gen_pt_NoNu<<" NoNu"<<std::endl;

     double ptSmeared = reco_pt;
     double ptSmearedUp = reco_pt;
     double ptSmearedDown = reco_pt;
     
     for(unsigned int j=0;j<=6;++j){
      if(reco_eta >= etaRange[j] && reco_eta < etaRange[j+1]){
       sf = scale[j];
       sfUp = scaleUp[j];
       sfDown = scaleDn[j];
      }
     }

     double gen_pt = gen_pt_NoNu;
     double gen_eta = gen_eta_NoNu;
     double gen_phi = gen_phi_NoNu;
     //double gen_pt = gen_pt_Nu;
     //double gen_eta = gen_eta_Nu;
     //double gen_phi = gen_phi_Nu;
     if(gen_pt > 0){
      ptSmeared     = std::max(0., gen_pt + sf     * (reco_pt - gen_pt) );
      ptSmearedUp   = std::max(0., gen_pt + sfUp   * (reco_pt - gen_pt) );
      ptSmearedDown = std::max(0., gen_pt + sfDown * (reco_pt - gen_pt) );
     } 
     else{
      ptSmeared     = std::max(0., reco_pt );
      ptSmearedUp   = std::max(0., reco_pt );
      ptSmearedDown = std::max(0., reco_pt );
     } 

    // std::cout<<"reco pT: "<<reco_pt<<" smeared pT: "<<ptSmeared<<" pt_gen: "<<gen_pt<<std::endl;
     math::PtEtaPhiMLorentzVector p4(ptSmeared,jet.eta(),jet.phi(),jet.mass());
     jet.setP4(p4);
     jet.addUserFloat("pt_gen_Nu",gen_pt);
     jet.addUserFloat("eta_gen_Nu",gen_eta);
     jet.addUserFloat("phi_gen_Nu",gen_phi);
     jet.addUserFloat("pt_gen_NoNu",gen_pt_NoNu);
     jet.addUserFloat("eta_gen_NoNu",gen_eta_NoNu);
     jet.addUserFloat("phi_gen_NoNu",gen_phi_NoNu);
     jet.addUserFloat("pt_uncorr",reco_pt);
     jet.addUserFloat("pt_smearedUp",ptSmearedUp);
     jet.addUserFloat("pt_smearedDown",ptSmearedDown);
     smeared->push_back(jet);
    }
   iEvent.put(smeared);
   }
  }
  // ----- member data -----
  edm::InputTag src_;
  edm::InputTag srcGenJet_;

};
