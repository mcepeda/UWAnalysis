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
// Original Author:  Michail Bachtis
//         Created:  Sun Jan 31 15:04:57 CST 2010
// $Id: GenJetCleaner.h,v 1.1 2011/06/18 00:41:29 bachtis Exp $
//
//
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

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/getRef.h" 

#include "UWAnalysis/DataFormats/interface/DressedLepton.h"




class GenJetCleaner : public edm::EDProducer {
   public:
  
  

  explicit GenJetCleaner(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src")),
  srcDressedLepton_(iConfig.getParameter<edm::InputTag>("srcDressedLepton")),
  onlyGoodLeptons_(iConfig.getParameter<bool>("onlyGoodLeptons")),
  excludeTaus_(iConfig.getUntrackedParameter<bool>("excludeTaus",true))
  {
    produces<pat::JetCollection>();
  }
  
  ~GenJetCleaner() {}
   private:

  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    using namespace edm;
    using namespace reco;
    using namespace std;
    std::auto_ptr<pat::JetCollection> jets(new pat::JetCollection);
    Handle<pat::JetCollection > cands;
    edm::Handle <DressedLeptonCollection> Leptons;
    if(!iEvent.getByLabel(srcDressedLepton_, Leptons)) edm::LogError("")<<"No Dressed Leptons!";
    if(!iEvent.getByLabel(src_,cands)) edm::LogError("")<<"No Jets!";

      for(unsigned int  i=0;i!=cands->size();++i){
	pat::Jet jet = cands->at(i);

          bool thisisnotajet=false;
          for( size_t i = 0; i < Leptons->size(); ++ i ) {
                DressedLepton  genpart = Leptons->at(i);
		if(genpart.status()!=1) continue;
                if(onlyGoodLeptons_ && genpart.Type()!=1) continue;
                if(excludeTaus_ && genpart.Type()==3) continue;
                double dR=deltaR(jet.eta(), jet.phi(), genpart.eta(), genpart.phi());
                if(dR<0.5)  thisisnotajet=true;
          }
          if(thisisnotajet) continue;

	jets->push_back(jet);

      }  
    
    iEvent.put(jets);

  }

      // ----------member data ---------------------------
      edm::InputTag src_;
      edm::InputTag srcDressedLepton_;
      bool onlyGoodLeptons_;
      bool excludeTaus_;

};

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(GenJetCleaner);
