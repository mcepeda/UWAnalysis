/**
* @file PATMuonIDer.cc
* @author T.M.Perry
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

/**
* @class PATMuonIDer
* @brief Produces a collection of muons with ID information saved as userFloats
*/
class PATMuonIDer : public edm::EDProducer
{

    public:
        PATMuonIDer(const edm::ParameterSet& pset);
        virtual ~PATMuonIDer () {}
        void produce( edm::Event& evt, const edm::EventSetup& es );

    private:
        typedef reco::Candidate::LorentzVector FourVec;
        edm::InputTag _src;
};



PATMuonIDer::PATMuonIDer( const edm::ParameterSet& pset)
{
    _src = pset.getParameter<edm::InputTag>("src");
    produces<pat::MuonCollection>();
}



/**
* Create a new collection of pat::Muons with looseID and
* push the collection to the event.
*/
void PATMuonIDer::produce( edm::Event& evt, const edm::EventSetup& es )
{
    std::auto_ptr<pat::MuonCollection> out(new pat::MuonCollection);

    edm::Handle<pat::MuonCollection> muons;
    evt.getByLabel( _src, muons );

    for ( size_t i = 0; i < muons->size(); ++i )
    {
       double looseID = 0.;
       pat::Muon muon = muons->at(i);
//       std::cout<<"isPF: "<<muon.isPFMuon()<<std::endl;
//       std::cout<<"isGlobal: "<<muon.isGlobalMuon()<<std::endl;
//       std::cout<<"isTracker: "<<muon.isTrackerMuon()<<std::endl;
       if(muon.isPFMuon()){
        if(muon.isGlobalMuon() || muon.isTrackerMuon()){
          looseID = 1.;
        }
       }
        
//       double tightIDtom = 0.;
//       if(muon.isGlobalMuon() && muon.isPFMuon()){
//        if(muon.globalTrack.normalizedChi2() < 10){
//         if(muon.track.hitPattern.trackerLayersWithMeasurement() > 5){
//          if(muon.globalTrack.hitPattern.numberOfValidMuonHits() > 0){
//           if(muon.innerTrack.hitPattern.numberOfValidPixelHits() > 0){
//            if(abs(muon.dB()) < 0.2){
//             if(muon.numberOfMatchedStations() > 1){
//              tightID = 1.;
//             }
//            }
//           }
//          }
//         }
//        }
//        muon.addUserFloat("tightIDtom",tightIDtom);
        muon.addUserFloat("looseID",looseID);
        out->push_back(muon);
    }

    evt.put( out );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMuonIDer);
