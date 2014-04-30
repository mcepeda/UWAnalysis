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
* Create a new collection of pat::Muons with the Rochester-corrected p4 and
* push the collection to the event.
*/
void PATMuonIDer::produce( edm::Event& evt, const edm::EventSetup& es )
{
    std::auto_ptr<pat::MuonCollection> out(new pat::MuonCollection);

    edm::Handle<pat::MuonCollection> muons;
    evt.getByLabel( _src, muons );

    for ( size_t i = 0; i < muons->size(); ++i )
    {
        pat::Muon muon = muons->at(i);
        
    }

    evt.put( out );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMuonIDer);