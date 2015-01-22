/**
* @file PATMuonVertexWeighter.cc
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

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "boost/filesystem.hpp"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

/**
* @class PATMuonVertexWeighter
* @brief Adds nrVtx, PU, lumiWeight to Muons as Userfloat
*/
class PATMuonVertexWeighter : public edm::EDProducer
{

    public:
        PATMuonVertexWeighter(const edm::ParameterSet& pset);
        virtual ~PATMuonVertexWeighter () {}
        void produce( edm::Event& evt, const edm::EventSetup& es );

    private:
        typedef reco::Candidate::LorentzVector FourVec;
        edm::InputTag srcMuon_;
        edm::InputTag srcVert_;
        edm::InputTag srcPU_;
};

PATMuonVertexWeighter::PATMuonVertexWeighter( const edm::ParameterSet& pset)
{
    srcMuon_ = pset.getParameter<edm::InputTag>("srcMuon");
    srcVert_ = pset.getParameter<edm::InputTag>("srcVert");
    srcPU_   = pset.getParameter<edm::InputTag>("srcPU");
    produces<pat::MuonCollection>();
}

void PATMuonVertexWeighter::produce( edm::Event& evt, const edm::EventSetup& es )
{
    std::auto_ptr<pat::MuonCollection> out(new pat::MuonCollection);

    edm::Handle<pat::MuonCollection> muons;
    evt.getByLabel( srcMuon_, muons );

    edm::Handle<reco::VertexCollection> vertices;
    evt.getByLabel(srcVert_, vertices);
    int nvtx = vertices->size();

    float nrPU = 0.;
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    if(evt.getByLabel(srcPU_, PupInfo)) {
      for(std::vector<PileupSummaryInfo>::const_iterator i = PupInfo->begin();
          i!=PupInfo->end();++i) {
            int BX = i->getBunchCrossing();
            if(BX==0){
              nrPU =i->getTrueNumInteractions();
            }
       }
      }
    float preWeight;
    std::string base = std::getenv("CMSSW_BASE");
    std::string fPUMCloc =   "/src/UWAnalysis/Configuration/data/MC_Summer12_PU_S10-600bins.root";
    std::string fPUDATAloc = "/src/UWAnalysis/Configuration/data/Data_Pileup_2012_ReReco-600bins.root";
    std::string fPUMCname =   base+fPUMCloc;
    std::string fPUDATAname = base+fPUDATAloc;
    bool fPUMCis   = boost::filesystem::exists( fPUMCname   );
    bool fPUDATAis = boost::filesystem::exists( fPUDATAname );
    edm::LumiReWeighting *LumiWeights;
    LumiWeights = new edm::LumiReWeighting(fPUMCname,fPUDATAname,"pileup","pileup");
    preWeight = LumiWeights->weight(nrPU);
    //~LumiWeights; 
    LumiWeights->clearHist();

    for ( size_t i = 0; i < muons->size(); ++i )
    {
        pat::Muon muon = muons->at(i);
        muon.addUserFloat("nvtx",nvtx);
        muon.addUserFloat("nrPU",nrPU);
        muon.addUserFloat("preWeight",preWeight); 
        out->push_back(muon);
    }

    evt.put( out );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMuonVertexWeighter);
