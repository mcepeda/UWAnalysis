/**
* @file PATMETVertexWeighter.cc
* @author T.M.Perry
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "boost/filesystem.hpp"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

/**
* @class PATMETVertexWeighter
* @brief Adds nrVtx, PU, lumiWeight to METs as Userfloat
*/
class PATMETVertexWeighter : public edm::EDProducer
{

    public:
        PATMETVertexWeighter(const edm::ParameterSet& pset);
        virtual ~PATMETVertexWeighter () {}
        void produce( edm::Event& evt, const edm::EventSetup& es );

    private:
        typedef reco::Candidate::LorentzVector FourVec;
        edm::InputTag srcMET_;
        edm::InputTag srcVert_;
        edm::InputTag srcPU_;
};

PATMETVertexWeighter::PATMETVertexWeighter( const edm::ParameterSet& pset)
{
    srcMET_ = pset.getParameter<edm::InputTag>("srcMET");
    srcVert_ = pset.getParameter<edm::InputTag>("srcVert");
    srcPU_   = pset.getParameter<edm::InputTag>("srcPU");
    produces<pat::METCollection>();
}

void PATMETVertexWeighter::produce( edm::Event& evt, const edm::EventSetup& es )
{
    std::auto_ptr<pat::METCollection> out(new pat::METCollection);

    edm::Handle<pat::METCollection> mets;
    evt.getByLabel( srcMET_, mets );

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
    float lumiWeightPU;
    std::string base = std::getenv("CMSSW_BASE");
    std::string fPUMCloc =   "/src/UWAnalysis/Configuration/data/MC_Summer12_PU_S10-600bins.root";
    std::string fPUDATAloc = "/src/UWAnalysis/Configuration/data/Data_Pileup_2012_ReReco-600bins.root";
    std::string fPUMCname =   base+fPUMCloc;
    std::string fPUDATAname = base+fPUDATAloc;
    bool fPUMCis   = boost::filesystem::exists( fPUMCname   );
    bool fPUDATAis = boost::filesystem::exists( fPUDATAname );
    edm::LumiReWeighting *LumiWeights;
    LumiWeights = new edm::LumiReWeighting(fPUMCname,fPUDATAname,"pileup","pileup");
    lumiWeightPU = LumiWeights->weight(nrPU);
    //std::cout<<"PATMETVertexWeighter.cc "<<lumiWeightPU<<std::endl;
    LumiWeights->clearHist();

    pat::MET met = mets->at(0);
    met.addUserFloat("nvtx",nvtx);
    met.addUserFloat("nrPU",nrPU);
    met.addUserFloat("SF_lumiWeightPU",lumiWeightPU); 
    out->push_back(met);

    evt.put( out );
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMETVertexWeighter);
