/**
* @file PATeleSFer.cc
* @author T.M.Perry
*/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "UWAnalysis/RecoTools/interface/table.h"
/**
* @class PATeleSFer
* @brief embeds trigger/id scalefactors as userFloats
*/
class PATeleSFer : public edm::EDProducer
{

    public:
        PATeleSFer(const edm::ParameterSet& pset);
        virtual ~PATeleSFer () {}
        void produce( edm::Event& iEvent, const edm::EventSetup& es );

    private:
        typedef reco::Candidate::LorentzVector FourVec;
        edm::InputTag src_;
};



PATeleSFer::PATeleSFer( const edm::ParameterSet& pset)
{
    src_ = pset.getParameter<edm::InputTag>("src");
    produces<pat::ElectronCollection>();
}

void PATeleSFer::produce( edm::Event& iEvent, const edm::EventSetup& es )
{
    table* SFtable_IDIso;
    table* SFtable_HLT;
    std::string base = std::getenv("CMSSW_BASE");
    SFtable_IDIso = new table(base+"/src/UWAnalysis/RecoTools/data/SF_ele_IDIso.txt");
    SFtable_HLT = new table(base+"/src/UWAnalysis/RecoTools/data/SF_ele_HLT.txt");

    std::auto_ptr<pat::ElectronCollection> out( new pat::ElectronCollection);
    edm::Handle<pat::ElectronCollection> electrons;
    if(iEvent.getByLabel(src_,electrons))
      for(unsigned int i=0;i<electrons->size();++i) {
        float SF_IDIso = 1.;
        float SF_IDIso_errUp = 0.;
        float SF_IDIso_errDn = 0.;
        float SF_HLT = 1.;
        float SF_HLT_errUp = 0.;
        float SF_HLT_errDn = 0.;
        float SF_IDIsoHLT = 1.;
        float SF_IDIsoHLT_errUp = 0.;
        float SF_IDIsoHLT_errDn = 0.;
        pat::Electron electron = electrons->at(i);

        float electron_pt  = electron.pt();
        float electron_eta = electron.eta();

        SF_IDIso = SFtable_IDIso->Val(electron_pt,electron_eta);
        SF_IDIso_errUp = SFtable_IDIso->ErrUp(electron_pt,electron_eta);
        SF_IDIso_errDn = SFtable_IDIso->ErrDn(electron_pt,electron_eta);

        SF_HLT = SFtable_HLT->Val(electron_pt,electron_eta);
        SF_HLT_errUp = SFtable_HLT->ErrUp(electron_pt,electron_eta);
        SF_HLT_errDn = SFtable_HLT->ErrDn(electron_pt,electron_eta);

        SF_IDIsoHLT = SF_IDIso * SF_HLT;
        SF_IDIsoHLT_errUp = sqrt( pow(SF_IDIso_errUp,2) + pow(SF_HLT_errUp,2) );
        SF_IDIsoHLT_errDn = sqrt( pow(SF_IDIso_errDn,2) + pow(SF_HLT_errDn,2) );

        //std::cout<<"RecoTools/plugins/PATeleSFer.cc "
        //<<electron_pt<<" "<<electron_eta<<"\n "<<
        //"IDIso: "<<SF_IDIso<<"+"<<SF_IDIso_errUp<<"-"<<SF_IDIso_errDn<<"\n "<<
        //"HLT:   "<<SF_HLT<<"+"<<SF_HLT_errUp<<"-"<<SF_HLT_errDn<<std::endl;

        electron.addUserFloat("SF_IDIso",SF_IDIso);
        electron.addUserFloat("SF_IDIso_errUp",SF_IDIso_errUp);
        electron.addUserFloat("SF_IDIso_errDn",SF_IDIso_errDn);
        //electron.addUserFloat("SF_HLT",SF_HLT);
        //electron.addUserFloat("SF_HLT_errUp",SF_HLT_errUp);
        //electron.addUserFloat("SF_HLT_errDn",SF_HLT_errDn);
        electron.addUserFloat("SF_IDIsoHLT",SF_IDIsoHLT);
        electron.addUserFloat("SF_IDIsoHLT_errUp",SF_IDIsoHLT_errUp);
        electron.addUserFloat("SF_IDIsoHLT_errDn",SF_IDIsoHLT_errDn);
        out->push_back(electron);
    }

    iEvent.put( out );
    delete SFtable_IDIso;
    delete SFtable_HLT;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATeleSFer);
