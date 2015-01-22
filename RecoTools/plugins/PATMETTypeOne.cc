#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// system include files
#include <memory>
#include <algorithm>
#include <TLorentzVector.h>

// user include files

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

class PATMETTypeOne : public edm::EDProducer {
    public:

        PATMETTypeOne(const edm::ParameterSet& pset);
        virtual ~PATMETTypeOne(){}
        void produce(edm::Event& iEvent, const edm::EventSetup& es);
        pat::Jet jetMatchL1(pat::Jet jetL123, edm::Handle<pat::JetCollection> jetsL1);
	math::PtEtaPhiMLorentzVector findCorrectionType1(
          edm::Handle<pat::JetCollection> l123Jets, bool isMC);
        math::PtEtaPhiMLorentzVector findCorrectionXY(int nrVtx, bool isMC);
        math::PtEtaPhiMLorentzVector antiCorrectionType1(
          edm::Handle<pat::JetCollection> l123Jets);
    private:
	edm::InputTag srcMET_;
        edm::InputTag srcJ123_;
        edm::InputTag srcVert_;
        bool mc_;
};

PATMETTypeOne::PATMETTypeOne(
        const edm::ParameterSet& pset):
     srcMET_(pset.getParameter<edm::InputTag>("srcMET")),
     srcJ123_(pset.getParameter<edm::InputTag>("srcJ123")),
     srcVert_(pset.getParameter<edm::InputTag>("srcVert")),
     mc_(pset.getParameter<bool>("mc"))
     {
      produces<pat::METCollection>();
     }

void PATMETTypeOne::produce(edm::Event& iEvent, const edm::EventSetup& es) {

    using namespace edm;
    using namespace reco;

    // get the collections and define the new output MET
    std::auto_ptr<pat::METCollection> outMET(new pat::METCollection);

    edm::Handle<pat::JetCollection> jetsL123;
    iEvent.getByLabel(srcJ123_,jetsL123);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(srcVert_, vertices);
    int Nvtx = vertices->size();
    //std::cout<<"NVTX: "<<Nvtx<<std::endl;
    //std::cout<<"isMC: "<<mc_<<std::endl;

    edm::Handle<pat::METCollection> inMET;

    math::PtEtaPhiMLorentzVector correctionType1;
    math::PtEtaPhiMLorentzVector correctionXY;
    math::PtEtaPhiMLorentzVector inP4;
    math::PtEtaPhiMLorentzVector outP4;
    math::PtEtaPhiMLorentzVector antiType1;

    if(iEvent.getByLabel(srcMET_,inMET)){
     pat::MET met=inMET->at(0);
     pat::MET newmet=inMET->at(0);
 

     correctionType1 = findCorrectionType1(jetsL123, mc_);
     antiType1 = antiCorrectionType1(jetsL123);
     correctionXY = findCorrectionXY(Nvtx, mc_);
     if(mc_){
      inP4 = met.p4();
      outP4 = inP4 - correctionXY;
     }
     //if(mc_){outP4 = inP4;}
     if(!mc_){
      inP4 = met.p4();
      //inP4 = met.userCand("ues+")->p4();
      outP4 = inP4 + antiType1;//- correctionXY;
     }
     //outP4 = inP4 + correctionType1 ;//- correctionXY;
     newmet.setP4(outP4);
//     std::cout<<"MET In   pt "<<inP4.pt()      <<" eta "<<inP4.eta()<<" phi "<<inP4.phi()<<std::endl;
//     std::cout<<"MET Corr pt "<<correctionType1.pt()<<" eta "<<correctionType1.eta()<<" phi "<<correctionType1.phi()<<std::endl;
//     std::cout<<"MET Out  pt "<<outP4.pt()     <<" eta "<<outP4.eta()<<" phi "<<outP4.phi()<<std::endl;
//     std::cout<<"MET New  pt "<<newmet.pt()    <<" eta "<<newmet.eta()<<" phi "<<newmet.phi()<<std::endl;
//     std::cout<<std::endl;
//     std::cout<<"MET In   px "<<inP4.px()      <<" py "<<inP4.py()<<" pz "<<inP4.pz()<<std::endl;
//     std::cout<<"MET Corr px "<<correctionType1.px()<<" py "<<correctionType1.py()<<" pz "<<correctionType1.pz()<<std::endl;
//     std::cout<<"MET Out  px "<<outP4.px()     <<" py "<<outP4.py()<<" pz "<<outP4.pz()<<std::endl;
//     std::cout<<"MET New  px "<<newmet.px()    <<" py "<<newmet.py()<<" pz "<<newmet.pz()<<std::endl;
//     std::cout<<std::endl;

     outMET->push_back(newmet);
    }
    iEvent.put(outMET);
}

   // matches jets with L1 correction to jets with L123 corrections
pat::Jet PATMETTypeOne::jetMatchL1(pat::Jet jetL123, edm::Handle<pat::JetCollection> jetsL1){
 pat::Jet retVal;
 //const pat::Jet* retVal = 0;
 //double PT=-20;
 double dRbestMatch = 0.5;
   for (unsigned int i=0; i<jetsL1->size(); i++){
   pat::Jet  jetL1=jetsL1->at(i);
     double dR = deltaR(jetL123.p4(), jetL1.p4());
     if ( dR < dRbestMatch ) {
     retVal = jetL1;
     dRbestMatch = dR;
   }
 }
 return retVal;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/METTypeOneType2Formulae
math::PtEtaPhiMLorentzVector PATMETTypeOne::antiCorrectionType1(edm::Handle<pat::JetCollection> l123Jets)
{
 math::PtEtaPhiMLorentzVector corr;
 //std::cout<<corr.pt()<<" "<<corr.eta()<<" "<<corr.phi()<<std::endl;
 for(unsigned int i=0;i!=l123Jets->size();++i){
  pat::Jet jet123 = l123Jets->at(i);
  if (jet123.pt() > 10.)
  {
    corr += jet123.correctedP4("L1FastJet") - jet123.correctedP4("L3Absolute");
    //else{corr += jet123.correctedP4("L1FastJet") - jet123.correctedP4("L2L3Residual");}
   //corr += jet123.correctedP4("L1FastJet") - jet123.p4();
//  std::cout<<"L1   px "<<jet123.correctedP4("L1FastJet").px()
//               <<" py "<<jet123.correctedP4("L1FastJet").px()
//               <<" pz "<<jet123.correctedP4("L1FastJet").px()<<std::endl;
//  std::cout<<"L123 px "<<jet123.p4().px()
//               <<" py "<<jet123.p4().px()
//               <<" pz "<<jet123.p4().px()<<std::endl<<std::endl;
  }
 }
 return corr;
}
// https://twiki.cern.ch/twiki/bin/view/CMS/METTypeOneType2Formulae
math::PtEtaPhiMLorentzVector PATMETTypeOne::findCorrectionType1(edm::Handle<pat::JetCollection> l123Jets, bool isMC)
{
 math::PtEtaPhiMLorentzVector corr;
 //std::cout<<corr.pt()<<" "<<corr.eta()<<" "<<corr.phi()<<std::endl;
 for(unsigned int i=0;i!=l123Jets->size();++i){
  pat::Jet jet123 = l123Jets->at(i);
  if (jet123.pt() > 10.)
  {
    if(isMC){corr += jet123.correctedP4("L1FastJet") - jet123.correctedP4("L3Absolute");}
    else{corr += jet123.correctedP4("L1FastJet") - jet123.correctedP4("L2L3Residual");}
   //corr += jet123.correctedP4("L1FastJet") - jet123.p4();
//  std::cout<<"L1   px "<<jet123.correctedP4("L1FastJet").px()
//               <<" py "<<jet123.correctedP4("L1FastJet").px()
//               <<" pz "<<jet123.correctedP4("L1FastJet").px()<<std::endl;
//  std::cout<<"L123 px "<<jet123.p4().px()
//               <<" py "<<jet123.p4().px()
//               <<" pz "<<jet123.p4().px()<<std::endl<<std::endl;
  }
 }
 return corr;
}

// xy corrections from JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py
// https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_14_patch1/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py
math::PtEtaPhiMLorentzVector PATMETTypeOne::findCorrectionXY(int nrVtx, bool isMC)
{
 //math::XYZTLorentzVector corrXYZ;
 //math::PtEtaPhiMLorentzVector corr;
 //std::cout<<corrXYZ.pt()<<" "<<corrXYZ.px()<<" "<<corrXYZ.py()<<std::endl;
 //std::cout<<corr.pt()<<" "<<corr.px()<<" "<<corr.py()<<std::endl;
 if(isMC){
  math::XYZTLorentzVector corrXYZ(1.62861e-01 - 2.38517e-02*nrVtx,3.60860e-01 - 1.30335e-01*nrVtx,0,0);
  math::PtEtaPhiMLorentzVector corr(corrXYZ.pt(),corrXYZ.eta(),corrXYZ.phi(),corrXYZ.M());
  //std::cout<<corrXYZ.pt()<<" "<<corrXYZ.px()<<" "<<corrXYZ.py()<<std::endl;
  //std::cout<<corr.pt()<<" "<<corr.px()<<" "<<corr.py()<<std::endl;
  return corr;
  //corrXYZ.setP4(1.62861e-01 - 2.38517e-02*nrVtx,3.60860e-01 - 1.30335e-01*nrVtx,0,0);
  //corr.px() = 1.62861e-01 - 2.38517e-02*nrVtx;
  //corr.py() = 3.60860e-01 - 1.30335e-01*nrVtx;
 }
 else{
  math::XYZTLorentzVector corrXYZ(4.83642e-02 + 2.48870e-01*nrVtx,-1.50135e-01 - 8.27917e-02*nrVtx,0,0);
  math::PtEtaPhiMLorentzVector corr(corrXYZ.pt(),corrXYZ.eta(),corrXYZ.phi(),corrXYZ.M());
  //std::cout<<corrXYZ.pt()<<" "<<corrXYZ.px()<<" "<<corrXYZ.py()<<std::endl;
  //std::cout<<corr.pt()<<" "<<corr.px()<<" "<<corr.py()<<std::endl;
  return corr;
  //corrXYZ.setP4(4.83642e-02 + 2.48870e-01*nrVtx,-1.50135e-01 - 8.27917e-02*nrVtx,0,0);
  //corr.px() = 4.83642e-02 + 2.48870e-01*nrVtx;
  //corr.py() = -1.50135e-01 - 8.27917e-02*nrVtx;
 }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMETTypeOne);
