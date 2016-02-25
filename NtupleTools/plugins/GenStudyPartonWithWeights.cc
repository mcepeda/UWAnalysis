#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "ZSV/BAnalysis/interface/SimBHadron.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/angle.h"
#include "UWAnalysis/DataFormats/interface/DressedLepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


class GenStudyPartonWithWeights : public edm::EDFilter {

 public:
  GenStudyPartonWithWeights (const edm::ParameterSet &);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();
 private:
  edm::InputTag LHEParticleTag_;
  edm::InputTag srcDressedLepton_;
  edm::InputTag srcJets_;
  edm::InputTag srcSIMB_;
  edm::InputTag srcNeutrinos_;

  int filterByNUP_;
  int selLeptonPDGID_;
  bool chooseSign_;

  double ptCutJets_;
  double etaCutJets_;

  double ptCutLepton_;
  double etaCutLepton_;

  int nCentralJetsForVeto_;
  int nJetsForVeto_;
  int nMinJets_;
  int nBJets_;

  double inputCrossSection_;

  edm::InputTag pdfInfoTag_;
  edm::InputTag pdfSetNames_;
  std::vector<std::string> pdfShortNames_;



  std::map<std::string,TH1D*> h1_;
  std::map<std::string,TH2D*> h2_;

  Double_t nall;
  Double_t nNUP;
  Double_t nLeptonCuts;
  Double_t nMuonOnly;
  Double_t n2BPartons;
  Double_t n2Jets;
  Double_t n2TaggedJets;
  Double_t n2TaggedJetsV2;
  Double_t nVetoExtraJetsCentral;
  Double_t nVetoExtraJets;
  Double_t nsel;

  Double_t nNoWeightall;
  Double_t nNoWeightNUP;
  Double_t nNoWeightLeptonCuts;
  Double_t nNoWeightMuonOnly;
  Double_t nNoWeight2BPartons;
  Double_t nNoWeight2Jets;
  Double_t nNoWeight2TaggedJets;
  Double_t nNoWeight2TaggedJetsV2;
  Double_t nNoWeightVetoExtraJetsCentral;
  Double_t nNoWeightVetoExtraJets;
  Double_t nNoWeightsel;


};
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <map>
#include <memory>

using namespace edm;
using namespace std;
using namespace reco;


GenStudyPartonWithWeights::GenStudyPartonWithWeights( const ParameterSet & cfg ) :
 LHEParticleTag_(cfg.getUntrackedParameter<edm::InputTag> ("LHETag", edm::InputTag("source"))),
 srcDressedLepton_(cfg.getUntrackedParameter<edm::InputTag> ("DressedLepton", edm::InputTag("dressedLeptons"))),
 srcJets_(cfg.getUntrackedParameter<edm::InputTag> ("Jets", edm::InputTag("cleanGenJets"))),
 srcSIMB_(cfg.getUntrackedParameter<edm::InputTag> ("SIMB", edm::InputTag("bPartons"))),
 srcNeutrinos_(cfg.getUntrackedParameter<edm::InputTag> ("Neutrinos", edm::InputTag("neutrinos"))),
 filterByNUP_(cfg.getUntrackedParameter<int>("filterByNUP",-1)),
 selLeptonPDGID_(cfg.getUntrackedParameter<int>("selLeptonPDGID",-1)),
 chooseSign_(cfg.getUntrackedParameter<bool>("chooseSign",false)),
 ptCutJets_(cfg.getUntrackedParameter<double>("ptCutJets", 25)),
 etaCutJets_(cfg.getUntrackedParameter<double>("etaCutJets", 2.4)),
 ptCutLepton_(cfg.getUntrackedParameter<double>("ptCutLepton", 30)),
 etaCutLepton_(cfg.getUntrackedParameter<double>("etaCutLepton", 2.1)),
 nCentralJetsForVeto_(cfg.getUntrackedParameter<int>("nCentralJetsForVeto",2)),
 nJetsForVeto_(cfg.getUntrackedParameter<int>("nJetsForVeto",2)),
 nMinJets_(cfg.getUntrackedParameter<int>("nMinJets",2)),
 nBJets_(cfg.getUntrackedParameter<int>("nBJets",2)),
 inputCrossSection_(cfg.getUntrackedParameter<double>("inputCrossSection",1)),
 pdfInfoTag_(cfg.getUntrackedParameter<edm::InputTag> ("PdfInfoTag", edm::InputTag("generator"))),
 pdfSetNames_(cfg.getUntrackedParameter<edm::InputTag> ("PdfSetNames"))
{
}


void GenStudyPartonWithWeights::beginJob() {
 nall=0;
 nNUP=0;
 nLeptonCuts=0;
 nMuonOnly=0;
 n2BPartons=0;
 n2Jets=0;
 n2TaggedJets=0;
 n2TaggedJetsV2=0;
 nVetoExtraJetsCentral=0;
 nVetoExtraJets=0;
 nsel=0;

 nNoWeightall=0;
 nNoWeightNUP=0;
 nNoWeightLeptonCuts=0;
 nNoWeightMuonOnly=0;
 nNoWeight2BPartons=0;
 nNoWeight2Jets=0;
 nNoWeight2TaggedJets=0;
 nNoWeight2TaggedJetsV2=0;
 nNoWeightVetoExtraJetsCentral=0;
 nNoWeightVetoExtraJets=0;
 nNoWeightsel=0;


 edm::Service<TFileService> fs;
 h1_["LHE_PARTONMULTIPLICITY_SEL"] =fs->make<TH1D>("LHE_PARTONMULTIPLICITY_SEL","hepeup().NUP",20,0.,20.);
 h1_["Selection"]=fs->make<TH1D>("Selection","",20,0,20);
 h1_["CrossSection"]=fs->make<TH1D>("CrossSection","",20,0,20);
 h1_["Selection_NoWeight"]=fs->make<TH1D>("Selection_NoWeight","",20,0,20);

 // Plots
 h1_["WMT"]=fs->make<TH1D>("WMT","",500,0,500);
 h1_["WPT"]=fs->make<TH1D>("WPT","",500,0,500);

 h1_["LeadJetPt"]=fs->make<TH1D>("LeadJetPt","",500,0,500);
 h1_["SecondJetPt"]=fs->make<TH1D>("SecondJetPt","",500,0,500);
 h1_["LeptonPt"]=fs->make<TH1D>("LeptonPt","",500,0,500);
 h1_["NeutrinoPt"]=fs->make<TH1D>("NeutrinoPt","",500,0,500);
 h1_["LeadJetBPt"]=fs->make<TH1D>("LeadJetBPt","",500,0,500);
 h1_["SecondJetBPt"]=fs->make<TH1D>("SecondJetBPt","",500,0,500);

 h1_["LeadJetEta"]=fs->make<TH1D>("LeadJetEta","",100,-5,5);
 h1_["SecondJetEta"]=fs->make<TH1D>("SecondJetEta","",100,-5,5);
 h1_["LeptonEta"]=fs->make<TH1D>("LeptonEta","",100,-5,5);
 h1_["LeadJetBEta"]=fs->make<TH1D>("LeadJetBEta","",100,-5,5);
 h1_["SecondJetBEta"]=fs->make<TH1D>("SecondJetBEta","",100,-5,5);

 h1_["hNJetsCentral"]=fs->make<TH1D>("hNJetsCentral","",10,-0.5,9.5);
 h1_["hNJets"]=fs->make<TH1D>("hNJets","",10,-0.5,9.5);
 h1_["hNBJets"]=fs->make<TH1D>("hNBJets","",10,-0.5,9.5);

 h1_["MassDiBJets"]=fs->make<TH1D>("MassDiBJets","",500,0,500);
 h1_["PtDiBJets"]=fs->make<TH1D>("PtDiBJets","",500,0,500);
 h1_["DRDiBJets"]=fs->make<TH1D>("DRDiBJets","",100,0,10);
 h1_["DPhiDiBJets"]=fs->make<TH1D>("DPhiDiBJets","",100,-6.3,6.3);
 h1_["DEtaDiBJets"]=fs->make<TH1D>("DEtaDiBJets","",100,-5,5);
 h1_["PtBalanceBJets"]=fs->make<TH1D>("PtBalanceBJets","",100,0,1);


 // Test PDFWeight
 h1_["Initial_PDF"]=fs->make<TH1D>("Initial_PDF","",500,-0.5,499.5);
 h1_["AfterSelectionCount_PDF"]=fs->make<TH1D>("AfterSelectionCount_PDF","",500,-0.5,499.5);

 h2_["WMT_PDF"]=fs->make<TH2D>("WMT_PDF","",500,0,500,500,-0.5,499.5);
 h2_["WPT_PDF"]=fs->make<TH2D>("WPT_PDF","",500,0,500,500,-0.5,499.5);

 h2_["LeadJetPt_PDF"]=fs->make<TH2D>("LeadJetPt_PDF","",500,0,500,500,-0.5,499.5);
 h2_["SecondJetPt_PDF"]=fs->make<TH2D>("SecondJetPt_PDF","",500,0,500,500,-0.5,499.5);
 h2_["LeptonPt_PDF"]=fs->make<TH2D>("LeptonPt_PDF","",500,0,500,500,-0.5,499.5);
 h2_["NeutrinoPt_PDF"]=fs->make<TH2D>("NeutrinoPt_PDF","",500,0,500,500,-0.5,499.5);
 h2_["LeadJetBPt_PDF"]=fs->make<TH2D>("LeadJetBPt_PDF","",500,0,500,500,-0.5,499.5);
 h2_["SecondJetBPt_PDF"]=fs->make<TH2D>("SecondJetBPt_PDF","",500,0,500,500,-0.5,499.5);

 h2_["LeadJetEta_PDF"]=fs->make<TH2D>("LeadJetEta_PDF","",100,-5,5,500,-0.5,499.5);
 h2_["SecondJetEta_PDF"]=fs->make<TH2D>("SecondJetEta_PDF","",100,-5,5,500,-0.5,499.5);
 h2_["LeptonEta_PDF"]=fs->make<TH2D>("LeptonEta_PDF","",100,-5,5,500,-0.5,499.5);
 h2_["LeadJetBEta_PDF"]=fs->make<TH2D>("LeadJetBEta_PDF","",100,-5,5,500,-0.5,499.5);
 h2_["SecondJetBEta_PDF"]=fs->make<TH2D>("SecondJetBEta_PDF","",100,-5,5,500,-0.5,499.5);

 h2_["hNJetsCentral_PDF"]=fs->make<TH2D>("hNJetsCentral_PDF","",10,-0.5,9.5,500,-0.5,499.5);
 h2_["hNJets_PDF"]=fs->make<TH2D>("hNJets_PDF","",10,-0.5,9.5,500,-0.5,499.5);
 h2_["hNBJets_PDF"]=fs->make<TH2D>("hNBJets_PDF","",10,-0.5,9.5,500,-0.5,499.5);

 h2_["MassDiBJets_PDF"]=fs->make<TH2D>("MassDiBJets_PDF","",500,0,500,500,-0.5,499.5);
 h2_["PtDiBJets_PDF"]=fs->make<TH2D>("PtDiBJets_PDF","",500,0,500,500,-0.5,499.5);
 h2_["DRDiBJets_PDF"]=fs->make<TH2D>("DRDiBJets_PDF","",100,0,10,500,-0.5,499.5);
 h2_["DPhiDiBJets_PDF"]=fs->make<TH2D>("DPhiDiBJets_PDF","",100,-6.3,6.3,500,-0.5,499.5);
 h2_["DEtaDiBJets_PDF"]=fs->make<TH2D>("DEtaDiBJets_PDF","",100,-5,5,500,-0.5,499.5);
 h2_["PtBalanceBJets_PDF"]=fs->make<TH2D>("PtBalanceBJets_PDF","",100,0,1,500,-0.5,499.5);

}

void GenStudyPartonWithWeights::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"GEN LEVEL FILTERING"<<endl<<endl;
 cout<<"Total Analyzed = "<<nall <<" For- > "<<filterByNUP_<<" -> "<<inputCrossSection_<<endl ;
 cout<<"LHE Selection = "<<nNUP <<" - > "<<1./nall*nNUP<<endl ;
 cout<<"Lepton Cuts = "<<nLeptonCuts <<" - "<<inputCrossSection_/nNUP*nLeptonCuts<<endl ;
 cout<<"Lepton PDGID = "<<nMuonOnly <<" - "<<inputCrossSection_/nNUP*nMuonOnly<<endl ;
 cout<<"2 BPartons = "<<n2BPartons <<" - "<<inputCrossSection_/nNUP*n2BPartons<<endl ;
 cout<<"2 Jets = "<<n2Jets <<" - "<<inputCrossSection_/nNUP*n2Jets<<endl ;
 cout<<"Veto Central = "<<nVetoExtraJetsCentral<<" - "<<inputCrossSection_/nNUP*nVetoExtraJetsCentral<<endl;
 cout<<"2 TaggedJets = "<<n2TaggedJets <<" - "<<inputCrossSection_/nNUP*n2TaggedJets<<endl ;
 cout<<"2 TaggedJetsV2 = "<<n2TaggedJetsV2 <<" - "<<inputCrossSection_/nNUP*n2TaggedJetsV2<<endl ;
 cout<<"Veto Extra = "<<nVetoExtraJets <<" - "<<inputCrossSection_/nNUP*nVetoExtraJets<<endl ;

 cout<<"No PDF Reiweight?"<<nNoWeightall<<"   -->  "<<nNoWeight2TaggedJetsV2<<endl;
 
 cout<<"********************************************************************"<<endl;

 h1_["Selection"]->SetBinContent(1,nall);
 h1_["Selection"]->SetBinContent(2,nNUP);
 h1_["Selection"]->SetBinContent(3,nLeptonCuts);
 h1_["Selection"]->SetBinContent(4,nMuonOnly);
 h1_["Selection"]->SetBinContent(5,n2BPartons);
 h1_["Selection"]->SetBinContent(6,n2Jets);
 h1_["Selection"]->SetBinContent(7,nVetoExtraJetsCentral);
 h1_["Selection"]->SetBinContent(8,n2TaggedJets);
 h1_["Selection"]->SetBinContent(9,n2TaggedJetsV2);
 h1_["Selection"]->SetBinContent(10,nVetoExtraJets);


 h1_["CrossSection"]->SetBinContent(1,inputCrossSection_);
 h1_["CrossSection"]->SetBinContent(2,inputCrossSection_);
 h1_["CrossSection"]->SetBinContent(3,inputCrossSection_/nNUP*nLeptonCuts);
 h1_["CrossSection"]->SetBinContent(4,inputCrossSection_/nNUP*nMuonOnly);
 h1_["CrossSection"]->SetBinContent(5,inputCrossSection_/nNUP*n2BPartons);
 h1_["CrossSection"]->SetBinContent(6,inputCrossSection_/nNUP*n2Jets);
 h1_["CrossSection"]->SetBinContent(7,inputCrossSection_/nNUP*nVetoExtraJetsCentral);
 h1_["CrossSection"]->SetBinContent(8,inputCrossSection_/nNUP*n2TaggedJets);
 h1_["CrossSection"]->SetBinContent(9,inputCrossSection_/nNUP*n2TaggedJetsV2);
 h1_["CrossSection"]->SetBinContent(10,inputCrossSection_/nNUP*nVetoExtraJets);


 h1_["Selection_NoWeight"]->SetBinContent(1,nNoWeightall);
 h1_["Selection_NoWeight"]->SetBinContent(2,nNoWeightNUP);
 h1_["Selection_NoWeight"]->SetBinContent(3,nNoWeightLeptonCuts);
 h1_["Selection_NoWeight"]->SetBinContent(4,nNoWeightMuonOnly);
 h1_["Selection_NoWeight"]->SetBinContent(5,nNoWeight2BPartons);
 h1_["Selection_NoWeight"]->SetBinContent(6,nNoWeight2Jets);
 h1_["Selection_NoWeight"]->SetBinContent(7,nNoWeightVetoExtraJetsCentral);
 h1_["Selection_NoWeight"]->SetBinContent(8,nNoWeight2TaggedJets);
 h1_["Selection_NoWeight"]->SetBinContent(9,nNoWeight2TaggedJetsV2);
 h1_["Selection_NoWeight"]->SetBinContent(10,nNoWeightVetoExtraJets);


}

bool GenStudyPartonWithWeights::filter (Event & ev, const EventSetup &) {

 // First of all get the PDF weight:

 double weight=1;
 double pdfmemberWeight[500];
 for (int nWeights=0; nWeights<500; nWeights++){ pdfmemberWeight[nWeights]=0;}

 double NUP=0;
 edm::Handle<LHEEventProduct> lheeventinfo;
 if(!ev.getByLabel(LHEParticleTag_, lheeventinfo)){
  LogDebug("") << ">>> LHE info not found!!";
 }
 else{
  NUP=lheeventinfo->hepeup().NUP;

/*  cout<<"NUP?"<< lheeventinfo->hepeup().NUP<<endl;
  cout<<"AQCDUP?"<< lheeventinfo->hepeup().AQCDUP<<endl;
  cout<<"XPDWUP?"<< lheeventinfo->hepeup().XPDWUP.first<<"   "<<lheeventinfo->hepeup().XPDWUP.second<<endl;
  cout<<"IDUP?"<< lheeventinfo->hepeup().IDUP.at(0)<<"  "<<lheeventinfo->hepeup().IDUP.at(1)<<endl;
*/
  //cout<<lheeventinfo->weights().size()<<endl;
  for (unsigned int i=0; i<lheeventinfo->weights().size(); i++){
  //                  cout<<i<<"   "<<lheeventinfo->weights().at(i).id<<"    "<<lheeventinfo->weights().at(i).wgt<<endl;
                    pdfmemberWeight[i]=lheeventinfo->weights().at(i).wgt; 
  }
  //cout<<endl;
 }

 for (int nWeights=0; nWeights<500; nWeights++){
    if(pdfmemberWeight[nWeights]==0) break;
    h1_["Initial_PDF"]->Fill(nWeights,pdfmemberWeight[nWeights]);

  }

 nNoWeightall++;
 nall+=weight;

 if(filterByNUP_!=-1 && NUP!=filterByNUP_) return false;

 nNUP+=weight;
 nNoWeightNUP++;

 Handle<pat::JetCollection > Jets;
 edm::Handle <DressedLeptonCollection> Leptons;
 edm::Handle <GenParticleCollection> Neutrinos;
 edm::Handle <GenParticleRefVector> BPartons;;


 if(!ev.getByLabel(srcDressedLepton_, Leptons)) edm::LogError("")<<"No Dressed Leptons!";
 if(!ev.getByLabel(srcJets_,Jets)) edm::LogError("")<<"No Jets!";
 if(!ev.getByLabel(srcSIMB_, BPartons)) edm::LogError("")<<"No B Partons!!";
 if(!ev.getByLabel(srcNeutrinos_, Neutrinos)) edm::LogError("")<<"No Neutrinos!!";


 int leptonIndex=0;
 int countLeptons=0, countHardLeptons=0;
 for( size_t i = 0; i < Leptons->size(); ++ i ) {
  DressedLepton genpart = Leptons->at(i);
  if(genpart.pt()<1) continue; 
  if(genpart.status()!=1) continue;
  countLeptons++;
  if(genpart.Type()!=1) continue;
  countHardLeptons++;
  if(genpart.pt()<Leptons->at(leptonIndex).pt()) continue;
  leptonIndex=i;
 }

 if(countHardLeptons<1) { return false; }
 else if (countHardLeptons>1) {edm::LogError("")<<"You cannot have 2 hard leptons, this is a W event! NL:"<<countHardLeptons;}

 DressedLepton lepton=Leptons->at(leptonIndex);

 if(lepton.pt()<ptCutLepton_ || fabs(lepton.eta())>=etaCutLepton_) return false;

 nLeptonCuts+=weight;
 nNoWeightLeptonCuts++;

 if(selLeptonPDGID_!=0 && abs(lepton.pdgId())!=selLeptonPDGID_) return false;
 if(chooseSign_ && lepton.pdgId()*selLeptonPDGID_<0) return false;
 nNoWeightMuonOnly++;
 nMuonOnly+=weight;

 int nBPartons=0;

 for (size_t j=0; j<BPartons->size(); j++){
    const GenParticle & genpart = *(BPartons->at(j).get());
    if(abs(genpart.pdgId())!=5)  continue; 
    nBPartons++;
 }


 if(nBPartons<2) return false;

 nNoWeight2BPartons++;
 n2BPartons+=weight;

 double met_px=0, met_py=0;


 for (size_t j=0; j<Neutrinos->size(); j++){
  const reco::GenParticle Nu=Neutrinos->at(j);
  met_px+=Nu.px();
  met_py+=Nu.py();
 }
 double MET=sqrt(met_px*met_px+met_py*met_py);


 double mt_px=lepton.px()+met_px;
 double mt_py=lepton.py()+met_py;
 double wpt=sqrt(mt_px*mt_px+mt_py*mt_py);
 double mt_et=lepton.pt()+MET;
 double wmt=mt_et*mt_et-mt_px*mt_px-mt_py*mt_py;
 wmt=(wmt>0)?sqrt(wmt):0;


 double countJetsPT=0, countJets=0, countBJets=0,countBJetsMinPt=0;
 int lead=-1, second=-1;
 double leadPt=0, secondPt=0;
 for (unsigned int i=0; i<Jets->size(); i++){
  pat::Jet jet = Jets->at(i);

  if(jet.pt()<ptCutJets_) continue;
  countJetsPT++;
  if( fabs(jet.eta())>=etaCutJets_) continue;
  countJets++;

  if(jet.userFloat("matchedBsSt2")>0) countBJets++;
  if(jet.userFloat("matchedBsSt2")>0 && jet.userFloat("leadJetSt2Pt")>5) countBJetsMinPt++;

  if(jet.pt()>leadPt) {secondPt=leadPt; leadPt=jet.pt(); lead=i;}
  else if(jet.pt()>secondPt) {secondPt=jet.pt(); second=i;}


 }



/*
        jet.addUserFloat("matchedBsSt2",matchedSt2);
        jet.addUserFloat("minDRSt2",dRSt2);
        jet.addUserFloat("leadJetSt2Pt",leadJetSt2Pt);
        jet.addUserFloat("secondJetSt2Pt",secondJetSt2Pt);
*/

 if(countJets<nMinJets_) return false;

 pat::Jet jet1 = Jets->at(lead);
 pat::Jet jet2 = Jets->at(second);

 math::XYZTLorentzVector diJet=jet1.p4()+jet2.p4();

 double dEta=(jet1.eta()-jet2.eta());
 double dPhi=deltaPhi(jet1.phi(),jet2.phi());
 double dR=deltaR(jet1,jet2);



 n2Jets+=weight;
 nNoWeight2Jets++;


 if(countJets>nCentralJetsForVeto_) return false;

 nVetoExtraJetsCentral+=weight;
 nNoWeightVetoExtraJetsCentral++;


 if(countBJets<nBJets_) return false;

 n2TaggedJets+=weight;
 nNoWeight2TaggedJets++;


 if(jet1.userFloat("matchedBsSt2")<0 && jet2.userFloat("matchedBsSt2")<0) return false;
 if(countBJetsMinPt<nBJets_) return false;  // added a min BHadron Pt CUt

  if(jet1.userFloat("leadJetSt2Pt") == jet2.userFloat("leadJetSt2Pt") ) {
        cout<<"ARGGGGGGGH Two jets matched to the same parton!!!!"<<endl; return false;;}
   
 n2TaggedJetsV2+=weight;
 nNoWeight2TaggedJetsV2++;


 if(countJetsPT<=nJetsForVeto_) { 
        nVetoExtraJets+=weight;
        nNoWeightVetoExtraJets++;}

 h1_["LHE_PARTONMULTIPLICITY_SEL"]->Fill(NUP,weight);

 h1_["LeadJetPt"]->Fill(jet1.pt(),weight);
 h1_["SecondJetPt"]->Fill(jet2.pt(),weight);
 h1_["LeptonPt"]->Fill(lepton.pt(),weight);
 h1_["NeutrinoPt"]->Fill(MET,weight);
 h1_["LeadJetBPt"]->Fill(jet1.userFloat("leadJetSt2Pt"),weight);
 h1_["SecondJetBPt"]->Fill(jet2.userFloat("leadJetSt2Pt"),weight);

 h1_["LeadJetEta"]->Fill(jet1.eta(),weight);
 h1_["SecondJetEta"]->Fill(jet2.eta(),weight);
 h1_["LeptonEta"]->Fill(lepton.eta(),weight);
 h1_["LeadJetBEta"]->Fill(jet1.userFloat("leadJetBHadronEta"),weight);
 h1_["SecondJetBEta"]->Fill(jet2.userFloat("leadJetBHadronEta"),weight);

 h1_["WMT"]->Fill(wmt,weight);
 h1_["WPT"]->Fill(wpt,weight);

 h1_["hNJetsCentral"]->Fill(countJets,weight);
 h1_["hNJets"]->Fill(countJetsPT,weight);
 h1_["hNBJets"]->Fill(countBJets,weight);

 h1_["MassDiBJets"]->Fill(diJet.mass(),weight);
 h1_["PtDiBJets"]->Fill(diJet.pt(),weight);
 h1_["DRDiBJets"]->Fill(dR,weight);
 h1_["DPhiDiBJets"]->Fill(dPhi,weight);
 h1_["DEtaDiBJets"]->Fill(dEta,weight);
 h1_["PtBalanceBJets"]->Fill(diJet.pt()/(jet1.pt()+jet2.pt()),weight);

 // For PDF Band

// cout<<"Adding weights!"<<endl;
 
 for (int nWeights=0; nWeights<500; nWeights++){

    if(pdfmemberWeight[nWeights]==0) break; 

    h1_["AfterSelectionCount_PDF"]->Fill(nWeights,pdfmemberWeight[nWeights]);

// cout<<nWeights<<"   "<<pdfmemberWeight[nWeights]<<endl;

 h2_["LeadJetPt_PDF"]->Fill(jet1.pt(),nWeights,pdfmemberWeight[nWeights]);
 h2_["SecondJetPt_PDF"]->Fill(jet2.pt(),nWeights,pdfmemberWeight[nWeights]);
 h2_["LeptonPt_PDF"]->Fill(lepton.pt(),nWeights,pdfmemberWeight[nWeights]);
 h2_["NeutrinoPt_PDF"]->Fill(MET,nWeights,pdfmemberWeight[nWeights]);
 h2_["LeadJetBPt_PDF"]->Fill(jet1.userFloat("leadJetSt2Pt"),nWeights,pdfmemberWeight[nWeights]);
 h2_["SecondJetBPt_PDF"]->Fill(jet2.userFloat("leadJetSt2Pt"),nWeights,pdfmemberWeight[nWeights]);

 h2_["LeadJetEta_PDF"]->Fill(jet1.eta(),nWeights,pdfmemberWeight[nWeights]);
 h2_["SecondJetEta_PDF"]->Fill(jet2.eta(),nWeights,pdfmemberWeight[nWeights]);
 h2_["LeptonEta_PDF"]->Fill(lepton.eta(),nWeights,pdfmemberWeight[nWeights]);
 h2_["LeadJetBEta_PDF"]->Fill(jet1.userFloat("leadJetBHadronEta"),nWeights,pdfmemberWeight[nWeights]);
 h2_["SecondJetBEta_PDF"]->Fill(jet2.userFloat("leadJetBHadronEta"),nWeights,pdfmemberWeight[nWeights]);

 h2_["WMT_PDF"]->Fill(wmt,nWeights,pdfmemberWeight[nWeights]);
 h2_["WPT_PDF"]->Fill(wpt,nWeights,pdfmemberWeight[nWeights]);

 h2_["hNJetsCentral_PDF"]->Fill(countJets,nWeights,pdfmemberWeight[nWeights]);
 h2_["hNJets_PDF"]->Fill(countJetsPT,nWeights,pdfmemberWeight[nWeights]);
 h2_["hNBJets_PDF"]->Fill(countBJets,nWeights,pdfmemberWeight[nWeights]);

 h2_["MassDiBJets_PDF"]->Fill(diJet.mass(),nWeights,pdfmemberWeight[nWeights]);
 h2_["PtDiBJets_PDF"]->Fill(diJet.pt(),nWeights,pdfmemberWeight[nWeights]);
 h2_["DRDiBJets_PDF"]->Fill(dR,nWeights,pdfmemberWeight[nWeights]);
 h2_["DPhiDiBJets_PDF"]->Fill(dPhi,nWeights,pdfmemberWeight[nWeights]);
 h2_["DEtaDiBJets_PDF"]->Fill(dEta,nWeights,pdfmemberWeight[nWeights]);
 h2_["PtBalanceBJets_PDF"]->Fill(diJet.pt()/(jet1.pt()+jet2.pt()),nWeights,pdfmemberWeight[nWeights]);


 }
 return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenStudyPartonWithWeights);
