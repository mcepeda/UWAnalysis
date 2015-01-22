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


class GenStudy : public edm::EDFilter {

public:
 GenStudy (const edm::ParameterSet &);
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

 std::map<std::string,TH1D*> h1_;
 std::map<std::string,TH2D*> h2_;

 Int_t nall;
 Int_t nNUP;
 Int_t nLeptonCuts;
 Int_t nMuonOnly;
 Int_t n2BHadrons;
 Int_t n2Jets;
 Int_t n2TaggedJets;
 Int_t n2TaggedJetsV2;
 Int_t nVetoExtraJetsCentral;
 Int_t nVetoExtraJets;
 Int_t nsel;
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


GenStudy::GenStudy( const ParameterSet & cfg ) :
 LHEParticleTag_(cfg.getUntrackedParameter<edm::InputTag> ("LHETag", edm::InputTag("source"))),
 srcDressedLepton_(cfg.getUntrackedParameter<edm::InputTag> ("DressedLepton", edm::InputTag("dressedLeptons"))),
 srcJets_(cfg.getUntrackedParameter<edm::InputTag> ("Jets", edm::InputTag("cleanGenJets"))),
 srcSIMB_(cfg.getUntrackedParameter<edm::InputTag> ("SIMB", edm::InputTag("bHadrons"))),
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
 inputCrossSection_(cfg.getUntrackedParameter<double>("inputCrossSection",1))
{
}


void GenStudy::beginJob() {
 nall=0;
 nNUP=0;
 nLeptonCuts=0;
 nMuonOnly=0;
 n2BHadrons=0;
 n2Jets=0;
 n2TaggedJets=0;
 n2TaggedJetsV2=0;
 nVetoExtraJetsCentral=0;
 nVetoExtraJets=0;
 nsel=0;
 edm::Service<TFileService> fs;
 h1_["LHE_PARTONMULTIPLICITY_SEL"] =fs->make<TH1D>("LHE_PARTONMULTIPLICITY_SEL","hepeup().NUP",20,0.,20.);
 h1_["Selection"]=fs->make<TH1D>("Selection","",20,0,20);
 h1_["CrossSection"]=fs->make<TH1D>("CrossSection","",20,0,20);


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

}

void GenStudy::endJob() {
 cout<<"********************************************************************"<<endl;
 cout<<"GEN LEVEL FILTERING"<<endl<<endl;
 cout<<"Total Analyzed = "<<nall <<" For- > "<<filterByNUP_<<" -> "<<inputCrossSection_<<endl ;
 cout<<"LHE Selection = "<<nNUP <<" - > "<<1./nall*nNUP<<endl ;
 cout<<"Lepton Cuts = "<<nLeptonCuts <<" - "<<inputCrossSection_/nNUP*nLeptonCuts<<endl ;
 cout<<"Lepton PDGID = "<<nMuonOnly <<" - "<<inputCrossSection_/nNUP*nMuonOnly<<endl ;
 cout<<"2 BHadrons = "<<n2BHadrons <<" - "<<inputCrossSection_/nNUP*n2BHadrons<<endl ;
 cout<<"2 Jets = "<<n2Jets <<" - "<<inputCrossSection_/nNUP*n2Jets<<endl ;
 cout<<"Veto Central = "<<nVetoExtraJetsCentral<<" - "<<inputCrossSection_/nNUP*nVetoExtraJetsCentral<<endl;
 cout<<"2 TaggedJets = "<<n2TaggedJets <<" - "<<inputCrossSection_/nNUP*n2TaggedJets<<endl ;
 cout<<"2 TaggedJetsV2 = "<<n2TaggedJetsV2 <<" - "<<inputCrossSection_/nNUP*n2TaggedJetsV2<<endl ;
 cout<<"Veto Extra = "<<nVetoExtraJets <<" - "<<inputCrossSection_/nNUP*nVetoExtraJets<<endl ;
 cout<<"********************************************************************"<<endl;

 h1_["Selection"]->SetBinContent(1,nall);
 h1_["Selection"]->SetBinContent(2,nNUP);
 h1_["Selection"]->SetBinContent(3,nLeptonCuts);
 h1_["Selection"]->SetBinContent(4,nMuonOnly);
 h1_["Selection"]->SetBinContent(5,n2BHadrons);
 h1_["Selection"]->SetBinContent(6,n2Jets);
 h1_["Selection"]->SetBinContent(7,nVetoExtraJetsCentral);
 h1_["Selection"]->SetBinContent(8,n2TaggedJets);
 h1_["Selection"]->SetBinContent(9,n2TaggedJetsV2);
 h1_["Selection"]->SetBinContent(10,nVetoExtraJets);


 h1_["CrossSection"]->SetBinContent(1,inputCrossSection_);
 h1_["CrossSection"]->SetBinContent(2,inputCrossSection_);
 h1_["CrossSection"]->SetBinContent(3,inputCrossSection_/nNUP*nLeptonCuts);
 h1_["CrossSection"]->SetBinContent(4,inputCrossSection_/nNUP*nMuonOnly);
 h1_["CrossSection"]->SetBinContent(5,inputCrossSection_/nNUP*n2BHadrons);
 h1_["CrossSection"]->SetBinContent(6,inputCrossSection_/nNUP*n2Jets);
 h1_["CrossSection"]->SetBinContent(7,inputCrossSection_/nNUP*nVetoExtraJetsCentral);
 h1_["CrossSection"]->SetBinContent(8,inputCrossSection_/nNUP*n2TaggedJets);
 h1_["CrossSection"]->SetBinContent(9,inputCrossSection_/nNUP*n2TaggedJetsV2);
 h1_["CrossSection"]->SetBinContent(10,inputCrossSection_/nNUP*nVetoExtraJets);

}

bool GenStudy::filter (Event & ev, const EventSetup &) {
nall++;

 edm::Handle<LHEEventProduct> lheeventinfo;
 if(!ev.getByLabel(LHEParticleTag_, lheeventinfo)){
 LogDebug("") << ">>> LHE info not found!!";
 return false;
 }
 double NUP=lheeventinfo->hepeup().NUP;

 if(filterByNUP_!=-1 && NUP!=filterByNUP_) return false;

 nNUP++;

 Handle<pat::JetCollection > Jets;
 edm::Handle <DressedLeptonCollection> Leptons;
 edm::Handle <SimBHadronCollection> BHadrons;
 edm::Handle <GenParticleCollection> Neutrinos;


 if(!ev.getByLabel(srcDressedLepton_, Leptons)) edm::LogError("")<<"No Dressed Leptons!";
 if(!ev.getByLabel(srcJets_,Jets)) edm::LogError("")<<"No Jets!";
 if(!ev.getByLabel(srcSIMB_, BHadrons)) edm::LogError("")<<"No B Hadrons!!";
 if(!ev.getByLabel(srcNeutrinos_, Neutrinos)) edm::LogError("")<<"No Neutrinos!!";


 int leptonIndex=0;
 int countLeptons=0, countHardLeptons=0;
 for( size_t i = 0; i < Leptons->size(); ++ i ) {
 DressedLepton genpart = Leptons->at(i);
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

 nLeptonCuts++;

 if(selLeptonPDGID_!=0 && abs(lepton.pdgId())!=selLeptonPDGID_) return false;
 if(chooseSign_ && lepton.pdgId()*selLeptonPDGID_<0) return false;
 nMuonOnly++;

 if(BHadrons->size()<2) return false;
 n2BHadrons++;

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

 if(jet.userFloat("matchedHadron")>0) countBJets++;
 if(jet.userFloat("matchedHadron")>0 && jet.userFloat("leadJetBHadronPt")>5) countBJetsMinPt++;

 if(jet.pt()>leadPt) {secondPt=leadPt; leadPt=jet.pt(); lead=i;}
 else if(jet.pt()>secondPt) {secondPt=jet.pt(); second=i;}

 }

 if(countJets<nMinJets_) return false;

 pat::Jet jet1 = Jets->at(lead);
 pat::Jet jet2 = Jets->at(second);

 math::XYZTLorentzVector diJet=jet1.p4()+jet2.p4();

 double dEta=(jet1.eta()-jet2.eta());
 double dPhi=deltaPhi(jet1.phi(),jet2.phi());
 double dR=deltaR(jet1,jet2);


 n2Jets++;

 if(countJets>nCentralJetsForVeto_) return false;

 nVetoExtraJetsCentral++;

 if(countBJets<nBJets_) return false;

 n2TaggedJets++;

 if(jet1.userFloat("matchedHadron")<0 && jet2.userFloat("matchedHadron")<0) return false;

 n2TaggedJetsV2++;

 if(countJetsPT<=nJetsForVeto_) { nVetoExtraJets++;}


 h1_["LHE_PARTONMULTIPLICITY_SEL"]->Fill(NUP);

 h1_["LeadJetPt"]->Fill(jet1.pt());
 h1_["SecondJetPt"]->Fill(jet2.pt());
 h1_["LeptonPt"]->Fill(lepton.pt());
 h1_["NeutrinoPt"]->Fill(MET);
 h1_["LeadJetBPt"]->Fill(jet1.userFloat("leadJetBHadronPt"));
 h1_["SecondJetBPt"]->Fill(jet2.userFloat("leadJetBHadronPt"));

 h1_["LeadJetEta"]->Fill(jet1.eta());
 h1_["SecondJetEta"]->Fill(jet2.eta());
 h1_["LeptonEta"]->Fill(lepton.eta());
 h1_["LeadJetBEta"]->Fill(jet1.userFloat("leadJetBHadronEta"));
 h1_["SecondJetBEta"]->Fill(jet2.userFloat("leadJetBHadronEta"));

 h1_["WMT"]->Fill(wmt);
 h1_["WPT"]->Fill(wpt);

 h1_["hNJetsCentral"]->Fill(countJets);
 h1_["hNJets"]->Fill(countJetsPT);
 h1_["hNBJets"]->Fill(countBJets);

 h1_["MassDiBJets"]->Fill(diJet.mass());
 h1_["PtDiBJets"]->Fill(diJet.pt());
 h1_["DRDiBJets"]->Fill(dR);
 h1_["DPhiDiBJets"]->Fill(dPhi);
 h1_["DEtaDiBJets"]->Fill(dEta);
 h1_["PtBalanceBJets"]->Fill(diJet.pt()/(jet1.pt()+jet2.pt()));


return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenStudy);
