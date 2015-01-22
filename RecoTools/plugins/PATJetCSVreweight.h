// system include files
#include <memory>
#include <algorithm>
#include <TLorentzVector.h>
#include <string> 

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

#include "TFile.h"
#include "TROOT.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"

class PATJetCSVreweight : public edm::EDProducer{
 public:

  explicit PATJetCSVreweight(const edm::ParameterSet& iConfig):
   src_(iConfig.getParameter<edm::InputTag>("src"))
  {
   produces<pat::JetCollection>();
  }
   ~PATJetCSVreweight() {}


 private:


  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    /* Fill vectors of histograms (systematics) */
    TString base = std::getenv("CMSSW_BASE");
    TString n_CSVwgt_HF = "/src/UWAnalysis/Configuration/data/csv_rwt_hf.root";
    TString n_CSVwgt_LF = "/src/UWAnalysis/Configuration/data/csv_rwt_lf.root";
    TFile *f_CSVwgt_HF = new TFile ( base + n_CSVwgt_HF );
    TFile *f_CSVwgt_LF = new TFile ( base + n_CSVwgt_LF );

    TH1D* h_csv_wgt_hf[9][6];
    TH1D* c_csv_wgt_hf[5][6];
    TH1D* h_csv_wgt_lf[9][4][3];

    for( int iSys=0; iSys<9; iSys++ ){
      std::string syst_csv_suffix_hf = "final";
      std::string syst_csv_suffix_c = "final";
      std::string syst_csv_suffix_lf = "final";

      switch( iSys ){
      case 0:  // this is the nominal case
        break;
      case 1:  // JESUp
        syst_csv_suffix_hf = "final_JESUp";
        syst_csv_suffix_lf = "final_JESUp";
        syst_csv_suffix_c  = "final_cErr1Up";
        break;
      case 2:  // JESDown
        syst_csv_suffix_hf = "final_JESDown";
        syst_csv_suffix_lf = "final_JESDown";
        syst_csv_suffix_c  = "final_cErr1Down";
        break;
      case 3:  // purity up
        syst_csv_suffix_hf = "final_LFUp";
        syst_csv_suffix_lf = "final_HFUp";
        syst_csv_suffix_c  = "final_cErr2Up";
        break;
      case 4:  // purity down
        syst_csv_suffix_hf = "final_LFDown";
        syst_csv_suffix_lf = "final_HFDown";
        syst_csv_suffix_c  = "final_cErr2Down";
        break;
      case 5:  // stats1 up
        syst_csv_suffix_hf = "final_Stats1Up";
        syst_csv_suffix_lf = "final_Stats1Up";
        break;
      case 6:  // stats1 down
        syst_csv_suffix_hf = "final_Stats1Down";
        syst_csv_suffix_lf = "final_Stats1Down";
        break;
      case 7:  // stats2 up
        syst_csv_suffix_hf = "final_Stats2Up";
        syst_csv_suffix_lf = "final_Stats2Up";
        break;
      case 8:  // stats2 down
        syst_csv_suffix_hf = "final_Stats2Down";
        syst_csv_suffix_lf = "final_Stats2Down";
        break;
      }

      std::string plot;
      for( int iPt=0; iPt<6; iPt++ ) {
       plot ="csv_ratio_Pt"+std::to_string(iPt)+"_Eta0_"+syst_csv_suffix_hf;
       h_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( (TString)plot );
      }

      if( iSys<5 ){
        for( int iPt=0; iPt<6; iPt++ ) {
          plot = "c_csv_ratio_Pt"+std::to_string(iPt)+"_Eta0_"+syst_csv_suffix_c;
          c_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( (TString)plot );
        }
      }
      
      for( int iPt=0; iPt<4; iPt++ ){
        for( int iEta=0; iEta<3; iEta++ ) {
          plot = "csv_ratio_Pt"+std::to_string(iPt)+"_Eta"+std::to_string(iEta)+"_"+syst_csv_suffix_lf;
          h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)f_CSVwgt_LF->Get( (TString)plot );
        }
      }
    }


    f_CSVwgt_HF->Close();
    f_CSVwgt_LF->Close();
    f_CSVwgt_HF->Delete();
    f_CSVwgt_LF->Delete();



    int iSysHF = 0;
//    switch(iSysType){
//    case sysType::JESup:           iSysHF=1; break;
//    case sysType::JESdown:         iSysHF=2; break;
//    case sysType::CSVLFup:         iSysHF=3; break;
//    case sysType::CSVLFdown:       iSysHF=4; break;
//    case sysType::CSVHFStats1up:   iSysHF=5; break;
//    case sysType::CSVHFStats1down: iSysHF=6; break;
//    case sysType::CSVHFStats2up:   iSysHF=7; break;
//    case sysType::CSVHFStats2down: iSysHF=8; break;
//    default : iSysHF = 0; break;
//    }
  
    int iSysC = 0;
//    switch(iSysType){
//    case sysType::CSVCErr1up:   iSysC=1; break;
//    case sysType::CSVCErr1down: iSysC=2; break;
//    case sysType::CSVCErr2up:   iSysC=3; break;
//    case sysType::CSVCErr2down: iSysC=4; break;
//    default : iSysC = 0; break;
//    }
  
    int iSysLF = 0;
//    switch(iSysType){
//    case sysType::JESup:           iSysLF=1; break;
//    case sysType::JESdown:         iSysLF=2; break;
//    case sysType::CSVHFup:         iSysLF=3; break;
//    case sysType::CSVHFdown:       iSysLF=4; break;
//    case sysType::CSVLFStats1up:   iSysLF=5; break;
//    case sysType::CSVLFStats1down: iSysLF=6; break;
//    case sysType::CSVLFStats2up:   iSysLF=7; break;
//    case sysType::CSVLFStats2down: iSysLF=8; break;
//    default : iSysLF = 0; break;
//    }
  
    double csvWgthf = 1.;
    double csvWgtC  = 1.;
    double csvWgtlf = 1.;

   std::auto_ptr<pat::JetCollection> outJets(new pat::JetCollection);
   edm::Handle<pat::JetCollection> recoJets;

   if(iEvent.getByLabel(src_,recoJets)){
    for(unsigned int i=0;i!=recoJets->size();++i){
    
    pat::Jet jet = recoJets->at(i);
    double csv = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
    double jetPt = jet.pt();
    double jetAbsEta = fabs( jet.eta() );
    int flavor = abs( jet.partonFlavour() );

    int iPt = -1; int iEta = -1;
    if (jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100 && jetPt<160) iPt = 4;
    else if (jetPt >=160 && jetPt<10000) iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6) iEta = 1;
    else if ( jetAbsEta>=1.6) iEta = 2;
    //else if ( jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
    }
    else {
      if (iPt >=3) iPt=3;       /// [>60] 
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
    }
  
     double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;
     jet.addUserFloat("CSVreweight",csvWgtTotal);
     outJets->push_back(jet);
    }
   iEvent.put(outJets);
   }
  }
  // ----- member data -----
  edm::InputTag src_;
};
