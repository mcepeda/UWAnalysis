#ifndef UWAnalysis_RecoTools_CompositePtrCandidateTMEtAlgorithm_h
#define UWAnalysis_RecoTools_CompositePtrCandidateTMEtAlgorithm_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

#include "UWAnalysis/DataFormats/interface/CompositePtrCandidateTMEt.h"

#include "UWAnalysis/RecoTools/interface/candidateAuxFunctions.h"
#include "UWAnalysis/RecoTools/interface/METCalibrator.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h" 
#include "DataFormats/Candidate/interface/Candidate.h" 
#include "TLorentzVector.h"
#include "ZSV/BAnalysis/interface/SimBHadron.h"

#include <TMath.h>
#include "RooHist.h"
#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"

template<typename T>
class CompositePtrCandidateTMEtAlgorithm 
{
  typedef edm::Ptr<T> TPtr;
  typedef edm::Ptr<pat::MET> METPtr;
  typedef edm::Ptr<pat::Jet> JetPtr;
  typedef std::vector<JetPtr> JetPtrVector;
  typedef std::vector<SimBHadron> bhadronPtr;

 public:

  CompositePtrCandidateTMEtAlgorithm()
  {
    verbosity_ = 0;
  }
  ~CompositePtrCandidateTMEtAlgorithm() {}

  void setMETCalibrator(METCalibrator * calibrator) {calibrator_ = calibrator;}

  CompositePtrCandidateTMEt<T> buildCompositePtrCandidate(const TPtr visDecayProducts, 
							  METPtr met,
							  JetPtrVector pfJets,
							  edm::View<T> visProductCollection,
							  const reco::GenParticleCollection* genParticles,
    							  bhadronPtr bhadrons)
  {
    CompositePtrCandidateTMEt<T> compositePtrCandidate(visDecayProducts, met);
  
    if ( visDecayProducts.isNull() ) {
      edm::LogError ("CompositePtrCandidateTMEtAlgorithm") << " Pointer to visible Decay products invalid !!";
      return compositePtrCandidate;
    }

    if ( met.isNull() ) {
      edm::LogError ("CompositePtrCandidateTMEtAlgorithm") << " Pointer to missing transverse momentum invalid !!";
      return compositePtrCandidate;
    }
    //METPtr met1= met->userCand("type1");
    //candidatePtr met1= met->userCand("type1");
    //check here for correct met
    reco::Candidate::LorentzVector recoil(-met->px()-visDecayProducts->px(),-met->py()-visDecayProducts->py(),0.0,sqrt(pow(met->px()+visDecayProducts->px(),2)+pow(met->py()+visDecayProducts->py(),2)));

    compositePtrCandidate.setRecoil(recoil);
    compositePtrCandidate.setRecoilDPhi(deltaPhi(recoil.phi(),visDecayProducts->phi()));
    compositePtrCandidate.setCharge(visDecayProducts->charge());
    compositePtrCandidate.setMt(compMt(visDecayProducts->p4(), met->userCand("type1")->px(), met->userCand("type1")->py()));
    compositePtrCandidate.setDPhi(TMath::Abs(normalizedPhi(visDecayProducts->phi() - met->phi())));
    compositePtrCandidate.setPx(visDecayProducts->px() + met->px());
    compositePtrCandidate.setPy(visDecayProducts->py() + met->py());
    compositePtrCandidate.setPt(sqrt( (visDecayProducts->py() + met->py())*(visDecayProducts->py() + met->py()) + (visDecayProducts->px() + met->px())*(visDecayProducts->px() + met->px())));
    computeMETPz(compositePtrCandidate, met->px(), met->py());
    reco::Candidate::LorentzVector correctedMET = met->userCand("type1")->p4();

    if(calibrator_!=0)
      correctedMET = calibrator_->calibrate(met->p4(),visDecayProducts->p4(),reco::Candidate::LorentzVector(0.,0.,0.000001,0.00000001),genParticles);

    compositePtrCandidate.setCalibratedMET(correctedMET);
    //printf("correctedMET: %f\n",correctedMET.pt());
    compositePtrCandidate.setCorMt(compMt(visDecayProducts->p4(), met->userCand("type1")->px(), met->userCand("type1")->py()));
    compositePtrCandidate.setCorPt(sqrt( (visDecayProducts->py() + correctedMET.py())*(visDecayProducts->py() + correctedMET.py()) + (visDecayProducts->px() + correctedMET.px())*(visDecayProducts->px() + correctedMET.px())));
    compositePtrCandidate.setCorPx(visDecayProducts->px() + correctedMET.px());
    compositePtrCandidate.setCorPy(visDecayProducts->py() + correctedMET.py());
    //Jets
    JetPtrVector cleanedJets;
    JetPtrVector cleanedJetsAll;
    JetPtrVector cleanedJets5;

    int nJets24Pt20 = 0;
    int nJets24Pt25 = 0;
    int nJets24Pt30 = 0;
//    double etaJ1_allEta = -777;
//    double etaJ2_allEta = -777;
//    double etaJ3_allEta = -777;
//    double etaJ4_allEta = -777;
//    double phiJ1_allEta = -777;
//    double phiJ2_allEta = -777;
//    double phiJ3_allEta = -777;
//    double phiJ4_allEta = -777;
//    double ptJ1_allEta = -777;
//    double ptJ2_allEta = -777;
//    double ptJ3_allEta = -777;
//    double ptJ4_allEta = -777;


    for(unsigned int i=0;i<pfJets.size();++i){
//      if((reco::deltaR(pfJets.at(i)->p4(),visDecayProducts->p4())>0.5)&&fabs(pfJets.at(i)->userFloat("idLoose"))>0&&fabs(pfJets.at(i)->eta())<4.5){
//	cleanedJetsAll.push_back(pfJets.at(i));
//        if (i == 0){
//         etaJ1_allEta = pfJets.at(i)->eta();
//         phiJ1_allEta = pfJets.at(i)->phi();
//         ptJ1_allEta  = pfJets.at(i)->pt() ;
//        } 
//        if (i == 1){
//         etaJ2_allEta = pfJets.at(i)->eta();
//         phiJ2_allEta = pfJets.at(i)->phi();
//         ptJ2_allEta  = pfJets.at(i)->pt() ;
//        } 
//        if (i == 2){
//         etaJ3_allEta = pfJets.at(i)->eta();
//         phiJ3_allEta = pfJets.at(i)->phi();
//         ptJ3_allEta  = pfJets.at(i)->pt() ;
//        } 
//        if (i == 3){
//         etaJ4_allEta = pfJets.at(i)->eta();
//         phiJ4_allEta = pfJets.at(i)->phi();
//         ptJ4_allEta  = pfJets.at(i)->pt() ;
//        } 
//
//      }
      if(2>1){
      //if((reco::deltaR(pfJets.at(i)->p4(),visDecayProducts->p4())>0.5)&&fabs(pfJets.at(i)->userFloat("idLoose"))>0&&fabs(pfJets.at(i)->eta())<4.5&&pfJets.at(i)->pt()>20){
      //if((reco::deltaR(pfJets.at(i)->p4(),visDecayProducts->p4())>0.5)&&fabs(pfJets.at(i)->userFloat("idLoose"))>0&&fabs(pfJets.at(i)->eta())<2.4){
	cleanedJets.push_back(pfJets.at(i));
      }
      if((reco::deltaR(pfJets.at(i)->p4(),visDecayProducts->p4())>0.5)&&fabs(pfJets.at(i)->userFloat("idLoose"))>0&&fabs(pfJets.at(i)->eta())>2.4&&fabs(pfJets.at(i)->eta())<4.5&&pfJets.at(i)->pt()>20){
	cleanedJets5.push_back(pfJets.at(i));
	nJets24Pt20+=1;
	if(pfJets.at(i)->pt()>25){
	  nJets24Pt25+=1;
	  if(pfJets.at(i)->pt()>30){
	    nJets24Pt30+=1 ;
	  }
	}
      }
    }
    compositePtrCandidate.setJetValuesNonTracker(nJets24Pt20,nJets24Pt25,nJets24Pt30);

    //sort them by Pt
    sortRefVectorByPt(cleanedJets);
    double SV1SV2M;
    double SV1SV2Pt;
    double massZ;
    if(cleanedJets.size()>1){
      if(cleanedJets.at(0)->userFloat("SV_pt")>0&&cleanedJets.at(1)->userFloat("SV_pt")>0){
	math::PtEtaPhiMLorentzVector SV1(cleanedJets.at(0)->userFloat("SV_pt"),
					 cleanedJets.at(0)->userFloat("SV_eta"),
					 cleanedJets.at(0)->userFloat("SV_phi"),
					 cleanedJets.at(0)->userFloat("SV_M")
					 );
	math::PtEtaPhiMLorentzVector SV2(cleanedJets.at(1)->userFloat("SV_pt"),
					 cleanedJets.at(1)->userFloat("SV_eta"),
					 cleanedJets.at(1)->userFloat("SV_phi"),
					 cleanedJets.at(1)->userFloat("SV_M")
					 );
	SV1SV2M = (SV1+SV2).M();
	SV1SV2Pt = (SV1+SV2).pt();
	compositePtrCandidate.setSV1SV2M(SV1SV2M);
	compositePtrCandidate.setSV1SV2Pt(SV1SV2Pt);
      }
    }

    if(cleanedJets.size()>1)
      if(cleanedJets.at(0)->userFloat("massBpm_charge")!=0&&cleanedJets.at(1)->userFloat("massBpm_charge")!=0&&cleanedJets.at(0)->userFloat("muon1ZMass")<60&&cleanedJets.at(1)->userFloat("muon1ZMass")<60){
	math::PtEtaPhiMLorentzVector B1(cleanedJets.at(0)->userFloat("massBpm_pt"),
					cleanedJets.at(0)->userFloat("massBpm_eta"),
					cleanedJets.at(0)->userFloat("massBpm_phi"),
					cleanedJets.at(0)->userFloat("massBpm")
					);

	math::PtEtaPhiMLorentzVector B2(cleanedJets.at(0)->userFloat("massBpm_pt"),
					cleanedJets.at(0)->userFloat("massBpm_eta"),
					cleanedJets.at(0)->userFloat("massBpm_phi"),
					cleanedJets.at(0)->userFloat("massBpm")
					);
	massZ = (B1+B2).M();
	compositePtrCandidate.setMassZ(massZ);///fix here
      }

    unsigned int nJets = cleanedJets.size();

    unsigned int nJets5 = cleanedJets5.size();

    if(nJets>3)
      compositePtrCandidate.setJJVariables((cleanedJets.at(0)->p4()+cleanedJets.at(1)->p4()).M(),(cleanedJets.at(0)->p4()+cleanedJets.at(1)->p4()).pt(),(cleanedJets.at(2)->p4()+cleanedJets.at(3)->p4()).M());
    else if(nJets>1) 
      compositePtrCandidate.setJJVariables((cleanedJets.at(0)->p4()+cleanedJets.at(1)->p4()).M(),(cleanedJets.at(0)->p4()+cleanedJets.at(1)->p4()).pt(),0.0);
    else
      compositePtrCandidate.setJJVariables(0.0,0.0,0.0);

    bool TOP = false;    
    //    printf("bhadrons size %i\n", (int)bhadrons.size());
    if(bhadrons.size()>0)    
      if(bhadrons.at(0).motherids.size()>0){
	//printf("mothers size %i\n", (int)bhadrons.at(0).motherids.size());
	float mother1 = bhadrons.at(0).motherids.at(0);
	float mother2 = 0;
	if(bhadrons.size()>1)
	  if(bhadrons.at(1).motherids.size()>0)
	    bhadrons.at(1).motherids.at(0);

	//printf("compo Mother ID: %f\n",mother1);
	if(abs(mother1)==6||abs(mother2)==6)
	  TOP = true;
      }
    
    CSVL(compositePtrCandidate,cleanedJets,TOP);
    CSVM(compositePtrCandidate,cleanedJets,TOP);
    CSVT(compositePtrCandidate,cleanedJets,TOP);
    SSVHE(compositePtrCandidate,cleanedJets,TOP);

   //find the nearest jet to leg1
    double ht = visDecayProducts->pt();
    for(unsigned int k=0;k<nJets;++k){
      if(cleanedJets.at(k)->pt()>25)
	ht+=cleanedJets.at(k)->pt();

    }
    compositePtrCandidate.setJetValues(cleanedJets,nJets5,ht); //sets cleaned jets as  default jets()

    efficiencies(compositePtrCandidate);

    float metjj=0;
    float leptonjj=0;

    if(cleanedJets.size()>1) {
      leptonjj=(visDecayProducts->p4()
		+cleanedJets.at(0)->p4()
		+cleanedJets.at(1)->p4()
		).M();
       metjj=(correctedMET
	      +cleanedJets.at(0)->p4()
	      +cleanedJets.at(1)->p4()
	      ).M();
    }

    compositePtrCandidate.setSSVar(metjj,leptonjj);

    //--- compute gen. level quantities
    if ( genParticles ) {
      compGenQuantities(compositePtrCandidate, genParticles);
      compTrueMet(compositePtrCandidate, genParticles);
      compositePtrCandidate.setGenMt(compMt(compositePtrCandidate.p4Leptongen(), compositePtrCandidate.trueMEx(),compositePtrCandidate.trueMEy()));
    }

    return compositePtrCandidate;
  }

 private: 
  void efficiencies(CompositePtrCandidateTMEt<T>& compositePtrCandidate){

    // Muon ID Eficiency   DATA_over_MC_Tight_eta_pt20-500  for each eta, mult by a factor
    edm::FileInPath fileID("UWAnalysis/Configuration/data/MuonEfficiencies_Run2012ReReco_53X.root");
    // Muon ISO Efficiency  DATA_over_MC_combRelIsoPF04dBeta<012_Tight_eta_pt20-500
    edm::FileInPath fileIS("UWAnalysis/Configuration/data/MuonEfficiencies_ISO_Run_2012ReReco_53X.root");
    // Muon Trigger Efficiency  IsoMu24_eta2p1_DATA_over_MC_TightID_ETA_pt25-500_2012ABCD
    edm::FileInPath fileTR("UWAnalysis/Configuration/data/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root");
    
    TFile *fID = new TFile(fileID.fullPath().c_str());
    TFile *fIS = new TFile(fileIS.fullPath().c_str());
    TFile *fTR = new TFile(fileTR.fullPath().c_str());

    double eta = compositePtrCandidate.lepton()->eta();
    float weight_etaID = 1.0;
    float weight_etaIS = 1.0;
    float weight_etaTR = 1.0;

    TGraphAsymmErrors *gID = (TGraphAsymmErrors*)fID->Get("DATA_over_MC_Tight_eta_pt20-500");
    TGraphAsymmErrors *gIS = (TGraphAsymmErrors*)fIS->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_eta_pt20-500");
    TGraphAsymmErrors *gTR = (TGraphAsymmErrors*)fTR->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_ETA_pt25-500_2012ABCD");

    weight_etaID = (float) gID->Eval(eta);
    weight_etaIS = (float) gIS->Eval(eta);
    weight_etaTR = (float) gTR->Eval(eta);

    compositePtrCandidate.setEffWEIGHTeta_ID(weight_etaID);
    compositePtrCandidate.setEffWEIGHTeta_IS(weight_etaIS);
    compositePtrCandidate.setEffWEIGHTeta_TR(weight_etaTR);

    fID->Close();
    fIS->Close();
    fTR->Close();
    return;      
  }
  

  void CSVL(CompositePtrCandidateTMEt<T>& compositePtrCandidate,const JetPtrVector& jets,bool TOP){
    float j1pt =0.0;
    float j2pt =0.0;
    float j1parton =0.0;
    float j2parton =0.0;
    float j1eta =0.0;
    float j2eta =0.0;
    float j1csv =0.0;
    float j2csv =0.0;
    float SFCSVL1Jet = 0.0;
    float SFCSVL2Jet = 0.0;
    float SFCSVL1JetT = 0.0;
    float SFCSVL2JetT = 0.0;
    float effCSV1 = 0.0;
    float effCSV2 = 0.0;
    float fakeCSV1 = 0.0;
    float fakeCSV2 = 0.0;
    float value_ = 0.242;
    if(jets.size()>1){
      j1pt = jets.at(0)->pt();
      j2pt = jets.at(1)->pt();
      //printf("J1pt %f J2pt %f\n",j1pt,j2pt);
      j1eta = jets.at(0)->eta();
      j2eta = jets.at(1)->eta();
      //printf("J1eta %f J2eta %f\n",j1eta,j2eta);	    
      j1csv = jets.at(0)->bDiscriminator("combinedSecondaryVertexBJetTags");
      j2csv = jets.at(1)->bDiscriminator("combinedSecondaryVertexBJetTags");
      //printf("J1csv %f J2csv %f\n",j1csv,j2csv);	    
      j1parton = jets.at(0)->partonFlavour();
      j2parton = jets.at(1)->partonFlavour();
      if(j1csv>value_){
	if(fabs(j1parton)==5||fabs(j1parton)==4){
	  effCSV1 = 0.997942*((1.+(0.00923753*j1pt))/(1.+(0.0096119*j1pt)));
	  SFCSVL1Jet = effCSV1;
	  SFCSVL1JetT = 1.01; // this will need to be updated 
	}
	else{
	  fakeCSV1 = LightLookupCSVL(j1pt, j1eta);
	  SFCSVL1Jet = fakeCSV1;
	  SFCSVL1JetT = fakeCSV1;
	}
      }
      if(j2csv>value_){
	if(fabs(j2parton)==5||fabs(j2parton)==4){
	  effCSV2 = 0.997942*((1.+(0.00923753*j2pt))/(1.+(0.0096119*j2pt))); 
	  if(effCSV1>0){SFCSVL1Jet = (effCSV1+effCSV2)/2; SFCSVL2Jet = effCSV1*effCSV2;  SFCSVL1JetT = 1.01;  SFCSVL2JetT = 1.01*1.01;}
	  else if(fakeCSV1>0){SFCSVL1Jet = fakeCSV1*effCSV2; SFCSVL2Jet = fakeCSV1*effCSV2;SFCSVL1JetT = 1.01*fakeCSV1;SFCSVL2JetT=SFCSVL1JetT;}
	  else {SFCSVL1Jet = effCSV2;  SFCSVL1JetT = 1.01;}
	}
	else { 
	  fakeCSV2 = LightLookupCSVL(j2pt, j2eta);
	  if(effCSV1>0){ SFCSVL1Jet =(effCSV1*fakeCSV2); SFCSVL2Jet = fakeCSV2*effCSV1;  SFCSVL1JetT = 1.01*fakeCSV2; SFCSVL2JetT = 1.01*fakeCSV2;}
	  else if(fakeCSV1>0){SFCSVL1Jet =(fakeCSV1+fakeCSV2)/2; SFCSVL2Jet = fakeCSV1*fakeCSV2; SFCSVL2JetT = SFCSVL2Jet; SFCSVL1JetT = SFCSVL1Jet;}
	  else {SFCSVL1Jet = fakeCSV2;  SFCSVL1JetT = fakeCSV2;}
	  
	}
	
      }
    }

    if(TOP){
      compositePtrCandidate.setSFCSVL1(SFCSVL1JetT);
      compositePtrCandidate.setSFCSVL2(SFCSVL2JetT);
    }
    else{
      compositePtrCandidate.setSFCSVL1(SFCSVL1Jet);
      compositePtrCandidate.setSFCSVL2(SFCSVL2Jet);
    }
  }

  void CSVM(CompositePtrCandidateTMEt<T>& compositePtrCandidate,const JetPtrVector& jets, bool TOP){
    float j1pt =0.0;
    float j2pt =0.0;
    float j1parton =0.0;
    float j2parton =0.0;
    float j1eta =0.0;
    float j2eta =0.0;
    float j1csv =0.0;
    float j2csv =0.0;
    float SFCSVM1Jet = 0.0;
    float SFCSVM2Jet = 0.0;
    float SFCSVM1JetT = 0.0;
    float SFCSVM2JetT = 0.0;
    float effCSV1 = 0.0;
    float effCSV2 = 0.0;
    float fakeCSV1 = 0.0;
    float fakeCSV2 = 0.0;
    float value_ = 0.679;
    if(jets.size()>1){
      j1pt = jets.at(0)->pt();
      j2pt = jets.at(1)->pt();
      j1eta = jets.at(0)->eta();
      j2eta = jets.at(1)->eta();
      j1csv = jets.at(0)->bDiscriminator("combinedSecondaryVertexBJetTags");
      j2csv = jets.at(1)->bDiscriminator("combinedSecondaryVertexBJetTags");
      j1parton = jets.at(0)->partonFlavour();
      j2parton = jets.at(1)->partonFlavour();

      effCSV1 = 0;effCSV2 = 0;
      fakeCSV1 = 0; fakeCSV2 = 0; 
      SFCSVM1Jet =1;
      SFCSVM2Jet =1;
      //do loop for 1 Btagged Loose x = 0.679;
      float x = 0.679;  //isn't this (threshold) the same as value_
      if(j1csv>x){
	if(fabs(j1parton)==5||fabs(j1parton)==4){
	  if(TOP){
	    effCSV1 = 0.97; // will need to be updated
	  }
	  else{
	    effCSV1 = (0.938887+(0.00017124*j1pt))+(-2.76366e-07*(j1pt*j1pt));
	  }
	}
	else{
	  fakeCSV1 = LightLookupCSVM(j1pt, j1eta);
	}
	  }
      if(j2csv>x){      
	if(fabs(j2parton)==5||fabs(j2parton)==4){
	  if(TOP){
	    effCSV2 = 0.97;
	    if(effCSV1>0){SFCSVM1Jet =(effCSV1+effCSV2)/2; SFCSVM2Jet = effCSV1*effCSV2;}
	    else if(fakeCSV1>0){SFCSVM1Jet =fakeCSV1*effCSV2; SFCSVM2Jet = fakeCSV1*effCSV2;}
	    else SFCSVM1Jet = effCSV2;
	  }
	  else{
	    effCSV2 = (0.938887+(0.00017124*j2pt))+(-2.76366e-07*(j2pt*j2pt)); 
	    if(effCSV1>0){SFCSVM1Jet = (effCSV1+effCSV2)/2; SFCSVM2Jet = effCSV1*effCSV2;}
	    else if(fakeCSV1>0){SFCSVM1Jet = fakeCSV1*effCSV2; SFCSVM2Jet = fakeCSV1*effCSV2;}
	    else SFCSVM1Jet = effCSV2;
	  }
	}
	else { 
	  fakeCSV2 = LightLookupCSVM(j2pt, j2eta);
	  if(effCSV1>0){ SFCSVM1Jet =(effCSV1*fakeCSV2); SFCSVM2Jet = fakeCSV2*effCSV1;}
	  else if(fakeCSV1>0){SFCSVM1Jet =(fakeCSV1+fakeCSV2)/2; SFCSVM2Jet = fakeCSV1*fakeCSV2;}
	  else {SFCSVM1Jet = fakeCSV2;}
	  
	}
      }
      else if(fakeCSV1>0) SFCSVM1Jet = fakeCSV1;
      else if(effCSV1>0) SFCSVM1Jet = effCSV1;
    }

    compositePtrCandidate.setSFCSVM1(SFCSVM1Jet);
    compositePtrCandidate.setSFCSVM2(SFCSVM2Jet);      
    
    //if(j1csv>value_&&j2csv>value_)
    //printf("SF CSVM2 Top: %f\n",SFCSVM2JetT);
    //compositePtrCandidate.setSFCSVM1(SFCSVM1Jet);
    //compositePtrCandidate.setSFCSVM2(SFCSVM2Jet);
    //compositePtrCandidate.setSFCSVM1T(SFCSVM1JetT);
    //compositePtrCandidate.setSFCSVM2T(SFCSVM2JetT);
  }


  void CSVT(CompositePtrCandidateTMEt<T>& compositePtrCandidate,const JetPtrVector& jets, bool TOP){
    float j1pt =0.0;
    float j2pt =0.0;
    float j1parton =0.0;
    float j2parton =0.0;
    float j1eta =0.0;
    float j2eta =0.0;
    float j1csv =0.0;
    float j2csv =0.0;
    float SFCSVT1Jet = 0.0;
    float SFCSVT2Jet = 0.0;
    float SFCSVT1JetT = 0.0;
    float SFCSVT2JetT = 0.0;
    float effCSV1 = 0.0;
    float effCSV2 = 0.0;
    float fakeCSV1 = 0.0;
    float fakeCSV2 = 0.0;
    float value_ = 0.898;
    if(jets.size()>1){
      j1pt = jets.at(0)->pt();
      j2pt = jets.at(1)->pt();
      //printf("J1pt %f J2pt %f\n",j1pt,j2pt);
      j1eta = jets.at(0)->eta();
      j2eta = jets.at(1)->eta();
      //printf("J1eta %f J2eta %f\n",j1eta,j2eta);	    
      j1csv = jets.at(0)->bDiscriminator("combinedSecondaryVertexBJetTags");
      j2csv = jets.at(1)->bDiscriminator("combinedSecondaryVertexBJetTags");
      //printf("J1csv %f J2csv %f\n",j1csv,j2csv);	    
      j1parton = jets.at(0)->partonFlavour();
      j2parton = jets.at(1)->partonFlavour();
      if(j1csv>value_){
	if(fabs(j1parton)==5||fabs(j1parton)==4){
	  effCSV1 = (0.927563+(1.55479e-05*j1pt))+(-1.90666e-07*(j1pt*j1pt)); 
	  SFCSVT1Jet = effCSV1;
	  SFCSVT1JetT = 0.96;
	}
	else{
	  fakeCSV1 = LightLookupCSVT(j1pt, j1eta);
	  SFCSVT1Jet = fakeCSV1;
	  SFCSVT1JetT = fakeCSV1;
	}
      }
      if(j2csv>value_){
	if(fabs(j2parton)==5||fabs(j2parton)==4){
	  effCSV2 = (0.927563+(1.55479e-05*j2pt))+(-1.90666e-07*(j2pt*j2pt)); 
	  if(effCSV1>0){SFCSVT1Jet =(effCSV1+effCSV2)/2;SFCSVT1JetT = 0.96; SFCSVT2Jet = effCSV1*effCSV2;SFCSVT2JetT = 0.96*0.96;}
	  else if(fakeCSV1>0){SFCSVT1Jet =fakeCSV1*effCSV2; SFCSVT2Jet = fakeCSV1*effCSV2;SFCSVT1JetT = 0.96*fakeCSV2;SFCSVT2JetT = 0.96*fakeCSV2;}
	  else {SFCSVT1Jet = effCSV2; SFCSVT1JetT = 0.96;}
	}
	else{
	  fakeCSV2 = LightLookupCSVT(j2pt, j2eta);
	  if(effCSV1>0){SFCSVT1Jet =(effCSV1*fakeCSV2); SFCSVT2Jet = fakeCSV2*effCSV1; SFCSVT1JetT = 0.96*fakeCSV2; SFCSVT2JetT = SFCSVT1JetT;}
	  else if(fakeCSV1>0){SFCSVT1Jet =(fakeCSV1+fakeCSV2)/2; SFCSVT2Jet = fakeCSV1*fakeCSV2; SFCSVT1JetT=SFCSVT1Jet; SFCSVT2JetT=SFCSVT2Jet;}
	  else {SFCSVT1Jet = fakeCSV2;  SFCSVT1JetT = fakeCSV2;}
	}
      }
    }
    if(TOP){
    compositePtrCandidate.setSFCSVT1(SFCSVT1JetT);
    compositePtrCandidate.setSFCSVT2(SFCSVT2JetT);
    }
    else{
      compositePtrCandidate.setSFCSVT1(SFCSVT1Jet);
      compositePtrCandidate.setSFCSVT2(SFCSVT2Jet);
    }
  }
  

  void SSVHE(CompositePtrCandidateTMEt<T>& compositePtrCandidate,const JetPtrVector& jets, bool TOP){
    float j1pt =0.0;
    float j2pt =0.0;
    float j1parton =0.0;
    float j2parton =0.0;
    float j1eta =0.0;
    float j2eta =0.0;
    float j1ssv =0.0;
    float j2ssv =0.0;
    float SFSSVHEM1Jet = 0.0;
    float SFSSVHEM2Jet = 0.0;
    float SFSSVHEM1JetT = 0.0;
    float SFSSVHEM2JetT = 0.0;
    float effCSV1 = 0.0;
    float effCSV2 = 0.0;
    float fakeCSV1 = 0.0;
    float fakeCSV2 = 0.0;
    float value_ =1.74;
    if(jets.size()>1){
      j1pt = jets.at(0)->pt();
      j2pt = jets.at(1)->pt();
      //printf("J1pt %f J2pt %f\n",j1pt,j2pt);
      j1eta = jets.at(0)->eta();
      j2eta = jets.at(1)->eta();
      //printf("J1eta %f J2eta %f\n",j1eta,j2eta);	    
      j1ssv = jets.at(0)->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      j2ssv = jets.at(1)->bDiscriminator("simpeSecondaryVertexHighEffBJetTags");
      //printf("J1csv %f J2csv %f\n",j1ssv,j2ssv);	    
      j1parton = jets.at(0)->partonFlavour();
      j2parton = jets.at(1)->partonFlavour();
      if(j1ssv>value_){
	if(fabs(j1parton)==5||fabs(j1parton)==4){
	  effCSV1 = 0.896462*((1.+(0.00957275*j1pt))/(1.+(0.00837582*j1pt)));
	  SFSSVHEM1Jet = effCSV1;
	  SFSSVHEM1JetT= 0.96;
	}
	else{
	  fakeCSV1 = LightLookupSSVHEM(j1pt, j1eta);
	  SFSSVHEM1Jet = fakeCSV1;
	  SFSSVHEM1JetT= fakeCSV1;
	}
      }
      if(j2ssv>value_){
	if(fabs(j2parton)==5||fabs(j2parton)==4){
	  effCSV2 =0.896462*((1.+(0.00957275*j2pt))/(1.+(0.00837582*j2pt)));
	  if(effCSV1>0){SFSSVHEM1Jet = (effCSV1+effCSV2)/2; SFSSVHEM2Jet = effCSV1*effCSV2;SFSSVHEM2JetT =0.96*0.96;SFSSVHEM1JetT =0.96;}
	  else if(fakeCSV1>0){SFSSVHEM1Jet = fakeCSV1*effCSV2; SFSSVHEM2Jet = fakeCSV1*effCSV2;SFSSVHEM2JetT =0.96*fakeCSV1;SFSSVHEM1JetT =0.96*fakeCSV2;}
	  else {SFSSVHEM1Jet = effCSV2;SFSSVHEM1JetT=0.96;}
	}
	else { 
	  fakeCSV2 = LightLookupSSVHEM(j2pt, j2eta);
	  if(effCSV1>0){ SFSSVHEM1Jet =(effCSV1*fakeCSV2); SFSSVHEM2Jet = fakeCSV2*effCSV1;SFSSVHEM2JetT=0.96*fakeCSV2;SFSSVHEM1JetT=SFSSVHEM2JetT;}
	  else if(fakeCSV1>0){SFSSVHEM1Jet =(fakeCSV1+fakeCSV2)/2; SFSSVHEM2Jet = fakeCSV1*fakeCSV2;SFSSVHEM2JetT=SFSSVHEM2Jet;SFSSVHEM1JetT=SFSSVHEM1Jet;}
	  else {SFSSVHEM1Jet = fakeCSV2;SFSSVHEM2JetT=fakeCSV2;}
	}
      } 
    }
    if(TOP){
      compositePtrCandidate.setSFSSVHE1(SFSSVHEM1JetT);
      compositePtrCandidate.setSFSSVHE2(SFSSVHEM2JetT);
    }
    else{
      compositePtrCandidate.setSFSSVHE1(SFSSVHEM1Jet);
      compositePtrCandidate.setSFSSVHE2(SFSSVHEM2Jet);
    }
  }

 void computeMETPz(CompositePtrCandidateTMEt<T>& compositePtrCandidate, double MetPx, double MetPy){
	// WARNING: This assumes that the lepton+MET is a W boson
	// For coherence, this should be moved out of this general class (I am lazy right now) 
        double M_W  = 80.4;
	
        double pxmu = compositePtrCandidate.lepton()->px();
        double pymu = compositePtrCandidate.lepton()->py();
        double pzmu = compositePtrCandidate.lepton()->pz();
        double emu = compositePtrCandidate.lepton()->energy();
        double pznu = 0.;
        double otherSol_ = 0.;
        bool isComplex_=false;

        double a = M_W*M_W + 2.0*pxmu*MetPx + 2.0*pymu*MetPy;
        double A = 4.0*(emu*emu - pzmu*pzmu);
        double B = -4.0*a*pzmu;
        double C = 4.0*emu*emu*(MetPx*MetPx+MetPy*MetPy) - a*a;

        double tmproot = B*B - 4.0*A*C;

        if (tmproot<0) {
                isComplex_= true;
                pznu = - B/(2*A); // take real part of complex roots
                otherSol_ = pznu;
                }
        else {
                        isComplex_ = false;
                        double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
                        double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);

                          if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2; otherSol_ = tmpsol1;}
                          else { pznu = tmpsol1; otherSol_ = tmpsol2; }
                                // if pznu is > 300 pick the most central root
                                if ( pznu > 300. ) {
                                  if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pznu = tmpsol1; otherSol_ = tmpsol2; }
                                  else { pznu = tmpsol2; otherSol_ = tmpsol1; }
                                }
                        }
   compositePtrCandidate.setMetPz(pznu);
   compositePtrCandidate.setMetPzB(otherSol_);
   compositePtrCandidate.setIsComplex(isComplex_);
   double PX_TOTAL=pxmu+MetPx;
   double PY_TOTAL=pymu+MetPy;
   double PZ_TOTAL=pzmu+pznu;
   double EN_TOTAL=sqrt(PX_TOTAL*PX_TOTAL+PY_TOTAL*PY_TOTAL+PZ_TOTAL*PZ_TOTAL);
   compositePtrCandidate.setWP4(math::XYZTLorentzVector(PX_TOTAL,PY_TOTAL,PZ_TOTAL,EN_TOTAL));
 } 

  void compGenQuantities(CompositePtrCandidateTMEt<T>& compositePtrCandidate, const reco::GenParticleCollection* genParticles)
  {
    const reco::GenParticle* genLepton = findGenParticle(compositePtrCandidate.lepton()->p4(), *genParticles, 0.5, -1);
    if ( genLepton ) {
      compositePtrCandidate.setP4Leptongen(genLepton->p4());
      compositePtrCandidate.setP4VisLeptongen(getVisMomentum(genLepton, genParticles));
      compositePtrCandidate.setPdg(genLepton->pdgId());
    }
  }     

  void compTrueMet(CompositePtrCandidateTMEt<T>& compositePtrCandidate, const reco::GenParticleCollection* genParticles){
    double trueMEY  = 0.0;
    double trueMEX  = 0.0;;
    double trueMEZ  = 0.0;;
 
    for( unsigned i = 0; i < genParticles->size(); i++ ) {
    const reco::GenParticle& genpart = (*genParticles)[i];
      if( genpart.status() == 1 && fabs(genpart.eta()) < 5.0 ) { 
        if( std::abs(genpart.pdgId()) == 12 ||
            std::abs(genpart.pdgId()) == 14 ||
            std::abs(genpart.pdgId()) == 16  
          //  ||std::abs(genpart.pdgId()) < 7   || 
          //  std::abs(genpart.pdgId()) == 21 
		) {
          trueMEX += genpart.px();
          trueMEY += genpart.py();
          trueMEZ += genpart.pz();

        } 
      }
    }
    double true_met = sqrt( trueMEX*trueMEX + trueMEY*trueMEY );
    double true_phi = atan2(trueMEY,trueMEX);	

   compositePtrCandidate.setTrueMet(true_met);
   compositePtrCandidate.setTrueMetPhi(true_phi);
   compositePtrCandidate.setTrueMEy(trueMEY);
   compositePtrCandidate.setTrueMEx(trueMEX);
   compositePtrCandidate.setTrueMEz(trueMEZ);
  }





  double compMt(const reco::Candidate::LorentzVector& visParticle, 
		double metPx, double metPy)
  {
    double px = visParticle.px() + metPx;
    double py = visParticle.py() + metPy;
    double et = visParticle.Et() + TMath::Sqrt(metPx*metPx + metPy*metPy);
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) {
//      edm::LogWarning ("compMt") << " mt2 = " << mt2 << " must not be negative !!";
      return 0.;
    }
    return TMath::Sqrt(mt2);
  }

  class refVectorPtSorter {
  public:
    refVectorPtSorter(const JetPtrVector vec)
      {
	vec_ = vec;
      }

    refVectorPtSorter()
      {
      }


    ~refVectorPtSorter()
      {}

    bool operator()(size_t a , size_t b) {
      return (vec_.at(a)->pt() > vec_.at(b)->pt());
    }
  private:
    JetPtrVector vec_;
  };



  void sortRefVectorByPt(JetPtrVector& vec)
  {
    std::vector<size_t> indices;
    indices.reserve(vec.size());
    for(unsigned int i=0;i<vec.size();++i)
      indices.push_back(i);
    
    refVectorPtSorter sorter(vec);
    std::sort(indices.begin(),indices.end(),sorter);
        
    JetPtrVector sorted;
    sorted.reserve(vec.size());
    
    for(unsigned int i=0;i<indices.size();++i)
      sorted.push_back(vec.at(indices.at(i)));

    vec = sorted;
  }


  float LightLookupCSVM(float pt, float eta){
    if(fabs(eta)>0&&fabs(eta)<0.8)
      {//mean
	return ((1.07541+(0.00231827*pt))+(-4.74249e-06*(pt*pt)))+(2.70862e-09*(pt*(pt*pt)));
      }
    else if(fabs(eta)>0.8&&fabs(eta)<1.6)
      {//mean
      return ((1.05613+(0.00114031*pt))+(-2.56066e-06*(pt*pt)))+(1.67792e-09*(pt*(pt*pt)));
      }
    else if( fabs(eta)>1.6 && fabs(eta)<2.4)
      {//mean
      return ((1.05625+(0.000487231*pt))+(-2.22792e-06*(pt*pt)))+(1.70262e-09*(pt*(pt*pt))); 
    }
    return 0;
  }
  
  float LightLookupCSVL(float pt, float eta){
    
    if(fabs(eta)>0.0&&fabs(eta)<0.5){//mean
      
    return (float) ((1.01177+(0.0023066*pt))+(-4.56052e-06*(pt*pt)))+(2.57917e-09*(pt*(pt*pt)));
    }
    else if( fabs(eta)>0.5 && fabs(eta)<1.0){//mean
      
      return (float) ((0.975966+(0.00196354*pt))+(-3.83768e-06*(pt*pt)))+(2.17466e-09*(pt*(pt*pt)));
    }
    else if( fabs(eta)>1.0 && fabs(eta)<1.5){//mean
      
      return (float) ((0.93821+(0.00180935*pt))+(-3.86937e-06*(pt*pt)))+(2.43222e-09*(pt*(pt*pt)));
    }
    else if( fabs(eta)>1.5 && fabs(eta)<2.4){//mean
      
      return (float) ((1.00022+(0.0010998*pt))+(-3.10672e-06*(pt*pt)))+(2.35006e-09*(pt*(pt*pt)));
    }
    
    return 0;
  }
  
  float LightLookupCSVT(float pt, float eta){
    if(fabs(eta) >0.0 && fabs(eta)<2.4){
      return (float) ((1.00462+(0.00325971*pt))+(-7.79184e-06*(pt*pt)))+(5.22506e-09*(pt*(pt*pt)));
    }
    return 0;
  }
  
  float LightLookupSSVHEM(float pt, float eta){
    if( fabs(eta)>0.0 && fabs(eta)<0.8){
      return (float)((0.86318+(0.000801639*pt))+(-1.64119e-06*(pt*pt)))+(2.59121e-10*(pt*(pt*pt)));
    }
    else if( fabs(eta)>0.8 && fabs(eta)<1.6){
      return (float)((0.958973+(-0.000269555*pt))+(1.381e-06*(pt*pt)))+(-1.87744e-09*(pt*(pt*pt)));
    }
    else if( fabs(eta)>1.6 && fabs(eta)<2.4){
      return (float)((0.923033+(-0.000898227*pt))+(4.74565e-06*(pt*pt)))+(-6.11053e-09*(pt*(pt*pt)));
    }
    return 0;   
  }
  
  float LightLookupSSVHPT(float pt, float eta){
    
    return (float) ((0.97409+(0.000646241*pt))+(-2.86294e-06*(pt*pt)))+(2.79484e-09*(pt*(pt*pt)));
  }


  METCalibrator *calibrator_;

  
  int verbosity_;
  };

#endif 

