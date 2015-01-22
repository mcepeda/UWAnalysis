/********************************
Instantiates histoFiller.C with
the correct chain and naming
for a definable set of parameters

Usage: 
root -l -b
.L callHistoFiller.cc
callHistoFiller()

root -l -b -q callHistoFiller.cc

Author: T.M.Perry UW-Madison
********************************/

#include <TStyle.h>
#include <vector>
#include <TChain.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include <iostream>

void callHistoFiller()
{
 int error = 0;
 gROOT->ProcessLine(".L histoFiller.C++", &error);
 if(error!=0){std::cout<<"ERROR LOADING histoFiller.C"<<std::endl;}
 histoFiller m;

 TChain *theChain        = new TChain("muEleEventTree/eventTree");
 TChain *theChainJESUp   = new TChain("muEleEventTreeJetUp/eventTree");
 TChain *theChainJESDown = new TChain("muEleEventTreeJetDown/eventTree");
 TChain *theChainLESUp   = new TChain("muEleEventTreeMuonUp/eventTree");
 TChain *theChainLESDown = new TChain("muEleEventTreeMuonDown/eventTree");
 TTree* tree;
 TString dirName;
 
 UInt_t lumiEle = 19757;
 UInt_t lumiMu  = 19778;

 // define the sample criteria by filling vectors
 std::vector<TString> Samples;
 std::vector<TString> EMu;
 std::vector<TString> MC;
 std::vector<TString> Shifts;
 Bool_t isMu_;
 Bool_t isMC_;
 UInt_t lumi;
 UInt_t nrEvents;
 Float_t crossSec;
 MC.push_back("MC");
 //MC.push_back("Data"); // don't activate this, save data for another loop
 //Samples.push_back("TTbar_full");
 Samples.push_back("TTbar_semi");
 //Samples.push_back("T_s");
 //Samples.push_back("T_t");
 //Samples.push_back("T_tW");
 //Samples.push_back("Tbar_s");
 //Samples.push_back("Tbar_t");
 //Samples.push_back("Tbar_tW");
 //Samples.push_back("WJets");
 //Samples.push_back("W1Jet");
 //Samples.push_back("W2Jet");
 //Samples.push_back("W3Jet");
 //Samples.push_back("W4Jet");
 //Samples.push_back("Wbb4F");
 //Samples.push_back("WW");
 //Samples.push_back("WZ");
 //Samples.push_back("ZZ"); 
 EMu.push_back("ele");
 //EMu.push_back("mu");
 //Shifts.push_back("aucun");
 Shifts.push_back("SFs");
 Shifts.push_back("UESUp");
 Shifts.push_back("UESDown");
 Shifts.push_back("JESUp");
 Shifts.push_back("JESDown");
 Shifts.push_back("LESUp");
 Shifts.push_back("LESDown");
 for (std::vector<TString>::iterator mc = MC.begin(); mc != MC.end(); ++mc){
  std::cout<<"xxxxxx Starting New Loop: MC.mc = "<<*mc<<" xxxxxx"<<std::endl;
  if(mc->EqualTo("MC")) {isMC_ = kTRUE;} 
  else if(mc->EqualTo("Data")) {isMC_ = kFALSE;} 
  else {std::cout<<*mc<<" is your value for 'mc'"<<std::endl; break;}

  for (std::vector<TString>::iterator sample = Samples.begin(); sample != Samples.end(); ++sample){
   std::cout<<" xxxxxx Starting New Loop: Samples.sample = "<<*sample<<" xxxxxx"<<std::endl;
   theChain->Reset();
   theChainJESUp   ->Reset();
   theChainJESDown ->Reset();
   theChainLESUp   ->Reset();
   theChainLESDown ->Reset();
   if(sample->EqualTo("TTbar_full")){
    nrEvents = 10429331;
    crossSec = 25.257; //240.6*0.324*0.324
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_full-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_full-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_full-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_full-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_full-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("TTbar_semi")){
    nrEvents = 23503785;
    crossSec = 105.394; //240.6*0.324*0.676*2
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/*.root");
    std::cout<<"added TTbar_semi to chain"<<std::endl;
    //theChain->Add("/hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/MuEle-PATMC-patTuple_cfg-00277FF2-7B84-E211-9475-782BCB27B958.root");
   }
   if(sample->EqualTo("T_s")){
    nrEvents = 259961;
    crossSec = 3.79;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_T_s-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_T_s-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_T_s-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_T_s-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_T_s-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("T_t")){
    nrEvents = 3915598;
    crossSec = 18.44; //56.4*0.327
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_T_t-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_T_t-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_T_t-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_T_t-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_T_t-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("T_tW")){
    nrEvents = 476173;
    crossSec = 11.1;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_T_tW-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_T_tW-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_T_tW-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_T_tW-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_T_tW-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("Tbar_s")){
    nrEvents = 139974;
    crossSec = 1.76;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_s-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_s-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_s-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_s-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_s-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("Tbar_t")){
    nrEvents = 1711403;
    crossSec = 10.04;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_t-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_t-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_t-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_t-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_t-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("Tbar_tW")){
    nrEvents = 484502;
    crossSec = 11.1;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_tW-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_tW-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_tW-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_tW-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Tbar_tW-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("W1Jet")){
    nrEvents = 1234567; // doesn't matter - hardcoded by jType
    crossSec = 36703.2;
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p1-MuEle-PATMC_V7A/*root");
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p2-MuEle-PATMC_V19/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p1-MuEle-PATMC_V7A/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p2-MuEle-PATMC_V19/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p1-MuEle-PATMC_V7A/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p2-MuEle-PATMC_V19/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p1-MuEle-PATMC_V7A/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p2-MuEle-PATMC_V19/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p1-MuEle-PATMC_V7A/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W1Jet_p2-MuEle-PATMC_V19/*root");
   }
   if(sample->EqualTo("W2Jet")){
    nrEvents = 1234567; // doesn't matter - hardcoded by jType
    crossSec = 36703.2;
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p1-MuEle-PATMC_V7A/*root");
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p2-MuEle-PATMC_V19/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p1-MuEle-PATMC_V7A/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p2-MuEle-PATMC_V19/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p1-MuEle-PATMC_V7A/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p2-MuEle-PATMC_V19/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p1-MuEle-PATMC_V7A/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p2-MuEle-PATMC_V19/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p1-MuEle-PATMC_V7A/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W2Jet_p2-MuEle-PATMC_V19/*root");
   }
   if(sample->EqualTo("W3Jet")){
    nrEvents = 1234567; // doesn't matter - hardcoded by jType
    crossSec = 36703.2;
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p1-MuEle-PATMC_V7A/*root");
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p2-MuEle-PATMC_V19/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p1-MuEle-PATMC_V7A/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p2-MuEle-PATMC_V19/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p1-MuEle-PATMC_V7A/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p2-MuEle-PATMC_V19/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p1-MuEle-PATMC_V7A/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p2-MuEle-PATMC_V19/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p1-MuEle-PATMC_V7A/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W3Jet_p2-MuEle-PATMC_V19/*root");
   }
   if(sample->EqualTo("W4Jet")){
    nrEvents = 1234567; // doesn't matter - hardcoded by jType
    crossSec = 36703.2;
    theChain->Add("/hdfs/store/user/tperry/NeuYear_W4Jet-MuEle-PATMC_V7A/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W4Jet-MuEle-PATMC_V7A/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_W4Jet-MuEle-PATMC_V7A/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_W4Jet-MuEle-PATMC_V7A/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_W4Jet-MuEle-PATMC_V7A/*root");
   }
   if(sample->EqualTo("WJets")){
    nrEvents = 1234567; // doesn't matter - hardcoded by jType
    crossSec = 36703.2;
    theChain->Add("/hdfs/store/user/tperry/NeuYear_WJets_p1-MuEle-PATMC_V7A/*root");
    theChain->Add("/hdfs/store/user/tperry/NeuYear_WJets_p2-MuEle-PATMC_V7A/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p1-MuEle-PATMC_V7A/*root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p2-MuEle-PATMC_V7A/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p1-MuEle-PATMC_V7A/*root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p2-MuEle-PATMC_V7A/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p1-MuEle-PATMC_V7A/*root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p2-MuEle-PATMC_V7A/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p1-MuEle-PATMC_V7A/*root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NeuYear_WJets_p2-MuEle-PATMC_V7A/*root");
   }
   if(sample->EqualTo("Wbb4F")){
    nrEvents = 20294311;
    crossSec = 138.9; //46.3*3
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("Wbb4F")){
    nrEvents = 20294311;
    crossSec = 138.9; //46.3*3
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("WW")){
    nrEvents = 10000431;
    crossSec = 56.75;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_WW-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_WW-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_WW-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_WW-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_WW-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("WZ")){
    nrEvents = 9979363;
    crossSec = 33.21;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
   }
   if(sample->EqualTo("ZZ")){
    nrEvents = 9799908;
    crossSec = 8.26;
    theChain->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainJESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainJESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainLESUp   ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
    theChainLESDown ->Add("/hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*.root");
   }

   for (std::vector<TString>::iterator emu = EMu.begin(); emu != EMu.end(); ++emu){
    std::cout<<"  xxxxxx Starting New Loop: EMu.emu = "<<*emu<<" xxxxxx"<<std::endl;
    if(emu->EqualTo("ele")){isMu_ = kFALSE; lumi = lumiEle;}
    else if(emu->EqualTo("mu")) {isMu_ = kTRUE; lumi = lumiMu;}
    else {std::cout<<*emu<<" is your value for 'emu'"<<std::endl; break;}

    for (std::vector<TString>::iterator shift = Shifts.begin(); shift != Shifts.end(); ++shift){
     std::cout<<"   xxxxxx Starting New Loop: Shifts.shift = "<<*shift<<" xxxxxx"<<std::endl;
     //if     (shift->EqualTo("UESUp"))   { dirName="muEleEventTree"; }
     //else if(shift->EqualTo("UESDown")) { dirName="muEleEventTree"; }
     //else if(shift->EqualTo("JESUp"))   { dirName="muEleEventTreeJetUp"; }
     //else if(shift->EqualTo("JESDown")) { dirName="muEleEventTreeJetDown"; }
     //else if(shift->EqualTo("LESUp"))   { dirName="muEleEventTreeMuonUp"; }
     //else if(shift->EqualTo("LESDown")) { dirName="muEleEventTreeMuonDown"; }
     //else if(shift->EqualTo("SFs"))     { dirName="muEleEventTree"; }
     //else                               { dirName="muEleEventTree"; }
     
     //std::cout<<"theChain->GetEntries() "<<theChain->GetEntries()<<std::endl;
     //std::cout<<"theChain->GetFile()->GetPath() "<<theChain->GetFile()->GetPath()<<std::endl;

     //TDirectory * dir = (TDirectory*)theChain->Get( theChain->GetFile()->GetPath()+dirName );
     ////TDirectory * dir = (TDirectory*)theChain->GetFile()->Get( theChain->GetFile()->GetPath()+dirName );
     ////TDirectory * dir = (TDirectory*)theChain->GetFile()->Get(dirName);
     //dir->GetObject("eventTree",tree);
     //std::cout<<dir->GetName()<<tree->GetName()<<std::endl;
     //m.Init(tree, isMC_);         

     //  std::vector<TString> MC;      mc  {MC, Data}
     //  std::vector<TString> Samples; sample
     //   {TTbar_full, TTbar_semi, 
     //    T_s, T_t, T_tW, 
     //    Tbar_s, Tbar_t, Tbar_tW, 
     //    WJets, W1Jet, W2Jet, W3jet, W4Jet
     //    Wbb4F, WW, WZ, ZZ}
     //  std::vector<TString> EMu;     emu {ele, mu}
     //  std::vector<TString> Shifts;  shift 
     //   {aucun, SFs, UESUp, UESDown, JetUp, JetDown, LepUp, LepDown}

     if     (shift->EqualTo("UESUp"))   { m.Init(theChain, isMC_); }
     else if(shift->EqualTo("UESDown")) { m.Init(theChain, isMC_); }
     else if(shift->EqualTo("JESUp"))   { m.Init(theChainJESUp, isMC_); }
     else if(shift->EqualTo("JESDown")) { m.Init(theChainJESDown, isMC_); }
     else if(shift->EqualTo("LESUp"))   { m.Init(theChainLESUp, isMC_); }
     else if(shift->EqualTo("LESDown")) { m.Init(theChainLESDown, isMC_); }
     else if(shift->EqualTo("SFs"))     { m.Init(theChain, isMC_); }
     else                               { m.Init(theChain, isMC_); }
     
     TString name = *sample+"_"+*emu+*mc+"_"+*shift;
     std::cout<<"     "<<name<<std::endl;
     if(sample->EqualTo("WJets")
      || sample->EqualTo("W1Jet")
      || sample->EqualTo("W2Jet")
      || sample->EqualTo("W3Jet") 
      || sample->EqualTo("W4Jet") ){
      std::cout<<"     Inside W+Jets loop"<<std::endl;
      m.Loop("../plots/"+name+"Wbb", isMu_, isMC_, *shift, 4, lumi, nrEvents, crossSec);
      m.Loop("../plots/"+name+"Wcc", isMu_, isMC_, *shift, 3, lumi, nrEvents, crossSec);
      m.Loop("../plots/"+name+"Wc",  isMu_, isMC_, *shift, 2, lumi, nrEvents, crossSec);
      m.Loop("../plots/"+name+"Wl",  isMu_, isMC_, *shift, 1, lumi, nrEvents, crossSec);
     }
     m.Loop("../plots/"+name, isMu_, isMC_, *shift, 0, lumi, nrEvents, crossSec);
    }
   }
  }
 }
}

//      m.Loop("../plots/"+name, isMu_, isMC_, isQCD_, *cut, *shift, jType, lumi, nrEvents, crossSec);
//
//      chain_TTbar_full->Reset();
//      chain_TTbar_full_JetUp->Reset(); 
//      chain_TTbar_full_JetDown->Reset(); 
//      chain_TTbar_full_LepUp->Reset(); 
//      chain_TTbar_full_LepDown->Reset(); 
//
//      chain_TTbar_semi->Reset();           
//      chain_TTbar_semi_JetUp->Reset();     
//      chain_TTbar_semi_JetDown->Reset();   
//      chain_TTbar_semi_LepUp->Reset();     
//      chain_TTbar_semi_LepDown->Reset();   
//
//      chain_T_s->Reset();
//      chain_T_s_JetUp->Reset(); 
//      chain_T_s_JetDown->Reset(); 
//      chain_T_s_LepUp->Reset(); 
//      chain_T_s_LepDown->Reset(); 
//
//      chain_T_t->Reset();
//      chain_T_t_JetUp->Reset(); 
//      chain_T_t_JetDown->Reset(); 
//      chain_T_t_LepUp->Reset(); 
//      chain_T_t_LepDown->Reset(); 
//
//      chain_T_tW->Reset();
//      chain_T_tW_JetUp->Reset(); 
//      chain_T_tW_JetDown->Reset(); 
//      chain_T_tW_LepUp->Reset(); 
//      chain_T_tW_LepDown->Reset(); 
//
//      chain_Tbar_s->Reset();
//      chain_Tbar_s_JetUp->Reset(); 
//      chain_Tbar_s_JetDown->Reset(); 
//      chain_Tbar_s_LepUp->Reset(); 
//      chain_Tbar_s_LepDown->Reset(); 
//
//      chain_Tbar_t->Reset();
//      chain_Tbar_t_JetUp->Reset(); 
//      chain_Tbar_t_JetDown->Reset(); 
//      chain_Tbar_t_LepUp->Reset(); 
//      chain_Tbar_t_LepDown->Reset(); 
//
//      chain_Tbar_tW->Reset();
//      chain_Tbar_tW_JetUp->Reset(); 
//      chain_Tbar_tW_JetDown->Reset(); 
//      chain_Tbar_tW_LepUp->Reset(); 
//      chain_Tbar_tW_LepDown->Reset(); 
//
//      chain_W1Jet->Reset();
//      chain_W1Jet_JetUp->Reset(); 
//      chain_W1Jet_JetDown->Reset(); 
//      chain_W1Jet_LepUp->Reset(); 
//      chain_W1Jet_LepDown->Reset(); 
//
//      chain_W2Jet->Reset();
//      chain_W2Jet_JetUp->Reset(); 
//      chain_W2Jet_JetDown->Reset(); 
//      chain_W2Jet_LepUp->Reset(); 
//      chain_W2Jet_LepDown->Reset(); 
//
//      chain_W3Jet->Reset();
//      chain_W3Jet_JetUp->Reset(); 
//      chain_W3Jet_JetDown->Reset(); 
//      chain_W3Jet_LepUp->Reset(); 
//      chain_W3Jet_LepDown->Reset(); 
//
//      chain_W4Jet->Reset();
//      chain_W4Jet_JetUp->Reset(); 
//      chain_W4Jet_JetDown->Reset(); 
//      chain_W4Jet_LepUp->Reset(); 
//      chain_W4Jet_LepDown->Reset(); 
//
//      chain_WJets->Reset();
//      chain_WJets_JetUp->Reset(); 
//      chain_WJets_JetDown->Reset(); 
//      chain_WJets_LepUp->Reset(); 
//      chain_WJets_LepDown->Reset(); 
//
//      chain_Wbb4F->Reset();
//      chain_Wbb4F_JetUp->Reset(); 
//      chain_Wbb4F_JetDown->Reset(); 
//      chain_Wbb4F_LepUp->Reset(); 
//      chain_Wbb4F_LepDown->Reset(); 
//
//      chain_WW->Reset();
//      chain_WW_JetUp->Reset(); 
//      chain_WW_JetDown->Reset(); 
//      chain_WW_LepUp->Reset(); 
//      chain_WW_LepDown->Reset(); 
//
//      chain_WZ->Reset();
//      chain_WZ_JetUp->Reset(); 
//      chain_WZ_JetDown->Reset(); 
//      chain_WZ_LepUp->Reset(); 
//      chain_WZ_LepDown->Reset(); 
//
//      chain_ZZ->Reset();
//      chain_ZZ_JetUp->Reset(); 
//      chain_ZZ_JetDown->Reset(); 
//      chain_ZZ_LepUp->Reset(); 
//      chain_ZZ_LepDown->Reset(); 

// Samples.clear();
// Samples.push_back("Data");// [32]
// MC.clear();
// MC.push_back("Data");
// Shifts.clear();
// Shifts.push_back("aucun");
// for (std::vector<TString>::iterator mc = MC.begin(); mc != MC.end(); ++mc){
//  std::cout<<"xxxxxx Starting New Loop: MC.mc = "<<*mc<<" xxxxxx"<<std::endl;
//  for (std::vector<TString>::iterator cut = Cuts.begin(); cut != Cuts.end(); ++cut){
//   std::cout<<"xxxxxx Starting New Loop: Cuts.cut = "<<*cut<<" xxxxxx"<<std::endl;
//    for (std::vector<TString>::iterator sample = Samples.begin(); sample != Samples.end(); ++sample){
//     std::cout<<"xxxxxx Starting New Loop: Samples.sample = "<<*sample<<" xxxxxx"<<std::endl;
//     for (std::vector<TString>::iterator emu = EMu.begin(); emu != EMu.end(); ++emu){
//      std::cout<<"xxxxxx Starting New Loop: EMu.emu = "<<*emu<<" xxxxxx"<<std::endl;
//      for (std::vector<TString>::iterator shift = Shifts.begin(); shift != Shifts.end(); ++shift){
//       std::cout<<"xxxxxx Starting New Loop: Shifts.shift = "<<*shift<<" xxxxxx"<<std::endl;
//
//        // Name and set isMU_, isMC_, isQCD_ 
//        TString name = *sample+"_"+*emu+*mc+*qcd+"_"+*cut+*shift;
//        std::cout<<name<<std::endl;
//        if(emu->EqualTo("ele")){isMu_ = kFALSE; lumi = 19757;}
//        else if(emu->EqualTo("mu")) {isMu_ = kTRUE; lumi = 19778;}
//        else {std::cout<<*emu<<" is your value for 'emu'"<<std::endl; break;}
//        if(mc->EqualTo("MC")) {isMC_ = kTRUE;} 
//        else if(mc->EqualTo("Data")) {isMC_ = kFALSE;} 
//        else {std::cout<<*mc<<" is your value for 'mc'"<<std::endl; break;}
//        if(qcd->EqualTo("good")) {isQCD_ = kFALSE;}
//        else if(qcd->EqualTo("qcd")) {isQCD_ = kTRUE;}
//        else {std::cout<<*qcd<<" is your value for 'qcd'"<<std::endl; break;}
//        jType = 0; // for the timebeing, only change for WJets samples
//
//        //  initialize the correct chain
//        if(emu->EqualTo("mu")){
//         nrEvents = 1;
//         crossSec = 1;
//         chain_DataMu->Add("/hdfs/store/user/tperry/NouvelleYear_DataA_8TeVMu-Mu-PATData/*root");
//         chain_DataMu->Add("/hdfs/store/user/tperry/NouvelleYear_DataB_8TeVMu-Mu-PATData/*root");
//         chain_DataMu->Add("/hdfs/store/user/tperry/NouvelleYear_DataC_8TeVMu-Mu-PATData/*root");
//         chain_DataMu->Add("/hdfs/store/user/tperry/NouvelleYear_DataD_8TeVMu-Mu-PATData/*root");
//         m.Init(chain_DataMu,isMC_); 
//        }
//
//        if(emu->EqualTo("ele")){
//         nrEvents = 1;
//         crossSec = 1;
//         chain_DataEle->Add("/hdfs/store/user/tperry/NouvelleYear_DataA_8TeVEle-Ele-PATData/*root");
//         chain_DataEle->Add("/hdfs/store/user/tperry/NouvelleYear_DataB_8TeVEle-Ele-PATData/*root");
//         chain_DataEle->Add("/hdfs/store/user/tperry/NouvelleYear_DataC_8TeVEle-Ele-PATData/*root");
//         chain_DataEle->Add("/hdfs/store/user/tperry/NouvelleYear_DataD_8TeVEle-Ele-PATData/*root");
//         m.Init(chain_DataEle,isMC_); 
//        }
//
//        m.Loop("../plots/"+name, isMu_, isMC_, isQCD_, *cut, *shift, jType, lumi, nrEvents, crossSec);
//
//        chain_DataEle->Reset();
//        chain_DataMu->Reset();
//      }
//     }
//    }
//   }
//  }
// }
//}
