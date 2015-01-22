#define histoFiller_cxx
#include "histoFiller.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <fstream>

void histoFiller::Loop(TString outfileName, Bool_t isMu, Bool_t isMC,
      TString shift, Int_t jType, UInt_t luminosity, UInt_t nrEvents,
      Float_t crossSec)
{
 if (fChain == 0) return;
 std::cout<<"about to fChain->GetEntries();"<<std::endl;
 //Long64_t nrEntries = fChain->GetEntries();
 Long64_t nrEntries = fChain->GetEntriesFast();
 Long64_t nb = 0;

 std::cout<<"initializing wbb histograms"<<std::endl;
 // wbb
 TH1F *h_wbb_good_goodJ1_pt   = new TH1F("h_wbb_good_goodJ1_pt","",12,25,205); h_wbb_good_goodJ1_pt->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_pt   = new TH1F("h_wbb_good_goodJ2_pt","",12,25,205); h_wbb_good_goodJ2_pt->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_eta  = new TH1F("h_wbb_good_goodJ1_eta","",20,-3,5); h_wbb_good_goodJ1_eta->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_eta  = new TH1F("h_wbb_good_goodJ2_eta","",20,-3,5); h_wbb_good_goodJ2_eta->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_phi  = new TH1F("h_wbb_good_goodJ1_phi","",35,-3.4033,5.7594); h_wbb_good_goodJ1_phi->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_phi  = new TH1F("h_wbb_good_goodJ2_phi","",35,-3.4033,5.7594); h_wbb_good_goodJ2_phi->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_CSV  = new TH1F("h_wbb_good_goodJ1_CSV","",30,0,1); h_wbb_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_CSV  = new TH1F("h_wbb_good_goodJ2_CSV","",30,0,1); h_wbb_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_wbb_good_dRgoodJ1J2  = new TH1F("h_wbb_good_dRgoodJ1J2","",20,0,6); h_wbb_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_wbb_good_mt          = new TH1F("h_wbb_good_mt","",40,0,200); h_wbb_good_mt->Sumw2(); 
 TH1F *h_wbb_good_goodLep_pt  = new TH1F("h_wbb_good_goodLep_pt","",12,25,205); h_wbb_good_goodLep_pt->Sumw2(); 
 TH1F *h_wbb_good_goodLep_eta = new TH1F("h_wbb_good_goodLep_eta","",20,-3,5); h_wbb_good_goodLep_eta->Sumw2(); 
 TH1F *h_wbb_good_goodLep_phi = new TH1F("h_wbb_good_goodLep_phi","",35,-3.4033,5.7594); h_wbb_good_goodLep_phi->Sumw2(); 

 TH1F *h_wbb_good_goodJ1_pt_CSVUp   = new TH1F("h_wbb_good_goodJ1_pt_CSVUp","",12,25,205); h_wbb_good_goodJ1_pt_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_pt_CSVUp   = new TH1F("h_wbb_good_goodJ2_pt_CSVUp","",12,25,205); h_wbb_good_goodJ2_pt_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_eta_CSVUp  = new TH1F("h_wbb_good_goodJ1_eta_CSVUp","",20,-3,5); h_wbb_good_goodJ1_eta_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_eta_CSVUp  = new TH1F("h_wbb_good_goodJ2_eta_CSVUp","",20,-3,5); h_wbb_good_goodJ2_eta_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_phi_CSVUp  = new TH1F("h_wbb_good_goodJ1_phi_CSVUp","",35,-3.4033,5.7594); h_wbb_good_goodJ1_phi_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_phi_CSVUp  = new TH1F("h_wbb_good_goodJ2_phi_CSVUp","",35,-3.4033,5.7594); h_wbb_good_goodJ2_phi_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_CSV_CSVUp  = new TH1F("h_wbb_good_goodJ1_CSV_CSVUp","",30,0,1); h_wbb_good_goodJ1_CSV_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_CSV_CSVUp  = new TH1F("h_wbb_good_goodJ2_CSV_CSVUp","",30,0,1); h_wbb_good_goodJ2_CSV_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_dRgoodJ1J2_CSVUp  = new TH1F("h_wbb_good_dRgoodJ1J2_CSVUp","",20,0,6); h_wbb_good_dRgoodJ1J2_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_mt_CSVUp          = new TH1F("h_wbb_good_mt_CSVUp","",40,0,200); h_wbb_good_mt_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodLep_pt_CSVUp  = new TH1F("h_wbb_good_goodLep_pt_CSVUp","",12,25,205); h_wbb_good_goodLep_pt_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodLep_eta_CSVUp = new TH1F("h_wbb_good_goodLep_eta_CSVUp","",20,-3,5); h_wbb_good_goodLep_eta_CSVUp->Sumw2(); 
 TH1F *h_wbb_good_goodLep_phi_CSVUp = new TH1F("h_wbb_good_goodLep_phi_CSVUp","",35,-3.4033,5.7594); h_wbb_good_goodLep_phi_CSVUp->Sumw2(); 

 TH1F *h_wbb_good_goodJ1_pt_CSVDown   = new TH1F("h_wbb_good_goodJ1_pt_CSVDown","",12,25,205); h_wbb_good_goodJ1_pt_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_pt_CSVDown   = new TH1F("h_wbb_good_goodJ2_pt_CSVDown","",12,25,205); h_wbb_good_goodJ2_pt_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_eta_CSVDown  = new TH1F("h_wbb_good_goodJ1_eta_CSVDown","",20,-3,5); h_wbb_good_goodJ1_eta_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_eta_CSVDown  = new TH1F("h_wbb_good_goodJ2_eta_CSVDown","",20,-3,5); h_wbb_good_goodJ2_eta_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_phi_CSVDown  = new TH1F("h_wbb_good_goodJ1_phi_CSVDown","",35,-3.4033,5.7594); h_wbb_good_goodJ1_phi_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_phi_CSVDown  = new TH1F("h_wbb_good_goodJ2_phi_CSVDown","",35,-3.4033,5.7594); h_wbb_good_goodJ2_phi_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_CSV_CSVDown  = new TH1F("h_wbb_good_goodJ1_CSV_CSVDown","",30,0,1); h_wbb_good_goodJ1_CSV_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_CSV_CSVDown  = new TH1F("h_wbb_good_goodJ2_CSV_CSVDown","",30,0,1); h_wbb_good_goodJ2_CSV_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_dRgoodJ1J2_CSVDown  = new TH1F("h_wbb_good_dRgoodJ1J2_CSVDown","",20,0,6); h_wbb_good_dRgoodJ1J2_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_mt_CSVDown          = new TH1F("h_wbb_good_mt_CSVDown","",40,0,200); h_wbb_good_mt_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodLep_pt_CSVDown  = new TH1F("h_wbb_good_goodLep_pt_CSVDown","",12,25,205); h_wbb_good_goodLep_pt_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodLep_eta_CSVDown = new TH1F("h_wbb_good_goodLep_eta_CSVDown","",20,-3,5); h_wbb_good_goodLep_eta_CSVDown->Sumw2(); 
 TH1F *h_wbb_good_goodLep_phi_CSVDown = new TH1F("h_wbb_good_goodLep_phi_CSVDown","",35,-3.4033,5.7594); h_wbb_good_goodLep_phi_CSVDown->Sumw2(); 

 TH1F *h_wbb_good_goodJ1_pt_EMuUp   = new TH1F("h_wbb_good_goodJ1_pt_EMuUp","",12,25,205); h_wbb_good_goodJ1_pt_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_pt_EMuUp   = new TH1F("h_wbb_good_goodJ2_pt_EMuUp","",12,25,205); h_wbb_good_goodJ2_pt_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_eta_EMuUp  = new TH1F("h_wbb_good_goodJ1_eta_EMuUp","",20,-3,5); h_wbb_good_goodJ1_eta_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_eta_EMuUp  = new TH1F("h_wbb_good_goodJ2_eta_EMuUp","",20,-3,5); h_wbb_good_goodJ2_eta_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_phi_EMuUp  = new TH1F("h_wbb_good_goodJ1_phi_EMuUp","",35,-3.4033,5.7594); h_wbb_good_goodJ1_phi_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_phi_EMuUp  = new TH1F("h_wbb_good_goodJ2_phi_EMuUp","",35,-3.4033,5.7594); h_wbb_good_goodJ2_phi_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_CSV_EMuUp  = new TH1F("h_wbb_good_goodJ1_CSV_EMuUp","",30,0,1); h_wbb_good_goodJ1_CSV_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_CSV_EMuUp  = new TH1F("h_wbb_good_goodJ2_CSV_EMuUp","",30,0,1); h_wbb_good_goodJ2_CSV_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_dRgoodJ1J2_EMuUp  = new TH1F("h_wbb_good_dRgoodJ1J2_EMuUp","",20,0,6); h_wbb_good_dRgoodJ1J2_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_mt_EMuUp          = new TH1F("h_wbb_good_mt_EMuUp","",40,0,200); h_wbb_good_mt_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodLep_pt_EMuUp  = new TH1F("h_wbb_good_goodLep_pt_EMuUp","",12,25,205); h_wbb_good_goodLep_pt_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodLep_eta_EMuUp = new TH1F("h_wbb_good_goodLep_eta_EMuUp","",20,-3,5); h_wbb_good_goodLep_eta_EMuUp->Sumw2(); 
 TH1F *h_wbb_good_goodLep_phi_EMuUp = new TH1F("h_wbb_good_goodLep_phi_EMuUp","",35,-3.4033,5.7594); h_wbb_good_goodLep_phi_EMuUp->Sumw2(); 

 TH1F *h_wbb_good_goodJ1_pt_EMuDown   = new TH1F("h_wbb_good_goodJ1_pt_EMuDown","",12,25,205); h_wbb_good_goodJ1_pt_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_pt_EMuDown   = new TH1F("h_wbb_good_goodJ2_pt_EMuDown","",12,25,205); h_wbb_good_goodJ2_pt_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_eta_EMuDown  = new TH1F("h_wbb_good_goodJ1_eta_EMuDown","",20,-3,5); h_wbb_good_goodJ1_eta_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_eta_EMuDown  = new TH1F("h_wbb_good_goodJ2_eta_EMuDown","",20,-3,5); h_wbb_good_goodJ2_eta_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_phi_EMuDown  = new TH1F("h_wbb_good_goodJ1_phi_EMuDown","",35,-3.4033,5.7594); h_wbb_good_goodJ1_phi_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_phi_EMuDown  = new TH1F("h_wbb_good_goodJ2_phi_EMuDown","",35,-3.4033,5.7594); h_wbb_good_goodJ2_phi_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ1_CSV_EMuDown  = new TH1F("h_wbb_good_goodJ1_CSV_EMuDown","",30,0,1); h_wbb_good_goodJ1_CSV_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodJ2_CSV_EMuDown  = new TH1F("h_wbb_good_goodJ2_CSV_EMuDown","",30,0,1); h_wbb_good_goodJ2_CSV_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_dRgoodJ1J2_EMuDown  = new TH1F("h_wbb_good_dRgoodJ1J2_EMuDown","",20,0,6); h_wbb_good_dRgoodJ1J2_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_mt_EMuDown          = new TH1F("h_wbb_good_mt_EMuDown","",40,0,200); h_wbb_good_mt_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodLep_pt_EMuDown  = new TH1F("h_wbb_good_goodLep_pt_EMuDown","",12,25,205); h_wbb_good_goodLep_pt_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodLep_eta_EMuDown = new TH1F("h_wbb_good_goodLep_eta_EMuDown","",20,-3,5); h_wbb_good_goodLep_eta_EMuDown->Sumw2(); 
 TH1F *h_wbb_good_goodLep_phi_EMuDown = new TH1F("h_wbb_good_goodLep_phi_EMuDown","",35,-3.4033,5.7594); h_wbb_good_goodLep_phi_EMuDown->Sumw2(); 

 TH1F *h_wbb_qcd_goodJ1_pt   = new TH1F("h_wbb_qcd_goodJ1_pt","",12,25,205); h_wbb_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_pt   = new TH1F("h_wbb_qcd_goodJ2_pt","",12,25,205); h_wbb_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_eta  = new TH1F("h_wbb_qcd_goodJ1_eta","",20,-3,5); h_wbb_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_eta  = new TH1F("h_wbb_qcd_goodJ2_eta","",20,-3,5); h_wbb_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_phi  = new TH1F("h_wbb_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_wbb_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_phi  = new TH1F("h_wbb_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_wbb_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_CSV  = new TH1F("h_wbb_qcd_goodJ1_CSV","",30,0,1); h_wbb_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_CSV  = new TH1F("h_wbb_qcd_goodJ2_CSV","",30,0,1); h_wbb_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_wbb_qcd_dRgoodJ1J2  = new TH1F("h_wbb_qcd_dRgoodJ1J2","",20,0,6); h_wbb_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_wbb_qcd_mt          = new TH1F("h_wbb_qcd_mt","",40,0,200); h_wbb_qcd_mt->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_pt  = new TH1F("h_wbb_qcd_qcdLep_pt","",12,25,205); h_wbb_qcd_qcdLep_pt->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_eta = new TH1F("h_wbb_qcd_qcdLep_eta","",20,-3,5); h_wbb_qcd_qcdLep_eta->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_phi = new TH1F("h_wbb_qcd_qcdLep_phi","",35,-3.4033,5.7594); h_wbb_qcd_qcdLep_phi->Sumw2(); 

 TH1F *h_wbb_qcd_goodJ1_pt_CSVUp   = new TH1F("h_wbb_qcd_goodJ1_pt_CSVUp","",12,25,205); h_wbb_qcd_goodJ1_pt_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_pt_CSVUp   = new TH1F("h_wbb_qcd_goodJ2_pt_CSVUp","",12,25,205); h_wbb_qcd_goodJ2_pt_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_eta_CSVUp  = new TH1F("h_wbb_qcd_goodJ1_eta_CSVUp","",20,-3,5); h_wbb_qcd_goodJ1_eta_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_eta_CSVUp  = new TH1F("h_wbb_qcd_goodJ2_eta_CSVUp","",20,-3,5); h_wbb_qcd_goodJ2_eta_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_phi_CSVUp  = new TH1F("h_wbb_qcd_goodJ1_phi_CSVUp","",35,-3.4033,5.7594); h_wbb_qcd_goodJ1_phi_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_phi_CSVUp  = new TH1F("h_wbb_qcd_goodJ2_phi_CSVUp","",35,-3.4033,5.7594); h_wbb_qcd_goodJ2_phi_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_CSV_CSVUp  = new TH1F("h_wbb_qcd_goodJ1_CSV_CSVUp","",30,0,1); h_wbb_qcd_goodJ1_CSV_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_CSV_CSVUp  = new TH1F("h_wbb_qcd_goodJ2_CSV_CSVUp","",30,0,1); h_wbb_qcd_goodJ2_CSV_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_dRgoodJ1J2_CSVUp  = new TH1F("h_wbb_qcd_dRgoodJ1J2_CSVUp","",20,0,6); h_wbb_qcd_dRgoodJ1J2_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_mt_CSVUp          = new TH1F("h_wbb_qcd_mt_CSVUp","",40,0,200); h_wbb_qcd_mt_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_pt_CSVUp  = new TH1F("h_wbb_qcd_qcdLep_pt_CSVUp","",12,25,205); h_wbb_qcd_qcdLep_pt_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_eta_CSVUp = new TH1F("h_wbb_qcd_qcdLep_eta_CSVUp","",20,-3,5); h_wbb_qcd_qcdLep_eta_CSVUp->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_phi_CSVUp = new TH1F("h_wbb_qcd_qcdLep_phi_CSVUp","",35,-3.4033,5.7594); h_wbb_qcd_qcdLep_phi_CSVUp->Sumw2(); 

 TH1F *h_wbb_qcd_goodJ1_pt_CSVDown   = new TH1F("h_wbb_qcd_goodJ1_pt_CSVDown","",12,25,205); h_wbb_qcd_goodJ1_pt_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_pt_CSVDown   = new TH1F("h_wbb_qcd_goodJ2_pt_CSVDown","",12,25,205); h_wbb_qcd_goodJ2_pt_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_eta_CSVDown  = new TH1F("h_wbb_qcd_goodJ1_eta_CSVDown","",20,-3,5); h_wbb_qcd_goodJ1_eta_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_eta_CSVDown  = new TH1F("h_wbb_qcd_goodJ2_eta_CSVDown","",20,-3,5); h_wbb_qcd_goodJ2_eta_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_phi_CSVDown  = new TH1F("h_wbb_qcd_goodJ1_phi_CSVDown","",35,-3.4033,5.7594); h_wbb_qcd_goodJ1_phi_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_phi_CSVDown  = new TH1F("h_wbb_qcd_goodJ2_phi_CSVDown","",35,-3.4033,5.7594); h_wbb_qcd_goodJ2_phi_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_CSV_CSVDown  = new TH1F("h_wbb_qcd_goodJ1_CSV_CSVDown","",30,0,1); h_wbb_qcd_goodJ1_CSV_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_CSV_CSVDown  = new TH1F("h_wbb_qcd_goodJ2_CSV_CSVDown","",30,0,1); h_wbb_qcd_goodJ2_CSV_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_dRgoodJ1J2_CSVDown  = new TH1F("h_wbb_qcd_dRgoodJ1J2_CSVDown","",20,0,6); h_wbb_qcd_dRgoodJ1J2_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_mt_CSVDown          = new TH1F("h_wbb_qcd_mt_CSVDown","",40,0,200); h_wbb_qcd_mt_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_pt_CSVDown  = new TH1F("h_wbb_qcd_qcdLep_pt_CSVDown","",12,25,205); h_wbb_qcd_qcdLep_pt_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_eta_CSVDown = new TH1F("h_wbb_qcd_qcdLep_eta_CSVDown","",20,-3,5); h_wbb_qcd_qcdLep_eta_CSVDown->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_phi_CSVDown = new TH1F("h_wbb_qcd_qcdLep_phi_CSVDown","",35,-3.4033,5.7594); h_wbb_qcd_qcdLep_phi_CSVDown->Sumw2(); 

 TH1F *h_wbb_qcd_goodJ1_pt_EMuUp   = new TH1F("h_wbb_qcd_goodJ1_pt_EMuUp","",12,25,205); h_wbb_qcd_goodJ1_pt_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_pt_EMuUp   = new TH1F("h_wbb_qcd_goodJ2_pt_EMuUp","",12,25,205); h_wbb_qcd_goodJ2_pt_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_eta_EMuUp  = new TH1F("h_wbb_qcd_goodJ1_eta_EMuUp","",20,-3,5); h_wbb_qcd_goodJ1_eta_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_eta_EMuUp  = new TH1F("h_wbb_qcd_goodJ2_eta_EMuUp","",20,-3,5); h_wbb_qcd_goodJ2_eta_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_phi_EMuUp  = new TH1F("h_wbb_qcd_goodJ1_phi_EMuUp","",35,-3.4033,5.7594); h_wbb_qcd_goodJ1_phi_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_phi_EMuUp  = new TH1F("h_wbb_qcd_goodJ2_phi_EMuUp","",35,-3.4033,5.7594); h_wbb_qcd_goodJ2_phi_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_CSV_EMuUp  = new TH1F("h_wbb_qcd_goodJ1_CSV_EMuUp","",30,0,1); h_wbb_qcd_goodJ1_CSV_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_CSV_EMuUp  = new TH1F("h_wbb_qcd_goodJ2_CSV_EMuUp","",30,0,1); h_wbb_qcd_goodJ2_CSV_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_dRgoodJ1J2_EMuUp  = new TH1F("h_wbb_qcd_dRgoodJ1J2_EMuUp","",20,0,6); h_wbb_qcd_dRgoodJ1J2_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_mt_EMuUp          = new TH1F("h_wbb_qcd_mt_EMuUp","",40,0,200); h_wbb_qcd_mt_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_pt_EMuUp  = new TH1F("h_wbb_qcd_qcdLep_pt_EMuUp","",12,25,205); h_wbb_qcd_qcdLep_pt_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_eta_EMuUp = new TH1F("h_wbb_qcd_qcdLep_eta_EMuUp","",20,-3,5); h_wbb_qcd_qcdLep_eta_EMuUp->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_phi_EMuUp = new TH1F("h_wbb_qcd_qcdLep_phi_EMuUp","",35,-3.4033,5.7594); h_wbb_qcd_qcdLep_phi_EMuUp->Sumw2(); 

 TH1F *h_wbb_qcd_goodJ1_pt_EMuDown   = new TH1F("h_wbb_qcd_goodJ1_pt_EMuDown","",12,25,205); h_wbb_qcd_goodJ1_pt_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_pt_EMuDown   = new TH1F("h_wbb_qcd_goodJ2_pt_EMuDown","",12,25,205); h_wbb_qcd_goodJ2_pt_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_eta_EMuDown  = new TH1F("h_wbb_qcd_goodJ1_eta_EMuDown","",20,-3,5); h_wbb_qcd_goodJ1_eta_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_eta_EMuDown  = new TH1F("h_wbb_qcd_goodJ2_eta_EMuDown","",20,-3,5); h_wbb_qcd_goodJ2_eta_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_phi_EMuDown  = new TH1F("h_wbb_qcd_goodJ1_phi_EMuDown","",35,-3.4033,5.7594); h_wbb_qcd_goodJ1_phi_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_phi_EMuDown  = new TH1F("h_wbb_qcd_goodJ2_phi_EMuDown","",35,-3.4033,5.7594); h_wbb_qcd_goodJ2_phi_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ1_CSV_EMuDown  = new TH1F("h_wbb_qcd_goodJ1_CSV_EMuDown","",30,0,1); h_wbb_qcd_goodJ1_CSV_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_goodJ2_CSV_EMuDown  = new TH1F("h_wbb_qcd_goodJ2_CSV_EMuDown","",30,0,1); h_wbb_qcd_goodJ2_CSV_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_dRgoodJ1J2_EMuDown  = new TH1F("h_wbb_qcd_dRgoodJ1J2_EMuDown","",20,0,6); h_wbb_qcd_dRgoodJ1J2_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_mt_EMuDown          = new TH1F("h_wbb_qcd_mt_EMuDown","",40,0,200); h_wbb_qcd_mt_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_pt_EMuDown  = new TH1F("h_wbb_qcd_qcdLep_pt_EMuDown","",12,25,205); h_wbb_qcd_qcdLep_pt_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_eta_EMuDown = new TH1F("h_wbb_qcd_qcdLep_eta_EMuDown","",20,-3,5); h_wbb_qcd_qcdLep_eta_EMuDown->Sumw2(); 
 TH1F *h_wbb_qcd_qcdLep_phi_EMuDown = new TH1F("h_wbb_qcd_qcdLep_phi_EMuDown","",35,-3.4033,5.7594); h_wbb_qcd_qcdLep_phi_EMuDown->Sumw2(); 
                                                
 std::cout<<"initializing tt1m1e histograms"<<std::endl;
 // tt1m1e
 TH1F *h_tt1m1e_good_goodJ1_pt   = new TH1F("h_tt1m1e_good_goodJ1_pt","",12,25,205); h_tt1m1e_good_goodJ1_pt->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_pt   = new TH1F("h_tt1m1e_good_goodJ2_pt","",12,25,205); h_tt1m1e_good_goodJ2_pt->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_eta  = new TH1F("h_tt1m1e_good_goodJ1_eta","",20,-3,5); h_tt1m1e_good_goodJ1_eta->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_eta  = new TH1F("h_tt1m1e_good_goodJ2_eta","",20,-3,5); h_tt1m1e_good_goodJ2_eta->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_phi  = new TH1F("h_tt1m1e_good_goodJ1_phi","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ1_phi->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_phi  = new TH1F("h_tt1m1e_good_goodJ2_phi","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ2_phi->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_CSV  = new TH1F("h_tt1m1e_good_goodJ1_CSV","",30,0,1); h_tt1m1e_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_CSV  = new TH1F("h_tt1m1e_good_goodJ2_CSV","",30,0,1); h_tt1m1e_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_tt1m1e_good_dRgoodJ1J2  = new TH1F("h_tt1m1e_good_dRgoodJ1J2","",20,0,6); h_tt1m1e_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_tt1m1e_good_mt          = new TH1F("h_tt1m1e_good_mt","",40,0,200); h_tt1m1e_good_mt->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_pt  = new TH1F("h_tt1m1e_good_goodLep_pt","",12,25,205); h_tt1m1e_good_goodLep_pt->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_eta = new TH1F("h_tt1m1e_good_goodLep_eta","",20,-3,5); h_tt1m1e_good_goodLep_eta->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_phi = new TH1F("h_tt1m1e_good_goodLep_phi","",35,-3.4033,5.7594); h_tt1m1e_good_goodLep_phi->Sumw2(); 

 TH1F *h_tt1m1e_good_goodJ1_pt_CSVUp   = new TH1F("h_tt1m1e_good_goodJ1_pt_CSVUp","",12,25,205); h_tt1m1e_good_goodJ1_pt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_pt_CSVUp   = new TH1F("h_tt1m1e_good_goodJ2_pt_CSVUp","",12,25,205); h_tt1m1e_good_goodJ2_pt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_eta_CSVUp  = new TH1F("h_tt1m1e_good_goodJ1_eta_CSVUp","",20,-3,5); h_tt1m1e_good_goodJ1_eta_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_eta_CSVUp  = new TH1F("h_tt1m1e_good_goodJ2_eta_CSVUp","",20,-3,5); h_tt1m1e_good_goodJ2_eta_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_phi_CSVUp  = new TH1F("h_tt1m1e_good_goodJ1_phi_CSVUp","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ1_phi_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_phi_CSVUp  = new TH1F("h_tt1m1e_good_goodJ2_phi_CSVUp","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ2_phi_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_CSV_CSVUp  = new TH1F("h_tt1m1e_good_goodJ1_CSV_CSVUp","",30,0,1); h_tt1m1e_good_goodJ1_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_CSV_CSVUp  = new TH1F("h_tt1m1e_good_goodJ2_CSV_CSVUp","",30,0,1); h_tt1m1e_good_goodJ2_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_dRgoodJ1J2_CSVUp  = new TH1F("h_tt1m1e_good_dRgoodJ1J2_CSVUp","",20,0,6); h_tt1m1e_good_dRgoodJ1J2_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_mt_CSVUp          = new TH1F("h_tt1m1e_good_mt_CSVUp","",40,0,200); h_tt1m1e_good_mt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_pt_CSVUp  = new TH1F("h_tt1m1e_good_goodLep_pt_CSVUp","",12,25,205); h_tt1m1e_good_goodLep_pt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_eta_CSVUp = new TH1F("h_tt1m1e_good_goodLep_eta_CSVUp","",20,-3,5); h_tt1m1e_good_goodLep_eta_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_phi_CSVUp = new TH1F("h_tt1m1e_good_goodLep_phi_CSVUp","",35,-3.4033,5.7594); h_tt1m1e_good_goodLep_phi_CSVUp->Sumw2(); 

 TH1F *h_tt1m1e_good_goodJ1_pt_CSVDown   = new TH1F("h_tt1m1e_good_goodJ1_pt_CSVDown","",12,25,205); h_tt1m1e_good_goodJ1_pt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_pt_CSVDown   = new TH1F("h_tt1m1e_good_goodJ2_pt_CSVDown","",12,25,205); h_tt1m1e_good_goodJ2_pt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_eta_CSVDown  = new TH1F("h_tt1m1e_good_goodJ1_eta_CSVDown","",20,-3,5); h_tt1m1e_good_goodJ1_eta_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_eta_CSVDown  = new TH1F("h_tt1m1e_good_goodJ2_eta_CSVDown","",20,-3,5); h_tt1m1e_good_goodJ2_eta_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_phi_CSVDown  = new TH1F("h_tt1m1e_good_goodJ1_phi_CSVDown","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ1_phi_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_phi_CSVDown  = new TH1F("h_tt1m1e_good_goodJ2_phi_CSVDown","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ2_phi_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_CSV_CSVDown  = new TH1F("h_tt1m1e_good_goodJ1_CSV_CSVDown","",30,0,1); h_tt1m1e_good_goodJ1_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_CSV_CSVDown  = new TH1F("h_tt1m1e_good_goodJ2_CSV_CSVDown","",30,0,1); h_tt1m1e_good_goodJ2_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_dRgoodJ1J2_CSVDown  = new TH1F("h_tt1m1e_good_dRgoodJ1J2_CSVDown","",20,0,6); h_tt1m1e_good_dRgoodJ1J2_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_mt_CSVDown          = new TH1F("h_tt1m1e_good_mt_CSVDown","",40,0,200); h_tt1m1e_good_mt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_pt_CSVDown  = new TH1F("h_tt1m1e_good_goodLep_pt_CSVDown","",12,25,205); h_tt1m1e_good_goodLep_pt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_eta_CSVDown = new TH1F("h_tt1m1e_good_goodLep_eta_CSVDown","",20,-3,5); h_tt1m1e_good_goodLep_eta_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_phi_CSVDown = new TH1F("h_tt1m1e_good_goodLep_phi_CSVDown","",35,-3.4033,5.7594); h_tt1m1e_good_goodLep_phi_CSVDown->Sumw2(); 

 TH1F *h_tt1m1e_good_goodJ1_pt_EMuUp   = new TH1F("h_tt1m1e_good_goodJ1_pt_EMuUp","",12,25,205); h_tt1m1e_good_goodJ1_pt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_pt_EMuUp   = new TH1F("h_tt1m1e_good_goodJ2_pt_EMuUp","",12,25,205); h_tt1m1e_good_goodJ2_pt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_eta_EMuUp  = new TH1F("h_tt1m1e_good_goodJ1_eta_EMuUp","",20,-3,5); h_tt1m1e_good_goodJ1_eta_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_eta_EMuUp  = new TH1F("h_tt1m1e_good_goodJ2_eta_EMuUp","",20,-3,5); h_tt1m1e_good_goodJ2_eta_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_phi_EMuUp  = new TH1F("h_tt1m1e_good_goodJ1_phi_EMuUp","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ1_phi_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_phi_EMuUp  = new TH1F("h_tt1m1e_good_goodJ2_phi_EMuUp","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ2_phi_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_CSV_EMuUp  = new TH1F("h_tt1m1e_good_goodJ1_CSV_EMuUp","",30,0,1); h_tt1m1e_good_goodJ1_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_CSV_EMuUp  = new TH1F("h_tt1m1e_good_goodJ2_CSV_EMuUp","",30,0,1); h_tt1m1e_good_goodJ2_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_dRgoodJ1J2_EMuUp  = new TH1F("h_tt1m1e_good_dRgoodJ1J2_EMuUp","",20,0,6); h_tt1m1e_good_dRgoodJ1J2_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_mt_EMuUp          = new TH1F("h_tt1m1e_good_mt_EMuUp","",40,0,200); h_tt1m1e_good_mt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_pt_EMuUp  = new TH1F("h_tt1m1e_good_goodLep_pt_EMuUp","",12,25,205); h_tt1m1e_good_goodLep_pt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_eta_EMuUp = new TH1F("h_tt1m1e_good_goodLep_eta_EMuUp","",20,-3,5); h_tt1m1e_good_goodLep_eta_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_phi_EMuUp = new TH1F("h_tt1m1e_good_goodLep_phi_EMuUp","",35,-3.4033,5.7594); h_tt1m1e_good_goodLep_phi_EMuUp->Sumw2(); 

 TH1F *h_tt1m1e_good_goodJ1_pt_EMuDown   = new TH1F("h_tt1m1e_good_goodJ1_pt_EMuDown","",12,25,205); h_tt1m1e_good_goodJ1_pt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_pt_EMuDown   = new TH1F("h_tt1m1e_good_goodJ2_pt_EMuDown","",12,25,205); h_tt1m1e_good_goodJ2_pt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_eta_EMuDown  = new TH1F("h_tt1m1e_good_goodJ1_eta_EMuDown","",20,-3,5); h_tt1m1e_good_goodJ1_eta_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_eta_EMuDown  = new TH1F("h_tt1m1e_good_goodJ2_eta_EMuDown","",20,-3,5); h_tt1m1e_good_goodJ2_eta_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_phi_EMuDown  = new TH1F("h_tt1m1e_good_goodJ1_phi_EMuDown","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ1_phi_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_phi_EMuDown  = new TH1F("h_tt1m1e_good_goodJ2_phi_EMuDown","",35,-3.4033,5.7594); h_tt1m1e_good_goodJ2_phi_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ1_CSV_EMuDown  = new TH1F("h_tt1m1e_good_goodJ1_CSV_EMuDown","",30,0,1); h_tt1m1e_good_goodJ1_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodJ2_CSV_EMuDown  = new TH1F("h_tt1m1e_good_goodJ2_CSV_EMuDown","",30,0,1); h_tt1m1e_good_goodJ2_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_dRgoodJ1J2_EMuDown  = new TH1F("h_tt1m1e_good_dRgoodJ1J2_EMuDown","",20,0,6); h_tt1m1e_good_dRgoodJ1J2_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_mt_EMuDown          = new TH1F("h_tt1m1e_good_mt_EMuDown","",40,0,200); h_tt1m1e_good_mt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_pt_EMuDown  = new TH1F("h_tt1m1e_good_goodLep_pt_EMuDown","",12,25,205); h_tt1m1e_good_goodLep_pt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_eta_EMuDown = new TH1F("h_tt1m1e_good_goodLep_eta_EMuDown","",20,-3,5); h_tt1m1e_good_goodLep_eta_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_good_goodLep_phi_EMuDown = new TH1F("h_tt1m1e_good_goodLep_phi_EMuDown","",35,-3.4033,5.7594); h_tt1m1e_good_goodLep_phi_EMuDown->Sumw2(); 

 TH1F *h_tt1m1e_qcd_goodJ1_pt   = new TH1F("h_tt1m1e_qcd_goodJ1_pt","",12,25,205); h_tt1m1e_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_pt   = new TH1F("h_tt1m1e_qcd_goodJ2_pt","",12,25,205); h_tt1m1e_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_eta  = new TH1F("h_tt1m1e_qcd_goodJ1_eta","",20,-3,5); h_tt1m1e_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_eta  = new TH1F("h_tt1m1e_qcd_goodJ2_eta","",20,-3,5); h_tt1m1e_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_phi  = new TH1F("h_tt1m1e_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_phi  = new TH1F("h_tt1m1e_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_CSV  = new TH1F("h_tt1m1e_qcd_goodJ1_CSV","",30,0,1); h_tt1m1e_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_CSV  = new TH1F("h_tt1m1e_qcd_goodJ2_CSV","",30,0,1); h_tt1m1e_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_tt1m1e_qcd_dRgoodJ1J2  = new TH1F("h_tt1m1e_qcd_dRgoodJ1J2","",20,0,6); h_tt1m1e_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_tt1m1e_qcd_mt          = new TH1F("h_tt1m1e_qcd_mt","",40,0,200); h_tt1m1e_qcd_mt->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_pt  = new TH1F("h_tt1m1e_qcd_qcdLep_pt","",12,25,205); h_tt1m1e_qcd_qcdLep_pt->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_eta = new TH1F("h_tt1m1e_qcd_qcdLep_eta","",20,-3,5); h_tt1m1e_qcd_qcdLep_eta->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_phi = new TH1F("h_tt1m1e_qcd_qcdLep_phi","",35,-3.4033,5.7594); h_tt1m1e_qcd_qcdLep_phi->Sumw2(); 

 TH1F *h_tt1m1e_qcd_goodJ1_pt_CSVUp   = new TH1F("h_tt1m1e_qcd_goodJ1_pt_CSVUp","",12,25,205); h_tt1m1e_qcd_goodJ1_pt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_pt_CSVUp   = new TH1F("h_tt1m1e_qcd_goodJ2_pt_CSVUp","",12,25,205); h_tt1m1e_qcd_goodJ2_pt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_eta_CSVUp  = new TH1F("h_tt1m1e_qcd_goodJ1_eta_CSVUp","",20,-3,5); h_tt1m1e_qcd_goodJ1_eta_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_eta_CSVUp  = new TH1F("h_tt1m1e_qcd_goodJ2_eta_CSVUp","",20,-3,5); h_tt1m1e_qcd_goodJ2_eta_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_phi_CSVUp  = new TH1F("h_tt1m1e_qcd_goodJ1_phi_CSVUp","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ1_phi_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_phi_CSVUp  = new TH1F("h_tt1m1e_qcd_goodJ2_phi_CSVUp","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ2_phi_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_CSV_CSVUp  = new TH1F("h_tt1m1e_qcd_goodJ1_CSV_CSVUp","",30,0,1); h_tt1m1e_qcd_goodJ1_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_CSV_CSVUp  = new TH1F("h_tt1m1e_qcd_goodJ2_CSV_CSVUp","",30,0,1); h_tt1m1e_qcd_goodJ2_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_dRgoodJ1J2_CSVUp  = new TH1F("h_tt1m1e_qcd_dRgoodJ1J2_CSVUp","",20,0,6); h_tt1m1e_qcd_dRgoodJ1J2_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_mt_CSVUp          = new TH1F("h_tt1m1e_qcd_mt_CSVUp","",40,0,200); h_tt1m1e_qcd_mt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_pt_CSVUp  = new TH1F("h_tt1m1e_qcd_qcdLep_pt_CSVUp","",12,25,205); h_tt1m1e_qcd_qcdLep_pt_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_eta_CSVUp = new TH1F("h_tt1m1e_qcd_qcdLep_eta_CSVUp","",20,-3,5); h_tt1m1e_qcd_qcdLep_eta_CSVUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_phi_CSVUp = new TH1F("h_tt1m1e_qcd_qcdLep_phi_CSVUp","",35,-3.4033,5.7594); h_tt1m1e_qcd_qcdLep_phi_CSVUp->Sumw2(); 

 TH1F *h_tt1m1e_qcd_goodJ1_pt_CSVDown   = new TH1F("h_tt1m1e_qcd_goodJ1_pt_CSVDown","",12,25,205); h_tt1m1e_qcd_goodJ1_pt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_pt_CSVDown   = new TH1F("h_tt1m1e_qcd_goodJ2_pt_CSVDown","",12,25,205); h_tt1m1e_qcd_goodJ2_pt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_eta_CSVDown  = new TH1F("h_tt1m1e_qcd_goodJ1_eta_CSVDown","",20,-3,5); h_tt1m1e_qcd_goodJ1_eta_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_eta_CSVDown  = new TH1F("h_tt1m1e_qcd_goodJ2_eta_CSVDown","",20,-3,5); h_tt1m1e_qcd_goodJ2_eta_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_phi_CSVDown  = new TH1F("h_tt1m1e_qcd_goodJ1_phi_CSVDown","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ1_phi_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_phi_CSVDown  = new TH1F("h_tt1m1e_qcd_goodJ2_phi_CSVDown","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ2_phi_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_CSV_CSVDown  = new TH1F("h_tt1m1e_qcd_goodJ1_CSV_CSVDown","",30,0,1); h_tt1m1e_qcd_goodJ1_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_CSV_CSVDown  = new TH1F("h_tt1m1e_qcd_goodJ2_CSV_CSVDown","",30,0,1); h_tt1m1e_qcd_goodJ2_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_dRgoodJ1J2_CSVDown  = new TH1F("h_tt1m1e_qcd_dRgoodJ1J2_CSVDown","",20,0,6); h_tt1m1e_qcd_dRgoodJ1J2_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_mt_CSVDown          = new TH1F("h_tt1m1e_qcd_mt_CSVDown","",40,0,200); h_tt1m1e_qcd_mt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_pt_CSVDown  = new TH1F("h_tt1m1e_qcd_qcdLep_pt_CSVDown","",12,25,205); h_tt1m1e_qcd_qcdLep_pt_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_eta_CSVDown = new TH1F("h_tt1m1e_qcd_qcdLep_eta_CSVDown","",20,-3,5); h_tt1m1e_qcd_qcdLep_eta_CSVDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_phi_CSVDown = new TH1F("h_tt1m1e_qcd_qcdLep_phi_CSVDown","",35,-3.4033,5.7594); h_tt1m1e_qcd_qcdLep_phi_CSVDown->Sumw2(); 

 TH1F *h_tt1m1e_qcd_goodJ1_pt_EMuUp   = new TH1F("h_tt1m1e_qcd_goodJ1_pt_EMuUp","",12,25,205); h_tt1m1e_qcd_goodJ1_pt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_pt_EMuUp   = new TH1F("h_tt1m1e_qcd_goodJ2_pt_EMuUp","",12,25,205); h_tt1m1e_qcd_goodJ2_pt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_eta_EMuUp  = new TH1F("h_tt1m1e_qcd_goodJ1_eta_EMuUp","",20,-3,5); h_tt1m1e_qcd_goodJ1_eta_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_eta_EMuUp  = new TH1F("h_tt1m1e_qcd_goodJ2_eta_EMuUp","",20,-3,5); h_tt1m1e_qcd_goodJ2_eta_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_phi_EMuUp  = new TH1F("h_tt1m1e_qcd_goodJ1_phi_EMuUp","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ1_phi_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_phi_EMuUp  = new TH1F("h_tt1m1e_qcd_goodJ2_phi_EMuUp","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ2_phi_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_CSV_EMuUp  = new TH1F("h_tt1m1e_qcd_goodJ1_CSV_EMuUp","",30,0,1); h_tt1m1e_qcd_goodJ1_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_CSV_EMuUp  = new TH1F("h_tt1m1e_qcd_goodJ2_CSV_EMuUp","",30,0,1); h_tt1m1e_qcd_goodJ2_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_dRgoodJ1J2_EMuUp  = new TH1F("h_tt1m1e_qcd_dRgoodJ1J2_EMuUp","",20,0,6); h_tt1m1e_qcd_dRgoodJ1J2_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_mt_EMuUp          = new TH1F("h_tt1m1e_qcd_mt_EMuUp","",40,0,200); h_tt1m1e_qcd_mt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_pt_EMuUp  = new TH1F("h_tt1m1e_qcd_qcdLep_pt_EMuUp","",12,25,205); h_tt1m1e_qcd_qcdLep_pt_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_eta_EMuUp = new TH1F("h_tt1m1e_qcd_qcdLep_eta_EMuUp","",20,-3,5); h_tt1m1e_qcd_qcdLep_eta_EMuUp->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_phi_EMuUp = new TH1F("h_tt1m1e_qcd_qcdLep_phi_EMuUp","",35,-3.4033,5.7594); h_tt1m1e_qcd_qcdLep_phi_EMuUp->Sumw2(); 

 TH1F *h_tt1m1e_qcd_goodJ1_pt_EMuDown   = new TH1F("h_tt1m1e_qcd_goodJ1_pt_EMuDown","",12,25,205); h_tt1m1e_qcd_goodJ1_pt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_pt_EMuDown   = new TH1F("h_tt1m1e_qcd_goodJ2_pt_EMuDown","",12,25,205); h_tt1m1e_qcd_goodJ2_pt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_eta_EMuDown  = new TH1F("h_tt1m1e_qcd_goodJ1_eta_EMuDown","",20,-3,5); h_tt1m1e_qcd_goodJ1_eta_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_eta_EMuDown  = new TH1F("h_tt1m1e_qcd_goodJ2_eta_EMuDown","",20,-3,5); h_tt1m1e_qcd_goodJ2_eta_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_phi_EMuDown  = new TH1F("h_tt1m1e_qcd_goodJ1_phi_EMuDown","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ1_phi_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_phi_EMuDown  = new TH1F("h_tt1m1e_qcd_goodJ2_phi_EMuDown","",35,-3.4033,5.7594); h_tt1m1e_qcd_goodJ2_phi_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ1_CSV_EMuDown  = new TH1F("h_tt1m1e_qcd_goodJ1_CSV_EMuDown","",30,0,1); h_tt1m1e_qcd_goodJ1_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_goodJ2_CSV_EMuDown  = new TH1F("h_tt1m1e_qcd_goodJ2_CSV_EMuDown","",30,0,1); h_tt1m1e_qcd_goodJ2_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_dRgoodJ1J2_EMuDown  = new TH1F("h_tt1m1e_qcd_dRgoodJ1J2_EMuDown","",20,0,6); h_tt1m1e_qcd_dRgoodJ1J2_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_mt_EMuDown          = new TH1F("h_tt1m1e_qcd_mt_EMuDown","",40,0,200); h_tt1m1e_qcd_mt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_pt_EMuDown  = new TH1F("h_tt1m1e_qcd_qcdLep_pt_EMuDown","",12,25,205); h_tt1m1e_qcd_qcdLep_pt_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_eta_EMuDown = new TH1F("h_tt1m1e_qcd_qcdLep_eta_EMuDown","",20,-3,5); h_tt1m1e_qcd_qcdLep_eta_EMuDown->Sumw2(); 
 TH1F *h_tt1m1e_qcd_qcdLep_phi_EMuDown = new TH1F("h_tt1m1e_qcd_qcdLep_phi_EMuDown","",35,-3.4033,5.7594); h_tt1m1e_qcd_qcdLep_phi_EMuDown->Sumw2(); 
                                                
 std::cout<<"initializing tt3j histograms"<<std::endl;
 // tt3j
 TH1F *h_tt3j_good_goodJ1_pt   = new TH1F("h_tt3j_good_goodJ1_pt","",12,25,205); h_tt3j_good_goodJ1_pt->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_pt   = new TH1F("h_tt3j_good_goodJ2_pt","",12,25,205); h_tt3j_good_goodJ2_pt->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_eta  = new TH1F("h_tt3j_good_goodJ1_eta","",20,-3,5); h_tt3j_good_goodJ1_eta->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_eta  = new TH1F("h_tt3j_good_goodJ2_eta","",20,-3,5); h_tt3j_good_goodJ2_eta->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_phi  = new TH1F("h_tt3j_good_goodJ1_phi","",35,-3.4033,5.7594); h_tt3j_good_goodJ1_phi->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_phi  = new TH1F("h_tt3j_good_goodJ2_phi","",35,-3.4033,5.7594); h_tt3j_good_goodJ2_phi->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_CSV  = new TH1F("h_tt3j_good_goodJ1_CSV","",30,0,1); h_tt3j_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_CSV  = new TH1F("h_tt3j_good_goodJ2_CSV","",30,0,1); h_tt3j_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_tt3j_good_dRgoodJ1J2  = new TH1F("h_tt3j_good_dRgoodJ1J2","",20,0,6); h_tt3j_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_tt3j_good_mt          = new TH1F("h_tt3j_good_mt","",40,0,200); h_tt3j_good_mt->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_pt  = new TH1F("h_tt3j_good_goodLep_pt","",12,25,205); h_tt3j_good_goodLep_pt->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_eta = new TH1F("h_tt3j_good_goodLep_eta","",20,-3,5); h_tt3j_good_goodLep_eta->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_phi = new TH1F("h_tt3j_good_goodLep_phi","",35,-3.4033,5.7594); h_tt3j_good_goodLep_phi->Sumw2(); 

 TH1F *h_tt3j_good_goodJ1_pt_CSVUp   = new TH1F("h_tt3j_good_goodJ1_pt_CSVUp","",12,25,205); h_tt3j_good_goodJ1_pt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_pt_CSVUp   = new TH1F("h_tt3j_good_goodJ2_pt_CSVUp","",12,25,205); h_tt3j_good_goodJ2_pt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_eta_CSVUp  = new TH1F("h_tt3j_good_goodJ1_eta_CSVUp","",20,-3,5); h_tt3j_good_goodJ1_eta_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_eta_CSVUp  = new TH1F("h_tt3j_good_goodJ2_eta_CSVUp","",20,-3,5); h_tt3j_good_goodJ2_eta_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_phi_CSVUp  = new TH1F("h_tt3j_good_goodJ1_phi_CSVUp","",35,-3.4033,5.7594); h_tt3j_good_goodJ1_phi_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_phi_CSVUp  = new TH1F("h_tt3j_good_goodJ2_phi_CSVUp","",35,-3.4033,5.7594); h_tt3j_good_goodJ2_phi_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_CSV_CSVUp  = new TH1F("h_tt3j_good_goodJ1_CSV_CSVUp","",30,0,1); h_tt3j_good_goodJ1_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_CSV_CSVUp  = new TH1F("h_tt3j_good_goodJ2_CSV_CSVUp","",30,0,1); h_tt3j_good_goodJ2_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_dRgoodJ1J2_CSVUp  = new TH1F("h_tt3j_good_dRgoodJ1J2_CSVUp","",20,0,6); h_tt3j_good_dRgoodJ1J2_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_mt_CSVUp          = new TH1F("h_tt3j_good_mt_CSVUp","",40,0,200); h_tt3j_good_mt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_pt_CSVUp  = new TH1F("h_tt3j_good_goodLep_pt_CSVUp","",12,25,205); h_tt3j_good_goodLep_pt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_eta_CSVUp = new TH1F("h_tt3j_good_goodLep_eta_CSVUp","",20,-3,5); h_tt3j_good_goodLep_eta_CSVUp->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_phi_CSVUp = new TH1F("h_tt3j_good_goodLep_phi_CSVUp","",35,-3.4033,5.7594); h_tt3j_good_goodLep_phi_CSVUp->Sumw2(); 

 TH1F *h_tt3j_good_goodJ1_pt_CSVDown   = new TH1F("h_tt3j_good_goodJ1_pt_CSVDown","",12,25,205); h_tt3j_good_goodJ1_pt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_pt_CSVDown   = new TH1F("h_tt3j_good_goodJ2_pt_CSVDown","",12,25,205); h_tt3j_good_goodJ2_pt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_eta_CSVDown  = new TH1F("h_tt3j_good_goodJ1_eta_CSVDown","",20,-3,5); h_tt3j_good_goodJ1_eta_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_eta_CSVDown  = new TH1F("h_tt3j_good_goodJ2_eta_CSVDown","",20,-3,5); h_tt3j_good_goodJ2_eta_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_phi_CSVDown  = new TH1F("h_tt3j_good_goodJ1_phi_CSVDown","",35,-3.4033,5.7594); h_tt3j_good_goodJ1_phi_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_phi_CSVDown  = new TH1F("h_tt3j_good_goodJ2_phi_CSVDown","",35,-3.4033,5.7594); h_tt3j_good_goodJ2_phi_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_CSV_CSVDown  = new TH1F("h_tt3j_good_goodJ1_CSV_CSVDown","",30,0,1); h_tt3j_good_goodJ1_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_CSV_CSVDown  = new TH1F("h_tt3j_good_goodJ2_CSV_CSVDown","",30,0,1); h_tt3j_good_goodJ2_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_dRgoodJ1J2_CSVDown  = new TH1F("h_tt3j_good_dRgoodJ1J2_CSVDown","",20,0,6); h_tt3j_good_dRgoodJ1J2_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_mt_CSVDown          = new TH1F("h_tt3j_good_mt_CSVDown","",40,0,200); h_tt3j_good_mt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_pt_CSVDown  = new TH1F("h_tt3j_good_goodLep_pt_CSVDown","",12,25,205); h_tt3j_good_goodLep_pt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_eta_CSVDown = new TH1F("h_tt3j_good_goodLep_eta_CSVDown","",20,-3,5); h_tt3j_good_goodLep_eta_CSVDown->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_phi_CSVDown = new TH1F("h_tt3j_good_goodLep_phi_CSVDown","",35,-3.4033,5.7594); h_tt3j_good_goodLep_phi_CSVDown->Sumw2(); 

 TH1F *h_tt3j_good_goodJ1_pt_EMuUp   = new TH1F("h_tt3j_good_goodJ1_pt_EMuUp","",12,25,205); h_tt3j_good_goodJ1_pt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_pt_EMuUp   = new TH1F("h_tt3j_good_goodJ2_pt_EMuUp","",12,25,205); h_tt3j_good_goodJ2_pt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_eta_EMuUp  = new TH1F("h_tt3j_good_goodJ1_eta_EMuUp","",20,-3,5); h_tt3j_good_goodJ1_eta_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_eta_EMuUp  = new TH1F("h_tt3j_good_goodJ2_eta_EMuUp","",20,-3,5); h_tt3j_good_goodJ2_eta_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_phi_EMuUp  = new TH1F("h_tt3j_good_goodJ1_phi_EMuUp","",35,-3.4033,5.7594); h_tt3j_good_goodJ1_phi_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_phi_EMuUp  = new TH1F("h_tt3j_good_goodJ2_phi_EMuUp","",35,-3.4033,5.7594); h_tt3j_good_goodJ2_phi_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_CSV_EMuUp  = new TH1F("h_tt3j_good_goodJ1_CSV_EMuUp","",30,0,1); h_tt3j_good_goodJ1_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_CSV_EMuUp  = new TH1F("h_tt3j_good_goodJ2_CSV_EMuUp","",30,0,1); h_tt3j_good_goodJ2_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_dRgoodJ1J2_EMuUp  = new TH1F("h_tt3j_good_dRgoodJ1J2_EMuUp","",20,0,6); h_tt3j_good_dRgoodJ1J2_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_mt_EMuUp          = new TH1F("h_tt3j_good_mt_EMuUp","",40,0,200); h_tt3j_good_mt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_pt_EMuUp  = new TH1F("h_tt3j_good_goodLep_pt_EMuUp","",12,25,205); h_tt3j_good_goodLep_pt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_eta_EMuUp = new TH1F("h_tt3j_good_goodLep_eta_EMuUp","",20,-3,5); h_tt3j_good_goodLep_eta_EMuUp->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_phi_EMuUp = new TH1F("h_tt3j_good_goodLep_phi_EMuUp","",35,-3.4033,5.7594); h_tt3j_good_goodLep_phi_EMuUp->Sumw2(); 

 TH1F *h_tt3j_good_goodJ1_pt_EMuDown   = new TH1F("h_tt3j_good_goodJ1_pt_EMuDown","",12,25,205); h_tt3j_good_goodJ1_pt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_pt_EMuDown   = new TH1F("h_tt3j_good_goodJ2_pt_EMuDown","",12,25,205); h_tt3j_good_goodJ2_pt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_eta_EMuDown  = new TH1F("h_tt3j_good_goodJ1_eta_EMuDown","",20,-3,5); h_tt3j_good_goodJ1_eta_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_eta_EMuDown  = new TH1F("h_tt3j_good_goodJ2_eta_EMuDown","",20,-3,5); h_tt3j_good_goodJ2_eta_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_phi_EMuDown  = new TH1F("h_tt3j_good_goodJ1_phi_EMuDown","",35,-3.4033,5.7594); h_tt3j_good_goodJ1_phi_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_phi_EMuDown  = new TH1F("h_tt3j_good_goodJ2_phi_EMuDown","",35,-3.4033,5.7594); h_tt3j_good_goodJ2_phi_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ1_CSV_EMuDown  = new TH1F("h_tt3j_good_goodJ1_CSV_EMuDown","",30,0,1); h_tt3j_good_goodJ1_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodJ2_CSV_EMuDown  = new TH1F("h_tt3j_good_goodJ2_CSV_EMuDown","",30,0,1); h_tt3j_good_goodJ2_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_dRgoodJ1J2_EMuDown  = new TH1F("h_tt3j_good_dRgoodJ1J2_EMuDown","",20,0,6); h_tt3j_good_dRgoodJ1J2_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_mt_EMuDown          = new TH1F("h_tt3j_good_mt_EMuDown","",40,0,200); h_tt3j_good_mt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_pt_EMuDown  = new TH1F("h_tt3j_good_goodLep_pt_EMuDown","",12,25,205); h_tt3j_good_goodLep_pt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_eta_EMuDown = new TH1F("h_tt3j_good_goodLep_eta_EMuDown","",20,-3,5); h_tt3j_good_goodLep_eta_EMuDown->Sumw2(); 
 TH1F *h_tt3j_good_goodLep_phi_EMuDown = new TH1F("h_tt3j_good_goodLep_phi_EMuDown","",35,-3.4033,5.7594); h_tt3j_good_goodLep_phi_EMuDown->Sumw2(); 

 TH1F *h_tt3j_qcd_goodJ1_pt   = new TH1F("h_tt3j_qcd_goodJ1_pt","",12,25,205); h_tt3j_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_pt   = new TH1F("h_tt3j_qcd_goodJ2_pt","",12,25,205); h_tt3j_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_eta  = new TH1F("h_tt3j_qcd_goodJ1_eta","",20,-3,5); h_tt3j_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_eta  = new TH1F("h_tt3j_qcd_goodJ2_eta","",20,-3,5); h_tt3j_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_phi  = new TH1F("h_tt3j_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_phi  = new TH1F("h_tt3j_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_CSV  = new TH1F("h_tt3j_qcd_goodJ1_CSV","",30,0,1); h_tt3j_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_CSV  = new TH1F("h_tt3j_qcd_goodJ2_CSV","",30,0,1); h_tt3j_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_tt3j_qcd_dRgoodJ1J2  = new TH1F("h_tt3j_qcd_dRgoodJ1J2","",20,0,6); h_tt3j_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_tt3j_qcd_mt          = new TH1F("h_tt3j_qcd_mt","",40,0,200); h_tt3j_qcd_mt->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_pt  = new TH1F("h_tt3j_qcd_qcdLep_pt","",12,25,205); h_tt3j_qcd_qcdLep_pt->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_eta = new TH1F("h_tt3j_qcd_qcdLep_eta","",20,-3,5); h_tt3j_qcd_qcdLep_eta->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_phi = new TH1F("h_tt3j_qcd_qcdLep_phi","",35,-3.4033,5.7594); h_tt3j_qcd_qcdLep_phi->Sumw2(); 

 TH1F *h_tt3j_qcd_goodJ1_pt_CSVUp   = new TH1F("h_tt3j_qcd_goodJ1_pt_CSVUp","",12,25,205); h_tt3j_qcd_goodJ1_pt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_pt_CSVUp   = new TH1F("h_tt3j_qcd_goodJ2_pt_CSVUp","",12,25,205); h_tt3j_qcd_goodJ2_pt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_eta_CSVUp  = new TH1F("h_tt3j_qcd_goodJ1_eta_CSVUp","",20,-3,5); h_tt3j_qcd_goodJ1_eta_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_eta_CSVUp  = new TH1F("h_tt3j_qcd_goodJ2_eta_CSVUp","",20,-3,5); h_tt3j_qcd_goodJ2_eta_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_phi_CSVUp  = new TH1F("h_tt3j_qcd_goodJ1_phi_CSVUp","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ1_phi_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_phi_CSVUp  = new TH1F("h_tt3j_qcd_goodJ2_phi_CSVUp","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ2_phi_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_CSV_CSVUp  = new TH1F("h_tt3j_qcd_goodJ1_CSV_CSVUp","",30,0,1); h_tt3j_qcd_goodJ1_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_CSV_CSVUp  = new TH1F("h_tt3j_qcd_goodJ2_CSV_CSVUp","",30,0,1); h_tt3j_qcd_goodJ2_CSV_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_dRgoodJ1J2_CSVUp  = new TH1F("h_tt3j_qcd_dRgoodJ1J2_CSVUp","",20,0,6); h_tt3j_qcd_dRgoodJ1J2_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_mt_CSVUp          = new TH1F("h_tt3j_qcd_mt_CSVUp","",40,0,200); h_tt3j_qcd_mt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_pt_CSVUp  = new TH1F("h_tt3j_qcd_qcdLep_pt_CSVUp","",12,25,205); h_tt3j_qcd_qcdLep_pt_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_eta_CSVUp = new TH1F("h_tt3j_qcd_qcdLep_eta_CSVUp","",20,-3,5); h_tt3j_qcd_qcdLep_eta_CSVUp->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_phi_CSVUp = new TH1F("h_tt3j_qcd_qcdLep_phi_CSVUp","",35,-3.4033,5.7594); h_tt3j_qcd_qcdLep_phi_CSVUp->Sumw2(); 

 TH1F *h_tt3j_qcd_goodJ1_pt_CSVDown   = new TH1F("h_tt3j_qcd_goodJ1_pt_CSVDown","",12,25,205); h_tt3j_qcd_goodJ1_pt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_pt_CSVDown   = new TH1F("h_tt3j_qcd_goodJ2_pt_CSVDown","",12,25,205); h_tt3j_qcd_goodJ2_pt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_eta_CSVDown  = new TH1F("h_tt3j_qcd_goodJ1_eta_CSVDown","",20,-3,5); h_tt3j_qcd_goodJ1_eta_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_eta_CSVDown  = new TH1F("h_tt3j_qcd_goodJ2_eta_CSVDown","",20,-3,5); h_tt3j_qcd_goodJ2_eta_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_phi_CSVDown  = new TH1F("h_tt3j_qcd_goodJ1_phi_CSVDown","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ1_phi_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_phi_CSVDown  = new TH1F("h_tt3j_qcd_goodJ2_phi_CSVDown","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ2_phi_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_CSV_CSVDown  = new TH1F("h_tt3j_qcd_goodJ1_CSV_CSVDown","",30,0,1); h_tt3j_qcd_goodJ1_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_CSV_CSVDown  = new TH1F("h_tt3j_qcd_goodJ2_CSV_CSVDown","",30,0,1); h_tt3j_qcd_goodJ2_CSV_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_dRgoodJ1J2_CSVDown  = new TH1F("h_tt3j_qcd_dRgoodJ1J2_CSVDown","",20,0,6); h_tt3j_qcd_dRgoodJ1J2_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_mt_CSVDown          = new TH1F("h_tt3j_qcd_mt_CSVDown","",40,0,200); h_tt3j_qcd_mt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_pt_CSVDown  = new TH1F("h_tt3j_qcd_qcdLep_pt_CSVDown","",12,25,205); h_tt3j_qcd_qcdLep_pt_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_eta_CSVDown = new TH1F("h_tt3j_qcd_qcdLep_eta_CSVDown","",20,-3,5); h_tt3j_qcd_qcdLep_eta_CSVDown->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_phi_CSVDown = new TH1F("h_tt3j_qcd_qcdLep_phi_CSVDown","",35,-3.4033,5.7594); h_tt3j_qcd_qcdLep_phi_CSVDown->Sumw2(); 

 TH1F *h_tt3j_qcd_goodJ1_pt_EMuUp   = new TH1F("h_tt3j_qcd_goodJ1_pt_EMuUp","",12,25,205); h_tt3j_qcd_goodJ1_pt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_pt_EMuUp   = new TH1F("h_tt3j_qcd_goodJ2_pt_EMuUp","",12,25,205); h_tt3j_qcd_goodJ2_pt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_eta_EMuUp  = new TH1F("h_tt3j_qcd_goodJ1_eta_EMuUp","",20,-3,5); h_tt3j_qcd_goodJ1_eta_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_eta_EMuUp  = new TH1F("h_tt3j_qcd_goodJ2_eta_EMuUp","",20,-3,5); h_tt3j_qcd_goodJ2_eta_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_phi_EMuUp  = new TH1F("h_tt3j_qcd_goodJ1_phi_EMuUp","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ1_phi_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_phi_EMuUp  = new TH1F("h_tt3j_qcd_goodJ2_phi_EMuUp","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ2_phi_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_CSV_EMuUp  = new TH1F("h_tt3j_qcd_goodJ1_CSV_EMuUp","",30,0,1); h_tt3j_qcd_goodJ1_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_CSV_EMuUp  = new TH1F("h_tt3j_qcd_goodJ2_CSV_EMuUp","",30,0,1); h_tt3j_qcd_goodJ2_CSV_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_dRgoodJ1J2_EMuUp  = new TH1F("h_tt3j_qcd_dRgoodJ1J2_EMuUp","",20,0,6); h_tt3j_qcd_dRgoodJ1J2_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_mt_EMuUp          = new TH1F("h_tt3j_qcd_mt_EMuUp","",40,0,200); h_tt3j_qcd_mt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_pt_EMuUp  = new TH1F("h_tt3j_qcd_qcdLep_pt_EMuUp","",12,25,205); h_tt3j_qcd_qcdLep_pt_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_eta_EMuUp = new TH1F("h_tt3j_qcd_qcdLep_eta_EMuUp","",20,-3,5); h_tt3j_qcd_qcdLep_eta_EMuUp->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_phi_EMuUp = new TH1F("h_tt3j_qcd_qcdLep_phi_EMuUp","",35,-3.4033,5.7594); h_tt3j_qcd_qcdLep_phi_EMuUp->Sumw2(); 

 TH1F *h_tt3j_qcd_goodJ1_pt_EMuDown   = new TH1F("h_tt3j_qcd_goodJ1_pt_EMuDown","",12,25,205); h_tt3j_qcd_goodJ1_pt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_pt_EMuDown   = new TH1F("h_tt3j_qcd_goodJ2_pt_EMuDown","",12,25,205); h_tt3j_qcd_goodJ2_pt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_eta_EMuDown  = new TH1F("h_tt3j_qcd_goodJ1_eta_EMuDown","",20,-3,5); h_tt3j_qcd_goodJ1_eta_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_eta_EMuDown  = new TH1F("h_tt3j_qcd_goodJ2_eta_EMuDown","",20,-3,5); h_tt3j_qcd_goodJ2_eta_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_phi_EMuDown  = new TH1F("h_tt3j_qcd_goodJ1_phi_EMuDown","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ1_phi_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_phi_EMuDown  = new TH1F("h_tt3j_qcd_goodJ2_phi_EMuDown","",35,-3.4033,5.7594); h_tt3j_qcd_goodJ2_phi_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ1_CSV_EMuDown  = new TH1F("h_tt3j_qcd_goodJ1_CSV_EMuDown","",30,0,1); h_tt3j_qcd_goodJ1_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_goodJ2_CSV_EMuDown  = new TH1F("h_tt3j_qcd_goodJ2_CSV_EMuDown","",30,0,1); h_tt3j_qcd_goodJ2_CSV_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_dRgoodJ1J2_EMuDown  = new TH1F("h_tt3j_qcd_dRgoodJ1J2_EMuDown","",20,0,6); h_tt3j_qcd_dRgoodJ1J2_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_mt_EMuDown          = new TH1F("h_tt3j_qcd_mt_EMuDown","",40,0,200); h_tt3j_qcd_mt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_pt_EMuDown  = new TH1F("h_tt3j_qcd_qcdLep_pt_EMuDown","",12,25,205); h_tt3j_qcd_qcdLep_pt_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_eta_EMuDown = new TH1F("h_tt3j_qcd_qcdLep_eta_EMuDown","",20,-3,5); h_tt3j_qcd_qcdLep_eta_EMuDown->Sumw2(); 
 TH1F *h_tt3j_qcd_qcdLep_phi_EMuDown = new TH1F("h_tt3j_qcd_qcdLep_phi_EMuDown","",35,-3.4033,5.7594); h_tt3j_qcd_qcdLep_phi_EMuDown->Sumw2(); 
                                                
 std::cout<<"initializing wjj histograms"<<std::endl;
 // wjj
 TH1F *h_wjj_good_goodJ1_pt   = new TH1F("h_wjj_good_goodJ1_pt","",12,25,205); h_wjj_good_goodJ1_pt->Sumw2(); 
 TH1F *h_wjj_good_goodJ2_pt   = new TH1F("h_wjj_good_goodJ2_pt","",12,25,205); h_wjj_good_goodJ2_pt->Sumw2(); 
 TH1F *h_wjj_good_goodJ1_eta  = new TH1F("h_wjj_good_goodJ1_eta","",20,-3,5); h_wjj_good_goodJ1_eta->Sumw2(); 
 TH1F *h_wjj_good_goodJ2_eta  = new TH1F("h_wjj_good_goodJ2_eta","",20,-3,5); h_wjj_good_goodJ2_eta->Sumw2(); 
 TH1F *h_wjj_good_goodJ1_phi  = new TH1F("h_wjj_good_goodJ1_phi","",35,-3.4033,5.7594); h_wjj_good_goodJ1_phi->Sumw2(); 
 TH1F *h_wjj_good_goodJ2_phi  = new TH1F("h_wjj_good_goodJ2_phi","",35,-3.4033,5.7594); h_wjj_good_goodJ2_phi->Sumw2(); 
 TH1F *h_wjj_good_goodJ1_CSV  = new TH1F("h_wjj_good_goodJ1_CSV","",30,0,1); h_wjj_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_wjj_good_goodJ2_CSV  = new TH1F("h_wjj_good_goodJ2_CSV","",30,0,1); h_wjj_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_wjj_good_dRgoodJ1J2  = new TH1F("h_wjj_good_dRgoodJ1J2","",20,0,6); h_wjj_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_wjj_good_mt          = new TH1F("h_wjj_good_mt","",40,0,200); h_wjj_good_mt->Sumw2(); 
 TH1F *h_wjj_good_goodLep_pt  = new TH1F("h_wjj_good_goodLep_pt","",12,25,205); h_wjj_good_goodLep_pt->Sumw2(); 
 TH1F *h_wjj_good_goodLep_eta = new TH1F("h_wjj_good_goodLep_eta","",20,-3,5); h_wjj_good_goodLep_eta->Sumw2(); 
 TH1F *h_wjj_good_goodLep_phi = new TH1F("h_wjj_good_goodLep_phi","",35,-3.4033,5.7594); h_wjj_good_goodLep_phi->Sumw2(); 

 TH1F *h_wjj_qcd_goodJ1_pt   = new TH1F("h_wjj_qcd_goodJ1_pt","",12,25,205); h_wjj_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ2_pt   = new TH1F("h_wjj_qcd_goodJ2_pt","",12,25,205); h_wjj_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ1_eta  = new TH1F("h_wjj_qcd_goodJ1_eta","",20,-3,5); h_wjj_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ2_eta  = new TH1F("h_wjj_qcd_goodJ2_eta","",20,-3,5); h_wjj_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ1_phi  = new TH1F("h_wjj_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_wjj_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ2_phi  = new TH1F("h_wjj_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_wjj_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ1_CSV  = new TH1F("h_wjj_qcd_goodJ1_CSV","",30,0,1); h_wjj_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_wjj_qcd_goodJ2_CSV  = new TH1F("h_wjj_qcd_goodJ2_CSV","",30,0,1); h_wjj_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_wjj_qcd_dRgoodJ1J2  = new TH1F("h_wjj_qcd_dRgoodJ1J2","",20,0,6); h_wjj_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_wjj_qcd_mt          = new TH1F("h_wjj_qcd_mt","",40,0,200); h_wjj_qcd_mt->Sumw2(); 
 TH1F *h_wjj_qcd_qcdLep_pt  = new TH1F("h_wjj_qcd_qcdLep_pt","",12,25,205); h_wjj_qcd_qcdLep_pt->Sumw2(); 
 TH1F *h_wjj_qcd_qcdLep_eta = new TH1F("h_wjj_qcd_qcdLep_eta","",20,-3,5); h_wjj_qcd_qcdLep_eta->Sumw2(); 
 TH1F *h_wjj_qcd_qcdLep_phi = new TH1F("h_wjj_qcd_qcdLep_phi","",35,-3.4033,5.7594); h_wjj_qcd_qcdLep_phi->Sumw2(); 
                                                
 std::cout<<"initializing stt histograms"<<std::endl;
 // stt
 TH1F *h_stt_good_goodJ1_pt   = new TH1F("h_stt_good_goodJ1_pt","",12,25,205); h_stt_good_goodJ1_pt->Sumw2(); 
 TH1F *h_stt_good_goodJ2_pt   = new TH1F("h_stt_good_goodJ2_pt","",12,25,205); h_stt_good_goodJ2_pt->Sumw2(); 
 TH1F *h_stt_good_goodJ1_eta  = new TH1F("h_stt_good_goodJ1_eta","",20,-3,5); h_stt_good_goodJ1_eta->Sumw2(); 
 TH1F *h_stt_good_goodJ2_eta  = new TH1F("h_stt_good_goodJ2_eta","",20,-3,5); h_stt_good_goodJ2_eta->Sumw2(); 
 TH1F *h_stt_good_goodJ1_phi  = new TH1F("h_stt_good_goodJ1_phi","",35,-3.4033,5.7594); h_stt_good_goodJ1_phi->Sumw2(); 
 TH1F *h_stt_good_goodJ2_phi  = new TH1F("h_stt_good_goodJ2_phi","",35,-3.4033,5.7594); h_stt_good_goodJ2_phi->Sumw2(); 
 TH1F *h_stt_good_goodJ1_CSV  = new TH1F("h_stt_good_goodJ1_CSV","",30,0,1); h_stt_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_stt_good_goodJ2_CSV  = new TH1F("h_stt_good_goodJ2_CSV","",30,0,1); h_stt_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_stt_good_dRgoodJ1J2  = new TH1F("h_stt_good_dRgoodJ1J2","",20,0,6); h_stt_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_stt_good_mt          = new TH1F("h_stt_good_mt","",40,0,200); h_stt_good_mt->Sumw2(); 
 TH1F *h_stt_good_goodLep_pt  = new TH1F("h_stt_good_goodLep_pt","",12,25,205); h_stt_good_goodLep_pt->Sumw2(); 
 TH1F *h_stt_good_goodLep_eta = new TH1F("h_stt_good_goodLep_eta","",20,-3,5); h_stt_good_goodLep_eta->Sumw2(); 
 TH1F *h_stt_good_goodLep_phi = new TH1F("h_stt_good_goodLep_phi","",35,-3.4033,5.7594); h_stt_good_goodLep_phi->Sumw2(); 

 TH1F *h_stt_qcd_goodJ1_pt   = new TH1F("h_stt_qcd_goodJ1_pt","",12,25,205); h_stt_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_stt_qcd_goodJ2_pt   = new TH1F("h_stt_qcd_goodJ2_pt","",12,25,205); h_stt_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_stt_qcd_goodJ1_eta  = new TH1F("h_stt_qcd_goodJ1_eta","",20,-3,5); h_stt_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_stt_qcd_goodJ2_eta  = new TH1F("h_stt_qcd_goodJ2_eta","",20,-3,5); h_stt_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_stt_qcd_goodJ1_phi  = new TH1F("h_stt_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_stt_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_stt_qcd_goodJ2_phi  = new TH1F("h_stt_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_stt_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_stt_qcd_goodJ1_CSV  = new TH1F("h_stt_qcd_goodJ1_CSV","",30,0,1); h_stt_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_stt_qcd_goodJ2_CSV  = new TH1F("h_stt_qcd_goodJ2_CSV","",30,0,1); h_stt_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_stt_qcd_dRgoodJ1J2  = new TH1F("h_stt_qcd_dRgoodJ1J2","",20,0,6); h_stt_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_stt_qcd_mt          = new TH1F("h_stt_qcd_mt","",40,0,200); h_stt_qcd_mt->Sumw2(); 
 TH1F *h_stt_qcd_qcdLep_pt  = new TH1F("h_stt_qcd_qcdLep_pt","",12,25,205); h_stt_qcd_qcdLep_pt->Sumw2(); 
 TH1F *h_stt_qcd_qcdLep_eta = new TH1F("h_stt_qcd_qcdLep_eta","",20,-3,5); h_stt_qcd_qcdLep_eta->Sumw2(); 
 TH1F *h_stt_qcd_qcdLep_phi = new TH1F("h_stt_qcd_qcdLep_phi","",35,-3.4033,5.7594); h_stt_qcd_qcdLep_phi->Sumw2(); 
                                                
 std::cout<<"initializing dyjj histograms"<<std::endl;
 // dyjj
 TH1F *h_dyjj_good_goodJ1_pt   = new TH1F("h_dyjj_good_goodJ1_pt","",12,25,205); h_dyjj_good_goodJ1_pt->Sumw2(); 
 TH1F *h_dyjj_good_goodJ2_pt   = new TH1F("h_dyjj_good_goodJ2_pt","",12,25,205); h_dyjj_good_goodJ2_pt->Sumw2(); 
 TH1F *h_dyjj_good_goodJ1_eta  = new TH1F("h_dyjj_good_goodJ1_eta","",20,-3,5); h_dyjj_good_goodJ1_eta->Sumw2(); 
 TH1F *h_dyjj_good_goodJ2_eta  = new TH1F("h_dyjj_good_goodJ2_eta","",20,-3,5); h_dyjj_good_goodJ2_eta->Sumw2(); 
 TH1F *h_dyjj_good_goodJ1_phi  = new TH1F("h_dyjj_good_goodJ1_phi","",35,-3.4033,5.7594); h_dyjj_good_goodJ1_phi->Sumw2(); 
 TH1F *h_dyjj_good_goodJ2_phi  = new TH1F("h_dyjj_good_goodJ2_phi","",35,-3.4033,5.7594); h_dyjj_good_goodJ2_phi->Sumw2(); 
 TH1F *h_dyjj_good_goodJ1_CSV  = new TH1F("h_dyjj_good_goodJ1_CSV","",30,0,1); h_dyjj_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_dyjj_good_goodJ2_CSV  = new TH1F("h_dyjj_good_goodJ2_CSV","",30,0,1); h_dyjj_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_dyjj_good_dRgoodJ1J2  = new TH1F("h_dyjj_good_dRgoodJ1J2","",20,0,6); h_dyjj_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_dyjj_good_goodLep1_pt  = new TH1F("h_dyjj_good_goodLep1_pt","",12,25,205); h_dyjj_good_goodLep1_pt->Sumw2(); 
 TH1F *h_dyjj_good_goodLep1_eta = new TH1F("h_dyjj_good_goodLep1_eta","",20,-3,5); h_dyjj_good_goodLep1_eta->Sumw2(); 
 TH1F *h_dyjj_good_goodLep1_phi = new TH1F("h_dyjj_good_goodLep1_phi","",35,-3.4033,5.7594); h_dyjj_good_goodLep1_phi->Sumw2(); 
 TH1F *h_dyjj_good_goodLep2_pt  = new TH1F("h_dyjj_good_goodLep2_pt","",12,25,205); h_dyjj_good_goodLep2_pt->Sumw2(); 
 TH1F *h_dyjj_good_goodLep2_eta = new TH1F("h_dyjj_good_goodLep2_eta","",20,-3,5); h_dyjj_good_goodLep2_eta->Sumw2(); 
 TH1F *h_dyjj_good_goodLep2_phi = new TH1F("h_dyjj_good_goodLep2_phi","",35,-3.4033,5.7594); h_dyjj_good_goodLep2_phi->Sumw2(); 

 TH1F *h_dyjj_qcd_goodJ1_pt   = new TH1F("h_dyjj_qcd_goodJ1_pt","",12,25,205); h_dyjj_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ2_pt   = new TH1F("h_dyjj_qcd_goodJ2_pt","",12,25,205); h_dyjj_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ1_eta  = new TH1F("h_dyjj_qcd_goodJ1_eta","",20,-3,5); h_dyjj_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ2_eta  = new TH1F("h_dyjj_qcd_goodJ2_eta","",20,-3,5); h_dyjj_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ1_phi  = new TH1F("h_dyjj_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_dyjj_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ2_phi  = new TH1F("h_dyjj_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_dyjj_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ1_CSV  = new TH1F("h_dyjj_qcd_goodJ1_CSV","",30,0,1); h_dyjj_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_dyjj_qcd_goodJ2_CSV  = new TH1F("h_dyjj_qcd_goodJ2_CSV","",30,0,1); h_dyjj_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_dyjj_qcd_dRgoodJ1J2  = new TH1F("h_dyjj_qcd_dRgoodJ1J2","",20,0,6); h_dyjj_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_dyjj_qcd_qcdLep1_pt  = new TH1F("h_dyjj_qcd_qcdLep1_pt","",12,25,205); h_dyjj_qcd_qcdLep1_pt->Sumw2(); 
 TH1F *h_dyjj_qcd_qcdLep1_eta = new TH1F("h_dyjj_qcd_qcdLep1_eta","",20,-3,5); h_dyjj_qcd_qcdLep1_eta->Sumw2(); 
 TH1F *h_dyjj_qcd_qcdLep1_phi = new TH1F("h_dyjj_qcd_qcdLep1_phi","",35,-3.4033,5.7594); h_dyjj_qcd_qcdLep1_phi->Sumw2(); 
 TH1F *h_dyjj_qcd_qcdLep2_pt  = new TH1F("h_dyjj_qcd_qcdLep2_pt","",12,25,205); h_dyjj_qcd_qcdLep2_pt->Sumw2(); 
 TH1F *h_dyjj_qcd_qcdLep2_eta = new TH1F("h_dyjj_qcd_qcdLep2_eta","",20,-3,5); h_dyjj_qcd_qcdLep2_eta->Sumw2(); 
 TH1F *h_dyjj_qcd_qcdLep2_phi = new TH1F("h_dyjj_qcd_qcdLep2_phi","",35,-3.4033,5.7594); h_dyjj_qcd_qcdLep2_phi->Sumw2(); 

 std::cout<<"initializing dybb histograms"<<std::endl;
 // dybb
 TH1F *h_dybb_good_goodJ1_pt   = new TH1F("h_dybb_good_goodJ1_pt","",12,25,205); h_dybb_good_goodJ1_pt->Sumw2(); 
 TH1F *h_dybb_good_goodJ2_pt   = new TH1F("h_dybb_good_goodJ2_pt","",12,25,205); h_dybb_good_goodJ2_pt->Sumw2(); 
 TH1F *h_dybb_good_goodJ1_eta  = new TH1F("h_dybb_good_goodJ1_eta","",20,-3,5); h_dybb_good_goodJ1_eta->Sumw2(); 
 TH1F *h_dybb_good_goodJ2_eta  = new TH1F("h_dybb_good_goodJ2_eta","",20,-3,5); h_dybb_good_goodJ2_eta->Sumw2(); 
 TH1F *h_dybb_good_goodJ1_phi  = new TH1F("h_dybb_good_goodJ1_phi","",35,-3.4033,5.7594); h_dybb_good_goodJ1_phi->Sumw2(); 
 TH1F *h_dybb_good_goodJ2_phi  = new TH1F("h_dybb_good_goodJ2_phi","",35,-3.4033,5.7594); h_dybb_good_goodJ2_phi->Sumw2(); 
 TH1F *h_dybb_good_goodJ1_CSV  = new TH1F("h_dybb_good_goodJ1_CSV","",30,0,1); h_dybb_good_goodJ1_CSV->Sumw2(); 
 TH1F *h_dybb_good_goodJ2_CSV  = new TH1F("h_dybb_good_goodJ2_CSV","",30,0,1); h_dybb_good_goodJ2_CSV->Sumw2(); 
 TH1F *h_dybb_good_dRgoodJ1J2  = new TH1F("h_dybb_good_dRgoodJ1J2","",20,0,6); h_dybb_good_dRgoodJ1J2->Sumw2(); 
 TH1F *h_dybb_good_goodLep1_pt  = new TH1F("h_dybb_good_goodLep1_pt","",12,25,205); h_dybb_good_goodLep1_pt->Sumw2(); 
 TH1F *h_dybb_good_goodLep1_eta = new TH1F("h_dybb_good_goodLep1_eta","",20,-3,5); h_dybb_good_goodLep1_eta->Sumw2(); 
 TH1F *h_dybb_good_goodLep1_phi = new TH1F("h_dybb_good_goodLep1_phi","",35,-3.4033,5.7594); h_dybb_good_goodLep1_phi->Sumw2(); 
 TH1F *h_dybb_good_goodLep2_pt  = new TH1F("h_dybb_good_goodLep2_pt","",12,25,205); h_dybb_good_goodLep2_pt->Sumw2(); 
 TH1F *h_dybb_good_goodLep2_eta = new TH1F("h_dybb_good_goodLep2_eta","",20,-3,5); h_dybb_good_goodLep2_eta->Sumw2(); 
 TH1F *h_dybb_good_goodLep2_phi = new TH1F("h_dybb_good_goodLep2_phi","",35,-3.4033,5.7594); h_dybb_good_goodLep2_phi->Sumw2(); 

 TH1F *h_dybb_qcd_goodJ1_pt   = new TH1F("h_dybb_qcd_goodJ1_pt","",12,25,205); h_dybb_qcd_goodJ1_pt->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ2_pt   = new TH1F("h_dybb_qcd_goodJ2_pt","",12,25,205); h_dybb_qcd_goodJ2_pt->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ1_eta  = new TH1F("h_dybb_qcd_goodJ1_eta","",20,-3,5); h_dybb_qcd_goodJ1_eta->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ2_eta  = new TH1F("h_dybb_qcd_goodJ2_eta","",20,-3,5); h_dybb_qcd_goodJ2_eta->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ1_phi  = new TH1F("h_dybb_qcd_goodJ1_phi","",35,-3.4033,5.7594); h_dybb_qcd_goodJ1_phi->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ2_phi  = new TH1F("h_dybb_qcd_goodJ2_phi","",35,-3.4033,5.7594); h_dybb_qcd_goodJ2_phi->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ1_CSV  = new TH1F("h_dybb_qcd_goodJ1_CSV","",30,0,1); h_dybb_qcd_goodJ1_CSV->Sumw2(); 
 TH1F *h_dybb_qcd_goodJ2_CSV  = new TH1F("h_dybb_qcd_goodJ2_CSV","",30,0,1); h_dybb_qcd_goodJ2_CSV->Sumw2(); 
 TH1F *h_dybb_qcd_dRgoodJ1J2  = new TH1F("h_dybb_qcd_dRgoodJ1J2","",20,0,6); h_dybb_qcd_dRgoodJ1J2->Sumw2(); 
 TH1F *h_dybb_qcd_qcdLep1_pt  = new TH1F("h_dybb_qcd_qcdLep1_pt","",12,25,205); h_dybb_qcd_qcdLep1_pt->Sumw2(); 
 TH1F *h_dybb_qcd_qcdLep1_eta = new TH1F("h_dybb_qcd_qcdLep1_eta","",20,-3,5); h_dybb_qcd_qcdLep1_eta->Sumw2(); 
 TH1F *h_dybb_qcd_qcdLep1_phi = new TH1F("h_dybb_qcd_qcdLep1_phi","",35,-3.4033,5.7594); h_dybb_qcd_qcdLep1_phi->Sumw2(); 
 TH1F *h_dybb_qcd_qcdLep2_pt  = new TH1F("h_dybb_qcd_qcdLep2_pt","",12,25,205); h_dybb_qcd_qcdLep2_pt->Sumw2(); 
 TH1F *h_dybb_qcd_qcdLep2_eta = new TH1F("h_dybb_qcd_qcdLep2_eta","",20,-3,5); h_dybb_qcd_qcdLep2_eta->Sumw2(); 
 TH1F *h_dybb_qcd_qcdLep2_phi = new TH1F("h_dybb_qcd_qcdLep2_phi","",35,-3.4033,5.7594); h_dybb_qcd_qcdLep2_phi->Sumw2(); 

 std::cout<<"Nr. Entries: "<<nrEntries<<std::endl;
 for (Long64_t jentry=0; jentry<nrEntries;jentry++) {
  std::cout<<" entry "<<jentry<<std::endl;
  //if(jentry%500000==0) std::cout<<" entry "<<jentry<<std::endl;
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);
  if(!nb){std::cout<<"histoFiller.C can't fChain->GetEntry(jentry)"<<std::endl;}

  // split up the W samples
  //0 = ignore
  //1 = W+light
  //2 = W+c
  //3 = W+cc
  //4 = W+bb

  twoCs  = (abs(goodJ1_partonFlavour)==4 && abs(goodJ2_partonFlavour)==4);
  someCs = (abs(goodJ1_partonFlavour)==4 || abs(goodJ2_partonFlavour)==4);
  someBs = (abs(goodJ1_partonFlavour)==5 || abs(goodJ2_partonFlavour)==5);

  isWl  = !someBs && !someCs;
  isWc  = !someBs &&  someCs && !twoCs;
  isWcc = !someBs &&  someCs &&  twoCs;
  isWbb =  someBs;

  if( isMC && jType==1 && !isWl){ continue; }
  if( isMC && jType==2 && !isWc){ continue; }
  if( isMC && jType==3 && !isWcc){ continue; }
  if( isMC && jType==4 && !isWbb){ continue; }

  // take care of scaling for inclusive/exclusive WJets samples
  if( isMC && jType!=0){
   if( LHEProduct<5 ){
    std::cout<<"LHEProduct: "<<LHEProduct<<" Event: "<<EVENT<<std::endl;
    break;
   }
   // these numbers are the total number of events in the samples
   //  output from UWAnalysis/CRAB/MuNu/countEvents.sh 
   //      WJets, W1Jet, W2Jet, W3Jet, W4Jet 
   Float_t nrWnJ, nrW1J, nrW2J, nrW3J, nrW4J;
   nrWnJ = 75865454;
   nrW1J = 52593689;
   nrW2J = 64409521;
   nrW3J = 30358906;
   nrW4J = 13042592;

   std::vector<Float_t> ev;
   ev.push_back(nrWnJ);
   ev.push_back(nrW1J/0.175999);
   ev.push_back(nrW2J/0.0562617);
   ev.push_back(nrW3J/0.0168926);
   ev.push_back(nrW4J/0.00692218); 
   
   if( LHEProduct==5 ){nrEvents=ev[0];}
   else if( LHEProduct==6 ){nrEvents=ev[0]+ev[1];}
   else if( LHEProduct==7 ){nrEvents=ev[0]+ev[2];}
   else if( LHEProduct==8 ){nrEvents=ev[0]+ev[3];}
   else {nrEvents=ev[0]+ev[4];}
  }

  // define some filter variables
  oneGoodMuon = 
   HLT_IsoMu24_eta2p1_v_fired 
   && nrGoodMu==1 && nrVetoMu<=1 && nrVetoEle==0;
  twoGoodMuons = 
   HLT_IsoMu24_eta2p1_v_fired
   && nrGoodMu==2 && nrVetoMu<=2 && nrVetoEle==0;
  oneQCDMuon = 
   nrQCDMu==1 && nrVetoMu==0 && nrVetoEle==0;
  twoQCDMuons = 
   HLT_IsoMu24_eta2p1_v_fired
   && nrGoodMu==1 && nrQCDMu==1 && nrVetoMu<=1 && nrVetoEle==0;
 
  oneGoodElectron = 
   HLT_Ele27_WP80_fired
   && nrGoodEle==1 && nrVetoEle<=1 && nrVetoMu==0;
  twoGoodElectrons = 
   HLT_Ele27_WP80_fired
   && nrGoodEle==2 && nrVetoEle<=2 && nrVetoMu==0;
  oneQCDElectron = 
   nrQCDEle==1 && nrVetoEle==0 && nrVetoMu==0;
  twoQCDElectrons =  
   HLT_Ele27_WP80_fired
   && nrGoodEle==1 && nrQCDEle==1 && nrVetoEle<=1 && nrVetoMu==0;
  
  twoGoodLeptons = HLT_IsoMu24_eta2p1_v_fired
   && nrGoodEle==1 && nrVetoEle<=1
   && nrGoodMu==1  && nrVetoMu<=1;
  twoQCDLeptons = 
   nrQCDMu==1   && nrVetoMu==0 &&
   nrGoodEle==1 && nrVetoEle<=1;

  min2goodJs = 
   nrGoodJets>=2 && nrFwdJets==0;
  exactly2goodJs = 
   nrGoodJets==2 && nrFwdJets==0;
  min2goodBJs = 
   nrGoodJets>=2 && nrFwdJets==0
   && goodJ1_CSV>0.898 && goodJ2_CSV>0.898;
  exactly2goodBJs = 
   nrGoodJets==2 && nrFwdJets==0
   && goodJ1_CSV>0.898 && goodJ2_CSV>0.898;
  aGoodBJaFwdJ = 
   nrGoodJets==1 && nrFwdJets==1 
   && goodJ1_CSV>0.898;

  // SFs for CSV
  SF_top2BJs = goodJ1_SF_CSVT * goodJ2_SF_CSVT; 
  SF_top2BJs_errUp = goodJ1_SF_CSVT_errUp * goodJ2_SF_CSVT_errUp;
  SF_top2BJs_errDn = goodJ1_SF_CSVT_errDn * goodJ2_SF_CSVT_errDn; 
  SF_topBJ = goodJ1_SF_CSVT;

  //std::cout<<"SF_top2BJs: "<<SF_top2BJs
  // <<" SF_top2BJs_errUp: "<<SF_top2BJs_errUp
  // <<" SF_top2BJs_errDn: "<<SF_top2BJs_errDn<<std::endl;
  SF_CSVreweight = goodJ1_CSVreweight * goodJ2_CSVreweight;

  // Decide which filters our event passes
  // wbb Selection
  Bool_t pass_wbb_good=kFALSE;
  Bool_t pass_wbb_qcd=kFALSE;

  if( isMu ){
   if( oneGoodMuon && exactly2goodBJs ){ 
    pass_wbb_good=kTRUE; 
    nrEntries_wbb_good_postcut++;
    mt_good = mt_goodMuon;
    if( shift=="JESUp" ) { mt_good = mt_goodMuon_jesUp; }
    if( shift=="JESDown" ) { mt_good = mt_goodMuon_jesDn; }
    if( shift=="LESUp" )  { mt_good = mt_goodMuon_eesUp; }
    if( shift=="LESDown" )  { mt_good = 2 * mt_goodMuon - mt_goodMuon_eesUp; };// { mt = mt_goodMuon_eesDn; }
    if( shift=="UESUp" )  { mt_good = mt_goodMuon_uesUp; }
    if( shift=="UESDown" )  { mt_good = mt_goodMuon_uesDn; }
    lep_good_pt = goodMu_pt_vec->at(0);
    lep_good_eta = goodMu_eta_vec->at(0);
    lep_good_phi = goodMu_phi_vec->at(0);
    if( isMC ){
     SF_wbb_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0);
     SF_wbb_good_IDIsoHLT_errUp = SF_wbb_good_IDIsoHLT + SF_goodMu_IDIsoHLT_errUp->at(0);
     SF_wbb_good_IDIsoHLT_errDn = SF_wbb_good_IDIsoHLT - SF_goodMu_IDIsoHLT_errDn->at(0);
    }
   }
   if( oneQCDMuon && exactly2goodBJs ){ 
    pass_wbb_qcd=kTRUE; 
    nrEntries_wbb_qcd_postcut++;
    mt_qcd = mt_qcdMuon;
    if( shift=="JESUp" ) { mt_qcd = mt_qcdMuon_jesUp; }
    if( shift=="JESDown" ) { mt_qcd = mt_qcdMuon_jesDn; }
    if( shift=="LESUp" )  { mt_qcd = mt_qcdMuon_eesUp; }
    if( shift=="LESDown" )  { mt_qcd = 2 * mt_goodMuon - mt_qcdMuon_eesUp; }// { mt = mt_qcdMuon_eesDn; }
    if( shift=="UESUp" )  { mt_qcd = mt_qcdMuon_uesUp; }
    if( shift=="UESDown" )  { mt_qcd = mt_qcdMuon_uesDn; }
    lep_qcd_pt = qcdMu_pt;
    lep_qcd_eta = qcdMu_eta;
    lep_qcd_phi = qcdMu_phi;
    if( isMC ){
     SF_wbb_qcd_IDIsoHLT = SF_qcdMu_IDIso;
     SF_wbb_qcd_IDIsoHLT_errUp = SF_wbb_qcd_IDIsoHLT + SF_qcdMu_IDIso_errUp;
     SF_wbb_qcd_IDIsoHLT_errDn = SF_wbb_qcd_IDIsoHLT - SF_qcdMu_IDIso_errDn;
    }
   }
  }
  else {
   if( oneGoodElectron && exactly2goodBJs ){ 
    pass_wbb_good=kTRUE; 
    nrEntries_wbb_good_postcut++;
    mt_good = mt_goodElectron;
    if( shift=="JESUp" ) { mt_good = mt_goodElectron_jesUp; }
    if( shift=="JESDown" ) { mt_good = mt_goodElectron_jesDn; }
    if( shift=="LESUp" )  { mt_good = mt_goodElectron_eesUp; }
    if( shift=="LESDown" )  { mt_good = 2 * mt_goodElectron - mt_goodElectron_eesUp; };// { mt = mt_goodElectron_eesDn; }
    if( shift=="UESUp" )  { mt_good = mt_goodElectron_uesUp; }
    if( shift=="UESDown" )  { mt_good = mt_goodElectron_uesDn; }
    lep_good_pt = goodEle_pt_vec->at(0);
    lep_good_eta = goodEle_eta_vec->at(0);
    lep_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_wbb_good_IDIsoHLT = SF_goodEle_IDIsoHLT->at(0);
     SF_wbb_good_IDIsoHLT_errUp = SF_wbb_good_IDIsoHLT + SF_goodEle_IDIsoHLT_errUp->at(0);
     SF_wbb_good_IDIsoHLT_errDn = SF_wbb_good_IDIsoHLT - SF_goodEle_IDIsoHLT_errDn->at(0);
    }
   }
   if( oneQCDElectron && exactly2goodBJs ){ 
    pass_wbb_qcd=kTRUE; 
    nrEntries_wbb_qcd_postcut++;
    mt_qcd = mt_qcdElectron;
    if( shift=="JESUp" ) { mt_qcd = mt_qcdElectron_jesUp; }
    if( shift=="JESDown" ) { mt_qcd = mt_qcdElectron_jesDn; }
    if( shift=="LESUp" )  { mt_qcd = mt_qcdElectron_eesUp; }
    if( shift=="LESDown" )  { mt_qcd = 2 * mt_goodElectron - mt_qcdElectron_eesUp; }// { mt = mt_qcdElectron_eesDn; }
    if( shift=="UESUp" )  { mt_qcd = mt_qcdElectron_uesUp; }
    if( shift=="UESDown" )  { mt_qcd = mt_qcdElectron_uesDn; }
    lep_qcd_pt = qcdEle_pt;
    lep_qcd_eta = qcdEle_eta;
    lep_qcd_phi = qcdEle_phi;
    if( isMC ){
     SF_wbb_qcd_IDIsoHLT = SF_qcdEle_IDIso;
     SF_wbb_qcd_IDIsoHLT_errUp = SF_wbb_qcd_IDIsoHLT + SF_qcdEle_IDIso_errUp;
     SF_wbb_qcd_IDIsoHLT_errDn = SF_wbb_qcd_IDIsoHLT - SF_qcdEle_IDIso_errDn;
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good =       SF_wbb_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_good_CSVUp = SF_wbb_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errUp / nrEvents;
   weight_good_CSVDn = SF_wbb_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errDn / nrEvents;
   weight_good_EMuUp = SF_wbb_good_IDIsoHLT_errUp * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_good_EMuDn = SF_wbb_good_IDIsoHLT_errDn * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd =       SF_wbb_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd_CSVUp = SF_wbb_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errUp / nrEvents;
   weight_qcd_CSVDn = SF_wbb_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errDn / nrEvents;
   weight_qcd_EMuUp = SF_wbb_qcd_IDIsoHLT_errUp * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd_EMuDn = SF_wbb_qcd_IDIsoHLT_errDn * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
  }
  else{
   weight_good =       1.; 
   weight_good_CSVUp = 1.;
   weight_good_CSVDn = 1.;
   weight_good_EMuUp = 1.;
   weight_good_EMuDn = 1.;
   weight_qcd =       1.; 
   weight_qcd_CSVUp = 1.;
   weight_qcd_CSVDn = 1.;
   weight_qcd_EMuUp = 1.;
   weight_qcd_EMuDn = 1.;
  }

  if(pass_wbb_good){
   h_wbb_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_wbb_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_wbb_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_wbb_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_wbb_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_wbb_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_wbb_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_wbb_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_wbb_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_wbb_good_mt->Fill(mt_good,weight_good);
   h_wbb_good_goodLep_pt->Fill(lep_good_pt,weight_good);
   h_wbb_good_goodLep_eta->Fill(lep_good_eta,weight_good);
   h_wbb_good_goodLep_phi->Fill(lep_good_phi,weight_good);
 
   if( shift=="SFs" ){
    h_wbb_good_goodJ1_pt_CSVUp->Fill(goodJ1_pt,weight_good_CSVUp);
    h_wbb_good_goodJ2_pt_CSVUp->Fill(goodJ2_pt,weight_good_CSVUp);
    h_wbb_good_goodJ1_eta_CSVUp->Fill(goodJ1_eta,weight_good_CSVUp);
    h_wbb_good_goodJ2_eta_CSVUp->Fill(goodJ2_eta,weight_good_CSVUp);
    h_wbb_good_goodJ1_phi_CSVUp->Fill(goodJ1_phi,weight_good_CSVUp);
    h_wbb_good_goodJ2_phi_CSVUp->Fill(goodJ2_phi,weight_good_CSVUp);
    h_wbb_good_goodJ1_CSV_CSVUp->Fill(goodJ1_CSV,weight_good_CSVUp);
    h_wbb_good_goodJ2_CSV_CSVUp->Fill(goodJ2_CSV,weight_good_CSVUp);
    h_wbb_good_dRgoodJ1J2_CSVUp->Fill(dRgoodJ1J2,weight_good_CSVUp);
    h_wbb_good_mt_CSVUp->Fill(mt_good,weight_good_CSVUp);
    h_wbb_good_goodLep_pt_CSVUp->Fill(lep_good_pt,weight_good_CSVUp);
    h_wbb_good_goodLep_eta_CSVUp->Fill(lep_good_eta,weight_good_CSVUp);
    h_wbb_good_goodLep_phi_CSVUp->Fill(lep_good_phi,weight_good_CSVUp);
 
    h_wbb_good_goodJ1_pt_CSVDown->Fill(goodJ1_pt,weight_good_CSVDn);
    h_wbb_good_goodJ2_pt_CSVDown->Fill(goodJ2_pt,weight_good_CSVDn);
    h_wbb_good_goodJ1_eta_CSVDown->Fill(goodJ1_eta,weight_good_CSVDn);
    h_wbb_good_goodJ2_eta_CSVDown->Fill(goodJ2_eta,weight_good_CSVDn);
    h_wbb_good_goodJ1_phi_CSVDown->Fill(goodJ1_phi,weight_good_CSVDn);
    h_wbb_good_goodJ2_phi_CSVDown->Fill(goodJ2_phi,weight_good_CSVDn);
    h_wbb_good_goodJ1_CSV_CSVDown->Fill(goodJ1_CSV,weight_good_CSVDn);
    h_wbb_good_goodJ2_CSV_CSVDown->Fill(goodJ2_CSV,weight_good_CSVDn);
    h_wbb_good_dRgoodJ1J2_CSVDown->Fill(dRgoodJ1J2,weight_good_CSVDn);
    h_wbb_good_mt_CSVDown->Fill(mt_good,weight_good_CSVDn);
    h_wbb_good_goodLep_pt_CSVDown->Fill(lep_good_pt,weight_good_CSVDn);
    h_wbb_good_goodLep_eta_CSVDown->Fill(lep_good_eta,weight_good_CSVDn);
    h_wbb_good_goodLep_phi_CSVDown->Fill(lep_good_phi,weight_good_CSVDn);
 
    h_wbb_good_goodJ1_pt_EMuUp->Fill(goodJ1_pt,weight_good_EMuUp);
    h_wbb_good_goodJ2_pt_EMuUp->Fill(goodJ2_pt,weight_good_EMuUp);
    h_wbb_good_goodJ1_eta_EMuUp->Fill(goodJ1_eta,weight_good_EMuUp);
    h_wbb_good_goodJ2_eta_EMuUp->Fill(goodJ2_eta,weight_good_EMuUp);
    h_wbb_good_goodJ1_phi_EMuUp->Fill(goodJ1_phi,weight_good_EMuUp);
    h_wbb_good_goodJ2_phi_EMuUp->Fill(goodJ2_phi,weight_good_EMuUp);
    h_wbb_good_goodJ1_CSV_EMuUp->Fill(goodJ1_CSV,weight_good_EMuUp);
    h_wbb_good_goodJ2_CSV_EMuUp->Fill(goodJ2_CSV,weight_good_EMuUp);
    h_wbb_good_dRgoodJ1J2_EMuUp->Fill(dRgoodJ1J2,weight_good_EMuUp);
    h_wbb_good_mt_EMuUp->Fill(mt_good,weight_good_EMuUp);
    h_wbb_good_goodLep_pt_EMuUp->Fill(lep_good_pt,weight_good_EMuUp);
    h_wbb_good_goodLep_eta_EMuUp->Fill(lep_good_eta,weight_good_EMuUp);
    h_wbb_good_goodLep_phi_EMuUp->Fill(lep_good_phi,weight_good_EMuUp);

    h_wbb_good_goodJ1_pt_EMuDown->Fill(goodJ1_pt,weight_good_EMuDn);
    h_wbb_good_goodJ2_pt_EMuDown->Fill(goodJ2_pt,weight_good_EMuDn);
    h_wbb_good_goodJ1_eta_EMuDown->Fill(goodJ1_eta,weight_good_EMuDn);
    h_wbb_good_goodJ2_eta_EMuDown->Fill(goodJ2_eta,weight_good_EMuDn);
    h_wbb_good_goodJ1_phi_EMuDown->Fill(goodJ1_phi,weight_good_EMuDn);
    h_wbb_good_goodJ2_phi_EMuDown->Fill(goodJ2_phi,weight_good_EMuDn);
    h_wbb_good_goodJ1_CSV_EMuDown->Fill(goodJ1_CSV,weight_good_EMuDn);
    h_wbb_good_goodJ2_CSV_EMuDown->Fill(goodJ2_CSV,weight_good_EMuDn);
    h_wbb_good_dRgoodJ1J2_EMuDown->Fill(dRgoodJ1J2,weight_good_EMuDn);
    h_wbb_good_mt_EMuDown->Fill(mt_good,weight_good_EMuDn);
    h_wbb_good_goodLep_pt_EMuDown->Fill(lep_good_pt,weight_good_EMuDn);
    h_wbb_good_goodLep_eta_EMuDown->Fill(lep_good_eta,weight_good_EMuDn);
    h_wbb_good_goodLep_phi_EMuDown->Fill(lep_good_phi,weight_good_EMuDn);
   }
  } 
  if(pass_wbb_qcd){
   h_wbb_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_wbb_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_wbb_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_wbb_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_wbb_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_wbb_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_wbb_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_wbb_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_wbb_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_wbb_qcd_mt->Fill(mt_qcd,weight_qcd);
   h_wbb_qcd_qcdLep_pt->Fill(lep_qcd_pt,weight_qcd);
   h_wbb_qcd_qcdLep_eta->Fill(lep_qcd_eta,weight_qcd);
   h_wbb_qcd_qcdLep_phi->Fill(lep_qcd_phi,weight_qcd);
 
   if( shift=="SFs" ){
    h_wbb_qcd_goodJ1_pt_CSVUp->Fill(goodJ1_pt,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ2_pt_CSVUp->Fill(goodJ2_pt,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ1_eta_CSVUp->Fill(goodJ1_eta,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ2_eta_CSVUp->Fill(goodJ2_eta,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ1_phi_CSVUp->Fill(goodJ1_phi,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ2_phi_CSVUp->Fill(goodJ2_phi,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ1_CSV_CSVUp->Fill(goodJ1_CSV,weight_qcd_CSVUp);
    h_wbb_qcd_goodJ2_CSV_CSVUp->Fill(goodJ2_CSV,weight_qcd_CSVUp);
    h_wbb_qcd_dRgoodJ1J2_CSVUp->Fill(dRgoodJ1J2,weight_qcd_CSVUp);
    h_wbb_qcd_mt_CSVUp->Fill(mt_qcd,weight_qcd_CSVUp);
    h_wbb_qcd_qcdLep_pt_CSVUp->Fill(lep_qcd_pt,weight_qcd_CSVUp);
    h_wbb_qcd_qcdLep_eta_CSVUp->Fill(lep_qcd_eta,weight_qcd_CSVUp);
    h_wbb_qcd_qcdLep_phi_CSVUp->Fill(lep_qcd_phi,weight_qcd_CSVUp);
 
    h_wbb_qcd_goodJ1_pt_CSVDown->Fill(goodJ1_pt,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ2_pt_CSVDown->Fill(goodJ2_pt,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ1_eta_CSVDown->Fill(goodJ1_eta,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ2_eta_CSVDown->Fill(goodJ2_eta,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ1_phi_CSVDown->Fill(goodJ1_phi,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ2_phi_CSVDown->Fill(goodJ2_phi,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ1_CSV_CSVDown->Fill(goodJ1_CSV,weight_qcd_CSVDn);
    h_wbb_qcd_goodJ2_CSV_CSVDown->Fill(goodJ2_CSV,weight_qcd_CSVDn);
    h_wbb_qcd_dRgoodJ1J2_CSVDown->Fill(dRgoodJ1J2,weight_qcd_CSVDn);
    h_wbb_qcd_mt_CSVDown->Fill(mt_qcd,weight_qcd_CSVDn);
    h_wbb_qcd_qcdLep_pt_CSVDown->Fill(lep_qcd_pt,weight_qcd_CSVDn);
    h_wbb_qcd_qcdLep_eta_CSVDown->Fill(lep_qcd_eta,weight_qcd_CSVDn);
    h_wbb_qcd_qcdLep_phi_CSVDown->Fill(lep_qcd_phi,weight_qcd_CSVDn);
 
    h_wbb_qcd_goodJ1_pt_EMuUp->Fill(goodJ1_pt,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ2_pt_EMuUp->Fill(goodJ2_pt,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ1_eta_EMuUp->Fill(goodJ1_eta,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ2_eta_EMuUp->Fill(goodJ2_eta,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ1_phi_EMuUp->Fill(goodJ1_phi,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ2_phi_EMuUp->Fill(goodJ2_phi,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ1_CSV_EMuUp->Fill(goodJ1_CSV,weight_qcd_EMuUp);
    h_wbb_qcd_goodJ2_CSV_EMuUp->Fill(goodJ2_CSV,weight_qcd_EMuUp);
    h_wbb_qcd_dRgoodJ1J2_EMuUp->Fill(dRgoodJ1J2,weight_qcd_EMuUp);
    h_wbb_qcd_mt_EMuUp->Fill(mt_qcd,weight_qcd_EMuUp);
    h_wbb_qcd_qcdLep_pt_EMuUp->Fill(lep_qcd_pt,weight_qcd_EMuUp);
    h_wbb_qcd_qcdLep_eta_EMuUp->Fill(lep_qcd_eta,weight_qcd_EMuUp);
    h_wbb_qcd_qcdLep_phi_EMuUp->Fill(lep_qcd_phi,weight_qcd_EMuUp);

    h_wbb_qcd_goodJ1_pt_EMuDown->Fill(goodJ1_pt,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ2_pt_EMuDown->Fill(goodJ2_pt,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ1_eta_EMuDown->Fill(goodJ1_eta,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ2_eta_EMuDown->Fill(goodJ2_eta,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ1_phi_EMuDown->Fill(goodJ1_phi,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ2_phi_EMuDown->Fill(goodJ2_phi,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ1_CSV_EMuDown->Fill(goodJ1_CSV,weight_qcd_EMuDn);
    h_wbb_qcd_goodJ2_CSV_EMuDown->Fill(goodJ2_CSV,weight_qcd_EMuDn);
    h_wbb_qcd_dRgoodJ1J2_EMuDown->Fill(dRgoodJ1J2,weight_qcd_EMuDn);
    h_wbb_qcd_mt_EMuDown->Fill(mt_qcd,weight_qcd_EMuDn);
    h_wbb_qcd_qcdLep_pt_EMuDown->Fill(lep_qcd_pt,weight_qcd_EMuDn);
    h_wbb_qcd_qcdLep_eta_EMuDown->Fill(lep_qcd_eta,weight_qcd_EMuDn);
    h_wbb_qcd_qcdLep_phi_EMuDown->Fill(lep_qcd_phi,weight_qcd_EMuDn);
   }
  } 
 
  // tt3j Selection
  Bool_t pass_tt3j_good=kFALSE;
  Bool_t pass_tt3j_qcd=kFALSE;

  if( isMu ){
   if( oneGoodMuon && min2goodBJs ){ 
    pass_tt3j_good=kTRUE; 
    nrEntries_tt3j_good_postcut++;
    mt_good = mt_goodMuon;
    if( shift=="JESUp" ) { mt_good = mt_goodMuon_jesUp; }
    if( shift=="JESDown" ) { mt_good = mt_goodMuon_jesDn; }
    if( shift=="LESUp" )  { mt_good = mt_goodMuon_eesUp; }
    if( shift=="LESDown" )  { mt_good = 2 * mt_goodMuon - mt_goodMuon_eesUp; };// { mt = mt_goodMuon_eesDn; }
    if( shift=="UESUp" )  { mt_good = mt_goodMuon_uesUp; }
    if( shift=="UESDown" )  { mt_good = mt_goodMuon_uesDn; }
    lep_good_pt = goodMu_pt_vec->at(0);
    lep_good_eta = goodMu_eta_vec->at(0);
    lep_good_phi = goodMu_phi_vec->at(0);
    if( isMC ){
     SF_tt3j_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0);
     SF_tt3j_good_IDIsoHLT_errUp = SF_tt3j_good_IDIsoHLT + SF_goodMu_IDIsoHLT_errUp->at(0);
     SF_tt3j_good_IDIsoHLT_errDn = SF_tt3j_good_IDIsoHLT - SF_goodMu_IDIsoHLT_errDn->at(0);
    }
   }
   if( oneQCDMuon && min2goodBJs ){ 
    pass_tt3j_qcd=kTRUE; 
    nrEntries_tt3j_qcd_postcut++;
    mt_qcd = mt_qcdMuon;
    if( shift=="JESUp" ) { mt_qcd = mt_qcdMuon_jesUp; }
    if( shift=="JESDown" ) { mt_qcd = mt_qcdMuon_jesDn; }
    if( shift=="LESUp" )  { mt_qcd = mt_qcdMuon_eesUp; }
    if( shift=="LESDown" )  { mt_qcd = 2 * mt_goodMuon - mt_qcdMuon_eesUp; }// { mt = mt_qcdMuon_eesDn; }
    if( shift=="UESUp" )  { mt_qcd = mt_qcdMuon_uesUp; }
    if( shift=="UESDown" )  { mt_qcd = mt_qcdMuon_uesDn; }
    lep_qcd_pt = qcdMu_pt;
    lep_qcd_eta = qcdMu_eta;
    lep_qcd_phi = qcdMu_phi;
    if( isMC ){
     SF_tt3j_qcd_IDIsoHLT = SF_qcdMu_IDIso;
     SF_tt3j_qcd_IDIsoHLT_errUp = SF_tt3j_qcd_IDIsoHLT + SF_qcdMu_IDIso_errUp;
     SF_tt3j_qcd_IDIsoHLT_errDn = SF_tt3j_qcd_IDIsoHLT - SF_qcdMu_IDIso_errDn;
    }
   }
  }
  else {
   if( oneGoodElectron && min2goodBJs ){ 
    pass_tt3j_good=kTRUE; 
    nrEntries_tt3j_good_postcut++;
    mt_good = mt_goodElectron;
    if( shift=="JESUp" ) { mt_good = mt_goodElectron_jesUp; }
    if( shift=="JESDown" ) { mt_good = mt_goodElectron_jesDn; }
    if( shift=="LESUp" )  { mt_good = mt_goodElectron_eesUp; }
    if( shift=="LESDown" )  { mt_good = 2 * mt_goodElectron - mt_goodElectron_eesUp; };// { mt = mt_goodElectron_eesDn; }
    if( shift=="UESUp" )  { mt_good = mt_goodElectron_uesUp; }
    if( shift=="UESDown" )  { mt_good = mt_goodElectron_uesDn; }
    lep_good_pt = goodEle_pt_vec->at(0);
    lep_good_eta = goodEle_eta_vec->at(0);
    lep_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_tt3j_good_IDIsoHLT = SF_goodEle_IDIsoHLT->at(0);
     SF_tt3j_good_IDIsoHLT_errUp = SF_tt3j_good_IDIsoHLT + SF_goodEle_IDIsoHLT_errUp->at(0);
     SF_tt3j_good_IDIsoHLT_errDn = SF_tt3j_good_IDIsoHLT - SF_goodEle_IDIsoHLT_errDn->at(0);
    }
   }
   if( oneQCDElectron && exactly2goodBJs ){ 
    pass_tt3j_qcd=kTRUE; 
    nrEntries_tt3j_qcd_postcut++;
    mt_qcd = mt_qcdElectron;
    if( shift=="JESUp" ) { mt_qcd = mt_qcdElectron_jesUp; }
    if( shift=="JESDown" ) { mt_qcd = mt_qcdElectron_jesDn; }
    if( shift=="LESUp" )  { mt_qcd = mt_qcdElectron_eesUp; }
    if( shift=="LESDown" )  { mt_qcd = 2 * mt_goodElectron - mt_qcdElectron_eesUp; }// { mt = mt_qcdElectron_eesDn; }
    if( shift=="UESUp" )  { mt_qcd = mt_qcdElectron_uesUp; }
    if( shift=="UESDown" )  { mt_qcd = mt_qcdElectron_uesDn; }
    lep_qcd_pt = qcdEle_pt;
    lep_qcd_eta = qcdEle_eta;
    lep_qcd_phi = qcdEle_phi;
    if( isMC ){
     SF_tt3j_qcd_IDIsoHLT = SF_qcdEle_IDIso;
     SF_tt3j_qcd_IDIsoHLT_errUp = SF_tt3j_qcd_IDIsoHLT + SF_qcdEle_IDIso_errUp;
     SF_tt3j_qcd_IDIsoHLT_errDn = SF_tt3j_qcd_IDIsoHLT - SF_qcdEle_IDIso_errDn;
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good =       SF_tt3j_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_good_CSVUp = SF_tt3j_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errUp / nrEvents;
   weight_good_CSVDn = SF_tt3j_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errDn / nrEvents;
   weight_good_EMuUp = SF_tt3j_good_IDIsoHLT_errUp * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_good_EMuDn = SF_tt3j_good_IDIsoHLT_errDn * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd =       SF_tt3j_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd_CSVUp = SF_tt3j_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errUp / nrEvents;
   weight_qcd_CSVDn = SF_tt3j_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errDn / nrEvents;
   weight_qcd_EMuUp = SF_tt3j_qcd_IDIsoHLT_errUp * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd_EMuDn = SF_tt3j_qcd_IDIsoHLT_errDn * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
  }
  else{
   weight_good =       1.; 
   weight_good_CSVUp = 1.;
   weight_good_CSVDn = 1.;
   weight_good_EMuUp = 1.;
   weight_good_EMuDn = 1.;
   weight_qcd =       1.; 
   weight_qcd_CSVUp = 1.;
   weight_qcd_CSVDn = 1.;
   weight_qcd_EMuUp = 1.;
   weight_qcd_EMuDn = 1.;
  }

  if(pass_tt3j_good){
   h_tt3j_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_tt3j_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_tt3j_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_tt3j_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_tt3j_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_tt3j_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_tt3j_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_tt3j_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_tt3j_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_tt3j_good_mt->Fill(mt_good,weight_good);
   h_tt3j_good_goodLep_pt->Fill(lep_good_pt,weight_good);
   h_tt3j_good_goodLep_eta->Fill(lep_good_eta,weight_good);
   h_tt3j_good_goodLep_phi->Fill(lep_good_phi,weight_good);
 
   if( shift=="SFs" ){
    h_tt3j_good_goodJ1_pt_CSVUp->Fill(goodJ1_pt,weight_good_CSVUp);
    h_tt3j_good_goodJ2_pt_CSVUp->Fill(goodJ2_pt,weight_good_CSVUp);
    h_tt3j_good_goodJ1_eta_CSVUp->Fill(goodJ1_eta,weight_good_CSVUp);
    h_tt3j_good_goodJ2_eta_CSVUp->Fill(goodJ2_eta,weight_good_CSVUp);
    h_tt3j_good_goodJ1_phi_CSVUp->Fill(goodJ1_phi,weight_good_CSVUp);
    h_tt3j_good_goodJ2_phi_CSVUp->Fill(goodJ2_phi,weight_good_CSVUp);
    h_tt3j_good_goodJ1_CSV_CSVUp->Fill(goodJ1_CSV,weight_good_CSVUp);
    h_tt3j_good_goodJ2_CSV_CSVUp->Fill(goodJ2_CSV,weight_good_CSVUp);
    h_tt3j_good_dRgoodJ1J2_CSVUp->Fill(dRgoodJ1J2,weight_good_CSVUp);
    h_tt3j_good_mt_CSVUp->Fill(mt_good,weight_good_CSVUp);
    h_tt3j_good_goodLep_pt_CSVUp->Fill(lep_good_pt,weight_good_CSVUp);
    h_tt3j_good_goodLep_eta_CSVUp->Fill(lep_good_eta,weight_good_CSVUp);
    h_tt3j_good_goodLep_phi_CSVUp->Fill(lep_good_phi,weight_good_CSVUp);
 
    h_tt3j_good_goodJ1_pt_CSVDown->Fill(goodJ1_pt,weight_good_CSVDn);
    h_tt3j_good_goodJ2_pt_CSVDown->Fill(goodJ2_pt,weight_good_CSVDn);
    h_tt3j_good_goodJ1_eta_CSVDown->Fill(goodJ1_eta,weight_good_CSVDn);
    h_tt3j_good_goodJ2_eta_CSVDown->Fill(goodJ2_eta,weight_good_CSVDn);
    h_tt3j_good_goodJ1_phi_CSVDown->Fill(goodJ1_phi,weight_good_CSVDn);
    h_tt3j_good_goodJ2_phi_CSVDown->Fill(goodJ2_phi,weight_good_CSVDn);
    h_tt3j_good_goodJ1_CSV_CSVDown->Fill(goodJ1_CSV,weight_good_CSVDn);
    h_tt3j_good_goodJ2_CSV_CSVDown->Fill(goodJ2_CSV,weight_good_CSVDn);
    h_tt3j_good_dRgoodJ1J2_CSVDown->Fill(dRgoodJ1J2,weight_good_CSVDn);
    h_tt3j_good_mt_CSVDown->Fill(mt_good,weight_good_CSVDn);
    h_tt3j_good_goodLep_pt_CSVDown->Fill(lep_good_pt,weight_good_CSVDn);
    h_tt3j_good_goodLep_eta_CSVDown->Fill(lep_good_eta,weight_good_CSVDn);
    h_tt3j_good_goodLep_phi_CSVDown->Fill(lep_good_phi,weight_good_CSVDn);
 
    h_tt3j_good_goodJ1_pt_EMuUp->Fill(goodJ1_pt,weight_good_EMuUp);
    h_tt3j_good_goodJ2_pt_EMuUp->Fill(goodJ2_pt,weight_good_EMuUp);
    h_tt3j_good_goodJ1_eta_EMuUp->Fill(goodJ1_eta,weight_good_EMuUp);
    h_tt3j_good_goodJ2_eta_EMuUp->Fill(goodJ2_eta,weight_good_EMuUp);
    h_tt3j_good_goodJ1_phi_EMuUp->Fill(goodJ1_phi,weight_good_EMuUp);
    h_tt3j_good_goodJ2_phi_EMuUp->Fill(goodJ2_phi,weight_good_EMuUp);
    h_tt3j_good_goodJ1_CSV_EMuUp->Fill(goodJ1_CSV,weight_good_EMuUp);
    h_tt3j_good_goodJ2_CSV_EMuUp->Fill(goodJ2_CSV,weight_good_EMuUp);
    h_tt3j_good_dRgoodJ1J2_EMuUp->Fill(dRgoodJ1J2,weight_good_EMuUp);
    h_tt3j_good_mt_EMuUp->Fill(mt_good,weight_good_EMuUp);
    h_tt3j_good_goodLep_pt_EMuUp->Fill(lep_good_pt,weight_good_EMuUp);
    h_tt3j_good_goodLep_eta_EMuUp->Fill(lep_good_eta,weight_good_EMuUp);
    h_tt3j_good_goodLep_phi_EMuUp->Fill(lep_good_phi,weight_good_EMuUp);

    h_tt3j_good_goodJ1_pt_EMuDown->Fill(goodJ1_pt,weight_good_EMuDn);
    h_tt3j_good_goodJ2_pt_EMuDown->Fill(goodJ2_pt,weight_good_EMuDn);
    h_tt3j_good_goodJ1_eta_EMuDown->Fill(goodJ1_eta,weight_good_EMuDn);
    h_tt3j_good_goodJ2_eta_EMuDown->Fill(goodJ2_eta,weight_good_EMuDn);
    h_tt3j_good_goodJ1_phi_EMuDown->Fill(goodJ1_phi,weight_good_EMuDn);
    h_tt3j_good_goodJ2_phi_EMuDown->Fill(goodJ2_phi,weight_good_EMuDn);
    h_tt3j_good_goodJ1_CSV_EMuDown->Fill(goodJ1_CSV,weight_good_EMuDn);
    h_tt3j_good_goodJ2_CSV_EMuDown->Fill(goodJ2_CSV,weight_good_EMuDn);
    h_tt3j_good_dRgoodJ1J2_EMuDown->Fill(dRgoodJ1J2,weight_good_EMuDn);
    h_tt3j_good_mt_EMuDown->Fill(mt_good,weight_good_EMuDn);
    h_tt3j_good_goodLep_pt_EMuDown->Fill(lep_good_pt,weight_good_EMuDn);
    h_tt3j_good_goodLep_eta_EMuDown->Fill(lep_good_eta,weight_good_EMuDn);
    h_tt3j_good_goodLep_phi_EMuDown->Fill(lep_good_phi,weight_good_EMuDn);
   }
  } 
  if(pass_tt3j_qcd){
   h_tt3j_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_tt3j_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_tt3j_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_tt3j_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_tt3j_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_tt3j_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_tt3j_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_tt3j_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_tt3j_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_tt3j_qcd_mt->Fill(mt_qcd,weight_qcd);
   h_tt3j_qcd_qcdLep_pt->Fill(lep_qcd_pt,weight_qcd);
   h_tt3j_qcd_qcdLep_eta->Fill(lep_qcd_eta,weight_qcd);
   h_tt3j_qcd_qcdLep_phi->Fill(lep_qcd_phi,weight_qcd);
 
   if( shift=="SFs" ){
    h_tt3j_qcd_goodJ1_pt_CSVUp->Fill(goodJ1_pt,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ2_pt_CSVUp->Fill(goodJ2_pt,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ1_eta_CSVUp->Fill(goodJ1_eta,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ2_eta_CSVUp->Fill(goodJ2_eta,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ1_phi_CSVUp->Fill(goodJ1_phi,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ2_phi_CSVUp->Fill(goodJ2_phi,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ1_CSV_CSVUp->Fill(goodJ1_CSV,weight_qcd_CSVUp);
    h_tt3j_qcd_goodJ2_CSV_CSVUp->Fill(goodJ2_CSV,weight_qcd_CSVUp);
    h_tt3j_qcd_dRgoodJ1J2_CSVUp->Fill(dRgoodJ1J2,weight_qcd_CSVUp);
    h_tt3j_qcd_mt_CSVUp->Fill(mt_qcd,weight_qcd_CSVUp);
    h_tt3j_qcd_qcdLep_pt_CSVUp->Fill(lep_qcd_pt,weight_qcd_CSVUp);
    h_tt3j_qcd_qcdLep_eta_CSVUp->Fill(lep_qcd_eta,weight_qcd_CSVUp);
    h_tt3j_qcd_qcdLep_phi_CSVUp->Fill(lep_qcd_phi,weight_qcd_CSVUp);
 
    h_tt3j_qcd_goodJ1_pt_CSVDown->Fill(goodJ1_pt,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ2_pt_CSVDown->Fill(goodJ2_pt,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ1_eta_CSVDown->Fill(goodJ1_eta,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ2_eta_CSVDown->Fill(goodJ2_eta,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ1_phi_CSVDown->Fill(goodJ1_phi,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ2_phi_CSVDown->Fill(goodJ2_phi,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ1_CSV_CSVDown->Fill(goodJ1_CSV,weight_qcd_CSVDn);
    h_tt3j_qcd_goodJ2_CSV_CSVDown->Fill(goodJ2_CSV,weight_qcd_CSVDn);
    h_tt3j_qcd_dRgoodJ1J2_CSVDown->Fill(dRgoodJ1J2,weight_qcd_CSVDn);
    h_tt3j_qcd_mt_CSVDown->Fill(mt_qcd,weight_qcd_CSVDn);
    h_tt3j_qcd_qcdLep_pt_CSVDown->Fill(lep_qcd_pt,weight_qcd_CSVDn);
    h_tt3j_qcd_qcdLep_eta_CSVDown->Fill(lep_qcd_eta,weight_qcd_CSVDn);
    h_tt3j_qcd_qcdLep_phi_CSVDown->Fill(lep_qcd_phi,weight_qcd_CSVDn);
 
    h_tt3j_qcd_goodJ1_pt_EMuUp->Fill(goodJ1_pt,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ2_pt_EMuUp->Fill(goodJ2_pt,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ1_eta_EMuUp->Fill(goodJ1_eta,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ2_eta_EMuUp->Fill(goodJ2_eta,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ1_phi_EMuUp->Fill(goodJ1_phi,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ2_phi_EMuUp->Fill(goodJ2_phi,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ1_CSV_EMuUp->Fill(goodJ1_CSV,weight_qcd_EMuUp);
    h_tt3j_qcd_goodJ2_CSV_EMuUp->Fill(goodJ2_CSV,weight_qcd_EMuUp);
    h_tt3j_qcd_dRgoodJ1J2_EMuUp->Fill(dRgoodJ1J2,weight_qcd_EMuUp);
    h_tt3j_qcd_mt_EMuUp->Fill(mt_qcd,weight_qcd_EMuUp);
    h_tt3j_qcd_qcdLep_pt_EMuUp->Fill(lep_qcd_pt,weight_qcd_EMuUp);
    h_tt3j_qcd_qcdLep_eta_EMuUp->Fill(lep_qcd_eta,weight_qcd_EMuUp);
    h_tt3j_qcd_qcdLep_phi_EMuUp->Fill(lep_qcd_phi,weight_qcd_EMuUp);

    h_tt3j_qcd_goodJ1_pt_EMuDown->Fill(goodJ1_pt,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ2_pt_EMuDown->Fill(goodJ2_pt,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ1_eta_EMuDown->Fill(goodJ1_eta,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ2_eta_EMuDown->Fill(goodJ2_eta,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ1_phi_EMuDown->Fill(goodJ1_phi,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ2_phi_EMuDown->Fill(goodJ2_phi,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ1_CSV_EMuDown->Fill(goodJ1_CSV,weight_qcd_EMuDn);
    h_tt3j_qcd_goodJ2_CSV_EMuDown->Fill(goodJ2_CSV,weight_qcd_EMuDn);
    h_tt3j_qcd_dRgoodJ1J2_EMuDown->Fill(dRgoodJ1J2,weight_qcd_EMuDn);
    h_tt3j_qcd_mt_EMuDown->Fill(mt_qcd,weight_qcd_EMuDn);
    h_tt3j_qcd_qcdLep_pt_EMuDown->Fill(lep_qcd_pt,weight_qcd_EMuDn);
    h_tt3j_qcd_qcdLep_eta_EMuDown->Fill(lep_qcd_eta,weight_qcd_EMuDn);
    h_tt3j_qcd_qcdLep_phi_EMuDown->Fill(lep_qcd_phi,weight_qcd_EMuDn);
   }
  } 

  // tt1m1e Selection
  Bool_t pass_tt1m1e_good=kFALSE;
  Bool_t pass_tt1m1e_qcd=kFALSE;

  if( isMu ){
   if( twoGoodLeptons && exactly2goodBJs ){ 
    pass_tt1m1e_good=kTRUE; 
    nrEntries_tt1m1e_good_postcut++;
    mt_good = mt_goodMuon;
    if( shift=="JESUp" ) { mt_good = mt_goodMuon_jesUp; }
    if( shift=="JESDown" ) { mt_good = mt_goodMuon_jesDn; }
    if( shift=="LESUp" )  { mt_good = mt_goodMuon_eesUp; }
    if( shift=="LESDown" )  { mt_good = 2 * mt_goodMuon - mt_goodMuon_eesUp; };// { mt = mt_goodMuon_eesDn; }
    if( shift=="UESUp" )  { mt_good = mt_goodMuon_uesUp; }
    if( shift=="UESDown" )  { mt_good = mt_goodMuon_uesDn; }
    mu_good_pt = goodMu_pt_vec->at(0);
    mu_good_eta = goodMu_eta_vec->at(0);
    mu_good_phi = goodMu_phi_vec->at(0);
    ele_good_pt = goodEle_pt_vec->at(0);
    ele_good_eta = goodEle_eta_vec->at(0);
    ele_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_tt1m1e_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0) * SF_goodEle_IDIsoHLT->at(0); // SF_goodEle_IDIso->at(0)
     SF_tt1m1e_good_IDIsoHLT_errUp = SF_tt1m1e_good_IDIsoHLT 
      + sqrt( pow(SF_goodMu_IDIsoHLT_errUp->at(0),2) + pow(SF_goodEle_IDIsoHLT_errUp->at(0),2) );
      // + sqrt(SF_goodMu_IDIsoHLT_errUp->at(0)**2 + SF_goodEle_IDIso_errUp->at(0)**2);
     SF_tt1m1e_good_IDIsoHLT_errDn = SF_tt1m1e_good_IDIsoHLT
      - sqrt( pow(SF_goodMu_IDIsoHLT_errDn->at(0),2) + pow(SF_goodEle_IDIsoHLT_errDn->at(0),2) );
      // - sqrt(SF_goodMu_IDIsoHLT_errDn->at(0)**2 + SF_goodEle_IDIso_errDn->at(0)**2);
    }
   }
   if( twoQCDLeptons && exactly2goodBJs ){ 
    pass_tt1m1e_qcd=kTRUE; 
    nrEntries_tt1m1e_qcd_postcut++;
    mt_qcd = mt_qcdMuon;
    if( shift=="JESUp" ) { mt_qcd = mt_qcdMuon_jesUp; }
    if( shift=="JESDown" ) { mt_qcd = mt_qcdMuon_jesDn; }
    if( shift=="LESUp" )  { mt_qcd = mt_qcdMuon_eesUp; }
    if( shift=="LESDown" )  { mt_qcd = 2 * mt_goodMuon - mt_qcdMuon_eesUp; }// { mt = mt_qcdMuon_eesDn; }
    if( shift=="UESUp" )  { mt_qcd = mt_qcdMuon_uesUp; }
    if( shift=="UESDown" )  { mt_qcd = mt_qcdMuon_uesDn; }
    mu_qcd_pt = qcdMu_pt;
    mu_qcd_eta = qcdMu_eta;
    mu_qcd_phi = qcdMu_phi;
    ele_good_pt = goodEle_pt_vec->at(0);
    ele_good_eta = goodEle_eta_vec->at(0);
    ele_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_tt1m1e_qcd_IDIsoHLT = SF_qcdMu_IDIso * SF_goodEle_IDIsoHLT->at(0); // SF_goodEle_IDIso->at(0);
     SF_tt1m1e_qcd_IDIsoHLT_errUp = SF_tt1m1e_qcd_IDIsoHLT 
      + sqrt( pow(SF_qcdMu_IDIso_errUp,2) + pow(SF_goodEle_IDIsoHLT_errUp->at(0),2) );
      // + sqrt(SF_qcdMu_IDIso_errUp**2 + SF_goodEle_IDIso_errUp->at(0)**2);
     SF_tt1m1e_qcd_IDIsoHLT_errDn = SF_tt1m1e_qcd_IDIsoHLT
      - sqrt( pow(SF_qcdMu_IDIso_errDn,2) + pow(SF_goodEle_IDIsoHLT_errDn->at(0),2));
      // - sqrt(SF_qcdMu_IDIso_errDn**2 + SF_goodEle_IDIso_errDn->at(0)**2);
    }
   }
  }
  else{
   if( twoGoodLeptons && exactly2goodBJs ){ 
    pass_tt1m1e_good=kTRUE; 
    nrEntries_tt1m1e_good_postcut++;
    mt_good = mt_goodElectron;
    if( shift=="JESUp" ) { mt_good = mt_goodElectron_jesUp; }
    if( shift=="JESDown" ) { mt_good = mt_goodElectron_jesDn; }
    if( shift=="LESUp" )  { mt_good = mt_goodElectron_eesUp; }
    if( shift=="LESDown" )  { mt_good = 2 * mt_goodElectron - mt_goodElectron_eesUp; };// { mt = mt_goodElectron_eesDn; }
    if( shift=="UESUp" )  { mt_good = mt_goodElectron_uesUp; }
    if( shift=="UESDown" )  { mt_good = mt_goodElectron_uesDn; }
    mu_good_pt = goodMu_pt_vec->at(0);
    mu_good_eta = goodMu_eta_vec->at(0);
    mu_good_phi = goodMu_phi_vec->at(0);
    ele_good_pt = goodEle_pt_vec->at(0);
    ele_good_eta = goodEle_eta_vec->at(0);
    ele_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_tt1m1e_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0) * SF_goodEle_IDIsoHLT->at(0); // SF_goodEle_IDIso->at(0)
     SF_tt1m1e_good_IDIsoHLT_errUp = SF_tt1m1e_good_IDIsoHLT 
      + sqrt( pow(SF_goodMu_IDIsoHLT_errUp->at(0),2) + pow(SF_goodEle_IDIsoHLT_errUp->at(0),2) );
      // + sqrt(SF_goodMu_IDIsoHLT_errUp->at(0)**2 + SF_goodEle_IDIso_errUp->at(0)**2);
     SF_tt1m1e_good_IDIsoHLT_errDn = SF_tt1m1e_good_IDIsoHLT
      - sqrt( pow(SF_goodMu_IDIsoHLT_errDn->at(0),2) + pow(SF_goodEle_IDIsoHLT_errDn->at(0),2) );
      // - sqrt(SF_goodMu_IDIsoHLT_errDn->at(0)**2 + SF_goodEle_IDIso_errDn->at(0)**2);
    }
   }
   if( twoQCDLeptons && exactly2goodBJs ){ 
    pass_tt1m1e_qcd=kTRUE; 
    nrEntries_tt1m1e_qcd_postcut++;
    mt_qcd = mt_qcdElectron;
    if( shift=="JESUp" ) { mt_qcd = mt_qcdElectron_jesUp; }
    if( shift=="JESDown" ) { mt_qcd = mt_qcdElectron_jesDn; }
    if( shift=="LESUp" )  { mt_qcd = mt_qcdElectron_eesUp; }
    if( shift=="LESDown" )  { mt_qcd = 2 * mt_goodElectron - mt_qcdElectron_eesUp; }// { mt = mt_qcdElectron_eesDn; }
    if( shift=="UESUp" )  { mt_qcd = mt_qcdElectron_uesUp; }
    if( shift=="UESDown" )  { mt_qcd = mt_qcdElectron_uesDn; }
    mu_qcd_pt = qcdMu_pt;
    mu_qcd_eta = qcdMu_eta;
    mu_qcd_phi = qcdMu_phi;
    ele_good_pt = goodEle_pt_vec->at(0);
    ele_good_eta = goodEle_eta_vec->at(0);
    ele_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_tt1m1e_qcd_IDIsoHLT = SF_qcdMu_IDIso * SF_goodEle_IDIsoHLT->at(0); // SF_goodEle_IDIso->at(0);
     SF_tt1m1e_qcd_IDIsoHLT_errUp = SF_tt1m1e_qcd_IDIsoHLT 
      + sqrt( pow(SF_qcdMu_IDIso_errUp,2) + pow(SF_goodEle_IDIsoHLT_errUp->at(0),2) );
      // + sqrt(SF_qcdMu_IDIso_errUp**2 + SF_goodEle_IDIso_errUp->at(0)**2);
     SF_tt1m1e_qcd_IDIsoHLT_errDn = SF_tt1m1e_qcd_IDIsoHLT
      - sqrt( pow(SF_qcdMu_IDIso_errDn,2) + pow(SF_goodEle_IDIsoHLT_errDn->at(0),2));
      // - sqrt(SF_qcdMu_IDIso_errDn**2 + SF_goodEle_IDIso_errDn->at(0)**2);
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good =       SF_tt1m1e_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_good_CSVUp = SF_tt1m1e_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errUp / nrEvents;
   weight_good_CSVDn = SF_tt1m1e_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errDn / nrEvents;
   weight_good_EMuUp = SF_tt1m1e_good_IDIsoHLT_errUp * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_good_EMuDn = SF_tt1m1e_good_IDIsoHLT_errDn * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd =       SF_tt1m1e_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd_CSVUp = SF_tt1m1e_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errUp / nrEvents;
   weight_qcd_CSVDn = SF_tt1m1e_qcd_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs_errDn / nrEvents;
   weight_qcd_EMuUp = SF_tt1m1e_qcd_IDIsoHLT_errUp * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd_EMuDn = SF_tt1m1e_qcd_IDIsoHLT_errDn * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
  }
  else{
   weight_good =       1.; 
   weight_good_CSVUp = 1.;
   weight_good_CSVDn = 1.;
   weight_good_EMuUp = 1.;
   weight_good_EMuDn = 1.;
   weight_qcd =       1.; 
   weight_qcd_CSVUp = 1.;
   weight_qcd_CSVDn = 1.;
   weight_qcd_EMuUp = 1.;
   weight_qcd_EMuDn = 1.;
  }

  if(pass_tt1m1e_good){
   h_tt1m1e_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_tt1m1e_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_tt1m1e_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_tt1m1e_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_tt1m1e_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_tt1m1e_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_tt1m1e_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_tt1m1e_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_tt1m1e_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_tt1m1e_good_mt->Fill(mt_good,weight_good);
   h_tt1m1e_good_goodLep_pt->Fill(lep_good_pt,weight_good);
   h_tt1m1e_good_goodLep_eta->Fill(lep_good_eta,weight_good);
   h_tt1m1e_good_goodLep_phi->Fill(lep_good_phi,weight_good);
 
   if( shift=="SFs" ){
    h_tt1m1e_good_goodJ1_pt_CSVUp->Fill(goodJ1_pt,weight_good_CSVUp);
    h_tt1m1e_good_goodJ2_pt_CSVUp->Fill(goodJ2_pt,weight_good_CSVUp);
    h_tt1m1e_good_goodJ1_eta_CSVUp->Fill(goodJ1_eta,weight_good_CSVUp);
    h_tt1m1e_good_goodJ2_eta_CSVUp->Fill(goodJ2_eta,weight_good_CSVUp);
    h_tt1m1e_good_goodJ1_phi_CSVUp->Fill(goodJ1_phi,weight_good_CSVUp);
    h_tt1m1e_good_goodJ2_phi_CSVUp->Fill(goodJ2_phi,weight_good_CSVUp);
    h_tt1m1e_good_goodJ1_CSV_CSVUp->Fill(goodJ1_CSV,weight_good_CSVUp);
    h_tt1m1e_good_goodJ2_CSV_CSVUp->Fill(goodJ2_CSV,weight_good_CSVUp);
    h_tt1m1e_good_dRgoodJ1J2_CSVUp->Fill(dRgoodJ1J2,weight_good_CSVUp);
    h_tt1m1e_good_mt_CSVUp->Fill(mt_good,weight_good_CSVUp);
    h_tt1m1e_good_goodLep_pt_CSVUp->Fill(lep_good_pt,weight_good_CSVUp);
    h_tt1m1e_good_goodLep_eta_CSVUp->Fill(lep_good_eta,weight_good_CSVUp);
    h_tt1m1e_good_goodLep_phi_CSVUp->Fill(lep_good_phi,weight_good_CSVUp);
 
    h_tt1m1e_good_goodJ1_pt_CSVDown->Fill(goodJ1_pt,weight_good_CSVDn);
    h_tt1m1e_good_goodJ2_pt_CSVDown->Fill(goodJ2_pt,weight_good_CSVDn);
    h_tt1m1e_good_goodJ1_eta_CSVDown->Fill(goodJ1_eta,weight_good_CSVDn);
    h_tt1m1e_good_goodJ2_eta_CSVDown->Fill(goodJ2_eta,weight_good_CSVDn);
    h_tt1m1e_good_goodJ1_phi_CSVDown->Fill(goodJ1_phi,weight_good_CSVDn);
    h_tt1m1e_good_goodJ2_phi_CSVDown->Fill(goodJ2_phi,weight_good_CSVDn);
    h_tt1m1e_good_goodJ1_CSV_CSVDown->Fill(goodJ1_CSV,weight_good_CSVDn);
    h_tt1m1e_good_goodJ2_CSV_CSVDown->Fill(goodJ2_CSV,weight_good_CSVDn);
    h_tt1m1e_good_dRgoodJ1J2_CSVDown->Fill(dRgoodJ1J2,weight_good_CSVDn);
    h_tt1m1e_good_mt_CSVDown->Fill(mt_good,weight_good_CSVDn);
    h_tt1m1e_good_goodLep_pt_CSVDown->Fill(lep_good_pt,weight_good_CSVDn);
    h_tt1m1e_good_goodLep_eta_CSVDown->Fill(lep_good_eta,weight_good_CSVDn);
    h_tt1m1e_good_goodLep_phi_CSVDown->Fill(lep_good_phi,weight_good_CSVDn);
 
    h_tt1m1e_good_goodJ1_pt_EMuUp->Fill(goodJ1_pt,weight_good_EMuUp);
    h_tt1m1e_good_goodJ2_pt_EMuUp->Fill(goodJ2_pt,weight_good_EMuUp);
    h_tt1m1e_good_goodJ1_eta_EMuUp->Fill(goodJ1_eta,weight_good_EMuUp);
    h_tt1m1e_good_goodJ2_eta_EMuUp->Fill(goodJ2_eta,weight_good_EMuUp);
    h_tt1m1e_good_goodJ1_phi_EMuUp->Fill(goodJ1_phi,weight_good_EMuUp);
    h_tt1m1e_good_goodJ2_phi_EMuUp->Fill(goodJ2_phi,weight_good_EMuUp);
    h_tt1m1e_good_goodJ1_CSV_EMuUp->Fill(goodJ1_CSV,weight_good_EMuUp);
    h_tt1m1e_good_goodJ2_CSV_EMuUp->Fill(goodJ2_CSV,weight_good_EMuUp);
    h_tt1m1e_good_dRgoodJ1J2_EMuUp->Fill(dRgoodJ1J2,weight_good_EMuUp);
    h_tt1m1e_good_mt_EMuUp->Fill(mt_good,weight_good_EMuUp);
    h_tt1m1e_good_goodLep_pt_EMuUp->Fill(lep_good_pt,weight_good_EMuUp);
    h_tt1m1e_good_goodLep_eta_EMuUp->Fill(lep_good_eta,weight_good_EMuUp);
    h_tt1m1e_good_goodLep_phi_EMuUp->Fill(lep_good_phi,weight_good_EMuUp);

    h_tt1m1e_good_goodJ1_pt_EMuDown->Fill(goodJ1_pt,weight_good_EMuDn);
    h_tt1m1e_good_goodJ2_pt_EMuDown->Fill(goodJ2_pt,weight_good_EMuDn);
    h_tt1m1e_good_goodJ1_eta_EMuDown->Fill(goodJ1_eta,weight_good_EMuDn);
    h_tt1m1e_good_goodJ2_eta_EMuDown->Fill(goodJ2_eta,weight_good_EMuDn);
    h_tt1m1e_good_goodJ1_phi_EMuDown->Fill(goodJ1_phi,weight_good_EMuDn);
    h_tt1m1e_good_goodJ2_phi_EMuDown->Fill(goodJ2_phi,weight_good_EMuDn);
    h_tt1m1e_good_goodJ1_CSV_EMuDown->Fill(goodJ1_CSV,weight_good_EMuDn);
    h_tt1m1e_good_goodJ2_CSV_EMuDown->Fill(goodJ2_CSV,weight_good_EMuDn);
    h_tt1m1e_good_dRgoodJ1J2_EMuDown->Fill(dRgoodJ1J2,weight_good_EMuDn);
    h_tt1m1e_good_mt_EMuDown->Fill(mt_good,weight_good_EMuDn);
    h_tt1m1e_good_goodLep_pt_EMuDown->Fill(lep_good_pt,weight_good_EMuDn);
    h_tt1m1e_good_goodLep_eta_EMuDown->Fill(lep_good_eta,weight_good_EMuDn);
    h_tt1m1e_good_goodLep_phi_EMuDown->Fill(lep_good_phi,weight_good_EMuDn);
   }
  } 
  if(pass_tt1m1e_qcd){
   h_tt1m1e_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_tt1m1e_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_tt1m1e_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_tt1m1e_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_tt1m1e_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_tt1m1e_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_tt1m1e_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_tt1m1e_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_tt1m1e_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_tt1m1e_qcd_mt->Fill(mt_qcd,weight_qcd);
   h_tt1m1e_qcd_qcdLep_pt->Fill(lep_qcd_pt,weight_qcd);
   h_tt1m1e_qcd_qcdLep_eta->Fill(lep_qcd_eta,weight_qcd);
   h_tt1m1e_qcd_qcdLep_phi->Fill(lep_qcd_phi,weight_qcd);
 
   if( shift=="SFs" ){
    h_tt1m1e_qcd_goodJ1_pt_CSVUp->Fill(goodJ1_pt,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ2_pt_CSVUp->Fill(goodJ2_pt,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ1_eta_CSVUp->Fill(goodJ1_eta,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ2_eta_CSVUp->Fill(goodJ2_eta,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ1_phi_CSVUp->Fill(goodJ1_phi,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ2_phi_CSVUp->Fill(goodJ2_phi,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ1_CSV_CSVUp->Fill(goodJ1_CSV,weight_qcd_CSVUp);
    h_tt1m1e_qcd_goodJ2_CSV_CSVUp->Fill(goodJ2_CSV,weight_qcd_CSVUp);
    h_tt1m1e_qcd_dRgoodJ1J2_CSVUp->Fill(dRgoodJ1J2,weight_qcd_CSVUp);
    h_tt1m1e_qcd_mt_CSVUp->Fill(mt_qcd,weight_qcd_CSVUp);
    h_tt1m1e_qcd_qcdLep_pt_CSVUp->Fill(lep_qcd_pt,weight_qcd_CSVUp);
    h_tt1m1e_qcd_qcdLep_eta_CSVUp->Fill(lep_qcd_eta,weight_qcd_CSVUp);
    h_tt1m1e_qcd_qcdLep_phi_CSVUp->Fill(lep_qcd_phi,weight_qcd_CSVUp);
 
    h_tt1m1e_qcd_goodJ1_pt_CSVDown->Fill(goodJ1_pt,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ2_pt_CSVDown->Fill(goodJ2_pt,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ1_eta_CSVDown->Fill(goodJ1_eta,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ2_eta_CSVDown->Fill(goodJ2_eta,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ1_phi_CSVDown->Fill(goodJ1_phi,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ2_phi_CSVDown->Fill(goodJ2_phi,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ1_CSV_CSVDown->Fill(goodJ1_CSV,weight_qcd_CSVDn);
    h_tt1m1e_qcd_goodJ2_CSV_CSVDown->Fill(goodJ2_CSV,weight_qcd_CSVDn);
    h_tt1m1e_qcd_dRgoodJ1J2_CSVDown->Fill(dRgoodJ1J2,weight_qcd_CSVDn);
    h_tt1m1e_qcd_mt_CSVDown->Fill(mt_qcd,weight_qcd_CSVDn);
    h_tt1m1e_qcd_qcdLep_pt_CSVDown->Fill(lep_qcd_pt,weight_qcd_CSVDn);
    h_tt1m1e_qcd_qcdLep_eta_CSVDown->Fill(lep_qcd_eta,weight_qcd_CSVDn);
    h_tt1m1e_qcd_qcdLep_phi_CSVDown->Fill(lep_qcd_phi,weight_qcd_CSVDn);
 
    h_tt1m1e_qcd_goodJ1_pt_EMuUp->Fill(goodJ1_pt,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ2_pt_EMuUp->Fill(goodJ2_pt,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ1_eta_EMuUp->Fill(goodJ1_eta,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ2_eta_EMuUp->Fill(goodJ2_eta,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ1_phi_EMuUp->Fill(goodJ1_phi,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ2_phi_EMuUp->Fill(goodJ2_phi,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ1_CSV_EMuUp->Fill(goodJ1_CSV,weight_qcd_EMuUp);
    h_tt1m1e_qcd_goodJ2_CSV_EMuUp->Fill(goodJ2_CSV,weight_qcd_EMuUp);
    h_tt1m1e_qcd_dRgoodJ1J2_EMuUp->Fill(dRgoodJ1J2,weight_qcd_EMuUp);
    h_tt1m1e_qcd_mt_EMuUp->Fill(mt_qcd,weight_qcd_EMuUp);
    h_tt1m1e_qcd_qcdLep_pt_EMuUp->Fill(lep_qcd_pt,weight_qcd_EMuUp);
    h_tt1m1e_qcd_qcdLep_eta_EMuUp->Fill(lep_qcd_eta,weight_qcd_EMuUp);
    h_tt1m1e_qcd_qcdLep_phi_EMuUp->Fill(lep_qcd_phi,weight_qcd_EMuUp);

    h_tt1m1e_qcd_goodJ1_pt_EMuDown->Fill(goodJ1_pt,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ2_pt_EMuDown->Fill(goodJ2_pt,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ1_eta_EMuDown->Fill(goodJ1_eta,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ2_eta_EMuDown->Fill(goodJ2_eta,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ1_phi_EMuDown->Fill(goodJ1_phi,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ2_phi_EMuDown->Fill(goodJ2_phi,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ1_CSV_EMuDown->Fill(goodJ1_CSV,weight_qcd_EMuDn);
    h_tt1m1e_qcd_goodJ2_CSV_EMuDown->Fill(goodJ2_CSV,weight_qcd_EMuDn);
    h_tt1m1e_qcd_dRgoodJ1J2_EMuDown->Fill(dRgoodJ1J2,weight_qcd_EMuDn);
    h_tt1m1e_qcd_mt_EMuDown->Fill(mt_qcd,weight_qcd_EMuDn);
    h_tt1m1e_qcd_qcdLep_pt_EMuDown->Fill(lep_qcd_pt,weight_qcd_EMuDn);
    h_tt1m1e_qcd_qcdLep_eta_EMuDown->Fill(lep_qcd_eta,weight_qcd_EMuDn);
    h_tt1m1e_qcd_qcdLep_phi_EMuDown->Fill(lep_qcd_phi,weight_qcd_EMuDn);
   }
  } 

  // wjj Selection
  Bool_t pass_wjj_good=kFALSE;
  Bool_t pass_wjj_qcd=kFALSE;

  if( isMu ){
   if( oneGoodMuon && exactly2goodJs ){ 
    pass_wjj_good=kTRUE; 
    nrEntries_wjj_good_postcut++;
    mt_good = mt_goodMuon;
    lep_good_pt = goodMu_pt_vec->at(0);
    lep_good_eta = goodMu_eta_vec->at(0);
    lep_good_phi = goodMu_phi_vec->at(0);
    if( isMC ){
     SF_wjj_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0);
    }
   }
   if( oneQCDMuon && exactly2goodJs ){ 
    pass_wjj_qcd=kTRUE; 
    nrEntries_wjj_qcd_postcut++;
    mt_qcd = mt_qcdMuon;
    lep_qcd_pt = qcdMu_pt;
    lep_qcd_eta = qcdMu_eta;
    lep_qcd_phi = qcdMu_phi;
    if( isMC ){
     SF_wjj_qcd_IDIsoHLT = SF_qcdMu_IDIso;
    }
   }
  }
  else {
   if( oneGoodElectron && exactly2goodJs ){ 
    pass_wjj_good=kTRUE; 
    nrEntries_wjj_good_postcut++;
    mt_good = mt_goodElectron;
    lep_good_pt = goodEle_pt_vec->at(0);
    lep_good_eta = goodEle_eta_vec->at(0);
    lep_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_wjj_good_IDIsoHLT = SF_goodEle_IDIsoHLT->at(0);
    }
   }
   if( oneQCDElectron && exactly2goodJs ){ 
    pass_wjj_qcd=kTRUE; 
    nrEntries_wjj_qcd_postcut++;
    mt_qcd = mt_qcdElectron;
    lep_qcd_pt = qcdEle_pt;
    lep_qcd_eta = qcdEle_eta;
    lep_qcd_phi = qcdEle_phi;
    if( isMC ){
     SF_wjj_qcd_IDIsoHLT = SF_qcdEle_IDIso;
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good = SF_wjj_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec / nrEvents;
   weight_qcd =  SF_wjj_qcd_IDIsoHLT  * SF_lumiWeightPU * luminosity * crossSec / nrEvents;
  }
  else{
   weight_good = 1.; 
   weight_qcd =  1.; 
  }

  if(pass_wjj_good){
   h_wjj_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_wjj_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_wjj_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_wjj_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_wjj_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_wjj_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_wjj_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_wjj_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_wjj_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_wjj_good_mt->Fill(mt_good,weight_good);
   h_wjj_good_goodLep_pt->Fill(lep_good_pt,weight_good);
   h_wjj_good_goodLep_eta->Fill(lep_good_eta,weight_good);
   h_wjj_good_goodLep_phi->Fill(lep_good_phi,weight_good);
  } 
  if(pass_wjj_qcd){
   h_wjj_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_wjj_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_wjj_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_wjj_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_wjj_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_wjj_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_wjj_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_wjj_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_wjj_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_wjj_qcd_mt->Fill(mt_qcd,weight_qcd);
   h_wjj_qcd_qcdLep_pt->Fill(lep_qcd_pt,weight_qcd);
   h_wjj_qcd_qcdLep_eta->Fill(lep_qcd_eta,weight_qcd);
   h_wjj_qcd_qcdLep_phi->Fill(lep_qcd_phi,weight_qcd);
  } 

  // stt Selection
  Bool_t pass_stt_good=kFALSE;
  Bool_t pass_stt_qcd=kFALSE;

  if( isMu ){
   if( oneGoodMuon && aGoodBJaFwdJ ){ 
    pass_stt_good=kTRUE; 
    nrEntries_stt_good_postcut++;
    mt_good = mt_goodMuon;
    lep_good_pt = goodMu_pt_vec->at(0);
    lep_good_eta = goodMu_eta_vec->at(0);
    lep_good_phi = goodMu_phi_vec->at(0);
    if( isMC ){
     SF_stt_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0);
    }
   }
   if( oneQCDMuon && aGoodBJaFwdJ ){ 
    pass_stt_qcd=kTRUE; 
    nrEntries_stt_qcd_postcut++;
    mt_qcd = mt_qcdMuon;
    lep_qcd_pt = qcdMu_pt;
    lep_qcd_eta = qcdMu_eta;
    lep_qcd_phi = qcdMu_phi;
    if( isMC ){
     SF_stt_qcd_IDIsoHLT = SF_qcdMu_IDIso;
    }
   }
  }
  else {
   if( oneGoodElectron && aGoodBJaFwdJ ){ 
    pass_stt_good=kTRUE; 
    nrEntries_stt_good_postcut++;
    mt_good = mt_goodElectron;
    lep_good_pt = goodEle_pt_vec->at(0);
    lep_good_eta = goodEle_eta_vec->at(0);
    lep_good_phi = goodEle_phi_vec->at(0);
    if( isMC ){
     SF_stt_good_IDIsoHLT = SF_goodEle_IDIsoHLT->at(0);
    }
   }
   if( oneQCDElectron && aGoodBJaFwdJ ){ 
    pass_stt_qcd=kTRUE; 
    nrEntries_stt_qcd_postcut++;
    mt_qcd = mt_qcdElectron;
    lep_qcd_pt = qcdEle_pt;
    lep_qcd_eta = qcdEle_eta;
    lep_qcd_phi = qcdEle_phi;
    if( isMC ){
     SF_stt_qcd_IDIsoHLT = SF_qcdMu_IDIso;
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good = SF_stt_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_topBJ / nrEvents;
   weight_qcd =  SF_stt_qcd_IDIsoHLT  * SF_lumiWeightPU * luminosity * crossSec * SF_topBJ / nrEvents;
  }
  else{
   weight_good = 1.; 
   weight_qcd =  1.; 
  }

  if(pass_stt_good){
   h_stt_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_stt_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_stt_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_stt_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_stt_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_stt_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_stt_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_stt_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_stt_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_stt_good_mt->Fill(mt_good,weight_good);
   h_stt_good_goodLep_pt->Fill(lep_good_pt,weight_good);
   h_stt_good_goodLep_eta->Fill(lep_good_eta,weight_good);
   h_stt_good_goodLep_phi->Fill(lep_good_phi,weight_good);
  } 
  if(pass_stt_qcd){
   h_stt_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_stt_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_stt_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_stt_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_stt_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_stt_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_stt_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_stt_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_stt_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_stt_qcd_mt->Fill(mt_qcd,weight_qcd);
   h_stt_qcd_qcdLep_pt->Fill(lep_qcd_pt,weight_qcd);
   h_stt_qcd_qcdLep_eta->Fill(lep_qcd_eta,weight_qcd);
   h_stt_qcd_qcdLep_phi->Fill(lep_qcd_phi,weight_qcd);
  } 

  // dyjj Selection
  Bool_t pass_dyjj_good=kFALSE;
  Bool_t pass_dyjj_qcd=kFALSE;
 
  if( isMu ){
   if( twoGoodMuons && exactly2goodJs ){ 
    pass_dyjj_good=kTRUE; 
    nrEntries_dyjj_good_postcut++;
    lep1_good_pt = goodMu_pt_vec->at(0);
    lep1_good_eta = goodMu_eta_vec->at(0);
    lep1_good_phi = goodMu_phi_vec->at(0);
    lep2_good_pt = goodMu_pt_vec->at(1);
    lep2_good_eta = goodMu_eta_vec->at(1);
    lep2_good_phi = goodMu_phi_vec->at(1);
    if( isMC ){
     SF_dyjj_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0) * SF_goodMu_IDIsoHLT->at(1);
    }
   }
   if( twoQCDMuons && exactly2goodJs ){ 
    pass_dyjj_qcd=kTRUE; 
    nrEntries_dyjj_qcd_postcut++;
    lep1_qcd_pt = qcdMu_pt;
    lep1_qcd_eta = qcdMu_eta;
    lep1_qcd_phi = qcdMu_phi;
    lep2_qcd_pt = goodMu_pt_vec->at(1);
    lep2_qcd_eta = goodMu_eta_vec->at(1);
    lep2_qcd_phi = goodMu_phi_vec->at(1);
    if( isMC ){
     SF_dyjj_qcd_IDIsoHLT = SF_qcdMu_IDIso * SF_goodMu_IDIsoHLT->at(0); // * SF_goodMu_IDIso->at(0);
    }
   }
  }
  else {
   if( twoGoodElectrons && exactly2goodJs ){ 
    pass_dyjj_good=kTRUE; 
    nrEntries_dyjj_good_postcut++;
    lep1_good_pt = goodEle_pt_vec->at(0);
    lep1_good_eta = goodEle_eta_vec->at(0);
    lep1_good_phi = goodEle_phi_vec->at(0);
    lep2_good_pt = goodEle_pt_vec->at(1);
    lep2_good_eta = goodEle_eta_vec->at(1);
    lep2_good_phi = goodEle_phi_vec->at(1);
    if( isMC ){
     SF_dyjj_good_IDIsoHLT = SF_goodEle_IDIsoHLT->at(0) * SF_goodEle_IDIsoHLT->at(1); // * SF_goodEle_IDIso->at(1);
    }
   }
   if( twoQCDElectrons && exactly2goodJs ){ 
    pass_dyjj_qcd=kTRUE; 
    nrEntries_dyjj_qcd_postcut++;
    lep1_qcd_pt = qcdEle_pt;
    lep1_qcd_eta = qcdEle_eta;
    lep1_qcd_phi = qcdEle_phi;
    lep2_qcd_pt = goodEle_pt_vec->at(1);
    lep2_qcd_eta = goodEle_eta_vec->at(1);
    lep2_qcd_phi = goodEle_phi_vec->at(1);
    if( isMC ){
     SF_dyjj_qcd_IDIsoHLT = SF_qcdEle_IDIso * SF_goodEle_IDIsoHLT->at(0); // * SF_goodEle_IDIso->at(0);
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good = SF_dyjj_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec / nrEvents;
   weight_qcd =  SF_dyjj_qcd_IDIsoHLT  * SF_lumiWeightPU * luminosity * crossSec / nrEvents;
  }
  else{
   weight_good = 1.; 
   weight_qcd =  1.; 
  }

  if(pass_dyjj_good){
   h_dyjj_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_dyjj_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_dyjj_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_dyjj_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_dyjj_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_dyjj_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_dyjj_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_dyjj_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_dyjj_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_dyjj_good_goodLep1_pt->Fill(lep1_good_pt,weight_good);
   h_dyjj_good_goodLep1_eta->Fill(lep1_good_eta,weight_good);
   h_dyjj_good_goodLep1_phi->Fill(lep1_good_phi,weight_good);
   h_dyjj_good_goodLep2_pt->Fill(lep2_good_pt,weight_good);
   h_dyjj_good_goodLep2_eta->Fill(lep2_good_eta,weight_good);
   h_dyjj_good_goodLep2_phi->Fill(lep2_good_phi,weight_good);
  } 
  if(pass_dyjj_qcd){
   h_dyjj_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_dyjj_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_dyjj_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_dyjj_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_dyjj_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_dyjj_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_dyjj_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_dyjj_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_dyjj_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_dyjj_qcd_qcdLep1_pt->Fill(lep1_qcd_pt,weight_qcd);
   h_dyjj_qcd_qcdLep1_eta->Fill(lep1_qcd_eta,weight_qcd);
   h_dyjj_qcd_qcdLep1_phi->Fill(lep1_qcd_phi,weight_qcd);
   h_dyjj_qcd_qcdLep2_pt->Fill(lep2_qcd_pt,weight_qcd);
   h_dyjj_qcd_qcdLep2_eta->Fill(lep2_qcd_eta,weight_qcd);
   h_dyjj_qcd_qcdLep2_phi->Fill(lep2_qcd_phi,weight_qcd);
  } 

  // dybb Selection
  Bool_t pass_dybb_good=kFALSE;
  Bool_t pass_dybb_qcd=kFALSE;
 
  if( isMu ){
   if( twoGoodMuons && exactly2goodBJs ){ 
    pass_dybb_good=kTRUE; 
    nrEntries_dybb_good_postcut++;
    lep1_good_pt = goodMu_pt_vec->at(0);
    lep1_good_eta = goodMu_eta_vec->at(0);
    lep1_good_phi = goodMu_phi_vec->at(0);
    lep2_good_pt = goodMu_pt_vec->at(1);
    lep2_good_eta = goodMu_eta_vec->at(1);
    lep2_good_phi = goodMu_phi_vec->at(1);
    if( isMC ){
     SF_dybb_good_IDIsoHLT = SF_goodMu_IDIsoHLT->at(0) * SF_goodMu_IDIsoHLT->at(1);
    }
   }
   if( twoQCDMuons && exactly2goodBJs ){ 
    pass_dybb_qcd=kTRUE; 
    nrEntries_dybb_qcd_postcut++;
    lep1_qcd_pt = qcdMu_pt;
    lep1_qcd_eta = qcdMu_eta;
    lep1_qcd_phi = qcdMu_phi;
    lep2_qcd_pt = goodMu_pt_vec->at(1);
    lep2_qcd_eta = goodMu_eta_vec->at(1);
    lep2_qcd_phi = goodMu_phi_vec->at(1);
    if( isMC ){
     SF_dybb_qcd_IDIsoHLT = SF_qcdMu_IDIso * SF_goodMu_IDIsoHLT->at(0); // * SF_goodMu_IDIso->at(0);
    }
   }
  }
  else {
   if( twoGoodElectrons && exactly2goodBJs ){ 
    pass_dybb_good=kTRUE; 
    nrEntries_dybb_good_postcut++;
    lep1_good_pt = goodEle_pt_vec->at(0);
    lep1_good_eta = goodEle_eta_vec->at(0);
    lep1_good_phi = goodEle_phi_vec->at(0);
    lep2_good_pt = goodEle_pt_vec->at(1);
    lep2_good_eta = goodEle_eta_vec->at(1);
    lep2_good_phi = goodEle_phi_vec->at(1);
    if( isMC ){
     SF_dybb_good_IDIsoHLT = SF_goodEle_IDIsoHLT->at(0) * SF_goodEle_IDIsoHLT->at(1); // * SF_goodEle_IDIso->at(1);
    }
   }
   if( twoQCDElectrons && exactly2goodBJs ){ 
    pass_dybb_qcd=kTRUE; 
    nrEntries_dybb_qcd_postcut++;
    lep1_qcd_pt = qcdEle_pt;
    lep1_qcd_eta = qcdEle_eta;
    lep1_qcd_phi = qcdEle_phi;
    lep2_qcd_pt = goodEle_pt_vec->at(1);
    lep2_qcd_eta = goodEle_eta_vec->at(1);
    lep2_qcd_phi = goodEle_phi_vec->at(1);
    if( isMC ){
     SF_dybb_qcd_IDIsoHLT = SF_qcdEle_IDIso * SF_goodEle_IDIsoHLT->at(0); // * SF_goodEle_IDIso->at(0);
    }
   }
  }
  // total weight
  if( isMC ){
   weight_good = SF_dybb_good_IDIsoHLT * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
   weight_qcd =  SF_dybb_qcd_IDIsoHLT  * SF_lumiWeightPU * luminosity * crossSec * SF_top2BJs / nrEvents;
  }
  else{
   weight_good = 1.; 
   weight_qcd =  1.; 
  }

  if(pass_dybb_good){
   h_dybb_good_goodJ1_pt->Fill(goodJ1_pt,weight_good);
   h_dybb_good_goodJ2_pt->Fill(goodJ2_pt,weight_good);
   h_dybb_good_goodJ1_eta->Fill(goodJ1_eta,weight_good);
   h_dybb_good_goodJ2_eta->Fill(goodJ2_eta,weight_good);
   h_dybb_good_goodJ1_phi->Fill(goodJ1_phi,weight_good);
   h_dybb_good_goodJ2_phi->Fill(goodJ2_phi,weight_good);
   h_dybb_good_goodJ1_CSV->Fill(goodJ1_CSV,weight_good);
   h_dybb_good_goodJ2_CSV->Fill(goodJ2_CSV,weight_good);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_dybb_good_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_good);
   h_dybb_good_goodLep1_pt->Fill(lep1_good_pt,weight_good);
   h_dybb_good_goodLep1_eta->Fill(lep1_good_eta,weight_good);
   h_dybb_good_goodLep1_phi->Fill(lep1_good_phi,weight_good);
   h_dybb_good_goodLep2_pt->Fill(lep2_good_pt,weight_good);
   h_dybb_good_goodLep2_eta->Fill(lep2_good_eta,weight_good);
   h_dybb_good_goodLep2_phi->Fill(lep2_good_phi,weight_good);
  } 
  if(pass_dybb_qcd){
   h_dybb_qcd_goodJ1_pt->Fill(goodJ1_pt,weight_qcd);
   h_dybb_qcd_goodJ2_pt->Fill(goodJ2_pt,weight_qcd);
   h_dybb_qcd_goodJ1_eta->Fill(goodJ1_eta,weight_qcd);
   h_dybb_qcd_goodJ2_eta->Fill(goodJ2_eta,weight_qcd);
   h_dybb_qcd_goodJ1_phi->Fill(goodJ1_phi,weight_qcd);
   h_dybb_qcd_goodJ2_phi->Fill(goodJ2_phi,weight_qcd);
   h_dybb_qcd_goodJ1_CSV->Fill(goodJ1_CSV,weight_qcd);
   h_dybb_qcd_goodJ2_CSV->Fill(goodJ2_CSV,weight_qcd);
   dRgoodJ1J2 = -1.;
    dphiJJ = abs(goodJ1_phi - goodJ2_phi);
    if( dphiJJ>3.1416 ){ dphiJJ-=6.2832; }
    dRgoodJ1J2 = sqrt( ((goodJ1_eta-goodJ2_eta)*(goodJ1_eta-goodJ2_eta)) + (dphiJJ*dphiJJ) );
   h_dybb_qcd_dRgoodJ1J2->Fill(dRgoodJ1J2,weight_qcd);
   h_dybb_qcd_qcdLep1_pt->Fill(lep1_qcd_pt,weight_qcd);
   h_dybb_qcd_qcdLep1_eta->Fill(lep1_qcd_eta,weight_qcd);
   h_dybb_qcd_qcdLep1_phi->Fill(lep1_qcd_phi,weight_qcd);
   h_dybb_qcd_qcdLep2_pt->Fill(lep2_qcd_pt,weight_qcd);
   h_dybb_qcd_qcdLep2_eta->Fill(lep2_qcd_eta,weight_qcd);
   h_dybb_qcd_qcdLep2_phi->Fill(lep2_qcd_phi,weight_qcd);
  } 
 } 

 Float_t bmin = h_wbb_good_mt->GetXaxis()->FindBin(0.);
 Float_t bmax = h_wbb_good_mt->GetXaxis()->FindBin(20.);
 // write a log
 std::ofstream logfile (outfileName+".log", std::ofstream::out);
 logfile<<"isMu = "<<isMu<<std::endl;
 logfile<<"isMC = "<<isMC<<std::endl;
 logfile<<"shift = "<<shift<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--                 W+bb Selection                  --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_wbb_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_wbb_good_goodLep_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_wbb_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_wbb_qcd_qcdLep_eta->Integral()<<std::endl;
  Float_t size_wbb_mt_good020 = h_wbb_good_mt->Integral(bmin,bmax);
  Float_t size_wbb_mt_qcd020 = h_wbb_qcd_mt->Integral(bmin,bmax);
 logfile<<"  Nr. Good Events in 0<mT<20: "<<size_wbb_mt_good020<<std::endl;
 logfile<<"  Nr. QCD Events in 0<mT<20: "<<size_wbb_mt_qcd020<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--                TT(1m1e) Selection               --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_tt1m1e_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_tt1m1e_good_goodLep_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_tt1m1e_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_tt1m1e_qcd_qcdLep_eta->Integral()<<std::endl;
  Float_t size_tt1m1e_mt_good020 = h_tt1m1e_good_mt->Integral(bmin,bmax);
  Float_t size_tt1m1e_mt_qcd020 = h_tt1m1e_qcd_mt->Integral(bmin,bmax);
 logfile<<"  Nr. Good Events in 0<mT<20: "<<size_tt1m1e_mt_good020<<std::endl;
 logfile<<"  Nr. QCD Events in 0<mT<20: "<<size_tt1m1e_mt_qcd020<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--                 TT(3j) Selection                --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_tt3j_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_tt3j_good_goodLep_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_tt3j_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_tt3j_qcd_qcdLep_eta->Integral()<<std::endl;
  Float_t size_tt3j_mt_good020 = h_tt3j_good_mt->Integral(bmin,bmax);
  Float_t size_tt3j_mt_qcd020 = h_tt3j_qcd_mt->Integral(bmin,bmax);
 logfile<<"  Nr. Good Events in 0<mT<20: "<<size_tt3j_mt_good020<<std::endl;
 logfile<<"  Nr. QCD Events in 0<mT<20: "<<size_tt3j_mt_qcd020<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--                  W+jj Selection                 --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_wjj_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_wjj_good_goodLep_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_wjj_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_wjj_qcd_qcdLep_eta->Integral()<<std::endl;
  Float_t size_wjj_mt_good020 = h_wjj_good_mt->Integral(bmin,bmax);
  Float_t size_wjj_mt_qcd020 = h_wjj_qcd_mt->Integral(bmin,bmax);
 logfile<<"  Nr. Good Events in 0<mT<20: "<<size_wjj_mt_good020<<std::endl;
 logfile<<"  Nr. QCD Events in 0<mT<20: "<<size_wjj_mt_qcd020<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--             Single Top (t) Selection            --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_stt_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_stt_good_goodLep_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_stt_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_stt_qcd_qcdLep_eta->Integral()<<std::endl;
  Float_t size_stt_mt_good020 = h_stt_good_mt->Integral(bmin,bmax);
  Float_t size_stt_mt_qcd020 = h_stt_qcd_mt->Integral(bmin,bmax);
 logfile<<"  Nr. Good Events in 0<mT<20: "<<size_stt_mt_good020<<std::endl;
 logfile<<"  Nr. QCD Events in 0<mT<20: "<<size_stt_mt_qcd020<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--              Drell-Yan+jj Selection             --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_dyjj_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_dyjj_good_goodLep1_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_dyjj_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_dyjj_qcd_qcdLep1_eta->Integral()<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"--              Drell-Yan+bb Selection             --"<<std::endl;
 logfile<<"-----------------------------------------------------"<<std::endl;
 logfile<<"  Nr. Initial Entries: "<<nrEntries<<std::endl;
 logfile<<"  Nr. Entries Passing Good Cut: "<<nrEntries_dybb_good_postcut<<std::endl;
 logfile<<"  Nr. Good Events After Weight: "<<h_dybb_good_goodLep1_eta->Integral()<<std::endl;
 logfile<<"  Nr. Entries Passing QCD Cut: "<<nrEntries_dybb_qcd_postcut<<std::endl;
 logfile<<"  Nr. QCD Events After Weight: "<<h_dybb_qcd_qcdLep1_eta->Integral()<<std::endl;
 logfile.close();
 
 // write these histograms to file
 TFile *outfile = new TFile(outfileName+".root","RECREATE");
 outfile->cd();
                                                                
 // wbb 
   h_wbb_good_goodJ1_pt->Write();            h_wbb_good_goodJ1_pt->Delete();
   h_wbb_good_goodJ2_pt->Write();            h_wbb_good_goodJ2_pt->Delete();
   h_wbb_good_goodJ1_eta->Write();           h_wbb_good_goodJ1_eta->Delete();
   h_wbb_good_goodJ2_eta->Write();           h_wbb_good_goodJ2_eta->Delete();
   h_wbb_good_goodJ1_phi->Write();           h_wbb_good_goodJ1_phi->Delete();
   h_wbb_good_goodJ2_phi->Write();           h_wbb_good_goodJ2_phi->Delete();
   h_wbb_good_goodJ1_CSV->Write();           h_wbb_good_goodJ1_CSV->Delete();
   h_wbb_good_goodJ2_CSV->Write();           h_wbb_good_goodJ2_CSV->Delete();
   h_wbb_good_dRgoodJ1J2->Write();           h_wbb_good_dRgoodJ1J2->Delete();
   h_wbb_good_mt->Write();                   h_wbb_good_mt->Delete();
   h_wbb_good_goodLep_pt->Write();           h_wbb_good_goodLep_pt->Delete();
   h_wbb_good_goodLep_eta->Write();          h_wbb_good_goodLep_eta->Delete();
   h_wbb_good_goodLep_phi->Write();          h_wbb_good_goodLep_phi->Delete();
                                            
   if( shift=="SFs" ){                      
    h_wbb_good_goodJ1_pt_CSVUp->Write();      
    h_wbb_good_goodJ2_pt_CSVUp->Write();      
    h_wbb_good_goodJ1_eta_CSVUp->Write();     
    h_wbb_good_goodJ2_eta_CSVUp->Write();     
    h_wbb_good_goodJ1_phi_CSVUp->Write();     
    h_wbb_good_goodJ2_phi_CSVUp->Write();     
    h_wbb_good_goodJ1_CSV_CSVUp->Write();     
    h_wbb_good_goodJ2_CSV_CSVUp->Write();     
    h_wbb_good_dRgoodJ1J2_CSVUp->Write();     
    h_wbb_good_mt_CSVUp->Write();             
    h_wbb_good_goodLep_pt_CSVUp->Write();     
    h_wbb_good_goodLep_eta_CSVUp->Write();    
    h_wbb_good_goodLep_phi_CSVUp->Write();    
                                            
    h_wbb_good_goodJ1_pt_CSVDown->Write();    
    h_wbb_good_goodJ2_pt_CSVDown->Write();    
    h_wbb_good_goodJ1_eta_CSVDown->Write();   
    h_wbb_good_goodJ2_eta_CSVDown->Write();   
    h_wbb_good_goodJ1_phi_CSVDown->Write();   
    h_wbb_good_goodJ2_phi_CSVDown->Write();   
    h_wbb_good_goodJ1_CSV_CSVDown->Write();   
    h_wbb_good_goodJ2_CSV_CSVDown->Write();   
    h_wbb_good_dRgoodJ1J2_CSVDown->Write();   
    h_wbb_good_mt_CSVDown->Write();           
    h_wbb_good_goodLep_pt_CSVDown->Write();   
    h_wbb_good_goodLep_eta_CSVDown->Write();  
    h_wbb_good_goodLep_phi_CSVDown->Write();  
                                            
    h_wbb_good_goodJ1_pt_EMuUp->Write();      
    h_wbb_good_goodJ2_pt_EMuUp->Write();      
    h_wbb_good_goodJ1_eta_EMuUp->Write();     
    h_wbb_good_goodJ2_eta_EMuUp->Write();     
    h_wbb_good_goodJ1_phi_EMuUp->Write();     
    h_wbb_good_goodJ2_phi_EMuUp->Write();     
    h_wbb_good_goodJ1_CSV_EMuUp->Write();     
    h_wbb_good_goodJ2_CSV_EMuUp->Write();     
    h_wbb_good_dRgoodJ1J2_EMuUp->Write();     
    h_wbb_good_mt_EMuUp->Write();             
    h_wbb_good_goodLep_pt_EMuUp->Write();     
    h_wbb_good_goodLep_eta_EMuUp->Write();    
    h_wbb_good_goodLep_phi_EMuUp->Write();    
                                            
    h_wbb_good_goodJ1_pt_EMuDown->Write();    
    h_wbb_good_goodJ2_pt_EMuDown->Write();    
    h_wbb_good_goodJ1_eta_EMuDown->Write();   
    h_wbb_good_goodJ2_eta_EMuDown->Write();   
    h_wbb_good_goodJ1_phi_EMuDown->Write();   
    h_wbb_good_goodJ2_phi_EMuDown->Write();   
    h_wbb_good_goodJ1_CSV_EMuDown->Write();   
    h_wbb_good_goodJ2_CSV_EMuDown->Write();   
    h_wbb_good_dRgoodJ1J2_EMuDown->Write();   
    h_wbb_good_mt_EMuDown->Write();           
    h_wbb_good_goodLep_pt_EMuDown->Write();   
    h_wbb_good_goodLep_eta_EMuDown->Write();  
    h_wbb_good_goodLep_phi_EMuDown->Write();  
   }                                         

   h_wbb_good_goodJ1_pt_CSVUp->Delete();
   h_wbb_good_goodJ2_pt_CSVUp->Delete();
   h_wbb_good_goodJ1_eta_CSVUp->Delete();
   h_wbb_good_goodJ2_eta_CSVUp->Delete();
   h_wbb_good_goodJ1_phi_CSVUp->Delete();
   h_wbb_good_goodJ2_phi_CSVUp->Delete();
   h_wbb_good_goodJ1_CSV_CSVUp->Delete();
   h_wbb_good_goodJ2_CSV_CSVUp->Delete();
   h_wbb_good_dRgoodJ1J2_CSVUp->Delete();
   h_wbb_good_mt_CSVUp->Delete();
   h_wbb_good_goodLep_pt_CSVUp->Delete();
   h_wbb_good_goodLep_eta_CSVUp->Delete();
   h_wbb_good_goodLep_phi_CSVUp->Delete();
                                             
   h_wbb_good_goodJ1_pt_CSVDown->Delete();
   h_wbb_good_goodJ2_pt_CSVDown->Delete();
   h_wbb_good_goodJ1_eta_CSVDown->Delete();
   h_wbb_good_goodJ2_eta_CSVDown->Delete();
   h_wbb_good_goodJ1_phi_CSVDown->Delete();
   h_wbb_good_goodJ2_phi_CSVDown->Delete();
   h_wbb_good_goodJ1_CSV_CSVDown->Delete();
   h_wbb_good_goodJ2_CSV_CSVDown->Delete();
   h_wbb_good_dRgoodJ1J2_CSVDown->Delete();
   h_wbb_good_mt_CSVDown->Delete();
   h_wbb_good_goodLep_pt_CSVDown->Delete();
   h_wbb_good_goodLep_eta_CSVDown->Delete();
   h_wbb_good_goodLep_phi_CSVDown->Delete();
                                             
   h_wbb_good_goodJ1_pt_EMuUp->Delete();
   h_wbb_good_goodJ2_pt_EMuUp->Delete();
   h_wbb_good_goodJ1_eta_EMuUp->Delete();
   h_wbb_good_goodJ2_eta_EMuUp->Delete();
   h_wbb_good_goodJ1_phi_EMuUp->Delete();
   h_wbb_good_goodJ2_phi_EMuUp->Delete();
   h_wbb_good_goodJ1_CSV_EMuUp->Delete();
   h_wbb_good_goodJ2_CSV_EMuUp->Delete();
   h_wbb_good_dRgoodJ1J2_EMuUp->Delete();
   h_wbb_good_mt_EMuUp->Delete();
   h_wbb_good_goodLep_pt_EMuUp->Delete();
   h_wbb_good_goodLep_eta_EMuUp->Delete();
   h_wbb_good_goodLep_phi_EMuUp->Delete();
                                             
   h_wbb_good_goodJ1_pt_EMuDown->Delete();
   h_wbb_good_goodJ2_pt_EMuDown->Delete();
   h_wbb_good_goodJ1_eta_EMuDown->Delete();
   h_wbb_good_goodJ2_eta_EMuDown->Delete();
   h_wbb_good_goodJ1_phi_EMuDown->Delete();
   h_wbb_good_goodJ2_phi_EMuDown->Delete();
   h_wbb_good_goodJ1_CSV_EMuDown->Delete();
   h_wbb_good_goodJ2_CSV_EMuDown->Delete();
   h_wbb_good_dRgoodJ1J2_EMuDown->Delete();
   h_wbb_good_mt_EMuDown->Delete();
   h_wbb_good_goodLep_pt_EMuDown->Delete();
   h_wbb_good_goodLep_eta_EMuDown->Delete();
   h_wbb_good_goodLep_phi_EMuDown->Delete();
                                            
   h_wbb_qcd_goodJ1_pt->Write();             h_wbb_qcd_goodJ1_pt->Delete();
   h_wbb_qcd_goodJ2_pt->Write();             h_wbb_qcd_goodJ2_pt->Delete();
   h_wbb_qcd_goodJ1_eta->Write();            h_wbb_qcd_goodJ1_eta->Delete();
   h_wbb_qcd_goodJ2_eta->Write();            h_wbb_qcd_goodJ2_eta->Delete();
   h_wbb_qcd_goodJ1_phi->Write();            h_wbb_qcd_goodJ1_phi->Delete();
   h_wbb_qcd_goodJ2_phi->Write();            h_wbb_qcd_goodJ2_phi->Delete();
   h_wbb_qcd_goodJ1_CSV->Write();            h_wbb_qcd_goodJ1_CSV->Delete();
   h_wbb_qcd_goodJ2_CSV->Write();            h_wbb_qcd_goodJ2_CSV->Delete();
   h_wbb_qcd_dRgoodJ1J2->Write();            h_wbb_qcd_dRgoodJ1J2->Delete();
   h_wbb_qcd_mt->Write();                    h_wbb_qcd_mt->Delete();
   h_wbb_qcd_qcdLep_pt->Write();             h_wbb_qcd_qcdLep_pt->Delete();
   h_wbb_qcd_qcdLep_eta->Write();            h_wbb_qcd_qcdLep_eta->Delete();
   h_wbb_qcd_qcdLep_phi->Write();            h_wbb_qcd_qcdLep_phi->Delete();
                                            
   if( shift=="SFs" ){                      
    h_wbb_qcd_goodJ1_pt_CSVUp->Write();       
    h_wbb_qcd_goodJ2_pt_CSVUp->Write();       
    h_wbb_qcd_goodJ1_eta_CSVUp->Write();      
    h_wbb_qcd_goodJ2_eta_CSVUp->Write();      
    h_wbb_qcd_goodJ1_phi_CSVUp->Write();      
    h_wbb_qcd_goodJ2_phi_CSVUp->Write();      
    h_wbb_qcd_goodJ1_CSV_CSVUp->Write();      
    h_wbb_qcd_goodJ2_CSV_CSVUp->Write();      
    h_wbb_qcd_dRgoodJ1J2_CSVUp->Write();      
    h_wbb_qcd_mt_CSVUp->Write();              
    h_wbb_qcd_qcdLep_pt_CSVUp->Write();       
    h_wbb_qcd_qcdLep_eta_CSVUp->Write();      
    h_wbb_qcd_qcdLep_phi_CSVUp->Write();      
                                            
    h_wbb_qcd_goodJ1_pt_CSVDown->Write();     
    h_wbb_qcd_goodJ2_pt_CSVDown->Write();     
    h_wbb_qcd_goodJ1_eta_CSVDown->Write();    
    h_wbb_qcd_goodJ2_eta_CSVDown->Write();    
    h_wbb_qcd_goodJ1_phi_CSVDown->Write();    
    h_wbb_qcd_goodJ2_phi_CSVDown->Write();    
    h_wbb_qcd_goodJ1_CSV_CSVDown->Write();    
    h_wbb_qcd_goodJ2_CSV_CSVDown->Write();    
    h_wbb_qcd_dRgoodJ1J2_CSVDown->Write();    
    h_wbb_qcd_mt_CSVDown->Write();            
    h_wbb_qcd_qcdLep_pt_CSVDown->Write();     
    h_wbb_qcd_qcdLep_eta_CSVDown->Write();    
    h_wbb_qcd_qcdLep_phi_CSVDown->Write();    
                                            
    h_wbb_qcd_goodJ1_pt_EMuUp->Write();       
    h_wbb_qcd_goodJ2_pt_EMuUp->Write();       
    h_wbb_qcd_goodJ1_eta_EMuUp->Write();      
    h_wbb_qcd_goodJ2_eta_EMuUp->Write();      
    h_wbb_qcd_goodJ1_phi_EMuUp->Write();      
    h_wbb_qcd_goodJ2_phi_EMuUp->Write();      
    h_wbb_qcd_goodJ1_CSV_EMuUp->Write();      
    h_wbb_qcd_goodJ2_CSV_EMuUp->Write();      
    h_wbb_qcd_dRgoodJ1J2_EMuUp->Write();      
    h_wbb_qcd_mt_EMuUp->Write();              
    h_wbb_qcd_qcdLep_pt_EMuUp->Write();       
    h_wbb_qcd_qcdLep_eta_EMuUp->Write();      
    h_wbb_qcd_qcdLep_phi_EMuUp->Write();      
                                            
    h_wbb_qcd_goodJ1_pt_EMuDown->Write();     
    h_wbb_qcd_goodJ2_pt_EMuDown->Write();     
    h_wbb_qcd_goodJ1_eta_EMuDown->Write();    
    h_wbb_qcd_goodJ2_eta_EMuDown->Write();    
    h_wbb_qcd_goodJ1_phi_EMuDown->Write();    
    h_wbb_qcd_goodJ2_phi_EMuDown->Write();    
    h_wbb_qcd_goodJ1_CSV_EMuDown->Write();    
    h_wbb_qcd_goodJ2_CSV_EMuDown->Write();    
    h_wbb_qcd_dRgoodJ1J2_EMuDown->Write();    
    h_wbb_qcd_mt_EMuDown->Write();            
    h_wbb_qcd_qcdLep_pt_EMuDown->Write();     
    h_wbb_qcd_qcdLep_eta_EMuDown->Write();    
    h_wbb_qcd_qcdLep_phi_EMuDown->Write();    
   }                                                             

   h_wbb_qcd_goodJ1_pt_CSVUp->Delete();
   h_wbb_qcd_goodJ2_pt_CSVUp->Delete();
   h_wbb_qcd_goodJ1_eta_CSVUp->Delete();
   h_wbb_qcd_goodJ2_eta_CSVUp->Delete();
   h_wbb_qcd_goodJ1_phi_CSVUp->Delete();
   h_wbb_qcd_goodJ2_phi_CSVUp->Delete();
   h_wbb_qcd_goodJ1_CSV_CSVUp->Delete();
   h_wbb_qcd_goodJ2_CSV_CSVUp->Delete();
   h_wbb_qcd_dRgoodJ1J2_CSVUp->Delete();
   h_wbb_qcd_mt_CSVUp->Delete();
   h_wbb_qcd_qcdLep_pt_CSVUp->Delete();
   h_wbb_qcd_qcdLep_eta_CSVUp->Delete();
   h_wbb_qcd_qcdLep_phi_CSVUp->Delete();
                                           
   h_wbb_qcd_goodJ1_pt_CSVDown->Delete();
   h_wbb_qcd_goodJ2_pt_CSVDown->Delete();
   h_wbb_qcd_goodJ1_eta_CSVDown->Delete();
   h_wbb_qcd_goodJ2_eta_CSVDown->Delete();
   h_wbb_qcd_goodJ1_phi_CSVDown->Delete();
   h_wbb_qcd_goodJ2_phi_CSVDown->Delete();
   h_wbb_qcd_goodJ1_CSV_CSVDown->Delete();
   h_wbb_qcd_goodJ2_CSV_CSVDown->Delete();
   h_wbb_qcd_dRgoodJ1J2_CSVDown->Delete();
   h_wbb_qcd_mt_CSVDown->Delete();
   h_wbb_qcd_qcdLep_pt_CSVDown->Delete();
   h_wbb_qcd_qcdLep_eta_CSVDown->Delete();
   h_wbb_qcd_qcdLep_phi_CSVDown->Delete();
                                           
   h_wbb_qcd_goodJ1_pt_EMuUp->Delete();
   h_wbb_qcd_goodJ2_pt_EMuUp->Delete();
   h_wbb_qcd_goodJ1_eta_EMuUp->Delete();
   h_wbb_qcd_goodJ2_eta_EMuUp->Delete();
   h_wbb_qcd_goodJ1_phi_EMuUp->Delete();
   h_wbb_qcd_goodJ2_phi_EMuUp->Delete();
   h_wbb_qcd_goodJ1_CSV_EMuUp->Delete();
   h_wbb_qcd_goodJ2_CSV_EMuUp->Delete();
   h_wbb_qcd_dRgoodJ1J2_EMuUp->Delete();
   h_wbb_qcd_mt_EMuUp->Delete();
   h_wbb_qcd_qcdLep_pt_EMuUp->Delete();
   h_wbb_qcd_qcdLep_eta_EMuUp->Delete();
   h_wbb_qcd_qcdLep_phi_EMuUp->Delete();
                                           
   h_wbb_qcd_goodJ1_pt_EMuDown->Delete();
   h_wbb_qcd_goodJ2_pt_EMuDown->Delete();
   h_wbb_qcd_goodJ1_eta_EMuDown->Delete();
   h_wbb_qcd_goodJ2_eta_EMuDown->Delete();
   h_wbb_qcd_goodJ1_phi_EMuDown->Delete();
   h_wbb_qcd_goodJ2_phi_EMuDown->Delete();
   h_wbb_qcd_goodJ1_CSV_EMuDown->Delete();
   h_wbb_qcd_goodJ2_CSV_EMuDown->Delete();
   h_wbb_qcd_dRgoodJ1J2_EMuDown->Delete();
   h_wbb_qcd_mt_EMuDown->Delete();
   h_wbb_qcd_qcdLep_pt_EMuDown->Delete();
   h_wbb_qcd_qcdLep_eta_EMuDown->Delete();
   h_wbb_qcd_qcdLep_phi_EMuDown->Delete();


 // tt3j 
   h_tt3j_good_goodJ1_pt->Write();            h_tt3j_good_goodJ1_pt->Delete();
   h_tt3j_good_goodJ2_pt->Write();            h_tt3j_good_goodJ2_pt->Delete();
   h_tt3j_good_goodJ1_eta->Write();           h_tt3j_good_goodJ1_eta->Delete();
   h_tt3j_good_goodJ2_eta->Write();           h_tt3j_good_goodJ2_eta->Delete();
   h_tt3j_good_goodJ1_phi->Write();           h_tt3j_good_goodJ1_phi->Delete();
   h_tt3j_good_goodJ2_phi->Write();           h_tt3j_good_goodJ2_phi->Delete();
   h_tt3j_good_goodJ1_CSV->Write();           h_tt3j_good_goodJ1_CSV->Delete();
   h_tt3j_good_goodJ2_CSV->Write();           h_tt3j_good_goodJ2_CSV->Delete();
   h_tt3j_good_dRgoodJ1J2->Write();           h_tt3j_good_dRgoodJ1J2->Delete();
   h_tt3j_good_mt->Write();                   h_tt3j_good_mt->Delete();
   h_tt3j_good_goodLep_pt->Write();           h_tt3j_good_goodLep_pt->Delete();
   h_tt3j_good_goodLep_eta->Write();          h_tt3j_good_goodLep_eta->Delete();
   h_tt3j_good_goodLep_phi->Write();          h_tt3j_good_goodLep_phi->Delete();
                                            
   if( shift=="SFs" ){                      
    h_tt3j_good_goodJ1_pt_CSVUp->Write();      
    h_tt3j_good_goodJ2_pt_CSVUp->Write();      
    h_tt3j_good_goodJ1_eta_CSVUp->Write();     
    h_tt3j_good_goodJ2_eta_CSVUp->Write();     
    h_tt3j_good_goodJ1_phi_CSVUp->Write();     
    h_tt3j_good_goodJ2_phi_CSVUp->Write();     
    h_tt3j_good_goodJ1_CSV_CSVUp->Write();     
    h_tt3j_good_goodJ2_CSV_CSVUp->Write();     
    h_tt3j_good_dRgoodJ1J2_CSVUp->Write();     
    h_tt3j_good_mt_CSVUp->Write();             
    h_tt3j_good_goodLep_pt_CSVUp->Write();     
    h_tt3j_good_goodLep_eta_CSVUp->Write();    
    h_tt3j_good_goodLep_phi_CSVUp->Write();    
                                            
    h_tt3j_good_goodJ1_pt_CSVDown->Write();    
    h_tt3j_good_goodJ2_pt_CSVDown->Write();    
    h_tt3j_good_goodJ1_eta_CSVDown->Write();   
    h_tt3j_good_goodJ2_eta_CSVDown->Write();   
    h_tt3j_good_goodJ1_phi_CSVDown->Write();   
    h_tt3j_good_goodJ2_phi_CSVDown->Write();   
    h_tt3j_good_goodJ1_CSV_CSVDown->Write();   
    h_tt3j_good_goodJ2_CSV_CSVDown->Write();   
    h_tt3j_good_dRgoodJ1J2_CSVDown->Write();   
    h_tt3j_good_mt_CSVDown->Write();           
    h_tt3j_good_goodLep_pt_CSVDown->Write();   
    h_tt3j_good_goodLep_eta_CSVDown->Write();  
    h_tt3j_good_goodLep_phi_CSVDown->Write();  
                                            
    h_tt3j_good_goodJ1_pt_EMuUp->Write();      
    h_tt3j_good_goodJ2_pt_EMuUp->Write();      
    h_tt3j_good_goodJ1_eta_EMuUp->Write();     
    h_tt3j_good_goodJ2_eta_EMuUp->Write();     
    h_tt3j_good_goodJ1_phi_EMuUp->Write();     
    h_tt3j_good_goodJ2_phi_EMuUp->Write();     
    h_tt3j_good_goodJ1_CSV_EMuUp->Write();     
    h_tt3j_good_goodJ2_CSV_EMuUp->Write();     
    h_tt3j_good_dRgoodJ1J2_EMuUp->Write();     
    h_tt3j_good_mt_EMuUp->Write();             
    h_tt3j_good_goodLep_pt_EMuUp->Write();     
    h_tt3j_good_goodLep_eta_EMuUp->Write();    
    h_tt3j_good_goodLep_phi_EMuUp->Write();    
                                            
    h_tt3j_good_goodJ1_pt_EMuDown->Write();    
    h_tt3j_good_goodJ2_pt_EMuDown->Write();    
    h_tt3j_good_goodJ1_eta_EMuDown->Write();   
    h_tt3j_good_goodJ2_eta_EMuDown->Write();   
    h_tt3j_good_goodJ1_phi_EMuDown->Write();   
    h_tt3j_good_goodJ2_phi_EMuDown->Write();   
    h_tt3j_good_goodJ1_CSV_EMuDown->Write();   
    h_tt3j_good_goodJ2_CSV_EMuDown->Write();   
    h_tt3j_good_dRgoodJ1J2_EMuDown->Write();   
    h_tt3j_good_mt_EMuDown->Write();           
    h_tt3j_good_goodLep_pt_EMuDown->Write();   
    h_tt3j_good_goodLep_eta_EMuDown->Write();  
    h_tt3j_good_goodLep_phi_EMuDown->Write();  
   }                                         

   h_tt3j_good_goodJ1_pt_CSVUp->Delete();
   h_tt3j_good_goodJ2_pt_CSVUp->Delete();
   h_tt3j_good_goodJ1_eta_CSVUp->Delete();
   h_tt3j_good_goodJ2_eta_CSVUp->Delete();
   h_tt3j_good_goodJ1_phi_CSVUp->Delete();
   h_tt3j_good_goodJ2_phi_CSVUp->Delete();
   h_tt3j_good_goodJ1_CSV_CSVUp->Delete();
   h_tt3j_good_goodJ2_CSV_CSVUp->Delete();
   h_tt3j_good_dRgoodJ1J2_CSVUp->Delete();
   h_tt3j_good_mt_CSVUp->Delete();
   h_tt3j_good_goodLep_pt_CSVUp->Delete();
   h_tt3j_good_goodLep_eta_CSVUp->Delete();
   h_tt3j_good_goodLep_phi_CSVUp->Delete();
                                              
   h_tt3j_good_goodJ1_pt_CSVDown->Delete();
   h_tt3j_good_goodJ2_pt_CSVDown->Delete();
   h_tt3j_good_goodJ1_eta_CSVDown->Delete();
   h_tt3j_good_goodJ2_eta_CSVDown->Delete();
   h_tt3j_good_goodJ1_phi_CSVDown->Delete();
   h_tt3j_good_goodJ2_phi_CSVDown->Delete();
   h_tt3j_good_goodJ1_CSV_CSVDown->Delete();
   h_tt3j_good_goodJ2_CSV_CSVDown->Delete();
   h_tt3j_good_dRgoodJ1J2_CSVDown->Delete();
   h_tt3j_good_mt_CSVDown->Delete();
   h_tt3j_good_goodLep_pt_CSVDown->Delete();
   h_tt3j_good_goodLep_eta_CSVDown->Delete();
   h_tt3j_good_goodLep_phi_CSVDown->Delete();
                                              
   h_tt3j_good_goodJ1_pt_EMuUp->Delete();
   h_tt3j_good_goodJ2_pt_EMuUp->Delete();
   h_tt3j_good_goodJ1_eta_EMuUp->Delete();
   h_tt3j_good_goodJ2_eta_EMuUp->Delete();
   h_tt3j_good_goodJ1_phi_EMuUp->Delete();
   h_tt3j_good_goodJ2_phi_EMuUp->Delete();
   h_tt3j_good_goodJ1_CSV_EMuUp->Delete();
   h_tt3j_good_goodJ2_CSV_EMuUp->Delete();
   h_tt3j_good_dRgoodJ1J2_EMuUp->Delete();
   h_tt3j_good_mt_EMuUp->Delete();
   h_tt3j_good_goodLep_pt_EMuUp->Delete();
   h_tt3j_good_goodLep_eta_EMuUp->Delete();
   h_tt3j_good_goodLep_phi_EMuUp->Delete();
                                              
   h_tt3j_good_goodJ1_pt_EMuDown->Delete();
   h_tt3j_good_goodJ2_pt_EMuDown->Delete();
   h_tt3j_good_goodJ1_eta_EMuDown->Delete();
   h_tt3j_good_goodJ2_eta_EMuDown->Delete();
   h_tt3j_good_goodJ1_phi_EMuDown->Delete();
   h_tt3j_good_goodJ2_phi_EMuDown->Delete();
   h_tt3j_good_goodJ1_CSV_EMuDown->Delete();
   h_tt3j_good_goodJ2_CSV_EMuDown->Delete();
   h_tt3j_good_dRgoodJ1J2_EMuDown->Delete();
   h_tt3j_good_mt_EMuDown->Delete();
   h_tt3j_good_goodLep_pt_EMuDown->Delete();
   h_tt3j_good_goodLep_eta_EMuDown->Delete();
   h_tt3j_good_goodLep_phi_EMuDown->Delete();

   h_tt3j_qcd_goodJ1_pt->Write();             h_tt3j_qcd_goodJ1_pt->Delete();
   h_tt3j_qcd_goodJ2_pt->Write();             h_tt3j_qcd_goodJ2_pt->Delete();
   h_tt3j_qcd_goodJ1_eta->Write();            h_tt3j_qcd_goodJ1_eta->Delete();
   h_tt3j_qcd_goodJ2_eta->Write();            h_tt3j_qcd_goodJ2_eta->Delete();
   h_tt3j_qcd_goodJ1_phi->Write();            h_tt3j_qcd_goodJ1_phi->Delete();
   h_tt3j_qcd_goodJ2_phi->Write();            h_tt3j_qcd_goodJ2_phi->Delete();
   h_tt3j_qcd_goodJ1_CSV->Write();            h_tt3j_qcd_goodJ1_CSV->Delete();
   h_tt3j_qcd_goodJ2_CSV->Write();            h_tt3j_qcd_goodJ2_CSV->Delete();
   h_tt3j_qcd_dRgoodJ1J2->Write();            h_tt3j_qcd_dRgoodJ1J2->Delete();
   h_tt3j_qcd_mt->Write();                    h_tt3j_qcd_mt->Delete();
   h_tt3j_qcd_qcdLep_pt->Write();             h_tt3j_qcd_qcdLep_pt->Delete();
   h_tt3j_qcd_qcdLep_eta->Write();            h_tt3j_qcd_qcdLep_eta->Delete();
   h_tt3j_qcd_qcdLep_phi->Write();            h_tt3j_qcd_qcdLep_phi->Delete();
                                            
   if( shift=="SFs" ){                      
    h_tt3j_qcd_goodJ1_pt_CSVUp->Write();       
    h_tt3j_qcd_goodJ2_pt_CSVUp->Write();       
    h_tt3j_qcd_goodJ1_eta_CSVUp->Write();      
    h_tt3j_qcd_goodJ2_eta_CSVUp->Write();      
    h_tt3j_qcd_goodJ1_phi_CSVUp->Write();      
    h_tt3j_qcd_goodJ2_phi_CSVUp->Write();      
    h_tt3j_qcd_goodJ1_CSV_CSVUp->Write();      
    h_tt3j_qcd_goodJ2_CSV_CSVUp->Write();      
    h_tt3j_qcd_dRgoodJ1J2_CSVUp->Write();      
    h_tt3j_qcd_mt_CSVUp->Write();              
    h_tt3j_qcd_qcdLep_pt_CSVUp->Write();       
    h_tt3j_qcd_qcdLep_eta_CSVUp->Write();      
    h_tt3j_qcd_qcdLep_phi_CSVUp->Write();      
                                            
    h_tt3j_qcd_goodJ1_pt_CSVDown->Write();     
    h_tt3j_qcd_goodJ2_pt_CSVDown->Write();     
    h_tt3j_qcd_goodJ1_eta_CSVDown->Write();    
    h_tt3j_qcd_goodJ2_eta_CSVDown->Write();    
    h_tt3j_qcd_goodJ1_phi_CSVDown->Write();    
    h_tt3j_qcd_goodJ2_phi_CSVDown->Write();    
    h_tt3j_qcd_goodJ1_CSV_CSVDown->Write();    
    h_tt3j_qcd_goodJ2_CSV_CSVDown->Write();    
    h_tt3j_qcd_dRgoodJ1J2_CSVDown->Write();    
    h_tt3j_qcd_mt_CSVDown->Write();            
    h_tt3j_qcd_qcdLep_pt_CSVDown->Write();     
    h_tt3j_qcd_qcdLep_eta_CSVDown->Write();    
    h_tt3j_qcd_qcdLep_phi_CSVDown->Write();    
                                            
    h_tt3j_qcd_goodJ1_pt_EMuUp->Write();       
    h_tt3j_qcd_goodJ2_pt_EMuUp->Write();       
    h_tt3j_qcd_goodJ1_eta_EMuUp->Write();      
    h_tt3j_qcd_goodJ2_eta_EMuUp->Write();      
    h_tt3j_qcd_goodJ1_phi_EMuUp->Write();      
    h_tt3j_qcd_goodJ2_phi_EMuUp->Write();      
    h_tt3j_qcd_goodJ1_CSV_EMuUp->Write();      
    h_tt3j_qcd_goodJ2_CSV_EMuUp->Write();      
    h_tt3j_qcd_dRgoodJ1J2_EMuUp->Write();      
    h_tt3j_qcd_mt_EMuUp->Write();              
    h_tt3j_qcd_qcdLep_pt_EMuUp->Write();       
    h_tt3j_qcd_qcdLep_eta_EMuUp->Write();      
    h_tt3j_qcd_qcdLep_phi_EMuUp->Write();      
                                            
    h_tt3j_qcd_goodJ1_pt_EMuDown->Write();     
    h_tt3j_qcd_goodJ2_pt_EMuDown->Write();     
    h_tt3j_qcd_goodJ1_eta_EMuDown->Write();    
    h_tt3j_qcd_goodJ2_eta_EMuDown->Write();    
    h_tt3j_qcd_goodJ1_phi_EMuDown->Write();    
    h_tt3j_qcd_goodJ2_phi_EMuDown->Write();    
    h_tt3j_qcd_goodJ1_CSV_EMuDown->Write();    
    h_tt3j_qcd_goodJ2_CSV_EMuDown->Write();    
    h_tt3j_qcd_dRgoodJ1J2_EMuDown->Write();    
    h_tt3j_qcd_mt_EMuDown->Write();            
    h_tt3j_qcd_qcdLep_pt_EMuDown->Write();     
    h_tt3j_qcd_qcdLep_eta_EMuDown->Write();    
    h_tt3j_qcd_qcdLep_phi_EMuDown->Write();    
   }                                                             

   h_tt3j_qcd_goodJ1_pt_CSVUp->Delete();       
   h_tt3j_qcd_goodJ2_pt_CSVUp->Delete();       
   h_tt3j_qcd_goodJ1_eta_CSVUp->Delete();      
   h_tt3j_qcd_goodJ2_eta_CSVUp->Delete();      
   h_tt3j_qcd_goodJ1_phi_CSVUp->Delete();      
   h_tt3j_qcd_goodJ2_phi_CSVUp->Delete();      
   h_tt3j_qcd_goodJ1_CSV_CSVUp->Delete();      
   h_tt3j_qcd_goodJ2_CSV_CSVUp->Delete();      
   h_tt3j_qcd_dRgoodJ1J2_CSVUp->Delete();      
   h_tt3j_qcd_mt_CSVUp->Delete();              
   h_tt3j_qcd_qcdLep_pt_CSVUp->Delete();       
   h_tt3j_qcd_qcdLep_eta_CSVUp->Delete();      
   h_tt3j_qcd_qcdLep_phi_CSVUp->Delete();      
                                           
   h_tt3j_qcd_goodJ1_pt_CSVDown->Delete();     
   h_tt3j_qcd_goodJ2_pt_CSVDown->Delete();     
   h_tt3j_qcd_goodJ1_eta_CSVDown->Delete();    
   h_tt3j_qcd_goodJ2_eta_CSVDown->Delete();    
   h_tt3j_qcd_goodJ1_phi_CSVDown->Delete();    
   h_tt3j_qcd_goodJ2_phi_CSVDown->Delete();    
   h_tt3j_qcd_goodJ1_CSV_CSVDown->Delete();    
   h_tt3j_qcd_goodJ2_CSV_CSVDown->Delete();    
   h_tt3j_qcd_dRgoodJ1J2_CSVDown->Delete();    
   h_tt3j_qcd_mt_CSVDown->Delete();            
   h_tt3j_qcd_qcdLep_pt_CSVDown->Delete();     
   h_tt3j_qcd_qcdLep_eta_CSVDown->Delete();    
   h_tt3j_qcd_qcdLep_phi_CSVDown->Delete();    
                                           
   h_tt3j_qcd_goodJ1_pt_EMuUp->Delete();       
   h_tt3j_qcd_goodJ2_pt_EMuUp->Delete();       
   h_tt3j_qcd_goodJ1_eta_EMuUp->Delete();      
   h_tt3j_qcd_goodJ2_eta_EMuUp->Delete();      
   h_tt3j_qcd_goodJ1_phi_EMuUp->Delete();      
   h_tt3j_qcd_goodJ2_phi_EMuUp->Delete();      
   h_tt3j_qcd_goodJ1_CSV_EMuUp->Delete();      
   h_tt3j_qcd_goodJ2_CSV_EMuUp->Delete();      
   h_tt3j_qcd_dRgoodJ1J2_EMuUp->Delete();      
   h_tt3j_qcd_mt_EMuUp->Delete();              
   h_tt3j_qcd_qcdLep_pt_EMuUp->Delete();       
   h_tt3j_qcd_qcdLep_eta_EMuUp->Delete();      
   h_tt3j_qcd_qcdLep_phi_EMuUp->Delete();      
                                           
   h_tt3j_qcd_goodJ1_pt_EMuDown->Delete();     
   h_tt3j_qcd_goodJ2_pt_EMuDown->Delete();     
   h_tt3j_qcd_goodJ1_eta_EMuDown->Delete();    
   h_tt3j_qcd_goodJ2_eta_EMuDown->Delete();    
   h_tt3j_qcd_goodJ1_phi_EMuDown->Delete();    
   h_tt3j_qcd_goodJ2_phi_EMuDown->Delete();    
   h_tt3j_qcd_goodJ1_CSV_EMuDown->Delete();    
   h_tt3j_qcd_goodJ2_CSV_EMuDown->Delete();    
   h_tt3j_qcd_dRgoodJ1J2_EMuDown->Delete();    
   h_tt3j_qcd_mt_EMuDown->Delete();            
   h_tt3j_qcd_qcdLep_pt_EMuDown->Delete();     
   h_tt3j_qcd_qcdLep_eta_EMuDown->Delete();    
   h_tt3j_qcd_qcdLep_phi_EMuDown->Delete();    

 // tt1m1e 
   h_tt1m1e_good_goodJ1_pt->Write();            h_tt1m1e_good_goodJ1_pt->Delete();
   h_tt1m1e_good_goodJ2_pt->Write();            h_tt1m1e_good_goodJ2_pt->Delete();
   h_tt1m1e_good_goodJ1_eta->Write();           h_tt1m1e_good_goodJ1_eta->Delete();
   h_tt1m1e_good_goodJ2_eta->Write();           h_tt1m1e_good_goodJ2_eta->Delete();
   h_tt1m1e_good_goodJ1_phi->Write();           h_tt1m1e_good_goodJ1_phi->Delete();
   h_tt1m1e_good_goodJ2_phi->Write();           h_tt1m1e_good_goodJ2_phi->Delete();
   h_tt1m1e_good_goodJ1_CSV->Write();           h_tt1m1e_good_goodJ1_CSV->Delete();
   h_tt1m1e_good_goodJ2_CSV->Write();           h_tt1m1e_good_goodJ2_CSV->Delete();
   h_tt1m1e_good_dRgoodJ1J2->Write();           h_tt1m1e_good_dRgoodJ1J2->Delete();
   h_tt1m1e_good_mt->Write();                   h_tt1m1e_good_mt->Delete();
   h_tt1m1e_good_goodLep_pt->Write();           h_tt1m1e_good_goodLep_pt->Delete();
   h_tt1m1e_good_goodLep_eta->Write();          h_tt1m1e_good_goodLep_eta->Delete();
   h_tt1m1e_good_goodLep_phi->Write();          h_tt1m1e_good_goodLep_phi->Delete();
                                            
   if( shift=="SFs" ){                      
    h_tt1m1e_good_goodJ1_pt_CSVUp->Write();      
    h_tt1m1e_good_goodJ2_pt_CSVUp->Write();      
    h_tt1m1e_good_goodJ1_eta_CSVUp->Write();     
    h_tt1m1e_good_goodJ2_eta_CSVUp->Write();     
    h_tt1m1e_good_goodJ1_phi_CSVUp->Write();     
    h_tt1m1e_good_goodJ2_phi_CSVUp->Write();     
    h_tt1m1e_good_goodJ1_CSV_CSVUp->Write();     
    h_tt1m1e_good_goodJ2_CSV_CSVUp->Write();     
    h_tt1m1e_good_dRgoodJ1J2_CSVUp->Write();     
    h_tt1m1e_good_mt_CSVUp->Write();             
    h_tt1m1e_good_goodLep_pt_CSVUp->Write();     
    h_tt1m1e_good_goodLep_eta_CSVUp->Write();    
    h_tt1m1e_good_goodLep_phi_CSVUp->Write();    
                                            
    h_tt1m1e_good_goodJ1_pt_CSVDown->Write();    
    h_tt1m1e_good_goodJ2_pt_CSVDown->Write();    
    h_tt1m1e_good_goodJ1_eta_CSVDown->Write();   
    h_tt1m1e_good_goodJ2_eta_CSVDown->Write();   
    h_tt1m1e_good_goodJ1_phi_CSVDown->Write();   
    h_tt1m1e_good_goodJ2_phi_CSVDown->Write();   
    h_tt1m1e_good_goodJ1_CSV_CSVDown->Write();   
    h_tt1m1e_good_goodJ2_CSV_CSVDown->Write();   
    h_tt1m1e_good_dRgoodJ1J2_CSVDown->Write();   
    h_tt1m1e_good_mt_CSVDown->Write();           
    h_tt1m1e_good_goodLep_pt_CSVDown->Write();   
    h_tt1m1e_good_goodLep_eta_CSVDown->Write();  
    h_tt1m1e_good_goodLep_phi_CSVDown->Write();  
                                            
    h_tt1m1e_good_goodJ1_pt_EMuUp->Write();      
    h_tt1m1e_good_goodJ2_pt_EMuUp->Write();      
    h_tt1m1e_good_goodJ1_eta_EMuUp->Write();     
    h_tt1m1e_good_goodJ2_eta_EMuUp->Write();     
    h_tt1m1e_good_goodJ1_phi_EMuUp->Write();     
    h_tt1m1e_good_goodJ2_phi_EMuUp->Write();     
    h_tt1m1e_good_goodJ1_CSV_EMuUp->Write();     
    h_tt1m1e_good_goodJ2_CSV_EMuUp->Write();     
    h_tt1m1e_good_dRgoodJ1J2_EMuUp->Write();     
    h_tt1m1e_good_mt_EMuUp->Write();             
    h_tt1m1e_good_goodLep_pt_EMuUp->Write();     
    h_tt1m1e_good_goodLep_eta_EMuUp->Write();    
    h_tt1m1e_good_goodLep_phi_EMuUp->Write();    
                                            
    h_tt1m1e_good_goodJ1_pt_EMuDown->Write();    
    h_tt1m1e_good_goodJ2_pt_EMuDown->Write();    
    h_tt1m1e_good_goodJ1_eta_EMuDown->Write();   
    h_tt1m1e_good_goodJ2_eta_EMuDown->Write();   
    h_tt1m1e_good_goodJ1_phi_EMuDown->Write();   
    h_tt1m1e_good_goodJ2_phi_EMuDown->Write();   
    h_tt1m1e_good_goodJ1_CSV_EMuDown->Write();   
    h_tt1m1e_good_goodJ2_CSV_EMuDown->Write();   
    h_tt1m1e_good_dRgoodJ1J2_EMuDown->Write();   
    h_tt1m1e_good_mt_EMuDown->Write();           
    h_tt1m1e_good_goodLep_pt_EMuDown->Write();   
    h_tt1m1e_good_goodLep_eta_EMuDown->Write();  
    h_tt1m1e_good_goodLep_phi_EMuDown->Write();  
   }                                         

   h_tt1m1e_good_goodJ1_pt_CSVUp->Delete();
   h_tt1m1e_good_goodJ2_pt_CSVUp->Delete();
   h_tt1m1e_good_goodJ1_eta_CSVUp->Delete();
   h_tt1m1e_good_goodJ2_eta_CSVUp->Delete();
   h_tt1m1e_good_goodJ1_phi_CSVUp->Delete();
   h_tt1m1e_good_goodJ2_phi_CSVUp->Delete();
   h_tt1m1e_good_goodJ1_CSV_CSVUp->Delete();
   h_tt1m1e_good_goodJ2_CSV_CSVUp->Delete();
   h_tt1m1e_good_dRgoodJ1J2_CSVUp->Delete();
   h_tt1m1e_good_mt_CSVUp->Delete();
   h_tt1m1e_good_goodLep_pt_CSVUp->Delete();
   h_tt1m1e_good_goodLep_eta_CSVUp->Delete();
   h_tt1m1e_good_goodLep_phi_CSVUp->Delete();
                                                
   h_tt1m1e_good_goodJ1_pt_CSVDown->Delete();
   h_tt1m1e_good_goodJ2_pt_CSVDown->Delete();
   h_tt1m1e_good_goodJ1_eta_CSVDown->Delete();
   h_tt1m1e_good_goodJ2_eta_CSVDown->Delete();
   h_tt1m1e_good_goodJ1_phi_CSVDown->Delete();
   h_tt1m1e_good_goodJ2_phi_CSVDown->Delete();
   h_tt1m1e_good_goodJ1_CSV_CSVDown->Delete();
   h_tt1m1e_good_goodJ2_CSV_CSVDown->Delete();
   h_tt1m1e_good_dRgoodJ1J2_CSVDown->Delete();
   h_tt1m1e_good_mt_CSVDown->Delete();
   h_tt1m1e_good_goodLep_pt_CSVDown->Delete();
   h_tt1m1e_good_goodLep_eta_CSVDown->Delete();
   h_tt1m1e_good_goodLep_phi_CSVDown->Delete();
                                                
   h_tt1m1e_good_goodJ1_pt_EMuUp->Delete();
   h_tt1m1e_good_goodJ2_pt_EMuUp->Delete();
   h_tt1m1e_good_goodJ1_eta_EMuUp->Delete();
   h_tt1m1e_good_goodJ2_eta_EMuUp->Delete();
   h_tt1m1e_good_goodJ1_phi_EMuUp->Delete();
   h_tt1m1e_good_goodJ2_phi_EMuUp->Delete();
   h_tt1m1e_good_goodJ1_CSV_EMuUp->Delete();
   h_tt1m1e_good_goodJ2_CSV_EMuUp->Delete();
   h_tt1m1e_good_dRgoodJ1J2_EMuUp->Delete();
   h_tt1m1e_good_mt_EMuUp->Delete();
   h_tt1m1e_good_goodLep_pt_EMuUp->Delete();
   h_tt1m1e_good_goodLep_eta_EMuUp->Delete();
   h_tt1m1e_good_goodLep_phi_EMuUp->Delete();
                                                
   h_tt1m1e_good_goodJ1_pt_EMuDown->Delete();
   h_tt1m1e_good_goodJ2_pt_EMuDown->Delete();
   h_tt1m1e_good_goodJ1_eta_EMuDown->Delete();
   h_tt1m1e_good_goodJ2_eta_EMuDown->Delete();
   h_tt1m1e_good_goodJ1_phi_EMuDown->Delete();
   h_tt1m1e_good_goodJ2_phi_EMuDown->Delete();
   h_tt1m1e_good_goodJ1_CSV_EMuDown->Delete();
   h_tt1m1e_good_goodJ2_CSV_EMuDown->Delete();
   h_tt1m1e_good_dRgoodJ1J2_EMuDown->Delete();
   h_tt1m1e_good_mt_EMuDown->Delete();
   h_tt1m1e_good_goodLep_pt_EMuDown->Delete();
   h_tt1m1e_good_goodLep_eta_EMuDown->Delete();
   h_tt1m1e_good_goodLep_phi_EMuDown->Delete();
                                            
   h_tt1m1e_qcd_goodJ1_pt->Write();             h_tt1m1e_qcd_goodJ1_pt->Delete();
   h_tt1m1e_qcd_goodJ2_pt->Write();             h_tt1m1e_qcd_goodJ2_pt->Delete();
   h_tt1m1e_qcd_goodJ1_eta->Write();            h_tt1m1e_qcd_goodJ1_eta->Delete();
   h_tt1m1e_qcd_goodJ2_eta->Write();            h_tt1m1e_qcd_goodJ2_eta->Delete();
   h_tt1m1e_qcd_goodJ1_phi->Write();            h_tt1m1e_qcd_goodJ1_phi->Delete();
   h_tt1m1e_qcd_goodJ2_phi->Write();            h_tt1m1e_qcd_goodJ2_phi->Delete();
   h_tt1m1e_qcd_goodJ1_CSV->Write();            h_tt1m1e_qcd_goodJ1_CSV->Delete();
   h_tt1m1e_qcd_goodJ2_CSV->Write();            h_tt1m1e_qcd_goodJ2_CSV->Delete();
   h_tt1m1e_qcd_dRgoodJ1J2->Write();            h_tt1m1e_qcd_dRgoodJ1J2->Delete();
   h_tt1m1e_qcd_mt->Write();                    h_tt1m1e_qcd_mt->Delete();
   h_tt1m1e_qcd_qcdLep_pt->Write();             h_tt1m1e_qcd_qcdLep_pt->Delete();
   h_tt1m1e_qcd_qcdLep_eta->Write();            h_tt1m1e_qcd_qcdLep_eta->Delete();
   h_tt1m1e_qcd_qcdLep_phi->Write();            h_tt1m1e_qcd_qcdLep_phi->Delete();
                                            
   if( shift=="SFs" ){                      
    h_tt1m1e_qcd_goodJ1_pt_CSVUp->Write();       
    h_tt1m1e_qcd_goodJ2_pt_CSVUp->Write();       
    h_tt1m1e_qcd_goodJ1_eta_CSVUp->Write();      
    h_tt1m1e_qcd_goodJ2_eta_CSVUp->Write();      
    h_tt1m1e_qcd_goodJ1_phi_CSVUp->Write();      
    h_tt1m1e_qcd_goodJ2_phi_CSVUp->Write();      
    h_tt1m1e_qcd_goodJ1_CSV_CSVUp->Write();      
    h_tt1m1e_qcd_goodJ2_CSV_CSVUp->Write();      
    h_tt1m1e_qcd_dRgoodJ1J2_CSVUp->Write();      
    h_tt1m1e_qcd_mt_CSVUp->Write();              
    h_tt1m1e_qcd_qcdLep_pt_CSVUp->Write();       
    h_tt1m1e_qcd_qcdLep_eta_CSVUp->Write();      
    h_tt1m1e_qcd_qcdLep_phi_CSVUp->Write();      
                                            
    h_tt1m1e_qcd_goodJ1_pt_CSVDown->Write();     
    h_tt1m1e_qcd_goodJ2_pt_CSVDown->Write();     
    h_tt1m1e_qcd_goodJ1_eta_CSVDown->Write();    
    h_tt1m1e_qcd_goodJ2_eta_CSVDown->Write();    
    h_tt1m1e_qcd_goodJ1_phi_CSVDown->Write();    
    h_tt1m1e_qcd_goodJ2_phi_CSVDown->Write();    
    h_tt1m1e_qcd_goodJ1_CSV_CSVDown->Write();    
    h_tt1m1e_qcd_goodJ2_CSV_CSVDown->Write();    
    h_tt1m1e_qcd_dRgoodJ1J2_CSVDown->Write();    
    h_tt1m1e_qcd_mt_CSVDown->Write();            
    h_tt1m1e_qcd_qcdLep_pt_CSVDown->Write();     
    h_tt1m1e_qcd_qcdLep_eta_CSVDown->Write();    
    h_tt1m1e_qcd_qcdLep_phi_CSVDown->Write();    
                                            
    h_tt1m1e_qcd_goodJ1_pt_EMuUp->Write();       
    h_tt1m1e_qcd_goodJ2_pt_EMuUp->Write();       
    h_tt1m1e_qcd_goodJ1_eta_EMuUp->Write();      
    h_tt1m1e_qcd_goodJ2_eta_EMuUp->Write();      
    h_tt1m1e_qcd_goodJ1_phi_EMuUp->Write();      
    h_tt1m1e_qcd_goodJ2_phi_EMuUp->Write();      
    h_tt1m1e_qcd_goodJ1_CSV_EMuUp->Write();      
    h_tt1m1e_qcd_goodJ2_CSV_EMuUp->Write();      
    h_tt1m1e_qcd_dRgoodJ1J2_EMuUp->Write();      
    h_tt1m1e_qcd_mt_EMuUp->Write();              
    h_tt1m1e_qcd_qcdLep_pt_EMuUp->Write();       
    h_tt1m1e_qcd_qcdLep_eta_EMuUp->Write();      
    h_tt1m1e_qcd_qcdLep_phi_EMuUp->Write();      
                                            
    h_tt1m1e_qcd_goodJ1_pt_EMuDown->Write();     
    h_tt1m1e_qcd_goodJ2_pt_EMuDown->Write();     
    h_tt1m1e_qcd_goodJ1_eta_EMuDown->Write();    
    h_tt1m1e_qcd_goodJ2_eta_EMuDown->Write();    
    h_tt1m1e_qcd_goodJ1_phi_EMuDown->Write();    
    h_tt1m1e_qcd_goodJ2_phi_EMuDown->Write();    
    h_tt1m1e_qcd_goodJ1_CSV_EMuDown->Write();    
    h_tt1m1e_qcd_goodJ2_CSV_EMuDown->Write();    
    h_tt1m1e_qcd_dRgoodJ1J2_EMuDown->Write();    
    h_tt1m1e_qcd_mt_EMuDown->Write();            
    h_tt1m1e_qcd_qcdLep_pt_EMuDown->Write();     
    h_tt1m1e_qcd_qcdLep_eta_EMuDown->Write();    
    h_tt1m1e_qcd_qcdLep_phi_EMuDown->Write();    
   }                                                             

   h_tt1m1e_qcd_goodJ1_pt_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ2_pt_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ1_eta_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ2_eta_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ1_phi_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ2_phi_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ1_CSV_CSVUp->Delete();
   h_tt1m1e_qcd_goodJ2_CSV_CSVUp->Delete();
   h_tt1m1e_qcd_dRgoodJ1J2_CSVUp->Delete();
   h_tt1m1e_qcd_mt_CSVUp->Delete();
   h_tt1m1e_qcd_qcdLep_pt_CSVUp->Delete();
   h_tt1m1e_qcd_qcdLep_eta_CSVUp->Delete();
   h_tt1m1e_qcd_qcdLep_phi_CSVUp->Delete();
                                              
   h_tt1m1e_qcd_goodJ1_pt_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ2_pt_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ1_eta_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ2_eta_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ1_phi_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ2_phi_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ1_CSV_CSVDown->Delete();
   h_tt1m1e_qcd_goodJ2_CSV_CSVDown->Delete();
   h_tt1m1e_qcd_dRgoodJ1J2_CSVDown->Delete();
   h_tt1m1e_qcd_mt_CSVDown->Delete();
   h_tt1m1e_qcd_qcdLep_pt_CSVDown->Delete();
   h_tt1m1e_qcd_qcdLep_eta_CSVDown->Delete();
   h_tt1m1e_qcd_qcdLep_phi_CSVDown->Delete();
                                              
   h_tt1m1e_qcd_goodJ1_pt_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ2_pt_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ1_eta_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ2_eta_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ1_phi_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ2_phi_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ1_CSV_EMuUp->Delete();
   h_tt1m1e_qcd_goodJ2_CSV_EMuUp->Delete();
   h_tt1m1e_qcd_dRgoodJ1J2_EMuUp->Delete();
   h_tt1m1e_qcd_mt_EMuUp->Delete();
   h_tt1m1e_qcd_qcdLep_pt_EMuUp->Delete();
   h_tt1m1e_qcd_qcdLep_eta_EMuUp->Delete();
   h_tt1m1e_qcd_qcdLep_phi_EMuUp->Delete();
                                              
   h_tt1m1e_qcd_goodJ1_pt_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ2_pt_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ1_eta_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ2_eta_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ1_phi_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ2_phi_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ1_CSV_EMuDown->Delete();
   h_tt1m1e_qcd_goodJ2_CSV_EMuDown->Delete();
   h_tt1m1e_qcd_dRgoodJ1J2_EMuDown->Delete();
   h_tt1m1e_qcd_mt_EMuDown->Delete();
   h_tt1m1e_qcd_qcdLep_pt_EMuDown->Delete();
   h_tt1m1e_qcd_qcdLep_eta_EMuDown->Delete();
   h_tt1m1e_qcd_qcdLep_phi_EMuDown->Delete();

   // wjj
   h_wjj_good_goodJ1_pt->Write();            h_wjj_good_goodJ1_pt->Delete();
   h_wjj_good_goodJ2_pt->Write();            h_wjj_good_goodJ2_pt->Delete();
   h_wjj_good_goodJ1_eta->Write();           h_wjj_good_goodJ1_eta->Delete();
   h_wjj_good_goodJ2_eta->Write();           h_wjj_good_goodJ2_eta->Delete();
   h_wjj_good_goodJ1_phi->Write();           h_wjj_good_goodJ1_phi->Delete();
   h_wjj_good_goodJ2_phi->Write();           h_wjj_good_goodJ2_phi->Delete();
   h_wjj_good_goodJ1_CSV->Write();           h_wjj_good_goodJ1_CSV->Delete();
   h_wjj_good_goodJ2_CSV->Write();           h_wjj_good_goodJ2_CSV->Delete();
   h_wjj_good_dRgoodJ1J2->Write();           h_wjj_good_dRgoodJ1J2->Delete();
   h_wjj_good_mt->Write();                   h_wjj_good_mt->Delete();
   h_wjj_good_goodLep_pt->Write();           h_wjj_good_goodLep_pt->Delete();
   h_wjj_good_goodLep_eta->Write();          h_wjj_good_goodLep_eta->Delete();
   h_wjj_good_goodLep_phi->Write();          h_wjj_good_goodLep_phi->Delete();
                                            
   h_wjj_qcd_goodJ1_pt->Write();             h_wjj_qcd_goodJ1_pt->Delete();
   h_wjj_qcd_goodJ2_pt->Write();             h_wjj_qcd_goodJ2_pt->Delete();
   h_wjj_qcd_goodJ1_eta->Write();            h_wjj_qcd_goodJ1_eta->Delete();
   h_wjj_qcd_goodJ2_eta->Write();            h_wjj_qcd_goodJ2_eta->Delete();
   h_wjj_qcd_goodJ1_phi->Write();            h_wjj_qcd_goodJ1_phi->Delete();
   h_wjj_qcd_goodJ2_phi->Write();            h_wjj_qcd_goodJ2_phi->Delete();
   h_wjj_qcd_goodJ1_CSV->Write();            h_wjj_qcd_goodJ1_CSV->Delete();
   h_wjj_qcd_goodJ2_CSV->Write();            h_wjj_qcd_goodJ2_CSV->Delete();
   h_wjj_qcd_dRgoodJ1J2->Write();            h_wjj_qcd_dRgoodJ1J2->Delete();
   h_wjj_qcd_mt->Write();                    h_wjj_qcd_mt->Delete();
   h_wjj_qcd_qcdLep_pt->Write();             h_wjj_qcd_qcdLep_pt->Delete();
   h_wjj_qcd_qcdLep_eta->Write();            h_wjj_qcd_qcdLep_eta->Delete();
   h_wjj_qcd_qcdLep_phi->Write();            h_wjj_qcd_qcdLep_phi->Delete();

   // stt
   h_stt_good_goodJ1_pt->Write();            h_stt_good_goodJ1_pt->Delete();
   h_stt_good_goodJ2_pt->Write();            h_stt_good_goodJ2_pt->Delete();
   h_stt_good_goodJ1_eta->Write();           h_stt_good_goodJ1_eta->Delete();
   h_stt_good_goodJ2_eta->Write();           h_stt_good_goodJ2_eta->Delete();
   h_stt_good_goodJ1_phi->Write();           h_stt_good_goodJ1_phi->Delete();
   h_stt_good_goodJ2_phi->Write();           h_stt_good_goodJ2_phi->Delete();
   h_stt_good_goodJ1_CSV->Write();           h_stt_good_goodJ1_CSV->Delete();
   h_stt_good_goodJ2_CSV->Write();           h_stt_good_goodJ2_CSV->Delete();
   h_stt_good_dRgoodJ1J2->Write();           h_stt_good_dRgoodJ1J2->Delete();
   h_stt_good_mt->Write();                   h_stt_good_mt->Delete();
   h_stt_good_goodLep_pt->Write();           h_stt_good_goodLep_pt->Delete();
   h_stt_good_goodLep_eta->Write();          h_stt_good_goodLep_eta->Delete();
   h_stt_good_goodLep_phi->Write();          h_stt_good_goodLep_phi->Delete();
                                            
   h_stt_qcd_goodJ1_pt->Write();             h_stt_qcd_goodJ1_pt->Delete();
   h_stt_qcd_goodJ2_pt->Write();             h_stt_qcd_goodJ2_pt->Delete();
   h_stt_qcd_goodJ1_eta->Write();            h_stt_qcd_goodJ1_eta->Delete();
   h_stt_qcd_goodJ2_eta->Write();            h_stt_qcd_goodJ2_eta->Delete();
   h_stt_qcd_goodJ1_phi->Write();            h_stt_qcd_goodJ1_phi->Delete();
   h_stt_qcd_goodJ2_phi->Write();            h_stt_qcd_goodJ2_phi->Delete();
   h_stt_qcd_goodJ1_CSV->Write();            h_stt_qcd_goodJ1_CSV->Delete();
   h_stt_qcd_goodJ2_CSV->Write();            h_stt_qcd_goodJ2_CSV->Delete();
   h_stt_qcd_dRgoodJ1J2->Write();            h_stt_qcd_dRgoodJ1J2->Delete();
   h_stt_qcd_mt->Write();                    h_stt_qcd_mt->Delete();
   h_stt_qcd_qcdLep_pt->Write();             h_stt_qcd_qcdLep_pt->Delete();
   h_stt_qcd_qcdLep_eta->Write();            h_stt_qcd_qcdLep_eta->Delete();
   h_stt_qcd_qcdLep_phi->Write();            h_stt_qcd_qcdLep_phi->Delete();

   // dyjj
   h_dyjj_good_goodJ1_pt->Write();            h_dyjj_good_goodJ1_pt->Delete();
   h_dyjj_good_goodJ2_pt->Write();            h_dyjj_good_goodJ2_pt->Delete();
   h_dyjj_good_goodJ1_eta->Write();           h_dyjj_good_goodJ1_eta->Delete();
   h_dyjj_good_goodJ2_eta->Write();           h_dyjj_good_goodJ2_eta->Delete();
   h_dyjj_good_goodJ1_phi->Write();           h_dyjj_good_goodJ1_phi->Delete();
   h_dyjj_good_goodJ2_phi->Write();           h_dyjj_good_goodJ2_phi->Delete();
   h_dyjj_good_goodJ1_CSV->Write();           h_dyjj_good_goodJ1_CSV->Delete();
   h_dyjj_good_goodJ2_CSV->Write();           h_dyjj_good_goodJ2_CSV->Delete();
   h_dyjj_good_dRgoodJ1J2->Write();           h_dyjj_good_dRgoodJ1J2->Delete();
   h_dyjj_good_goodLep1_pt->Write();          h_dyjj_good_goodLep1_pt->Delete();
   h_dyjj_good_goodLep1_eta->Write();         h_dyjj_good_goodLep1_eta->Delete();
   h_dyjj_good_goodLep1_phi->Write();         h_dyjj_good_goodLep1_phi->Delete();
   h_dyjj_good_goodLep2_pt->Write();          h_dyjj_good_goodLep2_pt->Delete();
   h_dyjj_good_goodLep2_eta->Write();         h_dyjj_good_goodLep2_eta->Delete();
   h_dyjj_good_goodLep2_phi->Write();         h_dyjj_good_goodLep2_phi->Delete();
                                            
   h_dyjj_qcd_goodJ1_pt->Write();             h_dyjj_qcd_goodJ1_pt->Delete();
   h_dyjj_qcd_goodJ2_pt->Write();             h_dyjj_qcd_goodJ2_pt->Delete();
   h_dyjj_qcd_goodJ1_eta->Write();            h_dyjj_qcd_goodJ1_eta->Delete();
   h_dyjj_qcd_goodJ2_eta->Write();            h_dyjj_qcd_goodJ2_eta->Delete();
   h_dyjj_qcd_goodJ1_phi->Write();            h_dyjj_qcd_goodJ1_phi->Delete();
   h_dyjj_qcd_goodJ2_phi->Write();            h_dyjj_qcd_goodJ2_phi->Delete();
   h_dyjj_qcd_goodJ1_CSV->Write();            h_dyjj_qcd_goodJ1_CSV->Delete();
   h_dyjj_qcd_goodJ2_CSV->Write();            h_dyjj_qcd_goodJ2_CSV->Delete();
   h_dyjj_qcd_dRgoodJ1J2->Write();            h_dyjj_qcd_dRgoodJ1J2->Delete();
   h_dyjj_qcd_qcdLep1_pt->Write();            h_dyjj_qcd_qcdLep1_pt->Delete();
   h_dyjj_qcd_qcdLep1_eta->Write();           h_dyjj_qcd_qcdLep1_eta->Delete();
   h_dyjj_qcd_qcdLep1_phi->Write();           h_dyjj_qcd_qcdLep1_phi->Delete();
   h_dyjj_qcd_qcdLep2_pt->Write();            h_dyjj_qcd_qcdLep2_pt->Delete();
   h_dyjj_qcd_qcdLep2_eta->Write();           h_dyjj_qcd_qcdLep2_eta->Delete();
   h_dyjj_qcd_qcdLep2_phi->Write();           h_dyjj_qcd_qcdLep2_phi->Delete();
                                            
   h_dybb_good_goodJ1_pt->Write();             h_dybb_good_goodJ1_pt->Delete();
   h_dybb_good_goodJ2_pt->Write();             h_dybb_good_goodJ2_pt->Delete();
   h_dybb_good_goodJ1_eta->Write();            h_dybb_good_goodJ1_eta->Delete();
   h_dybb_good_goodJ2_eta->Write();            h_dybb_good_goodJ2_eta->Delete();
   h_dybb_good_goodJ1_phi->Write();            h_dybb_good_goodJ1_phi->Delete();
   h_dybb_good_goodJ2_phi->Write();            h_dybb_good_goodJ2_phi->Delete();
   h_dybb_good_goodJ1_CSV->Write();            h_dybb_good_goodJ1_CSV->Delete();
   h_dybb_good_goodJ2_CSV->Write();            h_dybb_good_goodJ2_CSV->Delete();
   h_dybb_good_dRgoodJ1J2->Write();            h_dybb_good_dRgoodJ1J2->Delete();
   h_dybb_good_goodLep1_pt->Write();            h_dybb_good_goodLep1_pt->Delete();
   h_dybb_good_goodLep1_eta->Write();           h_dybb_good_goodLep1_eta->Delete();
   h_dybb_good_goodLep1_phi->Write();           h_dybb_good_goodLep1_phi->Delete();
   h_dybb_good_goodLep2_pt->Write();            h_dybb_good_goodLep2_pt->Delete();
   h_dybb_good_goodLep2_eta->Write();           h_dybb_good_goodLep2_eta->Delete();
   h_dybb_good_goodLep2_phi->Write();           h_dybb_good_goodLep2_phi->Delete();
                                            
   h_dybb_qcd_goodJ1_pt->Write();             h_dybb_qcd_goodJ1_pt->Delete();
   h_dybb_qcd_goodJ2_pt->Write();             h_dybb_qcd_goodJ2_pt->Delete();
   h_dybb_qcd_goodJ1_eta->Write();            h_dybb_qcd_goodJ1_eta->Delete();
   h_dybb_qcd_goodJ2_eta->Write();            h_dybb_qcd_goodJ2_eta->Delete();
   h_dybb_qcd_goodJ1_phi->Write();            h_dybb_qcd_goodJ1_phi->Delete();
   h_dybb_qcd_goodJ2_phi->Write();            h_dybb_qcd_goodJ2_phi->Delete();
   h_dybb_qcd_goodJ1_CSV->Write();            h_dybb_qcd_goodJ1_CSV->Delete();
   h_dybb_qcd_goodJ2_CSV->Write();            h_dybb_qcd_goodJ2_CSV->Delete();
   h_dybb_qcd_dRgoodJ1J2->Write();            h_dybb_qcd_dRgoodJ1J2->Delete();
   h_dybb_qcd_qcdLep1_pt->Write();            h_dybb_qcd_qcdLep1_pt->Delete();
   h_dybb_qcd_qcdLep1_eta->Write();           h_dybb_qcd_qcdLep1_eta->Delete();
   h_dybb_qcd_qcdLep1_phi->Write();           h_dybb_qcd_qcdLep1_phi->Delete();
   h_dybb_qcd_qcdLep2_pt->Write();            h_dybb_qcd_qcdLep2_pt->Delete();
   h_dybb_qcd_qcdLep2_eta->Write();           h_dybb_qcd_qcdLep2_eta->Delete();
   h_dybb_qcd_qcdLep2_phi->Write();           h_dybb_qcd_qcdLep2_phi->Delete();

 outfile->Close();
}
