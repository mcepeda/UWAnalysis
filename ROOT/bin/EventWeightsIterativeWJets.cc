#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TFile.h"
#include "TROOT.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFileMerger.h"
#include "boost/filesystem.hpp"
#include <cstdlib>

std::vector<float> data;
std::vector<float> mc;
edm::Lumi3DReWeighting *LumiWeights;
edm::LumiReWeighting *LumiWeightsOld;

void readdir(TDirectory *dir,optutl::CommandLineParser parser,std::vector<float> ev,int doPU,bool doRho,TH1F* puWeight,TH1F* rhoWeight); 

int main (int argc, char* argv[]) {
 optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
 parser.addOption("histoName",optutl::CommandLineParser::kString,"Counter Histogram Name","EventSummary");
 parser.addOption("weight",optutl::CommandLineParser::kDouble,"Weight to apply",1.0);
 parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__WEIGHT__");
 parser.addOption("doOneD",optutl::CommandLineParser::kInteger,"Do OneD",0);
 
 parser.parseArguments (argc, argv);

 TFile fPU,fPU2,fPU22,fPU3,fPU4;

 std::string base = std::getenv("CMSSW_BASE");
 std::string fPUloc =     "/src/UWAnalysis/Configuration/data/puInfo2012.root";
 std::string fPU3Dloc =   "/src/UWAnalysis/Configuration/data/puInfo3D.root";
 std::string fPUMC3Dloc = "/src/UWAnalysis/Configuration/data/puInfoMC3D.root";
 std::string fPUW3Dloc =  "/src/UWAnalysis/Configuration/data/Weight3D.root";
 std::string fPUW3Dbloc = "/src/UWAnalysis/Configuration/data/Weight3D.root";
 std::string fPUMCloc =   "/src/UWAnalysis/Configuration/data/MC_Summer12_PU_S10-600bins.root";
 std::string fPUDATAloc = "/src/UWAnalysis/Configuration/data/Data_Pileup_2012_ReReco-600bins.root";
 std::string fPUname =     base+fPUloc;
 std::string fPU3Dname =   base+fPU3Dloc ;
 std::string fPUMC3Dname = base+fPUMC3Dloc;
 std::string fPUW3Dname =  base+fPUW3Dloc;
 std::string fPUW3Dbname = base+fPUW3Dbloc ;
 std::string fPUMCname =   base+fPUMCloc;
 std::string fPUDATAname = base+fPUDATAloc;
 bool fPUis     = boost::filesystem::exists( fPUname     );
 bool fPU3Dis   = boost::filesystem::exists( fPU3Dname   );
 bool fPUMC3Dis = boost::filesystem::exists( fPUMC3Dname );
 bool fPUW3Dis  = boost::filesystem::exists( fPUW3Dname  );
 bool fPUW3Dbis = boost::filesystem::exists( fPUW3Dbname );
 bool fPUMCis   = boost::filesystem::exists( fPUMCname   );
 bool fPUDATAis = boost::filesystem::exists( fPUDATAname );

 //read PU info
 TH1F *puWeight=0;
 int doPU=0;
 if(fPUis){
  TFile *fPU = new TFile(TString(fPUname));
  puWeight = (TH1F*)fPU->Get("weight");
  doPU=1;
  printf("ENABLING PU WEIGHTING USING VERTICES\n");
 }
 else if(fPUW3Dis){
  TFile *fPU3 = new TFile(TString(fPUW3Dname));
  doPU=2;
  printf("ENABLING PU WEIGHTING USING 3D with ready distribution\n");
  fPU3->Close();
  LumiWeights = new edm::Lumi3DReWeighting(mc,data,"pileup");
  LumiWeights->weight3D_init(fPUW3Dname);
 }
 else if(fPUW3Dbis){
  TFile *fPU4 = new TFile(TString(fPUW3Dbname));
  doPU=2;
  printf("ENABLING PU WEIGHTING USING 3D with  distribution you just made\n");
  fPU4->Close();
  LumiWeights = new edm::Lumi3DReWeighting(mc,data,"pileup");
  LumiWeights->weight3D_init(fPUW3Dbname);
 }
 else if(fPU3Dis && fPUMC3Dis){
  TFile *fPU2 = new TFile( "$CMSSW_BASE/src/UWAnalysis/Configuration/data/puInfo3D.root");
  TFile *fPU22 = new TFile("$CMSSW_BASE/src/UWAnalysis/Configuration/data/puInfoMC3D.root");
  doPU=2;
  printf("ENABLING PU WEIGHTING USING 3D- I HAVE TO CALCULATE WEIGHTS SORRY\n");
  LumiWeights = new edm::Lumi3DReWeighting(
   "$CMSSW_BASE/src/UWAnalysis/Configuration/data/puInfoMC3D.root",
   "$CMSSW_BASE/src/UWAnalysis/Configuration/data/puInfo3D.root",
   "pileup","pileup","pileup");
  LumiWeights->weight3D_init(1.0);
 }
 if(parser.integerValue("doOneD")) {
  if(fPUMCis && fPUDATAis){
  doPU=3;
  LumiWeightsOld = new edm::LumiReWeighting(fPUMCname,fPUDATAname,"pileup","pileup");
  }
  else std::cout<<"Can't find MC and DATA PU distribution files"<<std::endl;
 }
 
 //rho reweighting
 TH1F *rhoWeight=0;
 bool doRho=false;
 std::string fRholoc = "/src/UWAnalysis/Configuration/data/rhoInfo.root";
 std::string fRhoname = base+fRholoc ;
 bool fRhois = boost::filesystem::exists(fRhoname);
 if(fRhois){
  TFile *fRho = new TFile("$CMSSW_BASE/src/UWAnalysis/Configuration/data/rhoInfo.root");
  rhoWeight = (TH1F*)fRho->Get("weight");
  doRho=true;
  printf("ENABLING Rho WEIGHTING\n");
 }
 
 TFile *w = new TFile("WJets.root","UPDATE");
 TH1F* evC  = (TH1F*)w->Get(parser.stringValue("histoName").c_str());
 float evW = evC->GetBinContent(1);
 w->Close();

 TFile *w1 = new TFile("W1Jet.root","UPDATE");
 TH1F* evC1  = (TH1F*)w1->Get(parser.stringValue("histoName").c_str());
 float evW1 = evC1->GetBinContent(1);
 w1->Close();   

 TFile *w2 = new TFile("W2Jet.root","UPDATE");
 TH1F* evC2  = (TH1F*)w2->Get(parser.stringValue("histoName").c_str());
 float evW2 = evC2->GetBinContent(1);
 w2->Close();

 TFile *w3 = new TFile("W3Jet.root","UPDATE");
 TH1F* evC3  = (TH1F*)w3->Get(parser.stringValue("histoName").c_str());
 float evW3 = evC3->GetBinContent(1);
 w3->Close();

 TFile *w4 = new TFile("W4Jet.root","UPDATE");
 TH1F* evC4  = (TH1F*)w4->Get(parser.stringValue("histoName").c_str());
 float evW4 = evC4->GetBinContent(1);
 w4->Close();
    
 printf("Found  %f W+nJet Events\n",evW);
 printf("Found  %f W+1Jet Events\n",evW1);
 printf("Found  %f W+2Jet Events\n",evW2);
 printf("Found  %f W+3Jet Events\n",evW3);
 printf("Found  %f W+4Jet Events\n",evW4);
 
 std::vector<float> ev;
 ev.push_back(evW);
 ev.push_back(evW1/0.175999);
 ev.push_back(evW2/0.0562617);
 ev.push_back(evW3/0.0168926);
 ev.push_back(evW4/0.00692218);
 
 TFile *f0 = new TFile("WJets.root","UPDATE");   
 readdir(f0,parser,ev,doPU,doRho,puWeight,rhoWeight);
 f0->Close();
 
 TFile *f1 = new TFile("W1Jet.root","UPDATE");   
 readdir(f1,parser,ev,doPU,doRho,puWeight,rhoWeight);
 f1->Close();
 
 TFile *f2 = new TFile("W2Jet.root","UPDATE");   
 readdir(f2,parser,ev,doPU,doRho,puWeight,rhoWeight);
 f2->Close();
 
 TFile *f3 = new TFile("W3Jet.root","UPDATE");   
 readdir(f3,parser,ev,doPU,doRho,puWeight,rhoWeight);
 f3->Close();
 
 TFile *f4 = new TFile("W4Jet.root","UPDATE");   
 readdir(f4,parser,ev,doPU,doRho,puWeight,rhoWeight);
 f4->Close();
 
// if(fPU!=0 && fPU->IsOpen())
//   fPU->Close();
//
// if(fPU2!=0 && fPU2->IsOpen())
//   fPU2->Close();
} 


void readdir(TDirectory *dir,optutl::CommandLineParser parser,std::vector<float> ev,int doPU,bool doRho,TH1F *puWeight,TH1F *rhoWeight){
  TDirectory *dirsav = gDirectory;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
   TString string = (TString) key->GetName();
   printf("Found key=%s \n",key->GetName());
  //to run faster, take out this line when doing systematics
  //if(string == "summary")  return;
  //if(string == "muNuEventTreeMuonUp")  return;

  TObject *obj = key->ReadObj();

  if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
   dir->cd(key->GetName());
   TDirectory *subdir = gDirectory;
   readdir(subdir,parser,ev,doPU,doRho,puWeight,rhoWeight);
   dirsav->cd();
  }
  else if(obj->IsA()->InheritsFrom(TTree::Class())) {
   TTree *t = (TTree*)obj;
   float weight;
 
   TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&weight,(parser.stringValue("branch")+"/F").c_str());
   int vertices;
   float bxm=0;
   float bx=0;
   float bxp=0;
   int LHEProduct=0;
   t->SetBranchAddress("LHEProduct",&LHEProduct);
 
   if(doPU==1)
    t->SetBranchAddress("vertices",&vertices);
   else if(doPU==2) {
    t->SetBranchAddress("puBXminus",&bxm);
    t->SetBranchAddress("puBX0",&bx);
    t->SetBranchAddress("puBXplus",&bxp);
   }
   else if( doPU==3 ){
    t->SetBranchAddress("puTruth",&bx);
   }
 
   float rho=0.0;
   if(doRho)
    t->SetBranchAddress("Rho",&rho);
 
   printf("Found tree -> weighting\n");
   for(Int_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);
    if(LHEProduct<5)
     printf("LHE Product is less then 5, check what is happeing!");
    else if(LHEProduct==5)
     weight = parser.doubleValue("weight")/(ev[0]);
    else if(LHEProduct==6)
     weight = parser.doubleValue("weight")/(ev[0]+ev[1]);
    else if(LHEProduct==7)
     weight = parser.doubleValue("weight")/(ev[0]+ev[2]);
    else if(LHEProduct==8)
     weight = parser.doubleValue("weight")/(ev[0]+ev[3]);
    else
     weight = parser.doubleValue("weight")/(ev[0]+ev[4]);
    	
    if(doPU==1) {
     int bin=puWeight->FindBin(vertices);
     if(bin>puWeight->GetNbinsX()){
      printf("Overflow using max bin\n");
      bin = puWeight->GetNbinsX();
     }
     weight*=puWeight->GetBinContent(bin);
     if(i==1)
      printf("PU WEIGHT = %f\n",puWeight->GetBinContent(puWeight->FindBin(vertices)));
    }
    else if(doPU==2) {
     float w = LumiWeights->weight3D( bxm,bx,bxp);
     if(i==1)
      printf("PU WEIGHT = %f\n",w);
     weight*=w;
    }
    else if(doPU==3) {
     weight*=LumiWeightsOld->weight(bx);
    }
    if(doRho) {
     weight*=rhoWeight->GetBinContent(rhoWeight->FindBin(rho));
     if(i==1)
      printf("RHO WEIGHT = %f\n",rhoWeight->GetBinContent(rhoWeight->FindBin(rho)));
    }
 
    newBranch->Fill();
   } 
   t->Write("",TObject::kOverwrite);
  }
 }
}
