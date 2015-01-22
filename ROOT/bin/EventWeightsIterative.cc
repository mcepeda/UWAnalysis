#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TFile.h"
#include "TROOT.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "boost/filesystem.hpp"
#include <cstdlib>

std::vector<float> data;
std::vector<float> mc;
edm::Lumi3DReWeighting *LumiWeights;
edm::LumiReWeighting *LumiWeightsOld;

void readdir(TDirectory *dir,optutl::CommandLineParser parser,float ev,int doPU,bool doRho,TH1F* puWeight,TH1F* rhoWeight); 

int main (int argc, char* argv[]){
 optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
 parser.addOption("histoName",optutl::CommandLineParser::kString,"Counter Histogram Name","EventSummary");
 parser.addOption("weight",optutl::CommandLineParser::kDouble,"Weight to apply",1.0);
 parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__WEIGHT__");
 parser.addOption("doOneD",optutl::CommandLineParser::kInteger,"Do OneD",0);
 parser.addOption("tHack",optutl::CommandLineParser::kInteger,"tHack",0);
 
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

 // 3D maybe not handled correctly - works for doOneD==1
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

 // FILE TO UNDERGO WEIGHTING
 TFile *f = new TFile(parser.stringValue("outputFile").c_str(),"UPDATE");
 // get total number of generated events
 TH1F* evC  = (TH1F*)f->Get(parser.stringValue("histoName").c_str());
 float ev = evC->GetBinContent(1);
 if(parser.integerValue("tHack")){
  ev = 2.497226e+06;
 }

 printf("Found  %f Events Counted\n",ev);
 readdir(f,parser,ev,doPU,doRho,puWeight,rhoWeight); //actually do the weighting
 f->Close();
} 

void readdir(TDirectory *dir,optutl::CommandLineParser parser,float ev,int doPU,bool doRho,TH1F *puWeight,TH1F *rhoWeight){

 TDirectory *dirsav = gDirectory;
 TIter next(dir->GetListOfKeys());
 TKey *key;

 while ((key = (TKey*)next())) {
  TString string = (TString) key->GetName();
  printf("Found key = %s \n",key->GetName());

  //to run faster, take out this line when doing systematics
//  if(string == "muNuEventTreeMuonUp")  return;

//  if(string == "summary")  return;

  TObject *obj = key->ReadObj();
  // if in a directory, recursively call function
  if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
   dir->cd(key->GetName());
   TDirectory *subdir = gDirectory;
   readdir(subdir,parser,ev,doPU,doRho,puWeight,rhoWeight);
   dirsav->cd();
  }
  // if at a tree, do the weighting
  else if(obj->IsA()->InheritsFrom(TTree::Class())) {
   TTree *t = (TTree*)obj;
   float weight = parser.doubleValue("weight")/(ev);
 
   TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&weight,(parser.stringValue("branch")+"/F").c_str());
   int vertices;
   float bxm=0;
   float bx=0;
   float bxp=0;

   if(doPU==1)
    t->SetBranchAddress("vertices",&vertices);
   else if(doPU==2) {
    t->SetBranchAddress("puBXminus",&bxm);
    t->SetBranchAddress("puBX0",&bx);
    t->SetBranchAddress("puBXplus",&bxp);
   }
   else if(doPU==3)
    t->SetBranchAddress("puTruth",&bx);

   float rho=0.0;
   if(doRho)
     t->SetBranchAddress("Rho",&rho);
   int x = 0;
   printf("Found tree -> weighting\n");
   for(Int_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);
    weight = parser.doubleValue("weight")/(ev);// = xsec / nr events
    if(doPU==1) {
     int bin=puWeight->FindBin(vertices); // bin w/ n vertices (in pu reweight histogram)
     if(bin > puWeight->GetNbinsX()){
      x +=1;
      printf("Overflow using max bin. This is time %i\n",x);
      bin = puWeight->GetNbinsX();
     }
     weight*=puWeight->GetBinContent(bin);//scale weight by entry in pu histo
     if(i==1){
      printf("PU WEIGHT = %f\n",puWeight->GetBinContent(puWeight->FindBin(vertices)));
     }
    }
    else if(doPU==2) {
     float w = LumiWeights->weight3D( bxm,bx,bxp);
      if(i==1){
       printf("PU WEIGHT = %f\n",w);
      }
      weight*=w;
    }
    else if(doPU==3) {
     weight*=LumiWeightsOld->weight(bx);
    }
    if(doRho) {
      weight*=rhoWeight->GetBinContent(rhoWeight->FindBin(rho));
      if(i==1){
       printf("RHO WEIGHT = %f\n",rhoWeight->GetBinContent(rhoWeight->FindBin(rho)));
      }
    }
    newBranch->Fill();
   }
   printf("Overflow bins used %i\n",x);
   t->Write("",TObject::kOverwrite);
   return;
  }
 }
}
