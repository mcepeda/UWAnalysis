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

void weightdir(TDirectory *dir,optutl::CommandLineParser parser,float ev); 

int main (int argc, char* argv[]){
 optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
 parser.addOption("histoName",optutl::CommandLineParser::kString,"Counter Histogram Name","EventSummary");
 parser.addOption("weight",optutl::CommandLineParser::kDouble,"Weight to apply",1.0);
 parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__WEIGHT__");
 
 parser.parseArguments (argc, argv);

 // FILE TO UNDERGO WEIGHTING
 TFile *f = new TFile(parser.stringValue("outputFile").c_str(),"UPDATE");
 TH1F* evC  = (TH1F*)f->Get(parser.stringValue("histoName").c_str());
 // get total number of generated events
 float ev = evC->GetBinContent(1);

 printf("Found  %f Events Counted\n",ev);
 weightdir(f,parser,ev); //actually do the weighting
 f->Close();
} 

void weightdir(TDirectory *dir,optutl::CommandLineParser parser,float ev){

 TDirectory *dirsav = gDirectory;
 TIter next(dir->GetListOfKeys());
 TKey *key;

 while ((key = (TKey*)next())) {
  TString string = (TString) key->GetName();
  printf("Found key = %s \n",key->GetName());

  //to run faster, take out this line when doing systematics
//  if(string == "muNuEventTreeMuonUp")  return;

  if(string == "summary")  return;

  TObject *obj = key->ReadObj();
  // if in a directory, recursively call function
  if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
   dir->cd(key->GetName());
   TDirectory *subdir = gDirectory;
   weightdir(subdir,parser,ev);
   dirsav->cd();
  }
  // if at a tree, do the weighting
  else if(obj->IsA()->InheritsFrom(TTree::Class())) {
   TTree *t = (TTree*)obj;
   float weight = 0.;
   float xOe = parser.doubleValue("weight")/(ev); // xsec/nr events
 
   TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&weight,(parser.stringValue("branch")+"/F").c_str());

   float pw = 0.;
   t->SetBranchAddress("preWeight",&pw);

   int x = 0;
   printf("Found tree -> weighting\n");
   for(Int_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);
    weight=pw*xOe;// = Data/MC * xsec / nr events
    newBranch->Fill();
   }
   printf("Overflow bins used %i\n",x);
   t->Write("",TObject::kOverwrite);
   return;
  }
 }
}

