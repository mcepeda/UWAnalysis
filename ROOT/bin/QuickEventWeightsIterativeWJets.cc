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

void weightdir(TDirectory *dir,optutl::CommandLineParser parser,std::vector<float> ev); 

int main (int argc, char* argv[]){
 optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
 parser.addOption("histoName",optutl::CommandLineParser::kString,"Counter Histogram Name","EventSummary");
 parser.addOption("weight",optutl::CommandLineParser::kDouble,"Weight to apply",1.0);
 parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__WEIGHT__");
 
 parser.parseArguments (argc, argv);


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
 weightdir(f0,parser,ev);
 f0->Close();

 TFile *f1 = new TFile("W1Jet.root","UPDATE");
 weightdir(f1,parser,ev);
 f1->Close();

 TFile *f2 = new TFile("W2Jet.root","UPDATE");
 weightdir(f2,parser,ev);
 f2->Close();

 TFile *f3 = new TFile("W3Jet.root","UPDATE");
 weightdir(f3,parser,ev);
 f3->Close();

 TFile *f4 = new TFile("W4Jet.root","UPDATE");
 weightdir(f4,parser,ev);
 f4->Close();

} 

void weightdir(TDirectory *dir,optutl::CommandLineParser parser,std::vector<float> ev){

 TDirectory *dirsav = gDirectory;
 TIter next(dir->GetListOfKeys());
 TKey *key;

 while ((key = (TKey*)next())) {
  TString string = (TString) key->GetName();
  printf("Found key = %s \n",key->GetName());

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
   TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&weight,(parser.stringValue("branch")+"/F").c_str());

   int LHEProduct=0;
   t->SetBranchAddress("LHEProduct",&LHEProduct);
   float pw = 0.;
   t->SetBranchAddress("preWeight",&pw);

   int x = 0;
   float xOe = 0;
   printf("Found tree -> weighting\n");
   for(Int_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);
    if(LHEProduct<5)
     printf("LHE Product is less then 5, check what is happeing!");
    else if(LHEProduct==5)
     xOe = parser.doubleValue("weight")/(ev[0]);
    else if(LHEProduct==6)
     xOe = parser.doubleValue("weight")/(ev[0]+ev[1]);
    else if(LHEProduct==7)
     xOe = parser.doubleValue("weight")/(ev[0]+ev[2]);
    else if(LHEProduct==8)
     xOe = parser.doubleValue("weight")/(ev[0]+ev[3]);
    else
     xOe = parser.doubleValue("weight")/(ev[0]+ev[4]);
    weight=pw*xOe;// = Data/MC * xsec / nr events
    newBranch->Fill();
   }
   printf("Overflow bins used %i\n",x);
   t->Write("",TObject::kOverwrite);
   return;
  }
 }
}

