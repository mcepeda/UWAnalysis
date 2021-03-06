#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
using namespace std;



void dumpContent(ofstream& f,TH1* h, int sf) {
  int m = h->GetNbinsX();
  for(int i = 0; i <= m+1; i++) {
    f << setw(9) << ((int) (sf * h->GetBinContent(i)));
    if(!((i+1)%10) || i == m+1) f << endl;
  }
}

void dumpError(ofstream& f,TH1* h, int sf) {
  int m = h->GetNbinsX();
  for(int i = 0; i <= m+1; i++) {

    float content=0.0;
    if(h->GetBinContent(i)>0.0)
      content = h->GetBinError(i)/h->GetBinContent(i);

    f << setw(9) << ((int) (sf * content));
    if(!((i+1)%10) || i == m+1) f << endl;
  }
}



void dumpNull(ofstream& f,TH1F *h) {
  int m = h->GetNbinsX();
  for(int i = 0; i <= m+1; i++) {
    f << setw(9) << 0;
    if(!((i+1)%10) ||i==m+1) f << endl;
  }
}


void createASCII(TString channel,TString mass,TString syst,TString histogram,float scaleFactor,float normalizationAbs = 0.0 ,float normalizationRel =1.0,TString dir="/",TString overrideH = "")
{
  gSystem->MakeDirectory(syst+dir);


  TFile f(channel+"-"+mass+"-"+syst+".root");

  if(f.IsZombie()) return;

  TH1F *h = (TH1F*)f.Get(histogram);

  //find the integral 
  if(normalizationAbs == 0.0) {
    h->Scale(normalizationRel);
  }
  else {
    h->Scale(1.0/h->Integral());
    h->Scale(normalizationAbs);
  }
  //create shape file 
  ofstream file;



  if(overrideH =="")
    file.open ((syst+"/"+dir+channel+"-"+histogram+".hst").Data());
  else
    file.open ((syst+"/"+dir+channel+"-"+overrideH+".hst").Data());
  file << setw(9) << 3;
  file << setw(9) << h->GetNbinsX()+2;
  file << setw(9) << ((int) (scaleFactor)); 
  file << setw(9) << ((int) (scaleFactor * h->Integral()));
  file << endl;


  if(channel=="ele-tau") dumpContent(file,h,scaleFactor);
  else dumpNull(file,h);
  if(channel=="mu-tau") dumpContent(file,h,scaleFactor);
  else dumpNull(file,h);
  if(channel=="ele-mu") dumpContent(file,h,scaleFactor);
  else dumpNull(file,h);
  file.close();	     



  //create error file 
  //if scale factor==1.0 assume Poisson weighted statistics
  //so dont create error file

  if(scaleFactor==1.0) { f.Close(); return;}


  ofstream file2;
  if(overrideH =="")
    file2.open ((syst+"/"+dir+channel+"-"+histogram+".unc").Data());
  else
    file2.open ((syst+"/"+dir+channel+"-"+overrideH+".unc").Data());

  file2 << setw(9) << 3;
  file2 << setw(9) << h->GetNbinsX()+2;
  file2 << setw(9) << ((int) (scaleFactor)); 
  file2 << setw(9) << ((int) (scaleFactor));
  file2 << endl;

  if(channel=="ele-tau") dumpError(file2,h,scaleFactor);
  else dumpNull(file2,h);
  if(channel=="mu-tau") dumpError(file2,h,scaleFactor);
  else dumpNull(file2,h);
  if(channel=="ele-mu") dumpError(file2,h,scaleFactor);
  else dumpNull(file2,h);
  file2.close();	     
  f.Close();
}


void createIterativeASCII(TString channel,TString mass,TString syst,std::vector<TString> histogramPrefix,std::vector<int> massPrefix,float scaleFactor)
{
  for(unsigned int i=0;i<histogramPrefix.size();++i)
    for(unsigned int j=0;j<massPrefix.size();++j) {
      TString histogram;
	histogram = histogramPrefix[i]+TString::Format("%d",massPrefix[j]);
      if(massPrefix[j]!=90) 
	createASCII(channel,mass,syst,histogram,scaleFactor,0.0,1.0,"/"+TString::Format("m%d",massPrefix[j])+"/",histogramPrefix[i]);
      else
	createASCII(channel,mass,syst,histogram,scaleFactor,0.0,1.0,"/"+TString::Format("m0%d",massPrefix[j])+"/",histogramPrefix[i]);



      }
}
 
void runFullChain(TString systematic,float scaleFactor,TString mass)
{

  //create Higgs list
  std::vector<TString> higgsPrefix;
  higgsPrefix.push_back("BBA");
  higgsPrefix.push_back("GGH");

  std::vector<int> higgsMass;
  higgsMass.push_back(90);
  higgsMass.push_back(100);
  higgsMass.push_back(120);
  higgsMass.push_back(130);
  higgsMass.push_back(140);
  higgsMass.push_back(160);
  higgsMass.push_back(180);
  higgsMass.push_back(200);
  higgsMass.push_back(250);
  higgsMass.push_back(300);
  higgsMass.push_back(350);
  higgsMass.push_back(400);
  higgsMass.push_back(450);
  higgsMass.push_back(500);
  




  //mu+tau
  createASCII("mu-tau",mass,systematic,"ZTT",scaleFactor,0.0 ,1.0);
  createASCII("mu-tau",mass,systematic,"QCDD",scaleFactor,131.0 ,1.0);
  createASCII("mu-tau",mass,systematic,"WMN",scaleFactor,54.8 ,1.0);
  createASCII("mu-tau",mass,systematic,"WTN",scaleFactor,14.7 ,1.0);
  createASCII("mu-tau",mass,systematic,"ZMM2",scaleFactor,6.4 ,1.0);
  createASCII("mu-tau",mass,systematic,"ZMMD",scaleFactor,13.2 ,1.0);
  createASCII("mu-tau",mass,systematic,"TTBar",scaleFactor,6.0 ,1.0);
  createASCII("mu-tau",mass,systematic,"DiBoson",scaleFactor,1.6 ,1.0);
  createASCII("mu-tau",mass,systematic,"DATA",1.0,0.0 ,1.0);

  createIterativeASCII("mu-tau",mass,systematic, higgsPrefix,higgsMass,scaleFactor);


  //e+tau  
  createASCII("ele-tau",mass,systematic,"ZTT",scaleFactor,0.0 ,1.0);
  createASCII("ele-tau",mass,systematic,"QCDD",scaleFactor,181. ,1.0);
  createASCII("ele-tau",mass,systematic,"WEN",scaleFactor,31.0 ,1.0);
  createASCII("ele-tau",mass,systematic,"WTN",scaleFactor,7.0 ,1.0);
  createASCII("ele-tau",mass,systematic,"ZEE2",scaleFactor,15.0 ,1.0);
  createASCII("ele-tau",mass,systematic,"ZEED",scaleFactor,109.0 ,1.0);
  createASCII("ele-tau",mass,systematic,"TTBar",scaleFactor,2.6 ,1.0);
  createASCII("ele-tau",mass,systematic,"DiBoson",scaleFactor,0.8 ,1.0);
  createASCII("ele-tau",mass,systematic,"DATA",1.0,0.0,1.0);

  createIterativeASCII("ele-tau",mass,systematic, higgsPrefix,higgsMass,scaleFactor);


  //e+mu  
  createASCII("ele-mu",mass,systematic,"ZTT",scaleFactor,0.0 ,1.0);
  createASCII("ele-mu",mass,systematic,"FAKES",scaleFactor,4.14 ,1.0);
  createASCII("ele-mu",mass,systematic,"EWK",scaleFactor,3.24 ,1.0);
  createASCII("ele-mu",mass,systematic,"TTBar",scaleFactor,7.5 ,1.0);
  createASCII("ele-mu",mass,systematic,"DATA",1.0,0.0 ,1.0);

  createIterativeASCII("ele-mu",mass,systematic, higgsPrefix,higgsMass,scaleFactor);
  



}

void ASCIIScripts() {
  TString mass =  "sv";
  float scaleFactor=100000.0;
  runFullChain("nominal",scaleFactor,mass);
  runFullChain("tauUp",scaleFactor,mass);
  runFullChain("tauDown",scaleFactor,mass);
  runFullChain("eleUp",scaleFactor,mass);
  runFullChain("eleDown",scaleFactor,mass);
  runFullChain("muUp",scaleFactor,mass);
  runFullChain("muDown",scaleFactor,mass);
  runFullChain("jetUp",scaleFactor,mass);
  runFullChain("jetDown",scaleFactor,mass);
  runFullChain("uncUp",scaleFactor,mass);
  runFullChain("uncDown",scaleFactor,mass);

 

}
