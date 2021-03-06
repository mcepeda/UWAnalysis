{
  gROOT->ProcessLine(".L UWAnalysis/ROOT/interactive/SimplePlotter.C+");
  gROOT->ProcessLine(".L UWAnalysis/ROOT/interactive/tdrstyle.C");
  setTDRStyle();
  
  SimplePlotter *plotter = new SimplePlotter();


    plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/W.root","W+jets","__WEIGHT__*(charge==0)",0,kRed+2,1);
    plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/TOP.root","t#bar{t}","__WEIGHT__*(charge==0)",0,kBlue-8,1);
    plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/DATA.root","SS Data","(charge!=0)/1130.",0,kMagenta-10,1);
    plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/Z.root","Z+jets","__WEIGHT__*(genTaus==0)*(charge==0)",0,kGreen-4,kGreen+2);
    plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/Z.root","Z#rightarrow #tau #tau","__WEIGHT__*(genTaus>0)*(charge==0)",0,kOrange-4,kBlack);

      plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/sm115.root","SM H(115) #rightarrow #tau #tau","1.38*__WEIGHT__",-1,kBlue,kBlue);
      plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/vbf115.root","","0.107*__WEIGHT__",-1,kBlue,kBlue);


  plotter->addFile("diTauEventTree/eventTree","sandbox/tt-latest/DATA.root","DATA","(charge==0)",1,kBlack,0);














}
