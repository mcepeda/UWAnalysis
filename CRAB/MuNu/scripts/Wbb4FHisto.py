#!/usr/bin/env python
'''
Makes a root file of Wbb 4 Flavor histograms
Author: T.M.Perry UW
'''
import ROOT
from ROOT import TGraph,THStack,TLegend,TCanvas,TFile,TH1F, TPad
from ROOT import TLatex
from ROOT import gROOT,gStyle
from ROOT import *
import os

import aHisto as h #function to make histograms
import TheCuts as ct  #function which makes cut strings
import histoRange as hr #manages range, lables for plots
import TheParameters as p

eventTreeLocations = [
'muNuEventTree/eventTree',
'muNuEventTreeMuonUp/eventTree',
'muNuEventTreeMuonDown/eventTree',
'muNuEventTreeJetUp/eventTree',
'muNuEventTreeJetDown/eventTree',
'muNuEventTreeUCEUp/eventTree',
'muNuEventTreeUCEDown/eventTree'
]
I=-1
F=-1

#get parameters (used in cutmaker)
lumi,bNr,btype,jNr,njetcut,jetcut,I,F,iso_value,antiIso_value,path,extraName,leafs,drawW,drawZ,drawQCD,drawData,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,eventTreeLocation,extraCut,Ctl_Andrea = p.arams()

CutsMCn, CutsMCnW, CutsMCi,CutsDatan,CutsDatai,CutsMCnwl,CutsMCiwl,CutsMCnwc,CutsMCiwc,CutsMCnwcc,CutsMCiwcc,CutsMCnwbb,CutsMCiwbb,CutsMCnT,CutsMCiT,CutsMCnTup,CutsMCiTup,CutsMCnTdn,CutsMCiTdn, cutMcNonIso_Bup, cutMcNonIsoW_Bup, cutMcIso_Bup, cutMcWlNonIso_Bup, cutMcWlIso_Bup, cutMcWcNonIso_Bup, cutMcWcIso_Bup, cutMcWccNonIso_Bup, cutMcWccIso_Bup, cutMcWbbNonIso_Bup, cutMcWbbIso_Bup, cutMcNonIso_Bdn, cutMcNonIsoW_Bdn, cutMcIso_Bdn, cutMcWlNonIso_Bdn, cutMcWlIso_Bdn, cutMcWcNonIso_Bdn, cutMcWcIso_Bdn, cutMcWccNonIso_Bdn, cutMcWccIso_Bdn, cutMcWbbNonIso_Bdn, cutMcWbbIso_Bdn = ct.cutmaker(
 leafs[0],iso_value,antiIso_value,lumi,bNr,btype,jNr,njetcut,jetcut,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,extraCut,Ctl_Andrea
)

for leaf in leafs:
 print("------------")
 print(leaf)
 print("------------")
 steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)
 version = os.environ.get('version')
 W4F_filename       = '/afs/hep.wisc.edu/cms/tperry/Wbb_CMSSW_5_3_14_patch1/src/UWAnalysis/CRAB/MuNu/data/%s/Wbb4F.root'%(version)
 W4F_file        = TFile( W4F_filename )
 
 outname = "%s%s_4F_%s"%(path,extraName,leaf)
 outFile=gROOT.FindObject('%s.root'%(outname))
 if outFile : outFile.Close()
 outFile = TFile('%s.root'%(outname),'RECREATE','4Flavor Histograms')
 
 for etl in eventTreeLocations:
  W4F_tree       =  W4F_file.Get(etl)
  if etl == 'muNuEventTreeMuonUp/eventTree' : ps ='_muonUp'
  if etl == 'muNuEventTreeMuonDown/eventTree' : ps ='_muonDown'
  if etl == 'muNuEventTreeJetUp/eventTree' : ps ='_jetUp'
  if etl == 'muNuEventTreeJetDown/eventTree' : ps ='_jetDown'
  if etl == 'muNuEventTreeUCEUp/eventTree': ps = '_uceUp'
  if etl == 'muNuEventTreeUCEDown/eventTree': ps = '_uceDown'
  if etl == 'muNuEventTree/eventTree' : ps = ''
  
  w4fh,w4fhSize,w4fhSizePart,w4fhEntries = h.gram(W4F_tree,leaf,xmin,xmax,steps,CutsMCiwbb,I,F)
  w4fh.SetName("Wbb4F%s"%(ps))
  print('  '+str(w4fhSize))

  if etl == 'muNuEventTree/eventTree':
   Wbb4F_BeffUp,wbbBupSize,wbbBupSizePart,wbbBupEntries = h.gram(W4F_tree,leaf,xmin,xmax,steps,cutMcWbbIso_Bup,I,F)
   Wbb4F_BeffUp.SetName('Wbb4F_BeffUp')
   Wbb4F_BeffDown,wbbBdnSize,wbbBdnSizePart,wbbBdnEntries = h.gram(W4F_tree,leaf,xmin,xmax,steps,cutMcWbbIso_Bdn,I,F)
   Wbb4F_BeffDown.SetName('Wbb4F_BeffDown')

  outFile.Write()
 outFile.Close()
 print('')
 print('Your File is here: '+outname+'.root')
 print('') 
