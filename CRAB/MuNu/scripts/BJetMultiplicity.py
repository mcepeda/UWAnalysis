#!/usr/bin/env python
'''
Makes a root file containing histograms for later plotting
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

#import scaleqcd as qs

qcdRescale = False
tRescale = False
ttRescale = False
vRescale = False
dRescale = False
sf_qcd = 1.
sf_drell = 1.
sf_st = 1.   
sf_st_t = 1. #(2.26033000000000000e+05)/2497226#1.#1.025431/2.497226   
sf_ttbar = 1.
sf_wjets = 1.
sf_vv = 1.

#get parameters (used in cutmaker)
lumi,bNr,btype,jNr,njetcut,jetcut,I,F,iso_value,antiIso_value,path,extraName,leafs,drawW,drawZ,drawQCD,drawData,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,eventTreeLocation,extraCut,Ctl_Andrea = p.arams() 

CutsMCn, CutsMCnW, CutsMCi,CutsDatan,CutsDatai,CutsMCnwl,CutsMCiwl,CutsMCnwc,CutsMCiwc,CutsMCnwcc,CutsMCiwcc,CutsMCnwbb,CutsMCiwbb,CutsMCnT,CutsMCiT,CutsMCnTup,CutsMCiTup,CutsMCnTdn,CutsMCiTdn, cutMcNonIso_Bup, cutMcNonIsoW_Bup, cutMcIso_Bup, cutMcWlNonIso_Bup, cutMcWlIso_Bup, cutMcWcNonIso_Bup, cutMcWcIso_Bup, cutMcWccNonIso_Bup, cutMcWccIso_Bup, cutMcWbbNonIso_Bup, cutMcWbbIso_Bup, cutMcNonIso_Bdn, cutMcNonIsoW_Bdn, cutMcIso_Bdn, cutMcWlNonIso_Bdn, cutMcWlIso_Bdn, cutMcWcNonIso_Bdn, cutMcWcIso_Bdn, cutMcWccNonIso_Bdn, cutMcWccIso_Bdn, cutMcWbbNonIso_Bdn, cutMcWbbIso_Bdn = ct.cutmaker(
 leafs[0],iso_value,antiIso_value,lumi,bNr,btype,jNr,njetcut,jetcut,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,extraCut,Ctl_Andrea
)

version = os.environ.get('version')
dpath = "/afs/hep.wisc.edu/cms/tperry/Wbb_CMSSW_5_3_14_patch1/src/UWAnalysis/CRAB/MuNu/data/%s/"%(version)
hpath = "/hdfs/store/user/tperry/%s/"%(version)
npath = "/nfs_scratch/tperry/%s/Filtered_"%(version)
data_filename     = npath+'Data.root'
t_t_filename      = dpath+'T_t.root'
t_s_filename      = dpath+'T_s.root'
t_tw_filename     = dpath+'T_tW.root'
tb_t_filename     = npath+'Tbar_t.root'
tb_s_filename     = npath+'Tbar_s.root'
tb_tw_filename    = npath+'Tbar_tW.root'
tt_semi_filename  = npath+'TTbar_semi.root'
tt_full_filename  = npath+'TTbar_full.root'
ww_filename       = dpath+'WW.root'
wz_filename       = dpath+'WZ.root'
zz_filename       = dpath+'ZZ.root'
wn_filename       = npath+'WJets.root'
w1_filename       = npath+'W1Jet.root'
w2_filename       = npath+'W2Jet.root'
w3_filename       = npath+'W3Jet.root'
w4_filename       = npath+'W4Jet.root'
z_filename        = dpath+'Drell.root'

#data_filename     = hpath+'Data.root'
#t_t_filename      = dpath+'T_t.root'
#t_s_filename      = dpath+'T_s.root'
#t_tw_filename     = dpath+'T_tW.root'
#tb_t_filename     = hpath+'Tbar_t.root'
#tb_s_filename     = hpath+'Tbar_s.root'
#tb_tw_filename    = hpath+'Tbar_tW.root'
#tt_semi_filename  = hpath+'TTbar_semi.root'
#tt_full_filename  = hpath+'TTbar_full.root'
#ww_filename       = dpath+'WW.root'
#wz_filename       = dpath+'WZ.root'
#zz_filename       = dpath+'ZZ.root'
#wn_filename       = hpath+'WJets.root'
#w1_filename       = hpath+'W1Jet.root'
#w2_filename       = hpath+'W2Jet.root'
#w3_filename       = hpath+'W3Jet.root'
#w4_filename       = hpath+'W4Jet.root'
#z_filename        = dpath+'Drell.root'

data_file      = TFile( data_filename )
t_t_file       = TFile( t_t_filename  )
t_s_file       = TFile( t_s_filename  )
t_tw_file      = TFile( t_tw_filename )
tb_t_file      = TFile( tb_t_filename )
tb_s_file      = TFile( tb_s_filename )
tb_tw_file     = TFile( tb_tw_filename)
tt_semi_file   = TFile( tt_semi_filename  )
tt_full_file   = TFile( tt_full_filename  )
ww_file        = TFile( ww_filename   )
wz_file        = TFile( wz_filename   )
zz_file        = TFile( zz_filename   )
wn_file        = TFile( wn_filename   )
w1_file        = TFile( w1_filename   )
w2_file        = TFile( w2_filename   )
w3_file        = TFile( w3_filename   )
w4_file        = TFile( w4_filename   )
z_file         = TFile( z_filename    )

eventTreeLocationNoShift = 'muNuEventTree/eventTree'
ETL_filtered   = "muNuEventTree_eventTree"
ETL_unfiltered = 'muNuEventTree/eventTree'

data_tree     =  data_file.Get(ETL_filtered)  # 
t_t_tree      =  t_t_file.Get(ETL_unfiltered)
t_s_tree      =  t_s_file.Get(ETL_unfiltered) 
t_tw_tree     =  t_tw_file.Get(ETL_unfiltered) 
tb_t_tree     =  tb_t_file.Get(ETL_filtered)  #
tb_s_tree     =  tb_s_file.Get(ETL_filtered)  #
tb_tw_tree    =  tb_tw_file.Get(ETL_filtered)  #
tt_semi_tree  =  tt_semi_file.Get(ETL_filtered)  #
tt_full_tree  =  tt_full_file.Get(ETL_filtered)  #
ww_tree       =  ww_file.Get(ETL_unfiltered)
wz_tree       =  wz_file.Get(ETL_unfiltered)
zz_tree       =  zz_file.Get(ETL_unfiltered)
wn_tree       =  wn_file.Get(ETL_filtered)  #
w1_tree       =  w1_file.Get(ETL_filtered)  #
w2_tree       =  w2_file.Get(ETL_filtered)  #
w3_tree       =  w3_file.Get(ETL_filtered)  #
w4_tree       =  w4_file.Get(ETL_filtered)  #
z_tree        =  z_file.Get(ETL_unfiltered)

#data_tree     =  data_file.Get(ETL_unfiltered) # 
#t_t_tree      =  t_t_file.Get(ETL_unfiltered)
#t_s_tree      =  t_s_file.Get(ETL_unfiltered) 
#t_tw_tree     =  t_tw_file.Get(ETL_unfiltered) 
#tb_t_tree     =  tb_t_file.Get(ETL_unfiltered) #
#tb_s_tree     =  tb_s_file.Get(ETL_unfiltered) #
#tb_tw_tree    =  tb_tw_file.Get(ETL_unfiltered) #
#tt_semi_tree  =  tt_semi_file.Get(ETL_unfiltered) #
#tt_full_tree  =  tt_full_file.Get(ETL_unfiltered) #
#ww_tree       =  ww_file.Get(ETL_unfiltered)
#wz_tree       =  wz_file.Get(ETL_unfiltered)
#zz_tree       =  zz_file.Get(ETL_unfiltered)
#wn_tree       =  wn_file.Get(ETL_unfiltered) #
#w1_tree       =  w1_file.Get(ETL_unfiltered) #
#w2_tree       =  w2_file.Get(ETL_unfiltered) #
#w3_tree       =  w3_file.Get(ETL_unfiltered) #
#w4_tree       =  w4_file.Get(ETL_unfiltered) #
#z_tree        =  z_file.Get(ETL_unfiltered)

#Start the Plotting Program
for leaf in leafs:

 steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)
 
 outFile=gROOT.FindObject(path+extraName+'_'+leaf+'.root')
 if outFile : outFile.Close()
 outFile = TFile(path+extraName+'_'+leaf+'.root','RECREATE','Demo ROOT file')
 
 log = open(path+extraName+'_'+leaf+'.log','w')

 log.write("Data    "+data_filename   +"\n") 
 log.write("t_t     "+t_t_filename    +"\n") 
 log.write("t_s     "+t_s_filename    +"\n") 
 log.write("t_tw    "+t_tw_filename   +"\n") 
 log.write("tb_t    "+tb_t_filename   +"\n") 
 log.write("tb_s    "+tb_s_filename   +"\n") 
 log.write("tb_tw   "+tb_tw_filename  +"\n") 
 log.write("tt_semi "+tt_semi_filename+"\n")
 log.write("tt_full "+tt_full_filename+"\n")
 log.write("ww      "+ww_filename     +"\n") 
 log.write("wz      "+wz_filename     +"\n") 
 log.write("zz      "+zz_filename     +"\n") 
 log.write("wn      "+wn_filename     +"\n") 
 log.write("w1      "+w1_filename     +"\n") 
 log.write("w2      "+w2_filename     +"\n") 
 log.write("w3      "+w3_filename     +"\n") 
 log.write("w4      "+w4_filename     +"\n") 
 log.write("z       "+z_filename      +"\n") 
 
 print('----------------------------')
 print('      --'+leaf+'--')
 print(extraName)
 print('----------------------------')
 
 if drawQCD:
  ###
  datanh,datanhSize,datanhSizePart,datanhEntries = h.gram(data_tree,leaf,xmin,xmax,steps,CutsDatan,I,F)
  print('For QCD Shape')
  print('data nonIso')
  datanh.SetName('datanh')
  print('  '+str(datanhSize))
  print('  '+str(datanhSizePart))
 ####  Drell
  print('  z nonIso')
  znh,znhSize,znhSizePart,znhEntries = h.gram(z_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  znh.SetName('znh')
  znh.Scale(sf_drell)
  print('    '+str(znhSize))
  print('    '+str(znhSizePart))
 ####  Dibosons
  print('  ww nonIso')
  wwnh,wwnhSize,wwnhSizePart,wwnhEntries = h.gram(ww_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  wwnh.SetName('wwnh')
  wwnh.Scale(sf_vv)
  print('    '+str(wwnhSize))
  print('    '+str(wwnhSizePart))
  ###
  print('  wz nonIso')
  wznh,wznhSize,wznhSizePart,wznhEntries = h.gram(wz_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  wznh.SetName('wznh')
  wznh.Scale(sf_vv)
  print('    '+str(wznhSize))
  print('    '+str(wznhSizePart))
  ###
  print('  ww nonIso')
  zznh,zznhSize,zznhSizePart,zznhEntries = h.gram(zz_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  zznh.SetName('zznh')
  zznh.Scale(sf_vv)
  print('    '+str(zznhSize))
  print('    '+str(zznhSizePart))
 ####  Single Top
  print('  t t nonIso')
  t_tnh,t_tnhSize,t_tnhSizePart,t_tnhEntries = h.gram(t_t_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  t_tnh.SetName('t_tnh')
  t_tnh.Scale(sf_st_t)
  print('    '+str(t_tnhSize))
  print('    '+str(t_tnhSizePart))
  ###
  print('  tb t nonIso')
  tb_tnh,tb_tnhSize,tb_tnhSizePart,tb_tnhEntries = h.gram(tb_t_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  tb_tnh.SetName('tb_tnh')
  tb_tnh.Scale(sf_st)
  print('    '+str(tb_tnhSize))
  print('    '+str(tb_tnhSizePart))
  ###
  print('  t s nonIso')
  t_snh,t_snhSize,t_snhSizePart,t_snhEntries = h.gram(t_s_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  t_snh.SetName('t_snh')
  t_snh.Scale(sf_st)
  print('    '+str(t_snhSize))
  print('    '+str(t_snhSizePart))
  ###
  print('  tb s nonIso')
  tb_snh,tb_snhSize,tb_snhSizePart,tb_snhEntries = h.gram(tb_s_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  tb_snh.SetName('tb_snh')
  tb_snh.Scale(sf_st)
  print('    '+str(tb_snhSize))
  print('    '+str(tb_snhSizePart))
  ###
  print('  t tw nonIso')
  t_twnh,t_twnhSize,t_twnhSizePart,t_twnhEntries = h.gram(t_tw_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  t_twnh.SetName('t_twnh')
  t_twnh.Scale(sf_st)
  print('    '+str(t_twnhSize))
  print('    '+str(t_twnhSizePart))
  ###
  print('  tb tw nonIso')
  tb_twnh,tb_twnhSize,tb_twnhSizePart,tb_twnhEntries = h.gram(tb_tw_tree,leaf,xmin,xmax,steps,CutsMCn,I,F)
  tb_twnh.SetName('tb_twnh')
  tb_twnh.Scale(sf_st)
  print('    '+str(tb_twnhSize))
  print('    '+str(tb_twnhSizePart))
 ####  TTbar
  print('  ttbar semi nonIso')
  tt_seminh,tt_seminhSize,tt_seminhSizePart,tt_seminhEntries = h.gram(tt_semi_tree,leaf,xmin,xmax,steps,CutsMCnT,I,F)
  tt_seminh.SetName('tt_seminh')
  tt_seminh.Scale(sf_ttbar)
  print('    '+str(tt_seminhSize))
  print('    '+str(tt_seminhSizePart))
  ###
  print('  ttbar full nonIso')
  tt_fullnh,tt_fullnhSize,tt_fullnhSizePart,tt_fullnhEntries = h.gram(tt_full_tree,leaf,xmin,xmax,steps,CutsMCnT,I,F)
  tt_fullnh.SetName('tt_fullnh')
  tt_fullnh.Scale(sf_ttbar)
  print('    '+str(tt_fullnhSize))
  print('    '+str(tt_fullnhSizePart))
 ####  W
  print('  wn nonIso')
  wnnh,wnnhSize,wnnhSizePart,wnhEntries = h.gram(wn_tree,leaf,xmin,xmax,steps,CutsMCnW,I,F)
  wnnh.SetName('wnnh')
  wnnh.Scale(sf_wjets)
  print('    '+str(wnnhSize))
  print('    '+str(wnnhSizePart))
  print('  w1 nonIso')
  ###
  w1nh,w1nhSize,w1nhSizePart,w1nhEntries = h.gram(w1_tree,leaf,xmin,xmax,steps,CutsMCnW,I,F)
  w1nh.SetName('w1nh')
  w1nh.Scale(sf_wjets)
  print('    '+str(w1nhSize))
  print('    '+str(w1nhSizePart))
  print('  w2 nonIso')
  ###
  w2nh,w2nhSize,w2nhSizePart,w2nhEntries = h.gram(w2_tree,leaf,xmin,xmax,steps,CutsMCnW,I,F)
  w2nh.SetName('w2nh')
  w2nh.Scale(sf_wjets)
  print('    '+str(w2nhSize))
  print('    '+str(w2nhSizePart))
  print('  w3 nonIso')
  ###
  w3nh,w3nhSize,w3nhSizePart,w3nhEntries = h.gram(w3_tree,leaf,xmin,xmax,steps,CutsMCnW,I,F)
  w3nh.SetName('w3nh')
  w3nh.Scale(sf_wjets)
  print('    '+str(w3nhSize))
  print('    '+str(w3nhSizePart))
  print('  w4 nonIso')
  ###
  w4nh,w4nhSize,w4nhSizePart,w4nhEntries = h.gram(w4_tree,leaf,xmin,xmax,steps,CutsMCnW,I,F)
  w4nh.SetName('w4nh')
  w4nh.Scale(sf_wjets)
  print('    '+str(w4nhSize))
  print('    '+str(w4nhSizePart))
 #
  print('subtracting non-iso MC from Data ')
  qh = datanh.Clone()
  qh.SetName('qh')
  qh.Add(znh,-1)    #Drell
  qh.Add(t_tnh,-1)  #Single Top
  qh.Add(tb_tnh,-1)
  qh.Add(t_snh,-1)
  qh.Add(tb_snh,-1)
  qh.Add(t_twnh,-1)
  qh.Add(tb_twnh,-1)
  qh.Add(wwnh,-1)   #Diboson
  qh.Add(wznh,-1)
  qh.Add(zznh,-1)
  qh.Add(tt_seminh,-1)  #TTbar
  qh.Add(tt_fullnh,-1)  #TTbar
  qh.Add(wnnh,-1)   #W
  qh.Add(w1nh,-1)
  qh.Add(w2nh,-1)
  qh.Add(w3nh,-1)
  qh.Add(w4nh,-1)
 
  print('qcd')
  qh.SetTitle('')
  bminq = qh.GetXaxis().FindBin(xmin)
  bmaxq = qh.GetXaxis().FindBin(xmax)
  bIq = qh.GetXaxis().FindBin(I)
  bFq = qh.GetXaxis().FindBin(F)
  qhSizePart = qh.Integral(bIq,bFq)
  qh.Scale(sf_qcd)
  qhSize = qh.Integral(bminq,bmaxq)
  qhEntries = qh.GetEntries()
  print('  '+str(qhSize))
  print('  '+str(qhSizePart))
 ####################################
 if drawData:
  ### Data
  print('')
  print('Isolated Variables Now')
  print("data Iso") 
  dataih,dataihSize,dataihSizePart,dataihEntries = h.gram(data_tree,leaf,xmin,xmax,steps,CutsDatai,I,F)
  dataih.SetName('dataih')
  print('  '+str(dataihSize))
  print('  '+str(dataihSizePart))
  dimax = dataih.GetMaximum()
#### Drell
 print("z Iso")
 zih,zihSize,zihSizePart,zihEntries = h.gram(z_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 zih.SetName('zih')
 zih.Scale(sf_drell)
 print('  '+str(zihSize))
 print('  '+str(zihSizePart))
#### Diboson
 print("ww Iso")
 wwih,wwihSize,wwihSizePart,wwihEntries = h.gram(ww_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 wwih.SetName('wwih')
 wwih.Scale(sf_vv)
 print('  '+str(wwihSize))
 print('  '+str(wwihSizePart))
 ###
 print("wz Iso")
 wzih,wzihSize,wzihSizePart,wzihEntries = h.gram(wz_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 wzih.SetName('wzih')
 wzih.Scale(sf_vv)
 print('  '+str(wzihSize))
 print('  '+str(wzihSizePart))
 ###
 print("zz Iso")
 zzih,zzihSize,zzihSizePart,zzihEntries = h.gram(zz_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 zzih.SetName('zzih')
 zzih.Scale(sf_vv)
 print('  '+str(zzihSize))
 print('  '+str(zzihSizePart))
#### Single Top
 print('t t Iso')
 t_tih,t_tihSize,t_tihSizePart,t_tihEntries = h.gram(t_t_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 t_tih.SetName('t_tih')
 t_tih.Scale(sf_st_t)
 print('  '+str(t_tihSize))
 print('  '+str(t_tihSizePart))
 ###
 print('tb t Iso')
 tb_tih,tb_tihSize,tb_tihSizePart,tb_tihEntries = h.gram(tb_t_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 tb_tih.SetName('tb_tih')
 tb_tih.Scale(sf_st)
 print('  '+str(tb_tihSize))
 print('  '+str(tb_tihSizePart))
 ###
 print('t s Iso')
 t_sih,t_sihSize,t_sihSizePart,t_sihEntries = h.gram(t_s_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 t_sih.SetName('t_sih')
 t_sih.Scale(sf_st)
 print('  '+str(t_sihSize))
 print('  '+str(t_sihSizePart))
 ###
 print('tb s Iso')
 tb_sih,tb_sihSize,tb_sihSizePart,tb_sihEntries = h.gram(tb_s_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 tb_sih.SetName('tb_sih')
 tb_sih.Scale(sf_st)
 print('  '+str(tb_sihSize))
 print('  '+str(tb_sihSizePart))
 ###
 print('t tw Iso')
 t_twih,t_twihSize,t_twihSizePart,t_twihEntries = h.gram(t_tw_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 t_twih.SetName('t_twih')
 t_twih.Scale(sf_st)
 print('  '+str(t_twihSize))
 print('  '+str(t_twihSizePart))
 ###
 print('tb tw Iso')
 tb_twih,tb_twihSize,tb_twihSizePart,tb_twihEntries = h.gram(tb_tw_tree,leaf,xmin,xmax,steps,CutsMCi,I,F)
 tb_twih.SetName('tb_twih')
 tb_twih.Scale(sf_st)
 print('  '+str(tb_twihSize))
 print('  '+str(tb_twihSizePart))
#### TTbar
 print('ttb semi Iso')
 tt_semiih,tt_semiihSize,tt_semiihSizePart,tt_semiihEntries = h.gram(tt_semi_tree,leaf,xmin,xmax,steps,CutsMCiT,I,F)
 tt_semiih.SetName('tt_semiih')
 tt_semiih.Scale(sf_ttbar)
 print('  '+str(tt_semiihSize))
 print('  '+str(tt_semiihSizePart))
 ###
 print('ttb full Iso')
 tt_fullih,tt_fullihSize,tt_fullihSizePart,tt_fullihEntries = h.gram(tt_full_tree,leaf,xmin,xmax,steps,CutsMCiT,I,F)
 tt_fullih.SetName('tt_fullih')
 tt_fullih.Scale(sf_ttbar)
 print('  '+str(tt_fullihSize))
 print('  '+str(tt_fullihSizePart))
#### W + Jets
 print('wl n Iso')
 wlnih,wlnihSize,wlnihSizePart,wlnihEntries = h.gram(wn_tree,leaf,xmin,xmax,steps,CutsMCiwl,I,F)
 wlnih.SetName('wlnih')
 wlnih.Scale(sf_wjets)
 print('  '+str(wlnihSize))
 print('  '+str(wlnihSizePart))
 ####
 print('wc n Iso')
 wcnih,wcnihSize,wcnihSizePart,wcnihEntries = h.gram(wn_tree,leaf,xmin,xmax,steps,CutsMCiwc,I,F)
 wcnih.SetName('wcnih')
 wcnih.Scale(sf_wjets)
 print('  '+str(wcnihSize))
 print('  '+str(wcnihSizePart))
 ####
 print('wcc n Iso')
 wccnih,wccnihSize,wccnihSizePart,wccnihEntries = h.gram(wn_tree,leaf,xmin,xmax,steps,CutsMCiwcc,I,F)
 wccnih.SetName('wccnih')
 wccnih.Scale(sf_wjets)
 print('  '+str(wccnihSize))
 print('  '+str(wccnihSizePart))
 ####
 print('wbb n Iso')
 wbbnih,wbbnihSize,wbbnihSizePart,wbbnihEntries = h.gram(wn_tree,leaf,xmin,xmax,steps,CutsMCiwbb,I,F)
 wbbnih.SetName('wbbnih')
 wbbnih.Scale(sf_wjets)
 print('  '+str(wbbnihSize))
 print('  '+str(wbbnihSizePart))
#####
 print('wl 1 Iso')
 wl1ih,wl1ihSize,wl1ihSizePart,wl1ihEntries = h.gram(w1_tree,leaf,xmin,xmax,steps,CutsMCiwl,I,F)
 wl1ih.SetName('wl1ih')
 wl1ih.Scale(sf_wjets)
 print('  '+str(wl1ihSize))
 print('  '+str(wl1ihSizePart))
 ####
 print('wc 1 Iso')
 wc1ih,wc1ihSize,wc1ihSizePart,wc1ihEntries = h.gram(w1_tree,leaf,xmin,xmax,steps,CutsMCiwc,I,F)
 wc1ih.SetName('wc1ih')
 wc1ih.Scale(sf_wjets)
 print('  '+str(wc1ihSize))
 print('  '+str(wc1ihSizePart))
 ####
 print('wcc 1 Iso')
 wcc1ih,wcc1ihSize,wcc1ihSizePart,wcc1ihEntries = h.gram(w1_tree,leaf,xmin,xmax,steps,CutsMCiwcc,I,F)
 wcc1ih.SetName('wcc1ih')
 wcc1ih.Scale(sf_wjets)
 print('  '+str(wcc1ihSize))
 print('  '+str(wcc1ihSizePart))
 ####
 print('wbb 1 Iso')
 wbb1ih,wbb1ihSize,wbb1ihSizePart,wbb1ihEntries = h.gram(w1_tree,leaf,xmin,xmax,steps,CutsMCiwbb,I,F)
 wbb1ih.SetName('wbb1ih')
 wbb1ih.Scale(sf_wjets)
 print('  '+str(wbb1ihSize))
 print('  '+str(wbb1ihSizePart))
#####
 print('wl 2 Iso')
 wl2ih,wl2ihSize,wl2ihSizePart,wl2ihEntries = h.gram(w2_tree,leaf,xmin,xmax,steps,CutsMCiwl,I,F)
 wl2ih.SetName('wl2ih')
 wl2ih.Scale(sf_wjets)
 print('  '+str(wl2ihSize))
 print('  '+str(wl2ihSizePart))
 ####
 print('wc 2 Iso')
 wc2ih,wc2ihSize,wc2ihSizePart,wc2ihEntries = h.gram(w2_tree,leaf,xmin,xmax,steps,CutsMCiwc,I,F)
 wc2ih.SetName('wc2ih')
 wc2ih.Scale(sf_wjets)
 print('  '+str(wc2ihSize))
 print('  '+str(wc2ihSizePart))
 ####
 print('wcc 2 Iso')
 wcc2ih,wcc2ihSize,wcc2ihSizePart,wcc2ihEntries = h.gram(w2_tree,leaf,xmin,xmax,steps,CutsMCiwcc,I,F)
 wcc2ih.SetName('wcc2ih')
 wcc2ih.Scale(sf_wjets)
 print('  '+str(wcc2ihSize))
 print('  '+str(wcc2ihSizePart))
 ####
 print('wbb 2 Iso')
 wbb2ih,wbb2ihSize,wbb2ihSizePart,wbb2ihEntries = h.gram(w2_tree,leaf,xmin,xmax,steps,CutsMCiwbb,I,F)
 wbb2ih.SetName('wbb2ih')
 wbb2ih.Scale(sf_wjets)
 print('  '+str(wbb2ihSize))
 print('  '+str(wbb2ihSizePart))
#####
 print('wl 3 Iso')
 wl3ih,wl3ihSize,wl3ihSizePart,wl3ihEntries = h.gram(w3_tree,leaf,xmin,xmax,steps,CutsMCiwl,I,F)
 wl3ih.SetName('wl3ih')
 wl3ih.Scale(sf_wjets)
 print('  '+str(wl3ihSize))
 print('  '+str(wl3ihSizePart))
 ####
 print('wc 3 Iso')
 wc3ih,wc3ihSize,wc3ihSizePart,wc3ihEntries = h.gram(w3_tree,leaf,xmin,xmax,steps,CutsMCiwc,I,F)
 wc3ih.SetName('wc3ih')
 wc3ih.Scale(sf_wjets)
 print('  '+str(wc3ihSize))
 print('  '+str(wc3ihSizePart))
 ####
 print('wcc 3 Iso')
 wcc3ih,wcc3ihSize,wcc3ihSizePart,wcc3ihEntries = h.gram(w3_tree,leaf,xmin,xmax,steps,CutsMCiwcc,I,F)
 wcc3ih.SetName('wcc3ih')
 wcc3ih.Scale(sf_wjets)
 print('  '+str(wcc3ihSize))
 print('  '+str(wcc3ihSizePart))
 ####
 print('wbb 3 Iso')
 wbb3ih,wbb3ihSize,wbb3ihSizePart,wbb3ihEntries = h.gram(w3_tree,leaf,xmin,xmax,steps,CutsMCiwbb,I,F)
 wbb3ih.SetName('wbb3ih')
 wbb3ih.Scale(sf_wjets)
 print('  '+str(wbb3ihSize))
 print('  '+str(wbb3ihSizePart))
#####
 print('wl 4 Iso')
 wl4ih,wl4ihSize,wl4ihSizePart,wl4ihEntries = h.gram(w4_tree,leaf,xmin,xmax,steps,CutsMCiwl,I,F)
 wl4ih.SetName('wl4ih')
 wl4ih.Scale(sf_wjets)
 print('  '+str(wl4ihSize))
 print('  '+str(wl4ihSizePart))
 ####
 print('wc 4 Iso')
 wc4ih,wc4ihSize,wc4ihSizePart,wc4ihEntries = h.gram(w4_tree,leaf,xmin,xmax,steps,CutsMCiwc,I,F)
 wc4ih.SetName('wc4ih')
 wc4ih.Scale(sf_wjets)
 print('  '+str(wc4ihSize))
 print('  '+str(wc4ihSizePart))
 ####
 print('wcc 4 Iso')
 wcc4ih,wcc4ihSize,wcc4ihSizePart,wcc4ihEntries = h.gram(w4_tree,leaf,xmin,xmax,steps,CutsMCiwcc,I,F)
 wcc4ih.SetName('wcc4ih')
 wcc4ih.Scale(sf_wjets)
 print('  '+str(wcc4ihSize))
 print('  '+str(wcc4ihSizePart))
 ####
 print('wbb 4 Iso')
 wbb4ih,wbb4ihSize,wbb4ihSizePart,wbb4ihEntries = h.gram(w4_tree,leaf,xmin,xmax,steps,CutsMCiwbb,I,F)
 wbb4ih.SetName('wbb4ih')
 wbb4ih.Scale(sf_wjets)
 print('  '+str(wbb4ihSize))
 print('  '+str(wbb4ihSizePart))
 
 outFile.Write()
 print('')
 print('Your File is here: '+path+extraName+'_'+leaf+'.root')
 print('')
 
 qhSize = max(qhSize,0.000001)
 qhSizePart = max(qhSizePart,0.000001)
 mcSize = \
  wlnihSize+wl1ihSize+wl2ihSize+wl3ihSize+wl4ihSize \
  +wcnihSize+wc1ihSize+wc2ihSize+wc3ihSize+wc4ihSize \
  +wccnihSize+wcc1ihSize+wcc2ihSize+wcc3ihSize+wcc4ihSize \
  +wbbnihSize+wbb1ihSize+wbb2ihSize+wbb3ihSize+wbb4ihSize \
  +tt_semiihSize+tt_fullihSize \
  +t_tihSize+tb_tihSize \
  +t_sihSize+tb_sihSize \
  +t_twihSize+tb_twihSize \
  +wwihSize+wzihSize+zzihSize \
  +zihSize 
 if drawQCD: mcSize = mcSize+qhSize
  #+qhSize

 mcSizePart = \
  wlnihSizePart+wl1ihSizePart+wl2ihSizePart+wl3ihSizePart+wl4ihSizePart+ \
  wcnihSizePart+wc1ihSizePart+wc2ihSizePart+wc3ihSizePart+wc4ihSizePart+ \
  wccnihSizePart+wcc1ihSizePart+wcc2ihSizePart+wcc3ihSizePart+wcc4ihSizePart+ \
  wbbnihSizePart+wbb1ihSizePart+wbb2ihSizePart+wbb3ihSizePart+wbb4ihSizePart+ \
  tt_semiihSizePart+tt_fullihSizePart+ \
  t_tihSizePart+tb_tihSizePart+ \
  t_sihSizePart+tb_sihSizePart+ \
  t_twihSizePart+tb_twihSizePart+ \
  wwihSizePart+wzihSizePart+zzihSizePart+ \
  zihSizePart 
 if drawQCD: mcSizePart = mcSizePart + qhSizePart
  #qhSizePart

 log.write('------------------------------------------------\n')
 if drawQCD:
  if leaf=='mt':
   log.write('You Probably Want to use this for QCD Scale\n')
   if qhSizePart!=0: SF = 1 + (dataihSizePart - mcSizePart)/qhSizePart
   else: SF = 1
   log.write('1 + (data - mc)/qcd = '+str(SF)+'\n')
  if leaf=='mt_new':
   log.write('You Probably Want to use this for QCD Scale\n')
   if qhSizePart!=0: SF = 1 + (dataihSizePart - mcSizePart)/qhSizePart
   else: SF = 1
   log.write('1 + (data - mc)/qcd = '+str(SF)+'\n')
 log.write('---------------------------\n')
 log.write('---------------------------\n')
 log.write('Non Isolated\n')
 log.write('---------------------------\n')
 log.write(' Cuts MC:      \n'+str(CutsMCn)+'\n\n')
 log.write(' Cuts MC Top:  \n'+str(CutsMCnT)+'\n\n')
 log.write(' Cuts Data:    \n'+str(CutsDatan)+'\n\n')
 log.write('Isolated\n')
 log.write('---------------------------\n')
 log.write(' Cuts MC:  \n'+str(CutsMCi)+'\n\n')
 log.write(' Cuts MC Top:  \n'+str(CutsMCiT)+'\n\n')
 log.write(' Cuts Data: \n'+str(CutsDatai)+'\n\n')
 log.write(' Cuts Wbb: \n'+str(CutsMCiwbb)+'\n\n')
 log.write(' Cuts Wcc: \n'+str(CutsMCiwcc)+'\n\n')
 log.write(' Cuts Wc: \n'+str(CutsMCiwc)+'\n\n')
 log.write(' Cuts Wl: \n'+str(CutsMCiwl)+'\n\n')
 if drawQCD:
  log.write('Anti-Isolated Sizes\n')
  log.write('---------------------------\n')
  log.write(' Wn Size:--------'+str(wnnhSize)+'\n')
  log.write(' W1 Size:--------'+str(w1nhSize)+'\n')
  log.write(' W2 Size:--------'+str(w2nhSize)+'\n')
  log.write(' W3 Size:--------'+str(w3nhSize)+'\n')
  log.write(' W4 Size:--------'+str(w4nhSize)+'\n')
  log.write(' tt semi Size:---'+str(tt_seminhSize)+'\n')
  log.write(' tt full Size:---'+str(tt_fullnhSize)+'\n')
  log.write(' t_t Size:-------'+str(t_tnhSize)+'\n')
  log.write(' tb_t Size:------'+str(tb_tnhSize)+'\n')
  log.write(' t_s Size:-------'+str(t_snhSize)+'\n')
  log.write(' tb_s Size:------'+str(tb_snhSize)+'\n')
  log.write(' t_tw Size:------'+str(t_twnhSize)+'\n')
  log.write(' tb_tw Size:-----'+str(tb_twnhSize)+'\n')
  log.write(' WW Size:--------'+str(wwnhSize)+'\n')
  log.write(' WZ Size:--------'+str(wznhSize)+'\n')
  log.write(' ZZ Size:--------'+str(zznhSize)+'\n')
  log.write(' Drell-Yan Size:-'+str(znhSize)+'\n')
  log.write(' QCD Size:-------'+str(qhSize)+'\n')
  log.write(' Data Size:      '+str(datanhSize)+'\n')
  log.write('---------------------------\n')
  log.write('Anti-Isolated Sizes from '+str(I)+' to '+str(F)+'\n')
  log.write('---------------------------\n')
  log.write(' Wn Size:--------'+str(wnnhSizePart)+'\n')
  log.write(' W1 Size:--------'+str(w1nhSizePart)+'\n')
  log.write(' W2 Size:--------'+str(w2nhSizePart)+'\n')
  log.write(' W3 Size:--------'+str(w3nhSizePart)+'\n')
  log.write(' W4 Size:--------'+str(w4nhSizePart)+'\n')
  log.write(' tt semi Size:---'+str(tt_seminhSizePart)+'\n')
  log.write(' tt full Size:---'+str(tt_fullnhSizePart)+'\n')
  log.write(' t_t Size:-------'+str(t_tnhSizePart)+'\n')
  log.write(' tb_t Size:------'+str(tb_tnhSizePart)+'\n')
  log.write(' t_s Size:-------'+str(t_snhSizePart)+'\n')
  log.write(' tb_s Size:------'+str(tb_snhSizePart)+'\n')
  log.write(' t_tw Size:------'+str(t_twnhSizePart)+'\n')
  log.write(' tb_tw Size:-----'+str(tb_twnhSizePart)+'\n')
  log.write(' WW Size:--------'+str(wwnhSizePart)+'\n')
  log.write(' WZ Size:--------'+str(wznhSizePart)+'\n')
  log.write(' ZZ Size:--------'+str(zznhSizePart)+'\n')
  log.write(' Drell-Yan Size:-'+str(znhSizePart)+'\n')
  log.write(' QCD Size:-------'+str(qhSizePart)+'\n')
  log.write(' Data Size:      '+str(datanhSizePart)+'\n')
  log.write('---------------------------------------------\n')
  log.write('QCD Scale:   '+str(sf_qcd)+'\n')
  log.write('Drell Scale: '+str(sf_drell)+'\n')
  log.write('t Scale:     '+str(sf_st)+'\n')
  log.write('t_t Scale:   '+str(sf_st_t)+'\n')
  log.write('tt Scale:    '+str(sf_ttbar)+'\n')
  log.write('W+x Scale:   '+str(sf_wjets)+'\n')
  log.write('---------------------------------------------\n\n')
 else:
  log.write('no qcd drawn\n\n')
 log.write('Isolated Sizes from '+str(I)+' to '+str(F)+'\n')
 log.write('---------------------------\n')
 log.write('W:\n')
 log.write('  W+light:--------'+str(wlnihSizePart)+'\n')
 log.write('  W+light:--------'+str(wl1ihSizePart)+'\n')
 log.write('  W+light:--------'+str(wl2ihSizePart)+'\n')
 log.write('  W+light:--------'+str(wl3ihSizePart)+'\n')
 log.write('  W+light:--------'+str(wl4ihSizePart)+'\n\n')
 log.write('  W+c:------------'+str(wcnihSizePart)+'\n')
 log.write('  W+c:------------'+str(wc1ihSizePart)+'\n')
 log.write('  W+c:------------'+str(wc2ihSizePart)+'\n')
 log.write('  W+c:------------'+str(wc3ihSizePart)+'\n')
 log.write('  W+c:------------'+str(wc4ihSizePart)+'\n\n')
 log.write('  W+cc:-----------'+str(wccnihSizePart)+'\n')
 log.write('  W+cc:-----------'+str(wcc1ihSizePart)+'\n')
 log.write('  W+cc:-----------'+str(wcc2ihSizePart)+'\n')
 log.write('  W+cc:-----------'+str(wcc3ihSizePart)+'\n')
 log.write('  W+cc:-----------'+str(wcc4ihSizePart)+'\n\n')
 log.write('  W+bb:-----------'+str(wbbnihSizePart)+'\n')
 log.write('  W+bb:-----------'+str(wbb1ihSizePart)+'\n')
 log.write('  W+bb:-----------'+str(wbb2ihSizePart)+'\n')
 log.write('  W+bb:-----------'+str(wbb3ihSizePart)+'\n')
 log.write('  W+bb:-----------'+str(wbb4ihSizePart)+'\n')
 log.write('tt semi Size:-----'+str(tt_semiihSizePart)+'\n')
 log.write('tt full Size:-----'+str(tt_fullihSizePart)+'\n')
 log.write('t_t Size:-------'+str(t_tihSizePart)+'\n')
 log.write('tb_t Size:------'+str(tb_tihSizePart)+'\n')
 log.write('t_s Size:-------'+str(t_sihSizePart)+'\n')
 log.write('tb_s Size:------'+str(tb_sihSizePart)+'\n')
 log.write('t_tw Size:------'+str(t_twihSizePart)+'\n')
 log.write('tb_tw Size:-----'+str(tb_twihSizePart)+'\n')
 log.write('WW Size:--------'+str(wwihSizePart)+'\n')
 log.write('WZ Size:--------'+str(wzihSizePart)+'\n')
 log.write('ZZ Size:--------'+str(zzihSizePart)+'\n')
 log.write('Drell-Yan Size:-'+str(zihSizePart)+'\n')
 if drawData:
  log.write('Data Size:          '+str(dataihSizePart)+'\n')
 else:
  log.write('no Data drawn\n')
 log.write('---------------------------------------------\n')
 log.write('Isolated Sizes\n')
 log.write('---------------------------\n')
 log.write('W:\n')
 log.write('  Wn+light:--------'+str(wlnihSize)+'\n')
 log.write('  W1+light:--------'+str(wl1ihSize)+'\n')
 log.write('  W2+light:--------'+str(wl2ihSize)+'\n')
 log.write('  W3+light:--------'+str(wl3ihSize)+'\n')
 log.write('  W4+light:--------'+str(wl4ihSize)+'\n\n')
 log.write('  Wn+c:------------'+str(wcnihSize)+'\n')
 log.write('  W1+c:------------'+str(wc1ihSize)+'\n')
 log.write('  W2+c:------------'+str(wc2ihSize)+'\n')
 log.write('  W3+c:------------'+str(wc3ihSize)+'\n')
 log.write('  W4+c:------------'+str(wc4ihSize)+'\n\n')
 log.write('  Wn+cc:-----------'+str(wccnihSize)+'\n')
 log.write('  W1+cc:-----------'+str(wcc1ihSize)+'\n')
 log.write('  W2+cc:-----------'+str(wcc2ihSize)+'\n')
 log.write('  W3+cc:-----------'+str(wcc3ihSize)+'\n')
 log.write('  W4+cc:-----------'+str(wcc4ihSize)+'\n\n')
 log.write('  Wn+bb:-----------'+str(wbbnihSize)+'\n')
 log.write('  W1+bb:-----------'+str(wbb1ihSize)+'\n')
 log.write('  W2+bb:-----------'+str(wbb2ihSize)+'\n')
 log.write('  W3+bb:-----------'+str(wbb3ihSize)+'\n')
 log.write('  W4+bb:-----------'+str(wbb4ihSize)+'\n')
 log.write('tt semi Size:---'+str(tt_semiihSize)+'\n')
 log.write('tt full Size:---'+str(tt_fullihSize)+'\n')
 log.write('t_t Size:-------'+str(t_tihSize)+'\n')
 log.write('tb_t Size:------'+str(tb_tihSize)+'\n')
 log.write('t_s Size:-------'+str(t_sihSize)+'\n')
 log.write('tb_s Size:------'+str(tb_sihSize)+'\n')
 log.write('t_tw Size:------'+str(t_twihSize)+'\n')
 log.write('tb_tw Size:-----'+str(tb_twihSize)+'\n')
 log.write('WW Size:--------'+str(wwihSize)+'\n')
 log.write('WZ Size:--------'+str(wzihSize)+'\n')
 log.write('ZZ Size:--------'+str(zzihSize)+'\n')
 log.write('Drell-Yan Size:-'+str(zihSize)+'\n')
 log.write('Total MC Size:      '+str(mcSize)+'\n')
 if drawData:
  log.write('Data Size:          '+str(dataihSize)+'\n')
 log.write('---------------------------------------------\n')
 log.write('---------------------------\n')
 log.write('Isolated Number of Entries\n')
 log.write('---------------------------\n')
 log.write('W:\n')
 log.write('  Wn+light:--------'+str(wlnihEntries)+'\n')
 log.write('  W1+light:--------'+str(wl1ihEntries)+'\n')
 log.write('  W2+light:--------'+str(wl2ihEntries)+'\n')
 log.write('  W3+light:--------'+str(wl3ihEntries)+'\n')
 log.write('  W4+light:--------'+str(wl4ihEntries)+'\n\n')
 log.write('  Wn+c:------------'+str(wcnihEntries)+'\n')
 log.write('  W1+c:------------'+str(wc1ihEntries)+'\n')
 log.write('  W2+c:------------'+str(wc2ihEntries)+'\n')
 log.write('  W3+c:------------'+str(wc3ihEntries)+'\n')
 log.write('  W4+c:------------'+str(wc4ihEntries)+'\n\n')
 log.write('  Wn+cc:-----------'+str(wccnihEntries)+'\n')
 log.write('  W1+cc:-----------'+str(wcc1ihEntries)+'\n')
 log.write('  W2+cc:-----------'+str(wcc2ihEntries)+'\n')
 log.write('  W3+cc:-----------'+str(wcc3ihEntries)+'\n')
 log.write('  W4+cc:-----------'+str(wcc4ihEntries)+'\n\n')
 log.write('  Wn+bb:-----------'+str(wbbnihEntries)+'\n')
 log.write('  W1+bb:-----------'+str(wbb1ihEntries)+'\n')
 log.write('  W2+bb:-----------'+str(wbb2ihEntries)+'\n')
 log.write('  W3+bb:-----------'+str(wbb3ihEntries)+'\n')
 log.write('  W4+bb:-----------'+str(wbb4ihEntries)+'\n')
 log.write('tt semi Entries:---'+str(tt_semiihEntries)+'\n')
 log.write('tt full Entries:---'+str(tt_fullihEntries)+'\n')
 log.write('t_t Entries:-------'+str(t_tihEntries)+'\n')
 log.write('tb_t Entries:------'+str(tb_tihEntries)+'\n')
 log.write('t_s Entries:-------'+str(t_sihEntries)+'\n')
 log.write('tb_s Entries:------'+str(tb_sihEntries)+'\n')
 log.write('t_tw Entries:------'+str(t_twihEntries)+'\n')
 log.write('tb_tw Entries:-----'+str(tb_twihEntries)+'\n')
 log.write('WW Entries:--------'+str(wwihEntries)+'\n')
 log.write('WZ Entries:--------'+str(wzihEntries)+'\n')
 log.write('ZZ Entries:--------'+str(zzihEntries)+'\n')
 log.write('Drell-Yan Entries:-'+str(zihEntries)+'\n')
 if drawData:
  log.write('Data Entries:          '+str(dataihEntries)+'\n')
 log.write('---------------------------------------------\n')
 
 print(".log written")
 log.close()
