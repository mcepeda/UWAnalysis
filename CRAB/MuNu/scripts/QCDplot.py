#!/usr/bin/env python
'''
Makes QCD plots with ratio from histograms made using makeHistos.py
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle
import os
import histoRange as hr
import cmsPrelim as cpr
import TheParameters as p
 
sf_qcd = 1

sf_Signal_Wbb    = 1. 
sf_Signal_Tbar   = 1. 
sf_Signal_tW     = 1. 
sf_Signal_TOP    = 1. 
sf_Signal_Z_jets = 1. 
sf_Signal_Wcc    = 1. 
sf_Signal_T      = 1. 
sf_Signal_WZ     = 1. 

sf_Top_T         = 1.#  (2.26033000000000000e+05)/2497226 
sf_Top_TOP       = 1. 
sf_Top_tW        = 1. 
sf_Top_Wbb       = 1. 
sf_Top_Tbar      = 1. 

rebin = 1
#canvas attributes
canx = 800 # for one plot on page
#canx = 550 # for two plots on page with text
#canx = 600 # for two plots on page just title
#canx = 500 # for three plots on page with text
#canx = 400 # for three plots on page with just title

cany = 900

#color scheme
d = 1

q =   ROOT.EColor.kRed+1
z =   ROOT.EColor.kOrange-3
d =   ROOT.EColor.kYellow-3
tt =  ROOT.EColor.kGreen+1
ts =  ROOT.EColor.kGreen-5
ttw = ROOT.EColor.kGreen+3 
ttb = ROOT.EColor.kGreen-9
wl =  ROOT.EColor.kAzure+10
wc =  ROOT.EColor.kBlue+1
wcc = ROOT.EColor.kAzure+2
wbb = 51#ROOT.EColor.kCyan

#7 TeV color scheme
#q =   ROOT.EColor.kGreen-3   
#z =   ROOT.EColor.kRed+2     
#d =   ROOT.EColor.kRed+3     
#tt =  ROOT.EColor.kOrange+7  
#ts =  ROOT.EColor.kOrange+6  
#ttw = ROOT.EColor.kOrange+5  
#ttb = ROOT.EColor.kOrange-2  
#wl =  ROOT.EColor.kAzure+1   
#wc =  ROOT.EColor.kMagenta-4 
#wcc = ROOT.EColor.kViolet-6  
#wbb = ROOT.EColor.kBlue-10   

tex = ROOT.TLatex()
tex.SetTextSize(0.07)
tex.SetTextAlign(13)
tex.SetNDC(True)
gStyle.SetOptStat('')

# get parameters to run on
lumi,bNr,btype,jNr,njetcut,jetcut,I,F,iso_value,antiIso_value,path,extraName,leafs,drawW,drawZ,drawQCD,drawData,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,eventTreeLocation,extraCut,Ctl_Andrea = p.arams() 

sf_Wbb    = 1. 
sf_Tbar   = 1. 
sf_tW     = 1. 
sf_TOP    = 1. 
sf_Z_jets = 1. 
sf_Wcc    = 1. 
sf_T      = 1. 
sf_WZ     = 1. 
sf_wl     = 1.
sf_wc     = 1.
#if Signal:
# sf_Wbb    = sf_Signal_Wbb   
# sf_Tbar   = sf_Signal_Tbar  
# sf_tW     = sf_Signal_tW    
# sf_TOP    = sf_Signal_TOP   
# sf_Z_jets = sf_Signal_Z_jets
# sf_Wcc    = sf_Signal_Wcc   
# sf_T      = sf_Signal_T     
# sf_WZ     = sf_Signal_WZ    
#if TT_m or TT_me:
# sf_Wbb    = sf_Top_Wbb  
# sf_Tbar   = sf_Top_Tbar 
# sf_tW     = sf_Top_tW   
# sf_TOP    = sf_Top_TOP  
# sf_T      = sf_Top_T    

for leaf in leafs:

 steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)
 
 if xunits is not None:
  xlabel = xtitle+' ['+xunits+']'
 else:
  xlabel = xtitle
 ylabel = 'Events/ %.0001f' %(float((xmax-xmin)*rebin)/(steps))
 title = xtitle #+' Data v MC'
 
 name = []
 name.append(extraName+'_'+leaf)
 
 for i in name:

  log = open(path+i+'.log','a')
  log.write('\n\nOn Plotting SF QCD: '+str(sf_qcd)+'\n') 
  log.write('On Plotting SF Drell: '+str(sf_Z_jets)+'\n')
  log.write('On Plotting SF Top: '+str(sf_T)+'\n')
  log.write('On Plotting SF tBar: '+str(sf_Tbar)+'\n')
  log.write('On Plotting SF tW: '+str(sf_tW)+'\n')
  log.write('On Plotting SF TTbar: '+str(sf_TOP)+'\n')
  log.write('On Plotting SF Wbb: '+str(sf_Wbb)+'\n')
  log.write('On Plotting SF Wcc: '+str(sf_Wcc)+'\n')
  log.write('On Plotting SF Diboson: '+str(sf_WZ)+'\n')
  log.write('----------------------------------------\n')
  log.write('----------------------------------------\n')
  log.close()

  c = TCanvas('c','Canvas Named c',canx,cany)
  version = os.environ.get('version')
  path = '../plots/%s_'%(version)
  print(path+i)
  theFile = TFile(path+i+'.root')
  
  datanh = theFile.Get('datanh')
  datanh.Rebin(rebin)
  datanh.SetMarkerStyle(22)
  datanh.SetMarkerSize(1.2)
  dmax = datanh.GetMaximum()
  
 #### Drell
  znh = theFile.Get('znh')
  znh.SetFillColor(z)
  znh.Rebin(rebin)
  znh.Scale(sf_Z_jets)
  znh.Draw()
 #### Diboson
  wwnh = theFile.Get('wwnh')
  wwnh.SetFillColor(d)
  wwnh.Rebin(rebin)
  wwnh.Scale(sf_WZ)
  wwnh.Draw()
  ###
  wznh = theFile.Get('wznh')
  wznh.SetFillColor(d)
  wznh.Rebin(rebin)
  wznh.Scale(sf_WZ)
  wznh.Draw()
  ###
  zznh = theFile.Get('zznh')
  zznh.SetFillColor(d)
  zznh.Rebin(rebin)
  zznh.Scale(sf_WZ)
  zznh.Draw()
  ####
  dnh = wwnh.Clone()
  dnh.SetName('dnh')
  dnh.Add(wznh)
  dnh.Add(zznh)
  dnh.Draw()
 #### Single Top
  # s
  t_snh = theFile.Get('t_snh')
  t_snh.SetFillColor(ts)
  t_snh.Rebin(rebin)
  t_snh.Scale(sf_T)
  t_snh.Draw()
  ###
  tb_snh = theFile.Get('tb_snh')
  tb_snh.SetFillColor(ts)
  tb_snh.Rebin(rebin)
  tb_snh.Scale(sf_T)
  tb_snh.Draw()
  #### 
  stsnh = t_snh.Clone()
  stsnh.SetName('stsnh')
  stsnh.Add(tb_snh)
  stsnh.Draw()
  ###
  # t
  t_tnh = theFile.Get('t_tnh')
  t_tnh.SetFillColor(tt)
  t_tnh.Rebin(rebin)
  t_tnh.Scale(sf_T)
  t_tnh.Draw()
  ###
  tb_tnh = theFile.Get('tb_tnh')
  tb_tnh.SetFillColor(tt)
  tb_tnh.Rebin(rebin)
  tb_tnh.Scale(sf_T)
  tb_tnh.Draw()
  #### 
  sttnh = t_tnh.Clone()
  sttnh.SetName('sttnh')
  sttnh.Add(tb_tnh)
  sttnh.Draw()
  ###
  # tw
  t_twnh = theFile.Get('t_twnh')
  t_twnh.SetFillColor(ttw)
  t_twnh.Rebin(rebin)
  t_twnh.Scale(sf_tW)
  t_twnh.Draw()
  ###
  tb_twnh = theFile.Get('tb_twnh')
  tb_twnh.SetFillColor(ttw)
  tb_twnh.Rebin(rebin)
  tb_twnh.Scale(sf_tW)
  tb_twnh.Draw()
  #### 
  sttwnh = t_twnh.Clone()
  sttwnh.SetName('sttwnh')
  sttwnh.Add(tb_twnh)
  sttwnh.Draw()
 #### TTbar
  tt_seminh = theFile.Get('tt_seminh')
  tt_seminh.SetFillColor(ttb)
  tt_seminh.Rebin(rebin)
  tt_seminh.Scale(sf_TOP)
  tt_seminh.Draw()
  ###
  tt_fullnh = theFile.Get('tt_fullnh')
  tt_fullnh.SetFillColor(ttb)
  tt_fullnh.Rebin(rebin)
  tt_fullnh.Scale(sf_TOP)
  tt_fullnh.Draw()
  ###
  ttbnh = tt_seminh.Clone()
  ttbnh.SetName('ttbnh')
  ttbnh.Add(tt_fullnh)
  ttbnh.Draw()

  bmin = ttbnh.GetXaxis().FindBin(xmin)
  bmax = ttbnh.GetXaxis().FindBin(xmax)
  print(ttbnh.Integral(bmin,bmax))
 #### W 
  w1nh = theFile.Get('w1nh')
  w1nh.SetFillColor(wbb)
  w1nh.Rebin(rebin)
  print(w1nh.Integral(bmin,bmax))
  ###
  w2nh = theFile.Get('w2nh')
  w2nh.SetFillColor(wbb)
  w2nh.Rebin(rebin)
  print(w2nh.Integral(bmin,bmax))
  ###
  w3nh = theFile.Get('w3nh')
  w3nh.SetFillColor(wbb)
  w3nh.Rebin(rebin)
  print(w3nh.Integral(bmin,bmax))
  ###
  w4nh = theFile.Get('w4nh')
  w4nh.SetFillColor(wbb)
  w4nh.Rebin(rebin)
  print(w4nh.Integral(bmin,bmax))
  ###
  wnnh = theFile.Get('wnnh')
  wnnh.SetFillColor(wbb)
  wnnh.Rebin(rebin)
  print(wnnh.Integral(bmin,bmax))
  ###
  wjnh = w1nh.Clone()
  wjnh.SetName("wjnh")
  wjnh.Add(w2nh)
  wjnh.Add(w3nh)
  wjnh.Add(w4nh)
  wjnh.Add(wnnh)
 
  hs = THStack('hs','')
  hs.SetTitle('')
  hs.Add(znh)
  hs.Add(dnh)
  hs.Add(sttnh)
  hs.Add(stsnh)
  hs.Add(sttwnh)
  hs.Add(ttbnh)
  hs.Add(wjnh)

  #### QCD
  qh = theFile.Get('qh')
  qh.SetLineColor(q)
  qh.SetLineStyle(2)
  qh.Rebin(rebin)
  qh.Scale(sf_qcd)
  qh.Draw()

  leg=TLegend(0.68,0.2,0.9,0.8)
  leg.AddEntry(datanh,'data')
  leg.AddEntry(wjnh,'W+Jets','f')
  leg.AddEntry(ttbnh,'t#bar{t}','f')
  leg.AddEntry(sttwnh,'t_tW','f')
  leg.AddEntry(stsnh,'t_s','f')
  leg.AddEntry(sttnh,'t_t','f')
  leg.AddEntry(dnh,'WW,WZ,ZZ','f')
  leg.AddEntry(znh,'Drell-Yan','f')
  leg.AddEntry(qh,'QCD','f')
  leg.SetFillColor(0)
  leg.SetBorderSize(0)
  
  hs.SetMaximum(1.2*dmax)
  c.cd()
  hs.Draw('hist')
  datanh.Draw('sames,E1')
  qh.Draw('sames,hist')
  leg.Draw('sames')
  cpr.prelim_alt(lumi)
  tex.SetTextAlign(11)#left, bottom
  tex.DrawLatex(0.17,0.9,title)
 #####
  c.cd()
  print('you just read '+leaf)
  #save2 = raw_input ('Press Enter to Continue (type save to save)\n')
  #if save2 == 'save':
  c.Print(path+"QCD"+i+'.png')
   #c.Print(path+i+'_unscaled.png')
  print('')
  c.Close()
