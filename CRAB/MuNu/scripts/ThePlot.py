#!/usr/bin/env python
'''
Makes nice plots with ratio from histograms made using makeHistos.py
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle
import histoRange as hr
import cmsPrelim as cpr
import TheParameters as p
 
# scale factors : sf_qcd = 1 + (data-allMC)/qcd in 0<mt<20

# Halloween
#sf_qcd = 0.00000000001
#sf_qcd = 0.0025316635942 # Wbb
sf_qcd = 0.00302154299378 # TT_3j
#sf_qcd = 1. # TT_1m1e
#sf_qcd = 0.0519607675029 # ST
#sf_qcd = 0.0253012710684 # WJets
#sf_qcd = 0.277049853807 # WJets + QCDwTRG

sf_Top_Wbb          = 1. #1.9553517118    / 1.57266658081
sf_Top_TTbar        = 1. #240.6/245.8 #1. #1366.90108757   / 1225.41492591
sf_Top_tW           = 1. #33.7173900878   / 34.2463617338
sf_Top_Tbar         = 1. #0.663758928599  / 0.6598701129
sf_Top_T            = 1. #0.957306964662  / 0.950968565694
sf_Top_QCD = sf_qcd * 1.#0.0724173672497 / 0.102551698685

sf_Signal_QCD = sf_qcd * 1. #(276.682527912 / 391.815720025) * sf_qcd
sf_Signal_Drell        = 1. #110.595743425 / 100.460683177
sf_Signal_VV           = 1. #115.84464714  / 114.132847588
sf_Signal_T            = 1. #415.707530564 / 389.801262543
sf_Signal_Tbar         = 1. #241.137591009 / 222.805572482
sf_Signal_tW           = 1. #89.1214765299 / 86.8481257456
sf_Signal_TTbar        = 1. #2738.7690035  / 2302.46459317
sf_Signal_Wl           = 1. #52.9669448168 / 74.6321442293
sf_Signal_Wcc          = 1. #130.726683911 / 144.534251764
sf_Signal_Wbb          = 1. #1160.15559794 / 946.523539255  
#
#sf_Top_Wbb       = 1.9553517118    / 1.57266658081
#sf_Top_TTbar     = 1366.90108757   / 1225.41492591
#sf_Top_tW        = 33.7173900878   / 34.2463617338
#sf_Top_Tbar      = 0.663758928599  / 0.6598701129
#sf_Top_T         = 0.957306964662  / 0.950968565694
#sf_Top_QCD       = 0.0724173672497 / 0.102551698685

ratioRange = 0.3
rebin = 1
errorBand = True
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

sf_QCD    = sf_qcd 
sf_Drell  = 1. 
sf_VV     = 1. 
sf_T      = 1. 
sf_Tbar   = 1. 
sf_tW     = 1. 
sf_TTbar  = 1. 
sf_Wl     = 1. 
sf_Wcc    = 1. 
sf_Wbb    = 1. 
if Signal:
 sf_QCD    = sf_Signal_QCD
 sf_Drell  = sf_Signal_Drell
 sf_VV     = sf_Signal_VV    
 sf_T      = sf_Signal_T     
 sf_Tbar   = sf_Signal_Tbar  
 sf_tW     = sf_Signal_tW    
 sf_TTbar  = sf_Signal_TTbar 
 sf_Wl     = sf_Signal_Wl
 sf_Wcc    = sf_Signal_Wcc   
 sf_Wbb    = sf_Signal_Wbb   
if TT_m or TT_me:
 sf_QCD    = sf_Top_QCD
 sf_Wbb    = sf_Top_Wbb  
 sf_Tbar   = sf_Top_Tbar 
 sf_tW     = sf_Top_tW   
 sf_TTbar  = sf_Top_TTbar
 sf_T      = sf_Top_T

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
  log.write('On Plotting SF Drell: '+str(sf_Drell)+'\n')
  log.write('On Plotting SF Top: '+str(sf_T)+'\n')
  log.write('On Plotting SF tBar: '+str(sf_Tbar)+'\n')
  log.write('On Plotting SF tW: '+str(sf_tW)+'\n')
  log.write('On Plotting SF TTbar: '+str(sf_TTbar)+'\n')
  log.write('On Plotting SF Wbb: '+str(sf_Wbb)+'\n')
  log.write('On Plotting SF Wcc: '+str(sf_Wcc)+'\n')
  log.write('On Plotting SF Diboson: '+str(sf_VV)+'\n')
  log.write('----------------------------------------\n')
  log.write('----------------------------------------\n')
  log.close()

  c = TCanvas('c','Canvas Named c',canx,cany)
  c.cd()
  p1 = TPad('p1','p1',0,0.3,1,1)
  p1.SetBottomMargin(0.08)
  p1.Draw()
  p1.cd()
  p1.SetLogy(setLogY)
  theFile = TFile(path+i+'.root')
  
  if drawData:
   dataih = theFile.Get('dataih')
   dataih.Rebin(rebin)
   dataih.SetMarkerStyle(22)
   dataih.SetMarkerSize(1.2)
   dmax = dataih.GetMaximum()
  
  if drawQCD:
  #### QCD
   qh = theFile.Get('qh')
   qh.SetFillColor(q)
   qh.Rebin(rebin)
   qh.Scale(sf_QCD)
   qh.Draw()
 #### Drell
  zih = theFile.Get('zih')
  zih.SetFillColor(z)
  zih.Rebin(rebin)
  zih.Scale(sf_Drell)
  zih.Draw()
 #### Diboson
  wwih = theFile.Get('wwih')
  wwih.SetFillColor(d)
  wwih.Rebin(rebin)
  wwih.Scale(sf_VV)
  wwih.Draw()
  ###
  wzih = theFile.Get('wzih')
  wzih.SetFillColor(d)
  wzih.Rebin(rebin)
  wzih.Scale(sf_VV)
  wzih.Draw()
  ###
  zzih = theFile.Get('zzih')
  zzih.SetFillColor(d)
  zzih.Rebin(rebin)
  zzih.Scale(sf_VV)
  zzih.Draw()
  ####
  dih = wwih.Clone()
  dih.SetName('dih')
  dih.Add(wzih)
  dih.Add(zzih)
  dih.Draw()
 #### Single Top
  # s
  t_sih = theFile.Get('t_sih')
  t_sih.SetFillColor(ts)
  t_sih.Rebin(rebin)
  t_sih.Scale(sf_T)
  t_sih.Draw()
  ###
  tb_sih = theFile.Get('tb_sih')
  tb_sih.SetFillColor(ts)
  tb_sih.Rebin(rebin)
  tb_sih.Scale(sf_Tbar)
  tb_sih.Draw()
  #### 
  stsih = t_sih.Clone()
  stsih.SetName('stsih')
  stsih.Add(tb_sih)
  stsih.Draw()
  ###
  # t
  t_tih = theFile.Get('t_tih')
  t_tih.SetFillColor(tt)
  t_tih.Rebin(rebin)
  t_tih.Scale(sf_T)
  t_tih.Draw()
  ###
  tb_tih = theFile.Get('tb_tih')
  tb_tih.SetFillColor(tt)
  tb_tih.Rebin(rebin)
  tb_tih.Scale(sf_Tbar)
  tb_tih.Draw()
  #### 
  sttih = t_tih.Clone()
  sttih.SetName('sttih')
  sttih.Add(tb_tih)
  sttih.Draw()
  ###
  # tw
  t_twih = theFile.Get('t_twih')
  t_twih.SetFillColor(ttw)
  t_twih.Rebin(rebin)
  t_twih.Scale(sf_tW)
  t_twih.Draw()
  ###
  tb_twih = theFile.Get('tb_twih')
  tb_twih.SetFillColor(ttw)
  tb_twih.Rebin(rebin)
  tb_twih.Scale(sf_tW)
  tb_twih.Draw()
  #### 
  sttwih = t_twih.Clone()
  sttwih.SetName('sttwih')
  sttwih.Add(tb_twih)
  sttwih.Draw()
 #### TTbar
  tt_semiih = theFile.Get('tt_semiih')
  tt_semiih.SetFillColor(ttb)
  tt_semiih.Rebin(rebin)
  tt_semiih.Scale(sf_TTbar)
  tt_semiih.Draw()
  ###
  tt_fullih = theFile.Get('tt_fullih')
  tt_fullih.SetFillColor(ttb)
  tt_fullih.Rebin(rebin)
  tt_fullih.Scale(sf_TTbar)
  tt_fullih.Draw()
  ###
  ttbih = tt_semiih.Clone()
  ttbih.SetName('ttbih')
  ttbih.Add(tt_fullih)
  ttbih.Draw()

  bmin = ttbih.GetXaxis().FindBin(xmin)
  bmax = ttbih.GetXaxis().FindBin(xmax)
  print(ttbih.Integral(bmin,bmax))
 #### W + Light
  print("W + Light")
  wlnih = theFile.Get('wlnih')
  wlnih.SetFillColor(wl)
  wlnih.Rebin(rebin)
  print(wlnih.Integral(bmin,bmax))
  ###
  wl1ih = theFile.Get('wl1ih')
  wl1ih.SetFillColor(wl)
  wl1ih.Rebin(rebin)
  print(wl1ih.Integral(bmin,bmax))
  ###
  wl2ih = theFile.Get('wl2ih')
  wl2ih.SetFillColor(wl)
  wl2ih.Rebin(rebin)
  print(wl2ih.Integral(bmin,bmax))
  ###
  wl3ih = theFile.Get('wl3ih')
  wl3ih.SetFillColor(wl)
  wl3ih.Rebin(rebin)
  print(wl3ih.Integral(bmin,bmax))
  ###
  wl4ih = theFile.Get('wl4ih')
  wl4ih.SetFillColor(wl)
  wl4ih.Rebin(rebin)
  print(wl4ih.Integral(bmin,bmax))
  ### 
  wlih = wlnih.Clone()
  wlih.SetName('wlih')
  wlih.Add(wl1ih)
  wlih.Add(wl2ih)
  wlih.Add(wl3ih)
  wlih.Add(wl4ih)
  wlih.Scale(sf_Wl)
  wlih.Draw()
  print(wlih.Integral(bmin,bmax))
 #### W + Charm
  print("W + Charm")
  wcnih = theFile.Get('wcnih')
  wcnih.SetFillColor(wc)
  wcnih.Rebin(rebin)
  print(wcnih.Integral(bmin,bmax))
  ###
  wc1ih = theFile.Get('wc1ih')
  wc1ih.SetFillColor(wc)
  wc1ih.Rebin(rebin)
  print(wc1ih.Integral(bmin,bmax))
  ###
  wc2ih = theFile.Get('wc2ih')
  wc2ih.SetFillColor(wc)
  wc2ih.Rebin(rebin)
  print(wc2ih.Integral(bmin,bmax))
  ###
  wc3ih = theFile.Get('wc3ih')
  wc3ih.SetFillColor(wc)
  wc3ih.Rebin(rebin)
  print(wc3ih.Integral(bmin,bmax))
  ###
  wc4ih = theFile.Get('wc4ih')
  wc4ih.SetFillColor(wc)
  wc4ih.Rebin(rebin)
  print(wc4ih.Integral(bmin,bmax))
  ###
  wcih = wcnih.Clone()
  wcih.SetName('wcih')
  wcih.Add(wc1ih)
  wcih.Add(wc2ih)
  wcih.Add(wc3ih)
  wcih.Add(wc4ih)
  wcih.Scale(sf_Wl)
  wcih.Draw()
  print(wcih.Integral(bmin,bmax))
 #### W + Charming
  print("W + Charming")
  wccnih = theFile.Get('wccnih')
  wccnih.SetFillColor(wcc)
  wccnih.Rebin(rebin)
  print(wccnih.Integral(bmin,bmax))
  ###
  wcc1ih = theFile.Get('wcc1ih')
  wcc1ih.SetFillColor(wcc)
  wcc1ih.Rebin(rebin)
  print(wcc1ih.Integral(bmin,bmax))
  ###
  wcc2ih = theFile.Get('wcc2ih')
  wcc2ih.SetFillColor(wcc)
  wcc2ih.Rebin(rebin)
  print(wcc2ih.Integral(bmin,bmax))
  ###
  wcc3ih = theFile.Get('wcc3ih')
  wcc3ih.SetFillColor(wcc)
  wcc3ih.Rebin(rebin)
  print(wcc3ih.Integral(bmin,bmax))
  ###
  wcc4ih = theFile.Get('wcc4ih')
  wcc4ih.SetFillColor(wcc)
  wcc4ih.Rebin(rebin)
  print(wcc4ih.Integral(bmin,bmax))
  ###
  wccih = wccnih.Clone()
  wccih.SetName('wccih')
  wccih.Add(wcc1ih)
  wccih.Add(wcc2ih)
  wccih.Add(wcc3ih)
  wccih.Add(wcc4ih)
  wccih.Scale(sf_Wcc)
  wccih.Draw()
  print(wccih.Integral(bmin,bmax))
 #### W + Beautiful
  print("W + Beautiful")
  wbbnih = theFile.Get('wbbnih')
  wbbnih.SetFillColor(wbb)
  wbbnih.Rebin(rebin)
  print(wbbnih.Integral(bmin,bmax))
  ###
  wbb1ih = theFile.Get('wbb1ih')
  wbb1ih.SetFillColor(wbb)
  wbb1ih.Rebin(rebin)
  print(wbb1ih.Integral(bmin,bmax))
  ###
  wbb2ih = theFile.Get('wbb2ih')
  wbb2ih.SetFillColor(wbb)
  wbb2ih.Rebin(rebin)
  print(wbb2ih.Integral(bmin,bmax))
  ###
  wbb3ih = theFile.Get('wbb3ih')
  wbb3ih.SetFillColor(wbb)
  wbb3ih.Rebin(rebin)
  print(wbb3ih.Integral(bmin,bmax))
  ###
  wbb4ih = theFile.Get('wbb4ih')
  wbb4ih.SetFillColor(wbb)
  wbb4ih.Rebin(rebin)
  print(wbb4ih.Integral(bmin,bmax))
  ###  
  wbbih = wbbnih.Clone()
  wbbih.SetName('wbbih')
  wbbih.Add(wbb1ih)
  wbbih.Add(wbb2ih)
  wbbih.Add(wbb3ih)
  wbbih.Add(wbb4ih)
  wbbih.Scale(sf_Wbb)
  wbbih.Draw()
  print(wbbih.Integral(bmin,bmax))
 
  hs = THStack('hs','')
  hs.SetTitle('')
  if drawQCD:
   hs.Add(qh)
  hs.Add(zih)
  hs.Add(dih)
  hs.Add(sttih)
  hs.Add(stsih)
  hs.Add(sttwih)
  hs.Add(ttbih)
  hs.Add(wlih)
  hs.Add(wcih)
  hs.Add(wccih)
  hs.Add(wbbih)

  hh = TH1F('hh','hh',steps,xmin,xmax)
  hh.Rebin(rebin)
  if drawQCD:
   hh.Add(qh)
  hh.Add(zih)
  hh.Add(dih)
  hh.Add(sttih)
  hh.Add(stsih)
  hh.Add(sttwih)
  hh.Add(ttbih)
  hh.Add(wbbih)
  hh.Add(wccih)
  hh.Add(wcih)
  hh.Add(wlih)
  if errorBand:
   hherr = hh.Clone()
   hherr.SetName('hherr')
   hherr.SetFillColor(ROOT.EColor.kBlue-6)
   hherr.SetFillStyle(3001)
 
  hs.Draw()
  hs.GetXaxis().SetTitle(xlabel)
  hs.GetXaxis().SetRangeUser(xmin,xmax)
#  if leaf=="mt" and not Legacy:
#   hs.GetXaxis().SetRangeUser(50,140)
  hs.GetYaxis().SetTitleOffset(1.5)
  hs.GetYaxis().SetTitle(ylabel)
  hsmax = hs.GetMaximum()
  
  leg=TLegend(0.68,0.2,0.9,0.8)
  if drawData:
   leg.AddEntry(dataih,'data')
 
  leg.AddEntry(wbbih,'W(#mu#nu)+b#bar{b}','f')
  leg.AddEntry(wccih,'W(#mu#nu)+c#bar{c}','f')
  leg.AddEntry(wcih,'W(#mu#nu)+c','f')
  leg.AddEntry(wlih,'W(#mu#nu)+light','f')
  leg.AddEntry(ttbih,'t#bar{t}','f')
  leg.AddEntry(sttwih,'t_tW','f')
  leg.AddEntry(stsih,'t_s','f')
  leg.AddEntry(sttih,'t_t','f')
  leg.AddEntry(dih,'WW,WZ,ZZ','f')
  leg.AddEntry(zih,'Drell-Yan','f')
  if drawQCD:
   leg.AddEntry(qh,'QCD','f')
  leg.SetFillColor(0)
  leg.SetBorderSize(0)
  
  theMax = hsmax 
  if drawData:
   theMax = max(dmax,hsmax) 
  hs.SetMaximum(1.2*theMax)
  c.cd()
  p1.cd()
  hs.Draw('hist')
  if errorBand:
   hherr.Draw('sames,E2')
  #hs.GetXaxis().SetRangeUser(70,110)
  if drawData:
   dataih.Draw('sames,E1')
  leg.Draw('sames')
  cpr.prelim_alt(lumi)
  tex.SetTextAlign(11)#left, bottom
  tex.DrawLatex(0.17,0.9,title)
  ## Draw Z Line
  if drawZ == True:
   zline = TLine(91.188,0,91.188,1.1*theMax)
   zline.SetLineWidth(2)
   zline.SetLineStyle(2)
   zline.Draw()
  ## Draw W Line
  if drawW == True:
   wline = TLine(80.385,0,80.385,1.1*theMax)
   wline.SetLineStyle(2)
   wline.SetLineWidth(2)
   wline.Draw()
  c.Update()
 #####
  c.cd()
  p2 = TPad('p2','p2',0,0,1,0.3)
  p2.SetTopMargin(0.1)
  p2.Draw()
  p2.cd()
 
  #datar = TH1F('datar','datar',steps,xmin,xmax)
  if drawData:
   datar = dataih.Clone()
  else:
   datar = hh.Clone()
  datar.SetName('datar')
  if errorBand:
   hNOerr = TH1F('hNOerr','hNOerr',steps,xmin,xmax)
   for j in range(steps+1):
    hNOerr.SetBinContent(j,hherr.GetBinContent(j))
    hNOerr.SetBinError(j,0)
   hrerr = hherr.Clone()
   hrerr.SetName("hrerr")
   hrerr.Divide(hNOerr)

  #if leaf =="mt" and not noMT:
  # datar.GetXaxis().SetRangeUser(50,140)
  datar.GetYaxis().SetRangeUser(1.-ratioRange,1.+ratioRange) 
  datar.GetYaxis().SetLabelSize(0.11)
  datar.Divide(hh)
  datar.Draw('ep')
  if errorBand: hrerr.Draw('sames,E2')
  #datar.GetXaxis().SetRangeUser(70,110)
  
  #if leaf=="mt" and not noMT:
  # l = TLine(50,1,140,1)
  # l.SetLineStyle(3)
  # l.Draw()
  #else:
  l = TLine(xmin,1,xmax,1)
  l.SetLineStyle(3)
  l.Draw()
  c.Update()
  print('you just read '+leaf)
  #save2 = raw_input ('Press Enter to Continue (type save to save)\n')
  #if save2 == 'save':
  c.Print(path+i+'.png')
   #c.Print(path+i+'_unscaled.png')
  print('')
  c.Close()
