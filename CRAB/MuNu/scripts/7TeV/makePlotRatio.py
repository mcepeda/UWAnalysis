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

leaf='mt'

lumi, steps, xmin, xmax, xtitle, xunits = hr.ranger(leaf)
lumiTitle = '#int L dt = %.1f fb^{-1}' %(lumi/1000.)
#xlabel = xtitle+' ['+xunits+']'
xlabel = 'Transverse Mass ['+xunits+']'
ylabel = 'Events/ %.001f' %((xmax-xmin)/(steps))
title = xtitle#+' Data v MC'
path = '../ForAN/'

jNr = 2
bNr = 0

name = []
name.append('sTopmt2j1bIso2unscaled_0_200')
#name.append('mt2j2bIso2unscaled_0_200')
#name.append('massDiMuon'+str(jNr)+'j'+str(bNr)+'bIso2rescaled_0_20')
#name.append('massDiMuon'+str(jNr)+'j'+str(bNr)+'bIso2rescaled_0_200')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso2rescaled_0_40')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso2rescaled_0_70')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso2rescaled_0_200')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso3rescaled_0_20')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso3rescaled_0_30')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso3rescaled_0_40')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso3rescaled_0_200')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso4rescaled_0_20')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso4rescaled_0_30')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso4rescaled_0_40')
#name.append('mt'+str(jNr)+'j'+str(bNr)+'bIso4rescaled_0_200')

#canvas attributes
canx = 775
cany = 900

tex = ROOT.TLatex()
tex.SetTextSize(0.07)
tex.SetTextAlign(13)
tex.SetNDC(True)
#color scheme
d = 1
q = ROOT.EColor.kRed+1
z = ROOT.EColor.kOrange-3
d = ROOT.EColor.kYellow-3
tt = ROOT.EColor.kGreen+1
ts = ROOT.EColor.kGreen-5
ttw = ROOT.EColor.kGreen+3 
ttb = ROOT.EColor.kGreen-9
wtau = ROOT.EColor.kAzure-7
wbx = ROOT.EColor.kCyan
wbb = 51#ROOT.EColor.kCyan
wcc = ROOT.EColor.kAzure+2
wc = ROOT.EColor.kBlue+1
wl = ROOT.EColor.kAzure+10

for i in name:
 c = TCanvas('c','Canvas Named c',canx,cany)
 p1 = TPad('p1','p1',0,0.3,1,1)
 p1.SetBottomMargin(0)
 p1.Draw()
 p1.cd()
 #p1.SetLogy()
 theFile = TFile(path+i+'.root')
 
 rebin = 2
 
 dataih = theFile.Get('dataih')
 dataih.Rebin(rebin)
 dataih.Sumw2()
 dataih.SetMarkerStyle(22)
 dataih.SetMarkerSize(1.2)
 dmax = dataih.GetMaximum()
 
 qh = theFile.Get('qh')
 qh.SetFillColor(q)
 qh.Rebin(rebin)
 zih = theFile.Get('zih')
 zih.SetFillColor(z)
 zih.Rebin(rebin)
 dih = theFile.Get('dih')
 dih.SetFillColor(d)
 dih.Rebin(rebin)
 ttih = theFile.Get('ttih')
 ttih.SetFillColor(tt)
 ttih.Rebin(rebin)
 tsih = theFile.Get('tsih')
 tsih.SetFillColor(ts)
 tsih.Rebin(rebin)
 ttwih = theFile.Get('ttwih')
 ttwih.SetFillColor(ttw)
 ttwih.Rebin(rebin)
 ttbih = theFile.Get('ttbih')
 ttbih.SetFillColor(ttb)
 ttbih.Rebin(rebin)
 wtih = theFile.Get('wtih')
 wtih.SetFillColor(wtau)
 wtih.Rebin(rebin)
 wbih = theFile.Get('wbih')
 wbih.SetFillColor(wbx)
 wbih.Rebin(rebin)
 wbbih = theFile.Get('wbbih')
 wbbih.SetFillColor(wbb)
 wbbih.Rebin(rebin)
 wccih = theFile.Get('wccih')
 wccih.SetFillColor(wcc)
 wccih.Rebin(rebin)
 wcih = theFile.Get('wcih')
 wcih.SetFillColor(wc)
 wcih.Rebin(rebin)
 wlih = theFile.Get('wlih')
 wlih.SetFillColor(wl)
 wlih.Rebin(rebin)
 
 hs = THStack('hs','')
 hs.SetTitle('')
 hs.Add(qh)
 hs.Add(zih)
 hs.Add(dih)
 hs.Add(ttih)
 hs.Add(tsih)
 hs.Add(ttwih)
 hs.Add(ttbih)
 hs.Add(wtih)
 hs.Add(wbih)
 hs.Add(wbbih)
 hs.Add(wccih)
 hs.Add(wcih)
 hs.Add(wlih)
 hs.Draw()
 hs.GetXaxis().SetTitle(xlabel)
 hs.GetXaxis().SetRangeUser(50,140)
 hs.GetYaxis().SetTitleOffset(1.5)
 hs.GetYaxis().SetTitle(ylabel)
 hsmax = hs.GetMaximum()
 
 leg=TLegend(0.68,0.2,0.9,0.8)
 leg.AddEntry(dataih,'data')
 leg.AddEntry(wlih,'W(#mu#nu)+light')
 leg.AddEntry(wcih,'W(#mu#nu)+c')
 leg.AddEntry(wccih,'W(#mu#nu)+cc')
 leg.AddEntry(wbbih,'W(#mu#nu)+bb')
 leg.AddEntry(wbih,'W(#mu#nu)+bx')
 leg.AddEntry(wtih,'W(#tau#nu)')
 leg.AddEntry(ttbih,'t#bar{t}')
 leg.AddEntry(ttwih,'t_tW')
 leg.AddEntry(tsih,'t_s')
 leg.AddEntry(ttih,'t_t')
 leg.AddEntry(dih,'WW,WZ')
 leg.AddEntry(zih,'Drell-Yan')
 leg.AddEntry(qh,'QCD')
 leg.SetFillColor(0)
 leg.SetBorderSize(0)
 
 if dmax>hsmax:
  hs.SetMaximum(1.2*dmax)
 else:
  hs.SetMaximum(1.2*hsmax)
# hs.SetMaximum(180)
 hs.Draw()
 dataih.Draw('sames')
 leg.Draw('sames')
 tex.DrawLatex(0,1,'Transverse Mass')
# tex.DrawLatex(0,1,title)
# tex.DrawLatex(0.15,0.88,str(jNr)+' Jets')
# tex.DrawLatex(0.15,0.8,str(bNr)+' bTags')
 tex.SetTextSize(0.06)
 tex.DrawLatex(0.26,0.88,'2 Jets')
 tex.DrawLatex(0.31,0.82,'one w/ b')
 tex.DrawLatex(0.31,0.76,'other w/ #eta > 2.8')
 tex.SetTextSize(0.03)
 tex.SetTextAlign(33)
 tex.DrawLatex(0.87,0.87,lumiTitle)
 tex.SetTextSize(0.07)
 tex.SetTextAlign(13)
# la = TLine(80.4,0,80.4,140)
# la = TLine(82,0,82,150)
# la.SetLineWidth(2)
# la.Draw('sames')
 c.Update()
#####
 c.cd()
 p2 = TPad('p2','p2',0,0,1,0.2)
 p2.SetTopMargin(0)
 p2.Draw()
 p2.cd()

 datar = TH1F('datar','datar',steps,xmin,xmax)
 hh = TH1F('hh','hh',steps,xmin,xmax)
 hh.Rebin(rebin)
 hh.Add(zih)
 hh.Add(dih)
 hh.Add(qh)
 hh.Add(ttih)
 hh.Add(tsih)
 hh.Add(ttwih)
 hh.Add(ttbih)
 hh.Add(wtih)
 hh.Add(wbih)
 hh.Add(wbbih)
 hh.Add(wccih)
 hh.Add(wcih)
 hh.Add(wlih)

 datar = dataih.Clone()
 datar.SetName('datar')
 datar.GetXaxis().SetRangeUser(50,140)
 datar.GetYaxis().SetRangeUser(0.5,1.5) 
 datar.GetYaxis().SetLabelSize(0.11)
 datar.Divide(hh)
 datar.Draw('ep')
 
 if 1 > 2:
  print('something is wrong')
 else:
#  l = TLine(xmin,1,xmax,1)
  l = TLine(50,1,144,1)
  l.SetLineStyle(3)
  l.Draw()
 c.Update()
 c.Print(path+i+'.png')
 c.Close()
