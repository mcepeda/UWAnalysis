#!/usr/bin/env python
'''
Compares the shapes of variables from two data.root files
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle
import os

import aHisto as h #function to make histograms
import histoRange as hr #manages range, lables for plots

version = os.environ.get('version')
path = '../plots/%s_4F5Fcomp_'%(version)

leafs = [
'J1_pt',
'J2_pt',
'muon_pt',
'met_pt',
'mt_new'
]

filenames5f = [
"../plots/%s_2j2525_2bt_J1_pt.root"%(version),
"../plots/%s_2j2525_2bt_J2_pt.root"%(version),
"../plots/%s_2j2525_2bt_muon_pt.root"%(version),
"../plots/%s_2j2525_2bt_met_pt.root"%(version),
"../plots/%s_2j2525_2bt_mt_new.root"%(version),
]
filenames4f = [
"../plots/%s_Wbb_4F_J1_pt.root"%(version), 
"../plots/%s_Wbb_4F_J2_pt.root"%(version),
"../plots/%s_Wbb_4F_muon_pt.root"%(version),
"../plots/%s_Wbb_4F_met_pt.root"%(version),
"../plots/%s_Wbb_4F_mt_new.root"%(version),
]

can = TCanvas('can','can',800,800)
for filename5f,filename4f,leaf in zip(filenames5f,filenames4f,leafs):
 
 file4f = TFile(filename4f)
 file5f = TFile(filename5f) 

 steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)

 wbb_n_plot = file5f.Get("wbbnih")
 wbb_n_plot.SetName("wbb_n_plot")
 wbb_n_plot.SetLineColor(ROOT.EColor.kGreen-3)
 wbb_n_plot.SetLineStyle(3)
 wbb_n_plot.SetLineWidth(1)
 wbb_1_plot = file5f.Get("wbb1ih")
 wbb_1_plot.SetName("wbb_1_plot")
 wbb_1_plot.SetLineColor(ROOT.EColor.kYellow-3)
 wbb_1_plot.SetLineStyle(3)
 wbb_1_plot.SetLineWidth(1)
 wbb_2_plot = file5f.Get("wbb2ih")
 wbb_2_plot.SetName("wbb_2_plot")
 wbb_2_plot.SetLineColor(ROOT.EColor.kRed-3)
 wbb_2_plot.SetLineStyle(3)
 wbb_2_plot.SetLineWidth(1)
 wbb_3_plot = file5f.Get("wbb3ih")
 wbb_3_plot.SetName("wbb_3_plot")
 wbb_3_plot.SetLineColor(ROOT.EColor.kCyan-3)
 wbb_3_plot.SetLineStyle(3)
 wbb_3_plot.SetLineWidth(1)
 wbb_4_plot = file5f.Get("wbb4ih")
 wbb_4_plot.SetName("wbb_4_plot")
 wbb_4_plot.SetLineColor(ROOT.EColor.kMagenta-3)
 wbb_4_plot.SetLineStyle(3)
 wbb_4_plot.SetLineWidth(1)

 wbb_5f_plot = wbb_n_plot.Clone()
 wbb_5f_plot.SetName("wbb_5f_plot")
 wbb_5f_plot.Add(wbb_1_plot)
 wbb_5f_plot.Add(wbb_2_plot)
 wbb_5f_plot.Add(wbb_3_plot)
 wbb_5f_plot.Add(wbb_4_plot)
 wbb_5f_plot.SetLineWidth(3)
 wbb_5f_plot.SetLineColor(ROOT.EColor.kBlue)
 wbb_5f_plot.SetLineStyle(1)

 wbb_4f_plot = file4f.Get("W4F_No")
 wbb_4f_plot.SetName("wbb_4f_plot")
 wbb_4f_plot.SetLineColor(ROOT.EColor.kRed)
 wbb_4f_plot.SetLineWidth(2)
 wbb_4f_plot.SetLineStyle(1)
 
 leg=TLegend(0.5,0.5,0.9,0.7)
 leg.SetFillColor(0)
 leg.SetBorderSize(0)
 leg.AddEntry(wbb_4f_plot,'4 Flavor')
 leg.AddEntry(wbb_5f_plot,'5 Flavor')
 #leg.AddEntry(wbb_n_plot,'5 Flavor nJets')
 #leg.AddEntry(wbb_1_plot,'5 Flavor 1Jets')
 #leg.AddEntry(wbb_2_plot,'5 Flavor 2Jets')
 #leg.AddEntry(wbb_3_plot,'5 Flavor 3Jets')
 #leg.AddEntry(wbb_4_plot,'5 Flavor 4Jets')

 wbb_n_plot.Scale(1./wbb_5f_plot.Integral(-100,10000))
 wbb_1_plot.Scale(1./wbb_5f_plot.Integral(-100,10000))
 wbb_2_plot.Scale(1./wbb_5f_plot.Integral(-100,10000))
 wbb_3_plot.Scale(1./wbb_5f_plot.Integral(-100,10000))
 wbb_4_plot.Scale(1./wbb_5f_plot.Integral(-100,10000))
 #wbb_n_plot.Scale(1./wbb_n_plot.Integral(-100,10000))
 #wbb_1_plot.Scale(1./wbb_1_plot.Integral(-100,10000))
 #wbb_2_plot.Scale(1./wbb_2_plot.Integral(-100,10000))
 #wbb_3_plot.Scale(1./wbb_3_plot.Integral(-100,10000))
 #wbb_4_plot.Scale(1./wbb_4_plot.Integral(-100,10000))
 wbb_4f_plot.Scale(1./wbb_4f_plot.Integral(-100,10000))
 wbb_5f_plot.Scale(1./wbb_5f_plot.Integral(-100,10000))

 gStyle.SetOptStat(0)
 can.cd()
 #hnwr.SetMaximum(0.07)
 wbb_5f_plot.Draw('hist')
 wbb_n_plot.SetTitle(leaf)
 #wbb_n_plot.Draw('hist,sames')
 #wbb_1_plot.Draw('hist,sames')
 #wbb_2_plot.Draw('hist,sames')
 #wbb_3_plot.Draw('hist,sames')
 #wbb_4_plot.Draw('hist,sames')
 wbb_5f_plot.SetTitle(leaf)
 wbb_5f_plot.Draw('hist,sames')
 wbb_4f_plot.Draw('hist,sames')
 leg.Draw('sames')

 can.Update()
 print('you just read '+leaf)
 #print('Nr. Events:  %s'%(hnwr.Integral()))
 #print('Nr. Events:  %s'%(hnwr.Integral(-10,1000)))
 #print('Nr. Entries: %s'%(hnwr.GetEntries()))
# save = raw_input ('Press Enter to Continue (type save to save)\n')
# if save == 'save':
 can.Print(path+leaf+'.png') 
 file4f.Close()
 file5f.Close()
