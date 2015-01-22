#!/usr/bin/env python
'''
Compares the shapes of variables from the three of us
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle

import aHisto as h #function to make histograms
import histoRange as hr #manages range, lables for plots

path = '../scripts/synch/Data_norm_comp_wbb'
#out_file = gROOT.FindObject("%s.root"%(path))
#if out_file: out_file.Close()
#out_file = TFile("%s.root"%(path),"RECREATE","synching")

tom_filename = '../scripts/synch/Data_norm_wbb.root' 
and_filename = '../scripts/synch/Data_norm_Andrea_SingleMuonB_v09_dumphistos.root'
tom_file = TFile(tom_filename)
and_file = TFile(and_filename)

tom_leafs = [
'J1_pt',
'J2_pt',
'muon_pt',
'met_pt',
'mt'
]
and_leafs = [
'h_first_bjet_pt',
'h_second_bjet_pt',
'h_first_muon_pt_bb',
'h_MET_bb',
'h_mt_wmnu_bb_wide'
]

jel_filenames = [
'../scripts/synch/Data_norm_Jelena_maxhJet_pt0hJet_ptSingleMuBJetNormalized.root',
'../scripts/synch/Data_norm_Jelena_minhJet_pt0hJet_ptSingleMuBJetNormalized.root',
'../scripts/synch/Data_norm_Jelena_vLepton_pt0SingleMuBOtherNormalized.root',
'../scripts/synch/Data_norm_Jelena_METtype1corretSingleMuBOtherNormalized.root',
'../scripts/synch/Data_norm_Jelena_VMtSingleMuBOtherNormalized.root'
]

I = -1
F = -1
can = TCanvas('can','can',800,800)

for leaf,a_leaf,jel_filename in zip(tom_leafs,and_leafs,jel_filenames):
 jel_file = TFile(jel_filename)
 steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)
 bmin = 1
 bmax = steps+1

 tomh = tom_file.Get(leaf)
 tomh.SetName('tomh')
 tomh.Scale(1./tomh.Integral())
 tomh.SetLineColor(ROOT.EColor.kGreen+3)
 tomh.SetLineWidth(3)
 tomh.Draw()
 tmax = tomh.GetMaximum()

 jelh = jel_file.Get('data_obs')
 jelh.SetName('jelh')
 jelh.Scale(1./jelh.Integral())
 jelh.SetLineColor(ROOT.EColor.kRed)
 jelh.SetLineWidth(2)
 jelh.Draw()
 jmax = jelh.GetMaximum()

 andh = and_file.Get(a_leaf)
 andh.SetName('andh')
 andh.Scale(1./andh.Integral())
 andh.SetLineColor(ROOT.EColor.kBlue)
 andh.SetLineWidth(1)
 andh.Draw()
 amax = andh.GetMaximum()

 leg=TLegend(0.5,0.5,0.9,0.7)
 leg.SetFillColor(0)
 leg.SetBorderSize(0)
 leg.AddEntry(tomh,'Tom W+bb','l')
 leg.AddEntry(jelh,'Jelena W+bb','l')
 leg.AddEntry(andh,'Andrea W+bb','l')

 can.cd()
 tomh.SetTitle(leaf)
 tomh.SetMaximum(1.1*max(tmax,jmax,amax))
 tomh.Draw('hist')
 jelh.Draw('hist,sames')
 andh.Draw('hist,sames')
 leg.Draw('sames')
 gStyle.SetOptStat('')

 can.Update()
 print('you just read '+leaf)
# save = raw_input ('Press Enter to Continue (type save to save)\n')
# if save == 'save':
 can.Print(path+'_'+leaf+'.png') 
 tomh.SetName('tom_'+leaf)
 jelh.SetName('jelena_'+leaf)
 andh.SetName('andrea_'+leaf)
 #out_file.Write()
 tomh.Delete()
 jelh.Delete()
 andh.Delete()
