#!/usr/bin/env python
'''
Compares the shapes of variables from two data.root files
Author: T.M.Perry UW
'''
import ROOT
from ROOT import THStack,TH1F,TFile
from ROOT import TLegend,TCanvas,TPad,TLatex,TLine
from ROOT import gROOT,gStyle,Double
from ROOT import *

import aHisto as h #function to make histograms
import histoRange as hr #manages range, lables for plots

file1_filename = '../data/post_synch_v6/Wbb4F.root' 

path = '../plots/21aout_W4F_easyMG'
out_file = TFile("%s.root"%(path),"RECREATE","heres a mg")
 
file1_file  = TFile( file1_filename)

eventTreeLocation='muNuEventTree/eventTree'
 
#leafs = [
##'J1_pt',
##'J2_pt',
##'muon_pt',
##'met_pt',
#'mt'
#]
leaf = 'mt'

cut='(((weightFactor * 19775 * weightEtaMuonID * weightEtaMuonIso * weightEtaMuonTrig) * ((((abs(J1_partonFlavour)==5 || abs(J1_partonFlavour)==4) * J1_CSVT_SFb)+(!(abs(J1_partonFlavour)==5 || abs(J1_partonFlavour)==4) * J1_CSVT_SFl)) * (((abs(J2_partonFlavour)==5 || abs(J2_partonFlavour)==4) * J2_CSVM_SFb)+(!(abs(J2_partonFlavour)==5 || abs(J2_partonFlavour)==4) * J2_CSVM_SFl)))) * ((muNuRelPFIsoDB_A<0.12 && HLT_IsoMu24_eta2p1_v_fired) && ((abs(muon_eta)<2.1 && muon_pt>30) && (nrMuLoose==1 && nrEle==0) && (J1_pt>25 && abs(J1_eta)<2.4 && J1_idLooseJet) && (J2_pt>25 && abs(J2_eta)<2.4 && J2_idLooseJet) && (J3_pt<25) && (nJets24Pt25==0) && (J1_CSVbtag>0.898) && (J2_CSVbtag>0.898) )))'
 #cut='(((weightFactor * 19775 * weightEtaMuonID * weightEtaMuonIso * weightEtaMuonTrig) * ((((abs(J1_partonFlavour)==5 || abs(J1_partonFlavour)==4) * J1_CSVT_SFb)+(!(abs(J1_partonFlavour)==5 || abs(J1_partonFlavour)==4) * J1_CSVT_SFl)) * (((abs(J2_partonFlavour)==5 || abs(J2_partonFlavour)==4) * J2_CSVM_SFb)+(!(abs(J2_partonFlavour)==5 || abs(J2_partonFlavour)==4) * J2_CSVM_SFl)))) * ((muNuRelPFIsoDB_A<0.12 && HLT_IsoMu24_eta2p1_v_fired) && ((abs(muon_eta)<2.1 && muon_pt>30) && (nrMuLoose==1 && nrEle==0) && (J1_pt>25 && abs(J1_eta)<2.4 && J1_idLooseJet) && (J2_pt>25 && abs(J2_eta)<2.4 && J2_idLooseJet) && (J3_pt<25) && (nJets24Pt25==0) && (J1_CSVbtag>0.898) && (J2_CSVbtag>0.898) && (mt>45))))'

I = -1
F = -1
#can = TCanvas('can','can',800,800)
#for eventTreeLocation in eventTreeLocations:
# print(eventTreeLocation)
# if eventTreeLocation == 'muNuEventTree/eventTree' : ps = 'No'
# if eventTreeLocation == 'muNuEventTreeMuonUp/eventTree' : ps ='muonUp'
# if eventTreeLocation == 'muNuEventTreeMuonDown/eventTree' : ps ='muonDown'
# if eventTreeLocation == 'muNuEventTreeJetUp/eventTree' : ps ='jetUp'
# if eventTreeLocation == 'muNuEventTreeJetDown/eventTree' : ps ='jetDown'
# if eventTreeLocation == 'muNuEventTreeUCEUp/eventTree': ps = 'UCEUp'
# if eventTreeLocation == 'muNuEventTreeUCEDown/eventTree': ps = 'UCEDown'
file1_tree  =  file1_file.Get(eventTreeLocation) 
#for leaf in leafs:
steps, xmin, xmax, xtitle, xunits, setLogY = hr.ranger(leaf)

hnwr,hnwrSize,hnwrSizePart,hnwrEvents = h.gram(file1_tree,leaf,xmin,xmax,steps,cut,I,F)
hnwr.SetName("hnwr")
out_file.Write()
##hnwr.Scale(1./hnwr.Integral())
#nwrMax = hnwr.GetMaximum()
#hnwr.SetLineColor(ROOT.EColor.kGreen+3)
#hnwr.SetLineWidth(3)
#
## leg=TLegend(0.5,0.5,0.9,0.7)
## leg.SetFillColor(0)
## leg.SetBorderSize(0)
## leg.AddEntry(hnwr,'Tom 2012B')
#
#can.cd()
##hnwr.SetMaximum(0.07)
#hnwr.SetTitle(leaf)
#hnwr.Draw('hist')
##leg.Draw('sames')
#gStyle.SetOptStat(1122)
#
##can.Update()
##err = Double(0.)
##print('you just read '+leaf)
##print('Nr. Events:  %s'%(hnwr.Integral()))
##print('Nr. Events:  %s'%(hnwr.IntegralAndError(-10,1000,err)))
##print('Error: %s'%(err))
##print('Nr. Entries: %s'%(hnwr.GetEntries()))
# save = raw_input ('Press Enter to Continue (type save to save)\n')
# if save == 'save':
#can.Print(path+'_'+leaf+'_'+ps+'.png') 
##hnwr.SetName(leaf)
##hnwr.SetName(leaf+ps)
##hnwr.Delete()
