#!/usr/bin/env python
'''
Makes a single histogram and gives its total size and over a range
Author: T.M.Perry UW
'''
import ROOT
from ROOT import *

def gram(tree,leaf,xmin=0,xmax=200,steps=20,cuts='(2>1)',I=-1,F=-1):
 
 histo = TH1F('histo','histo',steps,xmin,xmax)
 tree.Draw(leaf+'>>histo',cuts,'goff')

 histo.SetTitle('') 
 histo.Sumw2()
 gStyle.SetOptStat('')
 gROOT.SetStyle('Plain')

 bmin = histo.GetXaxis().FindBin(xmin)
 bmax = histo.GetXaxis().FindBin(xmax)
 if I == -1 and F == -1: 
  I = xmin
  F = xmax
 bminI = histo.GetXaxis().FindBin(I)
 bmaxF = histo.GetXaxis().FindBin(F)
 
 size = histo.Integral(bmin,bmax)
 sizePart = histo.Integral(bminI,bmaxF)
 entries = histo.GetEntries()
 return histo,size,sizePart,entries
