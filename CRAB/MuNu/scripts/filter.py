#CutNtuple.py
#Author: Aaron Levine, T.M.Perry UW Madison
#Produces reduced ntuples, to speed up the analysis

from sys import argv, stdout, stderr
import ROOT
import math
#Returns the number of lines in the files
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

counter = 0

# recursively loops through directories
# adding dirnames to an array
def tree_finder(key,dirs):
 obj = key.ReadObj()
 if(obj.IsA().InheritsFrom("TDirectory")):
  dirs.append(obj.GetName())
  for key2 in obj.GetListOfKeys():
   tree_finder(key2,dirs)
   return 
 else:
  dirs.append(obj.GetName())

#define the file(s) + cuts
#indir='/hdfs/store/user/tperry/Halloween/'
indir='/afs/hep.wisc.edu/cms/tperry/Wbb_CMSSW_5_3_14_patch1/src/UWAnalysis/CRAB/MuNu/data/Halloween/'
outdir='/nfs_scratch/tperry/Halloween/'
fname = "TTbar_full.txt"
cut_str = "( (J1_pt>25) && (J2_pt>25) && ((J1_CSVbtag>0.898) || (J2_CSVbtag>0.898)) )"
cut_str = "(((HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[0]<0.12 && muon_pt_vec[0]>30 && abs(muon_eta_vec[0])<2.1) && (nrMuLoose==1 && nrEleLoose==0) && (J1_pt>25 && abs(J1_eta)<2.4 && J1_idLooseJet) && (J2_pt>25 && abs(J2_eta)<2.4 && J2_idLooseJet) && ( !(J3_pt>25 && J3_idLooseJet) && !(J4_pt>25 && J4_idLooseJet)) && ( !(J3_pt>25 && abs(J3_eta)>2.4 && J3_idLooseJet) && !(J4_pt>25 && abs(J4_eta)>2.4 && J4_idLooseJet) && !(J5_pt>25 && abs(J5_eta)>2.4) && !(J6_pt>25 && abs(J6_eta)>2.4) && !(J7_pt>25 && abs(J7_eta)>2.4) && !(J8_pt>25 && abs(J8_eta)>2.4) ) && (J1_CSVbtag>0.898) && (J2_CSVbtag>0.898) ))"
numfiles = file_len(indir+fname)


## FROM SYNCH
#noFJ = "( !(J3_pt>25 && abs(J3_eta)>2.4 && J3_idLooseJet) && !(J4_pt>25 && abs(J4_eta)>2.4 && J4_idLooseJet) && !(J5_pt>25 && abs(J5_eta)>2.4) && !(J6_pt>25 && abs(J6_eta)>2.4) && !(J7_pt>25 && abs(J7_eta)>2.4) && !(J8_pt>25 && abs(J8_eta)>2.4) )"
#no3J = "( !(J3_pt>25 && J3_idLooseJet) && !(J4_pt>25 && J4_idLooseJet))"
##muonIso = "(muNuRelPFIsoDB_A_vec.size()>0)"
##muonIso = "(nrMu>0 && muNuRelPFIsoDB_A_vec[1]<0.12)"
##muonIso = "((muNuRelPFIsoDB_A_vec[0]<0.12))"
#
#muonIso = "(muNuRelPFIsoDB_A_vec[0]<0.12 && muon_pt_vec[0]>30 && abs(muon_eta_vec[0])<2.1 )"
##muonIso = "(muNuRelPFIsoDB_A_vec[1]<0.12 && muon_pt_vec[1]>30 && abs(muon_eta_vec[1])<2.1)"
##muonIso = "(muNuRelPFIsoDB_A_vec[2]<0.12 && muon_pt_vec[2]>30 && abs(muon_eta_vec[2])<2.1)"
#
##cut_str = "(( %s && HLT_IsoMu24_eta2p1_v_fired)&&((nrMuLoose==1 && nrEle==0) && (J1_pt>25 && abs(J1_eta)<2.4 && J1_idLooseJet) && (J2_pt>25 && abs(J2_eta)<2.4 && J2_idLooseJet) && %s && %s ))"%(muonIso,no3J,noFJ)
##cut_str = "(( %s && HLT_IsoMu24_eta2p1_v_fired)&&((nrMuLoose==1 && nrEle==0) && (J1_pt>25 && abs(J1_eta)<2.4 && J1_idLooseJet) && (J2_pt>25 && abs(J2_eta)<2.4 && J2_idLooseJet) && %s && %s && (mt_new>45)))"%(muonIso,no3J,noFJ)
#cut_str = "(( %s && HLT_IsoMu24_eta2p1_v_fired)&&((nrMuLoose==1 && nrEle==0) && (J1_pt>25 && abs(J1_eta)<2.4 && J1_idLooseJet) && (J2_pt>25 && abs(J2_eta)<2.4 && J2_idLooseJet) && %s && %s && (J1_CSVbtag>0.898) && (J2_CSVbtag>0.898) && (mt_new>45)))"%(muonIso,no3J,noFJ)



with open(indir+fname) as f:
 for x in f:
  print "Reading " +indir+x
  if counter % 5 == 0: 
   print "Processed " + str(counter) + " out of " + str(numfiles) + " files"
  #get file of current ntuples from input indirectory
  x = x.rstrip()
  new_fname = "FltdForJel_"+x
  print "Writing "+new_fname
  ntuple_file = ROOT.TFile(indir+x)

  Smalltree_file = ROOT.TFile(outdir+new_fname,"RECREATE")
  for key in ntuple_file.GetListOfKeys():
   ntuple_file.cd()
   dirnames = []
   morekeys = tree_finder(key,dirnames)
   if dirnames[0]=="summary": continue
   if dirnames[0]=="muNuEventTreeMuonUp": break
   if dirnames[0]=="LHEFilter": continue
   ntuple_file_spot = "/".join(dirnames)
   new_tree_name = "_".join(dirnames)
   print ntuple_file_spot
   tree = ntuple_file.Get(ntuple_file_spot)
   tree.SetName("tree")
   print("Before: %s"%tree.GetEntries())
   Smalltree_file.cd()
   new_ttree = ROOT.TTree(new_tree_name,new_tree_name)
   new_ttree = tree.CopyTree(cut_str)
   new_ttree.SetName(new_tree_name)
   print("After: %s"%new_ttree.GetEntries())
   tree.Delete()
   #new_ttree.Write("eventTree")
  Smalltree_file.Write()
  Smalltree_file.Close()
  counter = counter+1

