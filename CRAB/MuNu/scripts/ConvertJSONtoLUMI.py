#!/usr/bin/env python

# JSON files for 22Jan re-reco data are now available (covering runs from 190456 to 208686)
# Golden JSON
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
# Muon JSON
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt
#
# The integrated luminosity is 19.79 and 20.65 fb-1 for Golden and Muon scenario respectively (with lumiCalc2.py of lumitype:HF, datatag:v10, normtag:HFV2c, worktag:V04-02-03)

import os, re, sys, json

if len(sys.argv) != 3:
 print "Usage: ConvertJSONtoLUMI.py [inputFile.json] [outputFile.lumi]"
 sys.exit(0)

input_file = open(sys.argv[1],'r')

if os.path.exists(sys.argv[2]):
 print "Output file '%s' already exists.\nRename it or delete it first!!" % (sys.argv[2])
 sys.exit(-1)

output_file = open(sys.argv[2],'w')
line_to_write ="import FWCore.ParameterSet.Config as cms\n"
line_to_write +='source = cms.Source("PoolSource",\n'
line_to_write +='fileNames = cms.untracked.vstring(),\n'
line_to_write +=")\n"
line_to_write += "source.lumisToProcess = cms.untracked.VLuminosityBlockRange()\n"
output_file.write(line_to_write)

entry = json.load(input_file)
for run in sorted(entry):
 #print "Run %s >>>" % (run)
 for lumirange in entry[run]:
  #print "\tLuminosity section %s" % (lumirange)
  line_to_write = "source.lumisToProcess.append('%s:%s-%s:%s')\n" % (run, lumirange[0], run, lumirange[1])
  output_file.write(line_to_write)
