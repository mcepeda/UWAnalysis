#!/usr/bin/python

'''
adds up the total number of evnts run over MC in a list of files
pulling info from summary/results

usage:
python eventCounter.py list_of_file_locations.txt

Author: T.M.Perry UW-Madison
'''

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
total_events = 0

list_of_files = argv[1]

numfiles = file_len(list_of_files)
with open(list_of_files) as f:
 for filename in f:
  if counter % 100 == 0:
   print "Processed " + str(counter) + " out of " + str(numfiles) + " files"
  #get file of current ntuples from input directory
  filename = filename.rstrip()
  ntuple_file = ROOT.TFile(filename)
  summary = ntuple_file.Get("summary/results")
  nr_events = summary.GetBinContent(1)
  total_events+=nr_events
  counter = counter+1
print("Total Number of Events: %s"%(total_events))
