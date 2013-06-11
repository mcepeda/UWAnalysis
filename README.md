UWAnalysis
==========
UW Analysis for 53x (8Tev)


download this following the recipe below:

scram pro CMSSW CMSSW_5_3_7
cd CMSSW_5_3_7/src/

cmsenv

#set up cvs
#export CVSROOT=:gserver:<user_name>@cmssw.cvs.cern.ch:/local/reps/CMSSW
#kinit <user_name>@CERN.CH; aklog -cell cern.ch

git clone https://github.com/tmrhombus/UWAnalysis.git
./UWAnalysis/recipe

#check before compiling
#showtags #should give 38 (34 could mean cvs not up)

scram build -c
scramv1 build -j 8 
