#!/bin/bash

# WJets
echo "Making WJets Lists "
date +"%T"
find /hdfs/store/user/tperry/NouvelleYear_WJets_p1-MuEle-PATMC/*root >  data/NouvelleYear/list_WJets.txt
find /hdfs/store/user/tperry/NouvelleYear_WJets_p2-MuEle-PATMC/*root >> data/NouvelleYear/list_WJets.txt
find /hdfs/store/user/tperry/NouvelleYear_W1Jet_p1-MuEle-PATMC/*root >  data/NouvelleYear/list_W1Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_W1Jet_p2-MuEle-PATMC/*root >> data/NouvelleYear/list_W1Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_W2Jet_p1-MuEle-PATMC/*root >  data/NouvelleYear/list_W2Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_W2Jet_p2-MuEle-PATMC/*root >> data/NouvelleYear/list_W2Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_W3Jet_p1-MuEle-PATMC/*root >  data/NouvelleYear/list_W3Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_W3Jet_p2-MuEle-PATMC/*root >> data/NouvelleYear/list_W3Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_W4Jet-MuEle-PATMC/*root > data/NouvelleYear/list_W4Jet.txt
find /hdfs/store/user/tperry/NouvelleYear_Wbb4F-MuEle-PATMC/*root > data/NouvelleYear/list_Wbb4F.txt
#Top
echo "Making Top Lists"
date +"%T"
find /hdfs/store/user/tperry/NouvelleYear_TTbar_full-MuEle-PATMC/*root > data/NouvelleYear/list_TTbar_full.txt
find /hdfs/store/user/tperry/NouvelleYear_TTbar_semi-MuEle-PATMC/*root > data/NouvelleYear/list_TTbar_semi.txt
find /hdfs/store/user/tperry/NouvelleYear_T_s-MuEle-PATMC/*root >  data/NouvelleYear/list_T_s.txt
find /hdfs/store/user/tperry/NouvelleYear_T_t-MuEle-PATMC/*root >  data/NouvelleYear/list_T_t.txt
find /hdfs/store/user/tperry/NouvelleYear_T_tW-MuEle-PATMC/*root > data/NouvelleYear/list_T_tW.txt
find /hdfs/store/user/tperry/NouvelleYear_Tbar_s-MuEle-PATMC/*root >  data/NouvelleYear/list_Tbar_s.txt
find /hdfs/store/user/tperry/NouvelleYear_Tbar_t-MuEle-PATMC/*root >  data/NouvelleYear/list_Tbar_t.txt
find /hdfs/store/user/tperry/NouvelleYear_Tbar_tW-MuEle-PATMC/*root > data/NouvelleYear/list_Tbar_tW.txt
#DY/VV
echo "Making DY/VV Lists"
date +"%T"
find /hdfs/store/user/tperry/NouvelleYear_Drell-MuEle-PATMC/*root > data/NouvelleYear/list_Drell.txt
find /hdfs/store/user/tperry/NouvelleYear_WW-MuEle-PATMC/*root > data/NouvelleYear/list_WW.txt
find /hdfs/store/user/tperry/NouvelleYear_WZ-MuEle-PATMC/*root > data/NouvelleYear/list_WZ.txt
find /hdfs/store/user/tperry/NouvelleYear_ZZ-MuEle-PATMC/*root > data/NouvelleYear/list_ZZ.txt

touch data/NouvelleYear/initialEvents.txt

for list in \
"data/NouvelleYear/list_W1Jet.txt" \
"data/NouvelleYear/list_W2Jet.txt" \
"data/NouvelleYear/list_W3Jet.txt" \
"data/NouvelleYear/list_W4Jet.txt" \
"data/NouvelleYear/list_WJets.txt" \
"data/NouvelleYear/list_Drell.txt" \
"data/NouvelleYear/list_TTbar_full.txt" \
"data/NouvelleYear/list_TTbar_semi.txt" \
"data/NouvelleYear/list_T_s.txt" \
"data/NouvelleYear/list_T_t.txt" \
"data/NouvelleYear/list_T_tW.txt" \
"data/NouvelleYear/list_Tbar_s.txt" \
"data/NouvelleYear/list_Tbar_t.txt" \
"data/NouvelleYear/list_Tbar_tW.txt" \
"data/NouvelleYear/list_WW.txt" \
"data/NouvelleYear/list_WZ.txt" \
"data/NouvelleYear/list_Wbb4F.txt" \
"data/NouvelleYear/list_ZZ.txt" 

do
 echo "Reading Nr. Events found in " $list
 date +"%T"
 echo $list >> data/NouvelleYear/initialEvents.txt
 python scripts/eventCounter.py $list >> data/NouvelleYear/initialEvents.txt
done

