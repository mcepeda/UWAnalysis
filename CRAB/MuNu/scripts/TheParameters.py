#!/usr/bin/env python
'''
Parameters to be read by TheHistos.py and ThePlot.py
Author: T.M.Perry
'''
import os

def arams():
 lumi = 19783 # post_synch_v5,v6c

 btype = 't'
 njetcut = '25' #20,25,26,30,40
 jetcut = '25'
 I = 0
 F = 20
 iso_value = 0.12
 antiIso_value = 0.2
 extraCut = '(2>1)'
 #extraCut = '(J1SVMassb > 0 && J2SVMassb > 0)'
 #extraCut = '(vertixces > 25)'
 #extraCut = '(vertixces < 15)'

 #naming where output goes
 version = os.environ.get('version')
 path = '../plots/%s_'%(version)
 #extraName = 'QCDwTRG_'
 #extraName = 'MT45_rebin_'
 extraName = 'fullMT_BinsOfT_'
 #extraName = '1bSemiOff_J3FwdCut_'
 #extraName = '1bSemiOff_mu15_chargeOpp_'
 #extraName = 'bTag1M_noJetVeto_'

 drawQCD = True
 drawData = True

 noMT = False # control but no mT cut
 Control = False # 1 mu, 0b, 2+ jets
 Ctl_Andrea = False # 1 mu, 1+ b, 2+ jets 
 Signal = False # 1 mu, 2bs, 2 jets
 TT_m = True #1 mu, 2b, 3j
 TT_me = False # 1mu, 1e
 ST = False
 Z_Region = False # two muons, no mT cut
 Legacy = False #voids everything else and puts parametrs from 7Tev analysis
 Tomislav = False

 DataCard = False

 #eventTreeLocation = 'New_Tree'
 eventTreeLocation = 'muNuEventTree/eventTree'
 #eventTreeLocation = 'muNuEventTreeMuonUp/eventTree'
 #eventTreeLocation = 'muNuEventTreeMuonDown/eventTree'
 #eventTreeLocation = 'muNuEventTreeJetUp/eventTree'
 #eventTreeLocation = 'muNuEventTreeJetDown/eventTree'
 #eventTreeLocation = 'muNuEventTreeUCEUp/eventTree'
 #eventTreeLocation = 'muNuEventTreeUCEDown/eventTree'

 bNr = 0
 jNr = 2
 jetVeto = False

 if Signal or Legacy:
  #bNr = 0
  #bNr = 1
  bNr = 2
  jNr = 2
  jetVeto = True
 elif Control or noMT:
  bNr = 0 
  jNr = 2
  jetVeto = False
 elif Ctl_Andrea:
  bNr = 1 
  jNr = 2
  jetVeto = False
 elif TT_m:
  bNr = 2
  jNr = 3
  jetVeto = False
 elif TT_me:
  #bNr = 0
  #bNr = 1
  bNr = 2
  jNr = 2
  jetVeto = False
 elif Z_Region:
  bNr = 2
  jNr = 2
  jetVeto = False
 elif ST:
  bNr = 1
  jNr = 2
  jetVeto = True
 elif Tomislav:
  bNr = 0
  jNr = 2
  jetVeto = False
  njetcut = '20'
  jetcut = '20'

 drawW = False # line at W when plotting
 drawZ = False # line at Z when plotting

 if not DataCard:
  if eventTreeLocation == 'muNuEventTreeMuonUp/eventTree'  : extraName += 'mUP_'  
  if eventTreeLocation == 'muNuEventTreeMuonDown/eventTree': extraName += 'mDN_'  
  if eventTreeLocation == 'muNuEventTreeJetUp/eventTree'   : extraName += 'jUP_'  
  if eventTreeLocation == 'muNuEventTreeJetDown/eventTree' : extraName += 'jDN_'  
  if eventTreeLocation == 'muNuEventTreeUCEUp/eventTree'   : extraName += 'uUP_'
  if eventTreeLocation == 'muNuEventTreeUCEDown/eventTree' : extraName += 'uDN_'
 extraName += str(jNr)+'j'+njetcut+jetcut+'_'+str(bNr)+'b'+btype
 if Tomislav: extraName += '_tomi'
 if Z_Region: extraName += '_2mu'
 if noMT: extraName += '_noMT'
 if TT_me: extraName += '_1m1e'

 jet1_2dr2 = '(((J1_eta-J2_eta)*(J1_eta-J2_eta))+((J1_phi-J2_phi)*(J1_phi-J2_phi)))'

 leafs = [
'mt_new',
##'mt',
##'vertices',
#'J1_pt',
#'J1_eta',
#'muon_pt',
#'muon_eta',
#'met_pt',
#'met_phi',
#'ht',
#'mJ3J4',
#'J1_mass_SV_corrected',
#'J1_mass_SV_unweighted',
#'J1_mass_SV_weighted',
#'J2_pt',
#'J2_eta',
#'Wpt',
#'DiMuonMass',
#'DiMuonPt',
#'mJJ',
#'J1_mass_PV',
#'J1_mass_SSV',
#'J2_mass_PV',
#'J2_mass_SSV',
#'J2_mass_SV_corrected',
#'J2_mass_SV_unweighted',
#'J2_mass_SV_weighted',
#'J1_idBetaPu',
#'J1_idBetaClassicPu',
#'J1_idBetaStarPu',
#'J1_idBetaStarClassicPu',
#'J1_idBetaStarClassicModPu',
#'J1_PUID_fullDiscriminant',
#'J1_PUID_philv1Discriminant',
#'J1_fullIdLoosePu',
#'J2_fullIdLoosePu',
#'DiSVMass',
#'(J1_mass_SSV+J2_mass_SSV)',
#'J2_mass_SSV',
#'J1_nTracks_SV',
#'J2_nTracks_SV',
#'J1_mass_SSV_alt',
#'J1_mass_SSV',
#'J2_mass_SV_unweighted',
#'J2_mass_SV_weighted',
#'J1_normChi2_SV',
#'J2_normChi2_SV',
#'J1_flightDistance',
#'J1_flightDistanceError',
#'J1_normChi2_SV',
#'J1_pt_SSV',
#'J1_phi_SSV',
#'J1_eta_SSV',
#'J1_pt_PV',
#'J1_phi_PV',
#'J1_eta_PV',
#'J1_mass_PV',
#'(J1_mass_SV_weighted+J2_mass_SV_weighted)',
#'J1_phi',
#'muon_phi',
#'muonCharge',
#'J2_phi',
#'J2_mass_SSV_alt',
#'J2_mass_SSV',
#'J2_pt_SSV',
#'J2_phi_SSV',
#'J2_eta_SSV',
#'J2_pt_PV',
#'J2_phi_PV',
#'J2_eta_PV',
#'J2_mass_PV',
#'J3_pt',
#'J3_eta',
#'J3_phi',
#'J4_pt',
#'J4_eta',
#'J4_phi',
#'J1DR',
#'J2DR',
#'DiMuonMass',
#'nJetsPt25',
#'nJetsPt30',
#'J1_CSVbtag',
#'J2_CSVbtag',
#jet1_2dr2,
#'dz',
#'ptJJ',
#'mJJ',
]

 return lumi,bNr,btype,jNr,njetcut,jetcut,I,F,iso_value,antiIso_value,path,extraName,leafs,drawW,drawZ,drawQCD,drawData,jetVeto,Control,Z_Region,Legacy,noMT,TT_m,TT_me,ST,Signal,Tomislav,eventTreeLocation,extraCut,Ctl_Andrea
