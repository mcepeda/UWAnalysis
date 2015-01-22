#!/usr/bin/env python
'''
for Wbb analysis, input a few values and function outputs cut strings
Author: T.M.Perry UW-Madison
'''

def cutmaker(leaf='J1_pt',isolationValue=0.12,antiIsoValue=0.2,lumi=19759.,bnr=0,btype='t',jnr=2,njetPt='20',jetPt='20',jetVeto=False,Control=True,Z_Region=False,Legacy=False,noMT=False,TT_m=False,TT_me=False,ST=False,Signal=False,Tomislav=False,extraCut='(2>1)',Ctl_Andrea=False):

 # mT
 #mt='(2>1)'
 if leaf=="mt_new": mt='(2>1)'
 else: mt = '(mt_new>45)'

 # Jets
 j1Cent = '(J1_pt>%s && abs(J1_eta)<2.4 && J1_idLooseJet)'%(jetPt)
 j2Cent = '(J2_pt>%s && abs(J2_eta)<2.4 && J2_idLooseJet)'%(jetPt)
 j3Cent = '(J3_pt>%s && abs(J3_eta)<2.4 && J3_idLooseJet)'%(jetPt)
 j4Cent = '(J4_pt>%s && abs(J4_eta)<2.4 && J4_idLooseJet)'%(jetPt)
 #noFJ = '(2>1)'
 #noFJ = '(nJets24Pt25==0)'
 noFJ = "( !(J3_pt>25 && abs(J3_eta)>2.4 && J3_idLooseJet) && !(J4_pt>25 && abs(J4_eta)>2.4 && J4_idLooseJet) && !(J5_pt>25 && abs(J5_eta)>2.4) && !(J6_pt>25 && abs(J6_eta)>2.4) && !(J7_pt>25 && abs(J7_eta)>2.4) && !(J8_pt>25 && abs(J8_eta)>2.4) )"
 j1Fwd = '(J1_pt>%s && abs(J1_eta)>=2.4 && abs(J1_eta)<5.0 && J1_idLooseJet)'%(jetPt)
 j2Fwd = '(J2_pt>%s && abs(J2_eta)>=2.4 && abs(J2_eta)<5.0 && J2_idLooseJet)'%(jetPt)
 #j3Veto = '(J3_pt<%s)'%(njetPt)
 j3Veto = "( !(J3_pt>%s && J3_idLooseJet) && !(J4_pt>%s && J4_idLooseJet))"%(njetPt,njetPt)

  # bTagging
 realh_j1 = "(abs(J1_partonFlavour)==5 || abs(J1_partonFlavour)==4)"
 realh_j2 = "(abs(J2_partonFlavour)==5 || abs(J2_partonFlavour)==4)"
 realb_j1 = "(abs(J1_partonFlavour)==5)"
 realb_j2 = "(abs(J2_partonFlavour)==5)"
 realc_j1 = "(abs(J1_partonFlavour)==4)"
 realc_j2 = "(abs(J2_partonFlavour)==4)"

 if btype == 't':
  bcut = '0.898'
  j1btagSF = '((%s * J1_CSVT_SFb)+(!%s * J1_CSVT_SFl))'%(realh_j1,realh_j1) 
  j1btagSF_up  = '((%s * (J1_CSVT_SFb+J1_CSVT_errorb) ) + (%s * (J1_CSVT_SFb+2*J1_CSVT_errorb) ) + (!%s * J1_CSVT_SFl_up))'%(realb_j1,realc_j1,realh_j1) 
  j1btagSF_dn  = '((%s * (J1_CSVT_SFb-J1_CSVT_errorb) ) + (%s * (J1_CSVT_SFb-2*J1_CSVT_errorb) ) + (!%s * J1_CSVT_SFl_down))'%(realb_j1,realc_j1,realh_j1) 
  j2btagSF  = '((%s * J2_CSVT_SFb)+(!%s * J2_CSVT_SFl))'%(realh_j2,realh_j2) 
  j2btagSF_up  = '((%s * (J2_CSVT_SFb+J2_CSVT_errorb) ) + (%s * (J2_CSVT_SFb+2*J2_CSVT_errorb) ) + (!%s * J2_CSVT_SFl_up))'%(realb_j2,realc_j2,realh_j2) 
  j2btagSF_dn  = '((%s * (J2_CSVT_SFb-J2_CSVT_errorb) ) + (%s * (J2_CSVT_SFb-2*J2_CSVT_errorb) ) + (!%s * J2_CSVT_SFl_down))'%(realb_j2,realc_j2,realh_j2) 
 if btype == 'm':
  bcut = '0.679'
  j1btagSF  = '((%s * J1_CSVM_SFb)+(!%s * J1_CSVM_SFl))'%(realh_j1,realh_j1) 
  j2btagSF  = '((%s * J2_CSVM_SFb)+(!%s * J2_CSVM_SFl))'%(realh_j2,realh_j2) 
 if btype == 'l':
  bcut = '0.244'
  j1btagSF  = '((%s * J1_CSVL_SFb)+(!%s * J1_CSVL_SFl))'%(realh_j1,realh_j1) 
  j2btagSF  = '((%s * J2_CSVL_SFb)+(!%s * J2_CSVL_SFl))'%(realh_j2,realh_j2) 
  
 j1btag = '(J1_CSVbtag>%s)'%(bcut)
 j2btag = '(J2_CSVbtag>%s)'%(bcut)
 weight_2b = '('+j1btagSF+' * '+j2btagSF+')'
 weight_1b = j1btagSF
 weight_2bup = '('+j1btagSF_up+' * '+j2btagSF_up+')'
 weight_2bdn = '('+j1btagSF_dn+' * '+j2btagSF_dn+')'
 weight_1bup = j1btagSF_up
 weight_1bdn = j1btagSF_dn

 weightEff = 'weightEtaMuonID * weightEtaMuonIso * weightEtaMuonTrig'
 weight = '(QweightFactor * %s * %s)'%(lumi,weightEff)
 weightT = '(QweightFactor * %s * %s * weightTop)'%(lumi,weightEff)
 weightTup = '(QweightFactor * %s * %s * weightTop * weightTop)'%(lumi,weightEff)

 # Leptons
 oneMu_noEle = '(nrMuLoose==1 && nrEleLoose==0)'
 oneMu_oneEle = '(nrMuLoose==1 && nrEle==1 && ele_pt>30)'
 twoMu_noEle = '(nrMuLoose==2 && nrEleLoose==0)'

 oneMu_noEleQCD = '(nrMuLoose==0 && nrEleLoose==0)'
 oneMu_oneEleQCD = '(nrMuLoose==0 && nrEle==1 && ele_pt>30)'
 twoMu_noEleQCD = '(nrMuLoose==0 && nrEleLoose==0)'

 trgIsoMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[0]<0.12 && muon_pt_vec[0]>30 && abs(muon_eta_vec[0])<2.1)"
 #trgIsoMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[1]<0.12 && muon_pt_vec[1]>30 && abs(muon_eta_vec[1])<2.1)"
 #trgIsoMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[2]<0.12 && muon_pt_vec[2]>30 && abs(muon_eta_vec[2])<2.1)"
 qcdMu = "(muNuRelPFIsoDB_A_vec[0]>0.2 && muon_pt_vec[0]>30 && abs(muon_eta_vec[0])<2.1)"
 #qcdMu = "(muNuRelPFIsoDB_A_vec[1]>0.2 && muon_pt_vec[1]>30 && abs(muon_eta_vec[1])<2.1)"
 #qcdMu = "(muNuRelPFIsoDB_A_vec[2]>0.2 && muon_pt_vec[2]>30 && abs(muon_eta_vec[2])<2.1)"
 trgqcdMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[0]>0.2 && muon_pt_vec[0]>30 && abs(muon_eta_vec[0])<2.1)"
 #trgqcdMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[1]>0.2 && muon_pt_vec[1]>30 && abs(muon_eta_vec[1])<2.1)"
 #trgqcdMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[2]>0.2 && muon_pt_vec[2]>30 && abs(muon_eta_vec[2])<2.1)"


 if Control:
  skim= '('+trgIsoMu+' && '+oneMu_noEle+' && '+j1Cent+' && '+j2Cent+' && '+mt+')'
  #skimq='('+qcdMu+' && '+oneMu_noEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+mt+')'
  skimq='('+trgqcdMu+' && '+oneMu_noEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+mt+')'

 if Signal:
  skim= '('+trgIsoMu+' && '+oneMu_noEle+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
  skimq='('+qcdMu+' && '+oneMu_noEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
  weightBup = '('+weight+' * '+weight_2bup+')'
  weightBdn = '('+weight+' * '+weight_2bdn+')'
  weightT   = '('+weightT+' * '+weight_2b+')'
  weightTup = '('+weightTup+' * '+weight_2b+')'
  weight    = '('+weight+' * '+weight_2b+')'
 if TT_m:
  skim= '('+trgIsoMu+' && '+oneMu_noEle+' && '+j1Cent+' && '+j2Cent+' && '+j3Cent+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
  skimq='('+qcdMu+' && '+oneMu_noEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j3Cent+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
  weightBup = '('+weight+' * '+weight_2bup+')'
  weightBdn = '('+weight+' * '+weight_2bdn+')'
  weightT   = '('+weightT+' * '+weight_2b+')'
  weightTup = '('+weightTup+' * '+weight_2b+')'
  weight    = '('+weight+' * '+weight_2b+')'
 if Z_Region:
  skim= '('+trgIsoMu+' && '+twoMu_noEle+ ' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+')'
  skimq='('+qcdMu+' && '+twoMu_noEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+')'
  weightBup = '('+weight+' * '+weight_2bup+')'
  weightBdn = '('+weight+' * '+weight_2bdn+')'
  weightT   = '('+weightT+' * '+weight_2b+')'
  weightTup = '('+weightTup+' * '+weight_2b+')'
  weight    = '('+weight+' * '+weight_2b+')'
 if TT_me:
# # 2b tight
#  skim= '('+trgIsoMu+' && '+oneMu_oneEle+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
#  skimq='('+qcdMu+' && '+oneMu_oneEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
#  weightBup = '('+weight+' * '+weight_2bup+')'
#  weightBdn = '('+weight+' * '+weight_2bdn+')'
#  weightT   = '('+weightT+' * '+weight_2b+')'
#  weightTup = '('+weightTup+' * '+weight_2b+')'
#  weight    = '('+weight+' * '+weight_2b+')'
#
# # 2b 1T 1M
#  j1btagSF = '((%s * J1_CSVT_SFb)+(!%s * J1_CSVT_SFl))'%(realh_j1,realh_j1)
#  j1btagSF_up  = '((%s * (J1_CSVT_SFb+J1_CSVT_errorb) ) + (%s * (J1_CSVT_SFb+2*J1_CSVT_errorb) ) + (!%s * J1_CSVT_SFl_up))'%(realb_j1,realc_j1,realb_j1)
#  j1btagSF_dn  = '((%s * (J1_CSVT_SFb-J1_CSVT_errorb) ) + (%s * (J1_CSVT_SFb-2*J1_CSVT_errorb) ) + (!%s * J1_CSVT_SFl_down))'%(realb_j1,realc_j1,realb_j1)
#  j2btagSF  = '((%s * J2_CSVM_SFb)+(!%s * J2_CSVM_SFl))'%(realh_j2,realh_j2)
#  j2btagSF_up  = '((%s * (J2_CSVM_SFb+J2_CSVM_errorb) ) + (%s * (J2_CSVM_SFb+2*J2_CSVM_errorb) ) + (!%s * J2_CSVM_SFl_up))'%(realb_j2,realc_j2,realb_j2)
#  j2btagSF_dn  = '((%s * (J2_CSVM_SFb-J2_CSVM_errorb) ) + (%s * (J2_CSVM_SFb-2*J2_CSVM_errorb) ) + (!%s * J2_CSVM_SFl_down))'%(realb_j2,realc_j2,realb_j2)
#  j1btag = '(J1_CSVbtag>0.898)'
#  j2btag = '(J2_CSVbtag>0.679)'
#  weight_2b = '('+j1btagSF+' * '+j2btagSF+')'
#  weight_2bup = '('+j1btagSF_up+' * '+j2btagSF_up+')'
#  weight_2bdn = '('+j1btagSF_dn+' * '+j2btagSF_dn+')'
#
#  skim= '('+trgIsoMu+' && '+oneMu_oneEle+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
#  skimq='('+qcdMu+' && '+oneMu_oneEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+j2btag+' && '+mt+')'
#  weightBup = '('+weight+' * '+weight_2bup+')'
#  weightBdn = '('+weight+' * '+weight_2bdn+')'
#  weightT   = '('+weightT+' * '+weight_2b+')'
#  weightTup = '('+weightTup+' * '+weight_2b+')'
#  weight    = '('+weight+' * '+weight_2b+')'
#
# # 1bTag CSVM
#  j1btagSF = '((%s * J1_CSVM_SFb)+(!%s * J1_CSVM_SFl))'%(realh_j1,realh_j1)
#  j1btagSF_up  = '((%s * (J1_CSVM_SFb+J1_CSVM_errorb) ) + (%s * (J1_CSVM_SFb+2*J1_CSVM_errorb) ) + (!%s * J1_CSVM_SFl_up))'%(realb_j1,realc_j1,realb_j1)
#  j1btagSF_dn  = '((%s * (J1_CSVM_SFb-J1_CSVM_errorb) ) + (%s * (J1_CSVM_SFb-2*J1_CSVM_errorb) ) + (!%s * J1_CSVM_SFl_down))'%(realb_j1,realc_j1,realb_j1)
#  j1btag = '(J2_CSVbtag>0.679)'
#  weight_1b = '('+j1btagSF+')'
#  weight_1bup = '('+j1btagSF_up+')'
#  weight_1bdn = '('+j1btagSF_dn+')'

  skim= '('+trgIsoMu+' && '+oneMu_oneEle+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+mt+')'
  skimq='('+qcdMu+' && '+oneMu_oneEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j3Veto+' && '+noFJ+' && '+j1btag+' && '+mt+')'
  #skim= '('+trgIsoMu+' && '+oneMu_oneEle+' && '+j1Cent+' && '+j2Cent+' && '+j1btag+' && '+mt+')'
  #skimq='('+qcdMu+' && '+oneMu_oneEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j1btag+' && '+mt+')'
  weightBup = '('+weight+' * '+weight_1bup+')'
  weightBdn = '('+weight+' * '+weight_1bdn+')'
  weightT   = '('+weightT+' * '+weight_1b+')'
  weightTup = '('+weightTup+' * '+weight_1b+')'
  weight    = '('+weight+' * '+weight_1b+')'

#  # TTbar official analysis  http://arxiv.org/pdf/1312.7582v2.pdf
#  #trgIsoMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[0]<0.12 && muon_pt_vec[0]>20 && abs(muon_eta_vec[0])<2.1)"
#  trgIsoMu = "(HLT_IsoMu24_eta2p1_v_fired && muNuRelPFIsoDB_A_vec[0]<0.15 && muon_pt_vec[0]>20 && abs(muon_eta_vec[0])<2.1)"
#  qcdMu = "(muNuRelPFIsoDB_A_vec[0]>0.2 && muon_pt_vec[0]>20 && abs(muon_eta_vec[0])<2.1)"
#  j1Cent = '(J1_pt>30 && abs(J1_eta)<2.5 && J1_idLooseJet)'
#  j2Cent = '(J2_pt>30 && abs(J2_eta)<2.5 && J2_idLooseJet)'
#  oneMu_oneEle = '(nrMuLoose==1 && nrEle==1 && ele_pt>20)'
#  charge_opp = '(ele_charge * muon_charge)==-1'
#
#  skim= '('+trgIsoMu+' && '+oneMu_oneEle+' && '+j1Cent+' && '+j2Cent+' && '+j1btag+' && '+charge_opp+')'
#  skimq='('+qcdMu+' && '+oneMu_oneEleQCD+' && '+j1Cent+' && '+j2Cent+' && '+j1btag+' && '+charge_opp+')'
#  weightBup = '('+weight+' * '+weight_1bup+')'
#  weightBdn = '('+weight+' * '+weight_1bdn+')'
#  weightT   = '('+weightT+' * '+weight_1b+')'
#  weightTup = '('+weightTup+' * '+weight_1b+')'
#  weight    = '('+weight+' * '+weight_1b+')'

 if ST:
  skim= '('+trgIsoMu+' && '+oneMu_noEle+' && '+j1Cent+' && '+j2Fwd+' && '+j3Veto+' && '+j1btag+' && '+mt+')'
  skimq='('+qcdMu+' && '+oneMu_noEleQCD+' && '+j1Cent+' && '+j2Fwd+' && '+j3Veto+' && '+j1btag+' && '+mt+')'
  weightBup = '('+weight+' * '+weight_1bup+')'
  weightBdn = '('+weight+' * '+weight_1bdn+')'
  weightT   = '('+weightT+' * %s)'%(j1btagSF)
  weightTup = '('+weightTup+' * %s)'%(j1btagSF)
  weight    = '('+weight+' * %s)'%(j1btagSF)

 # for splitting up the W samples, we should probably match to jets
# genB = '(nbHadrons>0)' 
# genC = '((nrC+nrCbar)>0)'
# evenC = '(((nrC+nrCbar)%2)==0)'

 genB = '(abs(J1_partonFlavour)==5 || abs(J2_partonFlavour)==5)'
 genC = '(abs(J1_partonFlavour)==4 || abs(J2_partonFlavour)==4)'
 evenC = '(abs(J1_partonFlavour)==4 && abs(J2_partonFlavour)==4)'
 
 Wl = '((!'+genB+'&&!'+genC+'))'
 Wc = '((!'+genB+'&&'+genC+'&&!'+evenC+'))'
 Wcc = '((!'+genB+'&&'+genC+'&&'+evenC+'))'
 Wbb = '('+genB+')'

 cutDataNonIso   = '('+skimq+')' #Data Non Iso
 cutDataIso      = '('+skim+')'    #Data Iso

 cutMcNonIso     = '('+weight+' * ('+skimq+'))' #MC Non Iso
 cutMcNonIsoW    = cutMcNonIso #'('+weightW+'*('+NonIso+'&&'+theCut+'))' #MC Non IsoW
 cutMcIso        = '('+weight+' * ('+skim+'))' #MC Iso
 cutMcNonIsoTup  = '('+weightTup+' * ('+skimq+'))' #MC Non Iso
 cutMcIsoTup     = '('+weightTup+' * ('+skim+'))' #MC Iso
 cutMcNonIsoTdn  = cutMcNonIso
 cutMcIsoTdn     = cutMcIso
 cutMcNonIsoT    = '('+weightT+' * ('+skimq+'))' #MC Non Iso
 cutMcIsoT       = '('+weightT+' * ('+skim +'))' #MC Iso
 cutMcWlNonIso   = '('+weight+' *  ('+skimq+' && '+Wl+'))'
 cutMcWlIso      = '('+weight+' *  ('+skim+ ' && '+Wl+'))'
 cutMcWcNonIso   = '('+weight+' *  ('+skimq+' && '+Wc+'))'
 cutMcWcIso      = '('+weight+' *  ('+skim+ ' && '+Wc+'))'
 cutMcWccNonIso  = '('+weight+' *  ('+skimq+' && '+Wcc+'))'
 cutMcWccIso     = '('+weight+' *  ('+skim+ ' && '+Wcc+'))'
 cutMcWbbNonIso  = '('+weight+' *  ('+skimq+' && '+Wbb+'))'
 cutMcWbbIso     = '('+weight+' *  ('+skim+ ' && '+Wbb+'))'

 cutMcNonIso_Bup     = '1'
 cutMcNonIsoW_Bup    = '1' #MC Non IsoW
 cutMcIso_Bup        = '1'
 cutMcWlNonIso_Bup   = '1'
 cutMcWlIso_Bup      = '1'
 cutMcWcNonIso_Bup   = '1'
 cutMcWcIso_Bup      = '1'
 cutMcWccNonIso_Bup  = '1'
 cutMcWccIso_Bup     = '1'
 cutMcWbbNonIso_Bup  = '1'
 cutMcWbbIso_Bup     = '1'
 
 # B sf dn
 cutMcNonIso_Bdn     = '1'
 cutMcNonIsoW_Bdn    = '1' 
 cutMcIso_Bdn        = '1'
 cutMcWlNonIso_Bdn   = '1'
 cutMcWlIso_Bdn      = '1'
 cutMcWcNonIso_Bdn   = '1'
 cutMcWcIso_Bdn      = '1'
 cutMcWccNonIso_Bdn  = '1'
 cutMcWccIso_Bdn     = '1'
 cutMcWbbNonIso_Bdn  = '1'
 cutMcWbbIso_Bdn     = '1'

 if not Control:
  # B sf up
  cutMcNonIso_Bup     = '('+weightBup+' * ('+skimq+'))' #MC Non Iso
  cutMcNonIsoW_Bup    = cutMcNonIso_Bup #'('+weightBupW+'*('+NonIso+'&&'+theCut+'))' #MC Non IsoW
  cutMcIso_Bup        = '('+weightBup+' * ('+skim+'))' #MC Iso
  print("inside TheCut: Up: \n %s \n"%(cutMcIso_Bup))
  cutMcWlNonIso_Bup   = '('+weightBup+' * ('+skimq+' && '+Wl+'))'
  cutMcWlIso_Bup      = '('+weightBup+' * ('+skim+ ' && '+Wl+'))'
  cutMcWcNonIso_Bup   = '('+weightBup+' * ('+skimq+' && '+Wc+'))'
  cutMcWcIso_Bup      = '('+weightBup+' * ('+skim+ ' && '+Wc+'))'
  cutMcWccNonIso_Bup  = '('+weightBup+' * ('+skimq+' && '+Wcc+'))'
  cutMcWccIso_Bup     = '('+weightBup+' * ('+skim+ ' && '+Wcc+'))'
  cutMcWbbNonIso_Bup  = '('+weightBup+' * ('+skimq+' && '+Wbb+'))'
  cutMcWbbIso_Bup     = '('+weightBup+' * ('+skim+ ' && '+Wbb+'))'
 
  # B sf dn
  cutMcNonIso_Bdn     = '('+weightBdn+' * ('+skimq+'))' #MC Non Iso
  cutMcNonIsoW_Bdn    = cutMcNonIso_Bdn #'('+weightBdnW+'*('+NonIso+'&&'+theCut+'))' #MC Non IsoW
  cutMcIso_Bdn        = '('+weightBdn+' * ('+skim+'))' #MC Iso
  cutMcWlNonIso_Bdn   = '('+weightBdn+' * ('+skimq+' && '+Wl+'))'
  cutMcWlIso_Bdn      = '('+weightBdn+' * ('+skim+ ' && '+Wl+'))'
  cutMcWcNonIso_Bdn   = '('+weightBdn+' * ('+skimq+' && '+Wc+'))'
  cutMcWcIso_Bdn      = '('+weightBdn+' * ('+skim+ ' && '+Wc+'))'
  cutMcWccNonIso_Bdn  = '('+weightBdn+' * ('+skimq+' && '+Wcc+'))'
  cutMcWccIso_Bdn     = '('+weightBdn+' * ('+skim+ ' && '+Wcc+'))'
  cutMcWbbNonIso_Bdn  = '('+weightBdn+' * ('+skimq+' && '+Wbb+'))'
  cutMcWbbIso_Bdn     = '('+weightBdn+' * ('+skim+ ' && '+Wbb+'))'

 return cutMcNonIso, cutMcNonIsoW, cutMcIso, cutDataNonIso, cutDataIso, cutMcWlNonIso, cutMcWlIso, cutMcWcNonIso, cutMcWcIso, cutMcWccNonIso, cutMcWccIso, cutMcWbbNonIso, cutMcWbbIso, cutMcNonIsoT, cutMcIsoT, cutMcNonIsoTup, cutMcIsoTup, cutMcNonIsoTdn, cutMcIsoTdn, cutMcNonIso_Bup, cutMcNonIsoW_Bup, cutMcIso_Bup, cutMcWlNonIso_Bup, cutMcWlIso_Bup, cutMcWcNonIso_Bup, cutMcWcIso_Bup, cutMcWccNonIso_Bup, cutMcWccIso_Bup, cutMcWbbNonIso_Bup, cutMcWbbIso_Bup, cutMcNonIso_Bdn, cutMcNonIsoW_Bdn, cutMcIso_Bdn, cutMcWlNonIso_Bdn, cutMcWlIso_Bdn, cutMcWcNonIso_Bdn, cutMcWcIso_Bdn, cutMcWccNonIso_Bdn, cutMcWccIso_Bdn, cutMcWbbNonIso_Bdn, cutMcWbbIso_Bdn   
