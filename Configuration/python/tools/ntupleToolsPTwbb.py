import FWCore.ParameterSet.Config as cms

from UWAnalysis.Configuration.tools.analysisToolsPT import TriggerPaths
jetRanks = [0,1,2,3]
jetNames = ['J1_','J2_','J3_','J4_']

def makeJetUserFloat(floatName,xn='',source = 'wCandsJets'):
 PSet_List = []
 for rank,name in zip(jetRanks,jetNames):
  nameTag = name+floatName+xn
  PSet_List.append(cms.PSet(
        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
        src        = cms.InputTag(source),
        tag        = cms.string(nameTag),
        method     = cms.string('userFloat("'+floatName+'")'),
        rank       = cms.untracked.double(rank)
  ))
 return PSet_List

def makeJetUserInt(intName,xn='',source = 'wCandsJets'):
 PSet_List = []
 for rank,name in zip(jetRanks,jetNames):
  nameTag = name+intName+xn
  PSet_List.append(cms.PSet(
        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
        src        = cms.InputTag(source),
        tag        = cms.string(nameTag),
        method     = cms.string('userInt("'+intName+'")'),
        rank       = cms.untracked.double(rank)
  ))
 return PSet_List

def makeJetStringPar(strName,xn='',source='wCandsJets'):
 PSet_List = []
 for rank,name in zip(jetRanks,jetNames):
  nameTag = name+strName+xn
  PSet_List.append(cms.PSet(
        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
        src        = cms.InputTag(source),
        tag        = cms.string(nameTag),
        method     = cms.string(strName+'()'),
        rank       = cms.untracked.double(rank)
  ))
 return PSet_List

def makeJetString(strName,xn='',source='wCandsJets'):
 PSet_List = []
 for rank,name in zip(jetRanks,jetNames):
  nameTag = name+strName+xn
  PSet_List.append(cms.PSet(
        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
        src        = cms.InputTag(source),
        tag        = cms.string(nameTag),
        method     = cms.string(strName),
        rank       = cms.untracked.double(rank)
  ))
 return PSet_List

def makeJetStringName(strName,mthName,xn='',source='wCandsJets'):
 PSet_List = []
 for rank,name in zip(jetRanks,jetNames):
  nameTag = name+strName+xn
  PSet_List.append(cms.PSet(
        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
        src        = cms.InputTag(source),
        tag        = cms.string(nameTag),
        method     = cms.string(mthName),
        rank       = cms.untracked.double(rank)
  ))
 return PSet_List


def makeJetBTag(tagName,strName,source='wCandsJets'):
 PSet_List = []
 for rank,name in zip(jetRanks,jetNames):
  nameTag = name+tagName
  PSet_List.append(cms.PSet(
        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
        src        = cms.InputTag(source),
        tag        = cms.string(nameTag),
        method     = cms.string('bDiscriminator("'+strName+'")'),
        rank       = cms.untracked.double(rank)
  ))
 return PSet_List

def makeNJets(pt,source='wCandsJets'):
  PSet = cms.PSet(
        pluginType  = cms.string("PATMuonNuPairJetCountFiller"),
        src         = cms.InputTag(source),
        tag         = cms.string("nJetsPt"+str(pt)),
        method      = cms.string('pt()>'+str(pt)+'&&abs(eta())<2.4'),
        leadingOnly = cms.untracked.bool(True)
  )
  return PSet

def makeCollSize(srcName,tagName):
  PSet = cms.PSet(
        pluginType = cms.string("CollectionSizeFiller"),
        src        = cms.InputTag(srcName),
        tag        = cms.string(tagName)
  )
  return PSet

def makeZColl(tagName,methodName,sourceZ='diMuonsSorted'):
  PSet = cms.PSet(
        pluginType  = cms.string("PATMuPairFiller"),
        src         = cms.InputTag(sourceZ),
        tag         = cms.string(tagName),
        method      = cms.string(methodName),
        leadingOnly = cms.untracked.bool(True)
  )
  return PSet

def makeBasicEle(tagName,methodName,sourceEle='selectedPatElectrons'):
  PSet = cms.PSet(
        pluginType  = cms.string("PATElectronFiller"),
        src         = cms.InputTag(sourceEle),
        tag         = cms.string(tagName),
        method      = cms.string(methodName),
        leadingOnly = cms.untracked.bool(True)
  )
  return PSet

def makeMuNu(tagName,methodName,source='wCandsJets',lo=False):
  if lo:
   PSet = cms.PSet(
         pluginType  = cms.string("PATMuonNuPairFiller"),
         src         = cms.InputTag(source),
         tag         = cms.string(tagName),
         method      = cms.string(methodName),
         leadingOnly = cms.untracked.bool(True)
   )
  else:
   PSet = cms.PSet(
         pluginType = cms.string("PATMuonNuPairFiller"),
         src        = cms.InputTag(source),
         tag        = cms.string(tagName),
         method     = cms.string(methodName),
   )
  return PSet

def makeSimBHad(srcName,tagName,methodName):
  PSet = cms.PSet(
        pluginType  = cms.string("SimBHadronsFiller"),
        src         = cms.InputTag(srcName),
        tag         = cms.string(tagName),
        method      = cms.string(methodName+"()"),
        leadingOnly = cms.untracked.bool(False)
  )
  return PSet

def makeIVFBs(tagName,methodName):
  PSet = cms.PSet(
        pluginType  = cms.string("bCandidatesFiller"),
        src         = cms.InputTag('LCProducer','BCandFinalState'),
        tag         = cms.string(tagName),
        method      = cms.string(methodName+"()"),
        leadingOnly = cms.untracked.bool(False)
  )
  return PSet

def makeJetSVInfo(source = 'wCandsJets'):
    jetSVInfo = cms.PSet(
 # Jet b Tagging
    J1_TCHEbtag = makeJetBTag('TCHEbtag','trackCountingHighEffBJetTags',source)[0],
    J2_TCHEbtag = makeJetBTag('TCHEbtag','trackCountingHighEffBJetTags',source)[1],
    J3_TCHEbtag = makeJetBTag('TCHEbtag','trackCountingHighEffBJetTags',source)[2],
    J4_TCHEbtag = makeJetBTag('TCHEbtag','trackCountingHighEffBJetTags',source)[3],
    J1_TCHPbtag = makeJetBTag('TCHPbtag','trackCountingHighPurBJetTags',source)[0],
    J2_TCHPbtag = makeJetBTag('TCHPbtag','trackCountingHighPurBJetTags',source)[1],
    J3_TCHPbtag = makeJetBTag('TCHPbtag','trackCountingHighPurBJetTags',source)[2],
    J4_TCHPbtag = makeJetBTag('TCHPbtag','trackCountingHighPurBJetTags',source)[3],
    J1_CSVbtag = makeJetBTag('CSVbtag','combinedSecondaryVertexBJetTags',source)[0],
    J2_CSVbtag = makeJetBTag('CSVbtag','combinedSecondaryVertexBJetTags',source)[1],
    J3_CSVbtag = makeJetBTag('CSVbtag','combinedSecondaryVertexBJetTags',source)[2],
    J4_CSVbtag = makeJetBTag('CSVbtag','combinedSecondaryVertexBJetTags',source)[3],
    J1_CSVMVAbtag = makeJetBTag('CSVMVAbtag','combinedSecondaryVertexMVABJetTags',source)[0],
    J2_CSVMVAbtag = makeJetBTag('CSVMVAbtag','combinedSecondaryVertexMVABJetTags',source)[1],
    J3_CSVMVAbtag = makeJetBTag('CSVMVAbtag','combinedSecondaryVertexMVABJetTags',source)[2],
    J4_CSVMVAbtag = makeJetBTag('CSVMVAbtag','combinedSecondaryVertexMVABJetTags',source)[3],
 # from RecoTools/plugins/bTagSFer.h 
   J1_CSVM_SFb = makeJetUserFloat('CSVM_SFb','',source)[0],
   J2_CSVM_SFb = makeJetUserFloat('CSVM_SFb','',source)[1],
   J3_CSVM_SFb = makeJetUserFloat('CSVM_SFb','',source)[2],
   J4_CSVM_SFb = makeJetUserFloat('CSVM_SFb','',source)[3],
   J1_CSVM_errorb = makeJetUserFloat('CSVM_errorb','',source)[0],
   J2_CSVM_errorb = makeJetUserFloat('CSVM_errorb','',source)[1],
   J3_CSVM_errorb = makeJetUserFloat('CSVM_errorb','',source)[2],
   J4_CSVM_errorb = makeJetUserFloat('CSVM_errorb','',source)[3],
   J1_CSVM_SFl = makeJetUserFloat('CSVM_SFl','',source)[0],
   J2_CSVM_SFl = makeJetUserFloat('CSVM_SFl','',source)[1],
   J3_CSVM_SFl = makeJetUserFloat('CSVM_SFl','',source)[2],
   J4_CSVM_SFl = makeJetUserFloat('CSVM_SFl','',source)[3],
   J1_CSVM_SFl_down = makeJetUserFloat('CSVM_SFl_down','',source)[0],
   J2_CSVM_SFl_down = makeJetUserFloat('CSVM_SFl_down','',source)[1],
   J3_CSVM_SFl_down = makeJetUserFloat('CSVM_SFl_down','',source)[2],
   J4_CSVM_SFl_down = makeJetUserFloat('CSVM_SFl_down','',source)[3],
   J1_CSVM_SFl_up = makeJetUserFloat('CSVM_SFl_up','',source)[0],
   J2_CSVM_SFl_up = makeJetUserFloat('CSVM_SFl_up','',source)[1],
   J3_CSVM_SFl_up = makeJetUserFloat('CSVM_SFl_up','',source)[2],
   J4_CSVM_SFl_up = makeJetUserFloat('CSVM_SFl_up','',source)[3],

   J1_CSVT_SFb = makeJetUserFloat('CSVT_SFb','',source)[0],
   J2_CSVT_SFb = makeJetUserFloat('CSVT_SFb','',source)[1],
   J3_CSVT_SFb = makeJetUserFloat('CSVT_SFb','',source)[2],
   J4_CSVT_SFb = makeJetUserFloat('CSVT_SFb','',source)[3],
   J1_CSVT_errorb = makeJetUserFloat('CSVT_errorb','',source)[0],
   J2_CSVT_errorb = makeJetUserFloat('CSVT_errorb','',source)[1],
   J3_CSVT_errorb = makeJetUserFloat('CSVT_errorb','',source)[2],
   J4_CSVT_errorb = makeJetUserFloat('CSVT_errorb','',source)[3],
   J1_CSVT_SFl = makeJetUserFloat('CSVT_SFl','',source)[0],
   J2_CSVT_SFl = makeJetUserFloat('CSVT_SFl','',source)[1],
   J3_CSVT_SFl = makeJetUserFloat('CSVT_SFl','',source)[2],
   J4_CSVT_SFl = makeJetUserFloat('CSVT_SFl','',source)[3],
   J1_CSVT_SFl_down = makeJetUserFloat('CSVT_SFl_down','',source)[0],
   J2_CSVT_SFl_down = makeJetUserFloat('CSVT_SFl_down','',source)[1],
   J3_CSVT_SFl_down = makeJetUserFloat('CSVT_SFl_down','',source)[2],
   J4_CSVT_SFl_down = makeJetUserFloat('CSVT_SFl_down','',source)[3],
   J1_CSVT_SFl_up = makeJetUserFloat('CSVT_SFl_up','',source)[0],
   J2_CSVT_SFl_up = makeJetUserFloat('CSVT_SFl_up','',source)[1],
   J3_CSVT_SFl_up = makeJetUserFloat('CSVT_SFl_up','',source)[2],
   J4_CSVT_SFl_up = makeJetUserFloat('CSVT_SFl_up','',source)[3],

# from RecoTools/plugins/PATSSVJetEmbedder.h 
    J1_nTracks_SSV = makeJetUserFloat('nTracks_SSV','',source)[0],
    J2_nTracks_SSV = makeJetUserFloat('nTracks_SSV','',source)[1],
    J3_nTracks_SSV = makeJetUserFloat('nTracks_SSV','',source)[2],
    J4_nTracks_SSV = makeJetUserFloat('nTracks_SSV','',source)[3],
    J1_pt_SSV = makeJetUserFloat('pt_SSV','',source)[0],  
    J2_pt_SSV = makeJetUserFloat('pt_SSV','',source)[1],
    J3_pt_SSV = makeJetUserFloat('pt_SSV','',source)[2],
    J4_pt_SSV = makeJetUserFloat('pt_SSV','',source)[3],
    J1_eta_SSV = makeJetUserFloat('eta_SSV','',source)[0],  
    J2_eta_SSV = makeJetUserFloat('eta_SSV','',source)[1],
    J3_eta_SSV = makeJetUserFloat('eta_SSV','',source)[2],
    J4_eta_SSV = makeJetUserFloat('eta_SSV','',source)[3],
    J1_phi_SSV = makeJetUserFloat('phi_SSV','',source)[0],  
    J2_phi_SSV = makeJetUserFloat('phi_SSV','',source)[1],
    J3_phi_SSV = makeJetUserFloat('phi_SSV','',source)[2],
    J4_phi_SSV = makeJetUserFloat('phi_SSV','',source)[3],
    J1_mass_SSV = makeJetUserFloat('mass_SSV','',source)[0], # formerly J1SVMassb
    J2_mass_SSV = makeJetUserFloat('mass_SSV','',source)[1],
    J3_mass_SSV = makeJetUserFloat('mass_SSV','',source)[2],
    J4_mass_SSV = makeJetUserFloat('mass_SSV','',source)[3],
    J1_mass_SSV_alt = makeJetUserFloat('mass_SSV_alt','',source)[0],  
    J2_mass_SSV_alt = makeJetUserFloat('mass_SSV_alt','',source)[1],
    J3_mass_SSV_alt = makeJetUserFloat('mass_SSV_alt','',source)[2],
    J4_mass_SSV_alt = makeJetUserFloat('mass_SSV_alt','',source)[3],

# from RecoTools/plugins/PATCSVVertex.cc
    J1_pt_SV_weighted = makeJetUserFloat('pt_SV_weighted','',source)[0],
    J2_pt_SV_weighted = makeJetUserFloat('pt_SV_weighted','',source)[1],
    J3_pt_SV_weighted = makeJetUserFloat('pt_SV_weighted','',source)[2],
    J4_pt_SV_weighted = makeJetUserFloat('pt_SV_weighted','',source)[3],
    J1_pt_SV_unweighted = makeJetUserFloat('pt_SV_unweighted','',source)[0],
    J2_pt_SV_unweighted = makeJetUserFloat('pt_SV_unweighted','',source)[1],
    J3_pt_SV_unweighted = makeJetUserFloat('pt_SV_unweighted','',source)[2],
    J4_pt_SV_unweighted = makeJetUserFloat('pt_SV_unweighted','',source)[3],
    J1_mass_SV_weighted = makeJetUserFloat('mass_SV_weighted','',source)[0],
    J2_mass_SV_weighted = makeJetUserFloat('mass_SV_weighted','',source)[1],
    J3_mass_SV_weighted = makeJetUserFloat('mass_SV_weighted','',source)[2],
    J4_mass_SV_weighted = makeJetUserFloat('mass_SV_weighted','',source)[3],
    J1_mass_SV_unweighted = makeJetUserFloat('mass_SV_unweighted','',source)[0],
    J2_mass_SV_unweighted = makeJetUserFloat('mass_SV_unweighted','',source)[1],
    J3_mass_SV_unweighted = makeJetUserFloat('mass_SV_unweighted','',source)[2],
    J4_mass_SV_unweighted = makeJetUserFloat('mass_SV_unweighted','',source)[3],
    J1_normChi2_SV = makeJetUserFloat('normChi2_SV','',source)[0],
    J2_normChi2_SV = makeJetUserFloat('normChi2_SV','',source)[1],
    J3_normChi2_SV = makeJetUserFloat('normChi2_SV','',source)[2],
    J4_normChi2_SV = makeJetUserFloat('normChi2_SV','',source)[3],
    J1_nTracks_SV = makeJetUserFloat('nTracks_SV','',source)[0],
    J2_nTracks_SV = makeJetUserFloat('nTracks_SV','',source)[1],
    J3_nTracks_SV = makeJetUserFloat('nTracks_SV','',source)[2],
    J4_nTracks_SV = makeJetUserFloat('nTracks_SV','',source)[3],
    J1_mass_SV_corrected = makeJetUserFloat('mass_SV_corrected','',source)[0],
    J2_mass_SV_corrected = makeJetUserFloat('mass_SV_corrected','',source)[1],
    J3_mass_SV_corrected = makeJetUserFloat('mass_SV_corrected','',source)[2],
    J4_mass_SV_corrected = makeJetUserFloat('mass_SV_corrected','',source)[3],
    )
    return jetSVInfo


def makeJetIDInfo(source = 'wCandsJets'):
   jetIDInfo = cms.PSet( 
    
     # variables going in to Jet ID  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Documentation
    #J1_jetID_neutHadEnFr = makeJetStringName("jetID_neutHadEnFr","neutralHadronEnergyFraction()",xn='',source='wCandsJets')[0],
    #J1_jetID_neutEmEnFr  = makeJetStringName("jetID_neutEmEnFr","neutralEmEnergyFraction()",xn='',source='wCandsJets')[0],
    #J1_jetID_nrPFConstit = makeJetStringName("jetID_nrPFConstit","getPFConstituents().size()",xn='',source='wCandsJets')[0],
    #J1_jetID_chHadEnFr   = makeJetStringName("jetID_chHadEnFr","chargedHadronEnergyFraction()",xn='',source='wCandsJets')[0],
    #J1_jetID_chHadMult   = makeJetStringName("jetID_chHadMult","chargedHadronMultiplicity()",xn='',source='wCandsJets')[0],
    #J1_jetID_chMult      = makeJetStringName("jetID_chMult","chargedMultiplicity()",xn='',source='wCandsJets')[0],
    #J1_jetID_chEmEnFr    = makeJetStringName("jetID_chEmEnFr","chargedEmEnergyFraction()",xn='',source='wCandsJets')[0],
    #J2_jetID_neutHadEnFr = makeJetStringName("jetID_neutHadEnFr","neutralHadronEnergyFraction()",xn='',source='wCandsJets')[1],
    #J2_jetID_neutEmEnFr  = makeJetStringName("jetID_neutEmEnFr","neutralEmEnergyFraction()",xn='',source='wCandsJets')[1],
    #J2_jetID_nrPFConstit = makeJetStringName("jetID_nrPFConstit","getPFConstituents().size()",xn='',source='wCandsJets')[1],
    #J2_jetID_chHadEnFr   = makeJetStringName("jetID_chHadEnFr","chargedHadronEnergyFraction()",xn='',source='wCandsJets')[1],
    #J2_jetID_chHadMult   = makeJetStringName("jetID_chHadMult","chargedHadronMultiplicity()",xn='',source='wCandsJets')[1],
    #J2_jetID_chMult      = makeJetStringName("jetID_chMult","chargedMultiplicity()",xn='',source='wCandsJets')[1],
    #J2_jetID_chEmEnFr    = makeJetStringName("jetID_chEmEnFr","chargedEmEnergyFraction()",xn='',source='wCandsJets')[1],
    #J3_jetID_neutHadEnFr = makeJetStringName("jetID_neutHadEnFr","neutralHadronEnergyFraction()",xn='',source='wCandsJets')[2],
    #J3_jetID_neutEmEnFr  = makeJetStringName("jetID_neutEmEnFr","neutralEmEnergyFraction()",xn='',source='wCandsJets')[2],
    #J3_jetID_nrPFConstit = makeJetStringName("jetID_nrPFConstit","getPFConstituents().size()",xn='',source='wCandsJets')[2],
    #J3_jetID_chHadEnFr   = makeJetStringName("jetID_chHadEnFr","chargedHadronEnergyFraction()",xn='',source='wCandsJets')[2],
    #J3_jetID_chHadMult   = makeJetStringName("jetID_chHadMult","chargedHadronMultiplicity()",xn='',source='wCandsJets')[2],
    #J3_jetID_chMult      = makeJetStringName("jetID_chMult","chargedMultiplicity()",xn='',source='wCandsJets')[2],
    #J3_jetID_chEmEnFr    = makeJetStringName("jetID_chEmEnFr","chargedEmEnergyFraction()",xn='',source='wCandsJets')[2],
    #J4_jetID_neutHadEnFr = makeJetStringName("jetID_neutHadEnFr","neutralHadronEnergyFraction()",xn='',source='wCandsJets')[3],
    #J4_jetID_neutEmEnFr  = makeJetStringName("jetID_neutEmEnFr","neutralEmEnergyFraction()",xn='',source='wCandsJets')[3],
    #J4_jetID_nrPFConstit = makeJetStringName("jetID_nrPFConstit","getPFConstituents().size()",xn='',source='wCandsJets')[3],
    #J4_jetID_chHadEnFr   = makeJetStringName("jetID_chHadEnFr","chargedHadronEnergyFraction()",xn='',source='wCandsJets')[3],
    #J4_jetID_chHadMult   = makeJetStringName("jetID_chHadMult","chargedHadronMultiplicity()",xn='',source='wCandsJets')[3],
    #J4_jetID_chMult      = makeJetStringName("jetID_chMult","chargedMultiplicity()",xn='',source='wCandsJets')[3],
    #J4_jetID_chEmEnFr    = makeJetStringName("jetID_chEmEnFr","chargedEmEnergyFraction()",xn='',source='wCandsJets')[3],

    # Jet ID and Pu ID
    J1_idTightJet = makeJetUserFloat('idTight','Jet',source)[0], 
    J2_idTightJet = makeJetUserFloat('idTight','Jet',source)[1], 
    J3_idTightJet = makeJetUserFloat('idTight','Jet',source)[2], 
    J4_idTightJet = makeJetUserFloat('idTight','Jet',source)[3], 
    J1_idLooseJet = makeJetUserFloat('idLoose','Jet',source)[0], 
    J2_idLooseJet = makeJetUserFloat('idLoose','Jet',source)[1], 
    J3_idLooseJet = makeJetUserFloat('idLoose','Jet',source)[2], 
    J4_idLooseJet = makeJetUserFloat('idLoose','Jet',source)[3], 
    J1_fullIdTightPu = makeJetUserInt('fullIdTight','Pu',source)[0], 
    J2_fullIdTightPu = makeJetUserInt('fullIdTight','Pu',source)[1], 
    J3_fullIdTightPu = makeJetUserInt('fullIdTight','Pu',source)[2], 
    J4_fullIdTightPu = makeJetUserInt('fullIdTight','Pu',source)[3], 
    J1_fullIdLoosePu = makeJetUserInt('fullIdLoose','Pu',source)[0], 
    J2_fullIdLoosePu = makeJetUserInt('fullIdLoose','Pu',source)[1], 
    J3_fullIdLoosePu = makeJetUserInt('fullIdLoose','Pu',source)[2], 
    J4_fullIdLoosePu = makeJetUserInt('fullIdLoose','Pu',source)[3], 
    J1_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[0], 
    J2_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[1], 
    J3_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[2], 
    J4_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[3], 
    J1_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[0], 
    J2_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[1], 
    J3_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[2], 
    J4_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[3], 
    J1_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[0], 
    J2_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[1], 
    J3_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[2], 
    J4_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[3], 
    J1_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[0],
    J2_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[1],
    J3_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[2],
    J4_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[3],
    J1_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[0],
    J2_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[1],
    J3_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[2],
    J4_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[3],
    J1_idBetaStarPu = makeJetUserFloat('idBetaStar','Pu',source)[0],
    J2_idBetaStarPu = makeJetUserFloat('idBetaStar','Pu',source)[1],
    J3_idBetaStarPu = makeJetUserFloat('idBetaStar','Pu',source)[2],
    J4_idBetaStarPu = makeJetUserFloat('idBetaStar','Pu',source)[3],
    J1_idBetaStarClassicPu = makeJetUserFloat('idBetaStarClassic','Pu',source)[0],
    J2_idBetaStarClassicPu = makeJetUserFloat('idBetaStarClassic','Pu',source)[1],
    J3_idBetaStarClassicPu = makeJetUserFloat('idBetaStarClassic','Pu',source)[2],
    J4_idBetaStarClassicPu = makeJetUserFloat('idBetaStarClassic','Pu',source)[3],
    J1_idBetaStarClassicModPu = makeJetUserFloat('idBetaStarClassicMod','Pu',source)[0],
    J2_idBetaStarClassicModPu = makeJetUserFloat('idBetaStarClassicMod','Pu',source)[1],
    J3_idBetaStarClassicModPu = makeJetUserFloat('idBetaStarClassicMod','Pu',source)[2],
    J4_idBetaStarClassicModPu = makeJetUserFloat('idBetaStarClassicMod','Pu',source)[3],
    J1_id_nTrackPu = makeJetUserFloat('id_nTrack','Pu',source)[0],
    J2_id_nTrackPu = makeJetUserFloat('id_nTrack','Pu',source)[1],
    J3_id_nTrackPu = makeJetUserFloat('id_nTrack','Pu',source)[2],
    J4_id_nTrackPu = makeJetUserFloat('id_nTrack','Pu',source)[3],
    J1_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[0],
    J2_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[1],
    J3_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[2],
    J4_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[3],
    J1_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[0],
    J2_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[1],
    J3_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[2],
    J4_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[3],
   )
   return jetIDInfo

def makeJetCorrectionInfo():
    # Jets at Various Levels of Correction from RecoTools/plugins/PATJetOverloader.h
   jetCorrectionInfo = cms.PSet(
    J1_pt_L1 = makeJetUserFloat('pt_L1')[0],
    J1_eta_L1 = makeJetUserFloat('eta_L1')[0],
    J1_phi_L1 = makeJetUserFloat('phi_L1')[0],
    J1_pt_L2 = makeJetUserFloat('pt_L2')[0],
    J1_eta_L2 = makeJetUserFloat('eta_L2')[0],
    J1_phi_L2 = makeJetUserFloat('phi_L2')[0],
    J1_pt_L3 = makeJetUserFloat('pt_L3')[0],
    J1_eta_L3 = makeJetUserFloat('eta_L3')[0],
    J1_phi_L3 = makeJetUserFloat('phi_L3')[0],
    J1_pt_L23 = makeJetUserFloat('pt_L23')[0],
    J2_pt_L1 = makeJetUserFloat('pt_L1')[1],
    J2_eta_L1 = makeJetUserFloat('eta_L1')[1],
    J2_phi_L1 = makeJetUserFloat('phi_L1')[1],
    J2_pt_L2 = makeJetUserFloat('pt_L2')[1],
    J2_eta_L2 = makeJetUserFloat('eta_L2')[1],
    J2_phi_L2 = makeJetUserFloat('phi_L2')[1],
    J2_pt_L3 = makeJetUserFloat('pt_L3')[1],
    J2_eta_L3 = makeJetUserFloat('eta_L3')[1],
    J2_phi_L3 = makeJetUserFloat('phi_L3')[1],
    J2_pt_L23 = makeJetUserFloat('pt_L23')[1],
    J3_pt_L1 = makeJetUserFloat('pt_L1')[2],
    J3_eta_L1 = makeJetUserFloat('eta_L1')[2],
    J3_phi_L1 = makeJetUserFloat('phi_L1')[2],
    J3_pt_L2 = makeJetUserFloat('pt_L2')[2],
    J3_eta_L2 = makeJetUserFloat('eta_L2')[2],
    J3_phi_L2 = makeJetUserFloat('phi_L2')[2],
    J3_pt_L3 = makeJetUserFloat('pt_L3')[2],
    J3_eta_L3 = makeJetUserFloat('eta_L3')[2],
    J3_phi_L3 = makeJetUserFloat('phi_L3')[2],
    J3_pt_L23 = makeJetUserFloat('pt_L23')[2],
    J4_pt_L1 = makeJetUserFloat('pt_L1')[3],
    J4_eta_L1 = makeJetUserFloat('eta_L1')[3],
    J4_phi_L1 = makeJetUserFloat('phi_L1')[3],
    J4_pt_L2 = makeJetUserFloat('pt_L2')[3],
    J4_eta_L2 = makeJetUserFloat('eta_L2')[3],
    J4_phi_L2 = makeJetUserFloat('phi_L2')[3],
    J4_pt_L3 = makeJetUserFloat('pt_L3')[3],
    J4_eta_L3 = makeJetUserFloat('eta_L3')[3],
    J4_phi_L3 = makeJetUserFloat('phi_L3')[3],
    J4_pt_L23 = makeJetUserFloat('pt_L23')[3],
   )
   return jetCorrectionInfo


def makeCollections(source = 'wCandsJets', sourceZ = 'diMuonsSorted',sourceE = 'weCandsJets', sourceEle="selectedPatElectrons"):
 commonCollections = cms.PSet(         
    #electronPt = makeMuNu("electronPt","lepton.pt()",sourceE),
    PVs = cms.PSet(
        pluginType = cms.string("VertexSizeFiller"),
        src = cms.InputTag("primaryVertexFilter"),
        tag = cms.string("vertices")
    ), 
    pu = cms.PSet(
        pluginType = cms.string("PUFiller"),
        src = cms.InputTag("addPileupInfo"),
        tag = cms.string("pu"),
        ),
    trigger = cms.PSet(
        pluginType = cms.string("TriggerFiller"),
        src = cms.InputTag("patTrigger"),
        paths = cms.vstring(TriggerPaths)
    ),


    J1_flightDistance = makeJetUserFloat('flightDistance','',source)[0],
    J2_flightDistance = makeJetUserFloat('flightDistance','',source)[1],
    J1_errorFlightDistance = makeJetUserFloat('errorFlightDistance','',source)[0],
    J2_errorFlightDistance = makeJetUserFloat('errorFlightDistance','',source)[1],
    J1_DR = makeJetUserFloat('DR','',source)[0],
    J2_DR = makeJetUserFloat('DR','',source)[1],
    J1_ptRMS = makeJetUserFloat('ptRMS','',source)[0],
    J2_ptRMS = makeJetUserFloat('ptRMS','',source)[1],
    J1_dxy = makeJetUserFloat('dxy','_track',source)[0],
    J2_dxy = makeJetUserFloat('dxy','_track',source)[1],
    J1_dz = makeJetUserFloat('dz','_track',source)[0],
    J2_dz = makeJetUserFloat('dz','_track',source)[1],
    J1_pt_PV = makeJetUserFloat('pt_PV','',source)[0],  
    J2_pt_PV = makeJetUserFloat('pt_PV','',source)[1],
    J3_pt_PV = makeJetUserFloat('pt_PV','',source)[2],
    J4_pt_PV = makeJetUserFloat('pt_PV','',source)[3],
    J1_eta_PV = makeJetUserFloat('eta_PV','',source)[0],  
    J2_eta_PV = makeJetUserFloat('eta_PV','',source)[1],
    J3_eta_PV = makeJetUserFloat('eta_PV','',source)[2],
    J4_eta_PV = makeJetUserFloat('eta_PV','',source)[3],
    J1_phi_PV = makeJetUserFloat('phi_PV','',source)[0],  
    J2_phi_PV = makeJetUserFloat('phi_PV','',source)[1],
    J3_phi_PV = makeJetUserFloat('phi_PV','',source)[2],
    J4_phi_PV = makeJetUserFloat('phi_PV','',source)[3],
    J1_mass_PV = makeJetUserFloat('mass_PV','',source)[0], 
    J2_mass_PV = makeJetUserFloat('mass_PV','',source)[1],
    J3_mass_PV = makeJetUserFloat('mass_PV','',source)[2],
    J4_mass_PV = makeJetUserFloat('mass_PV','',source)[3],
    J1_partonFlavour = makeJetStringPar('partonFlavour','',source)[0],
    J2_partonFlavour = makeJetStringPar('partonFlavour','',source)[1],
    J1_muonMultiplicity = makeJetStringPar('muonMultiplicity','',source)[0],
    J2_muonMultiplicity = makeJetStringPar('muonMultiplicity','',source)[1],
    J1_chargedMultiplicity = makeJetStringPar('chargedMultiplicity','',source)[0],
    J2_chargedMultiplicity = makeJetStringPar('chargedMultiplicity','',source)[1],
    J1_electronMultiplicity = makeJetStringPar('electronMultiplicity','',source)[0],
    J2_electronMultiplicity = makeJetStringPar('electronMultiplicity','',source)[1],
    J1_photonMultiplicity = makeJetStringPar('photonMultiplicity','',source)[0],
    J2_photonMultiplicity = makeJetStringPar('photonMultiplicity','',source)[1],
    J1_jetCharge = makeJetString('jetCharge','',source)[0],
    J2_jetCharge = makeJetString('jetCharge','',source)[1],
    J1_mass = makeJetString('mass','',source)[0],
    J2_mass = makeJetString('mass','',source)[1],
    J1_photonEnergy = makeJetString('photonEnergy','',source)[0],
    J2_photonEnergy = makeJetString('photonEnergy','',source)[1],
    J1_photonEnergyFraction = makeJetString('photonEnergyFraction','',source)[0],
    J2_photonEnergyFraction = makeJetString('photonEnergyFraction','',source)[1],
    J1_electronEnergy = makeJetString('electronEnergy','',source)[0],
    J2_electronEnergy = makeJetString('electronEnergy','',source)[1],
    J1_chargedEmEnergyFraction = makeJetString('chargedEmEnergyFraction','',source)[0],
    J2_chargedEmEnergyFraction = makeJetString('chargedEmEnergyFraction','',source)[1],
    J1_chargedMuEnergy = makeJetString('chargedMuEnergy','',source)[0],
    J2_chargedMuEnergy = makeJetString('chargedMuEnergy','',source)[1],
    J1_muonEnergyFraction = makeJetString('muonEnergyFraction','',source)[0],
    J2_muonEnergyFraction = makeJetString('muonEnergyFraction','',source)[1],
    J1_chargedHadronMultiplicity = makeJetStringPar('chargedHadronMultiplicity','',source)[0],
    J2_chargedHadronMultiplicity = makeJetStringPar('chargedHadronMultiplicity','',source)[1],
# Jet Kinematic Variables
    J1_pt = makeJetStringPar('pt','',source)[0], # formerly highestJetPt
    J2_pt = makeJetStringPar('pt','',source)[1],
    J3_pt = makeJetStringPar('pt','',source)[2],
    J4_pt = makeJetStringPar('pt','',source)[3],
    J1_phi = makeJetStringPar('phi','',source)[0],
    J2_phi = makeJetStringPar('phi','',source)[1],
    J3_phi = makeJetStringPar('phi','',source)[2],
    J4_phi = makeJetStringPar('phi','',source)[3],
    J1_eta = makeJetStringPar('eta','',source)[0],
    J2_eta = makeJetStringPar('eta','',source)[1],
    J3_eta = makeJetStringPar('eta','',source)[2],
    J4_eta = makeJetStringPar('eta','',source)[3],

    nJetsPt20 = makeNJets(20,source),
    nJetsPt25 = makeNJets(25,source),
    nJetsPt30 = makeNJets(30,source),
    nJetsPt35 = makeNJets(35,source),
    nJetsPt40 = makeNJets(40,source),

    nJets24to5 = makeMuNu("nJets24to5","nJets24to5",source),
    nJets24Pt20 = makeMuNu("nJets24Pt20","nJets24Pt20",source),
    nJets24Pt25 = makeMuNu("nJets24Pt25","nJets24Pt25",source),
    nJets24Pt30 = makeMuNu("nJets24Pt30","nJets24Pt30",source),

    nrCbar = makeCollSize('cbarCands','nrCbar'), 
    nrC = makeCollSize('cCands','nrC'), 
    nrW = makeCollSize(source,'nrW'), 
    nrMu = makeCollSize('selectedPatMuons','nrMu'),
    nrEle = makeCollSize('selectedPatElectrons','nrEle'),
    nrMuLoose = makeCollSize('preselectedPatMuons','nrMuLoose'),
    nrEleLoose = makeCollSize('preselectedPatElectrons','nrEleLoose'),

# Few Electron Variables for MuEle control region
# BIG WARNING: This ID is super outdated. We need to go over the new EGamma ID and actualize it. WP80 is similar to the new "medium" one
    ptEle = makeBasicEle("ptEle","pt",sourceEle),
    phiEle = makeBasicEle("phiEle","phi",sourceEle),
    etaEle = makeBasicEle("etaEle","eta",sourceEle),
    wp80Ele = makeBasicEle("wp80Ele","userFloat('wp80')",sourceEle),
    chargeEle = makeBasicEle("chargeEle","charge",sourceEle),
    isoEleDB = makeBasicEle("isoEleDB","(userIso(0)+max(userIso(1)+neutralHadronIso()-0.5*userIso(2),0.0))/pt",sourceEle),

# Z Variables
    DiMuonMass = makeZColl("DiMuonMass","mass",sourceZ),
    DiMuonPt = makeZColl("DiMuonPt","mass",sourceZ),
    mu1_pt = makeZColl("mu1_pt","leg1.pt()",sourceZ),
    mu2_pt = makeZColl("mu2_pt","leg2.pt()",sourceZ),
    mu1_phi = makeZColl("mu1_phi","leg1.phi()",sourceZ),
    mu2_phi = makeZColl("mu2_phi","leg2.phi()",sourceZ),
    mu1_eta = makeZColl("mu1_eta","leg1.eta()",sourceZ),
    mu2_eta = makeZColl("mu2_eta","leg2.eta()",sourceZ),
    l1StdRelIso = makeZColl("l1StdRelIso",
     "(leg1.isolationR03.sumPt+leg1.isolationR03.emEt+leg1.isolationR03.hadEt)/leg1.pt()",sourceZ),
    l1PfIsoDB = makeZColl("l1PfIsoDB",
     "(leg1.userIso(0)+max(leg1.photonIso()+leg1.neutralHadronIso()-0.5*leg1.puChargedHadronIso,0.0))/leg1.pt()",sourceZ),
    l2StdRelIso = makeZColl("l2StdRelIso",
     "(leg2.isolationR03.sumPt+leg2.isolationR03.emEt+leg2.isolationR03.hadEt)/leg2.pt()",sourceZ),
    l2PfIsoDB = makeZColl("l2PfIsoDB",
     "(leg2.userIso(0)+max(leg2.photonIso()+leg2.neutralHadronIso()-0.5*leg2.puChargedHadronIso,0.0))/leg2.pt()",sourceZ),

    mJJ = makeMuNu("mJJ","mJJ",source,True),
    mJ3J4 = makeMuNu("mJ3J4","mJJ2",source,True),
    ptJJ = makeMuNu("ptJJ","ptJJ",source,True),

    muon1_pt = makeMuNu("muon1_pt","lepton().pt",source),
    muon1_eta = makeMuNu("muon1_eta","lepton.eta",source),
    muon1_phi = makeMuNu("muon1_phi","lepton.phi",source),
    muon1_charge = makeMuNu("muon1_charge","lepton.charge()",source),
    Wpt = makeMuNu("Wpt","corPt()",source),

    muon_globalMuon = makeMuNu("gloablMuon","isGlobalMuon()",source="wCandsJets"),
    #muon_PFMuon = makeMuNu("PFMuon","isPFMuon()",source="wCandsJets"),
    #muon_nChi2 = makeMuNu("nChi2","globalTrack()->normalizedChi2()",source="wCandsJets"),
    #muon_nValidHits = makeMuNu("nValidHits","globalTrack()->hitPattern().numberOfValidMuonHits()",source="wCandsJets"),
    #muon_nMatchedStations = makeMuNu("nMatchedStations","numberOfMatchedStations()",source="wCandsJets"),
    #muon_dxy = makeMuNu("dxy","muonBestTrack()->dxy(vertex->position())",source="wCandsJets"),
    #muon_dz = makeMuNu("dz","muonBestTrack()->dz(vertex->position())",source="wCandsJets"),
    #muon_nPixelHits = makeMuNu("nPixelHits","innerTrack()->hitPattern().numberOfValidPixelHits()",source="wCandsJets"),
    #muon_trackLayers = makeMuNu("trackLayers","track()->hitPattern().trackerLayersWithMeasurement()",source="wCandsJets"),

    met_pt = makeMuNu("met_pt","met().pt",source,True),
    met_eta = makeMuNu("met_eta","met().eta",source,True),
    met_phi = makeMuNu("met_phi","met().phi",source,True),
    mt = makeMuNu("mt","mt",source),
    
    metJJ = makeMuNu("metjj","metjj",source),
    leptonjj = makeMuNu("leptonjj","leptonjj",source),
    muNuDPhi = makeMuNu("muNuDPhi","dPhi",source),
    muNuRecoil = makeMuNu("muNuRecoil","recoil().pt()",source),
    muNuRelPFIso = makeMuNu("muNuRelPFIso",
     "(lepton.chargedHadronIso()+lepton.photonIso()+lepton.neutralHadronIso())/lepton.pt()",source),
#    PFIsoVeto = makeMuNu("PFIsoVeto","lepton.userIso(0)",source),
    PFIsoRho = makeMuNu("PFIsoRho","lepton.userFloat('rho')",source),
    muNuRelStdIso03 = makeMuNu("muNuRelStdIso03",
     "(lepton.isolationR03.sumPt+lepton.isolationR03.emEt+lepton.isolationR03.hadEt)/lepton.pt()",source),
    muNuRelPFIsoDB = makeMuNu("muNuRelPFIsoDB",
     "(lepton.userIso(0)+max(lepton.photonIso()+lepton.neutralHadronIso()-0.5*lepton.puChargedHadronIso,0.0))/lepton.pt()",
     source,True),
    ipDXY = makeMuNu("ipDXY","lepton.userFloat('ipDXY')",source,True),
     #"(lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt +lepton.pfIsolationR04().sumPhotonEt - 0.5*lepton.pfIsolationR04().sumPUPt,0.0))/lepton.pt()<0.12" ,
    muNuRelPFIsoDBAndrea = makeMuNu("muNuRelPFIsoDBAndrea",
      "((lepton.pfIsolationR04().sumChargedHadronPt + max(lepton.pfIsolationR04().sumNeutralHadronEt +lepton.pfIsolationR04().sumPhotonEt- 0.5*lepton.pfIsolationR04().sumPUPt,0))/lepton.pt)",
      source,True),
    muNuRelPFIsoDBAndrea_A= makeMuNu("muNuRelPFIsoDBAndrea_A",
      "lepton.pfIsolationR04().sumChargedHadronPt",source,True),
    muNuRelPFIsoDB_A = makeMuNu("muNuRelPFIsoDB_A","lepton.userIso(0)",source, True),
    muNuRelPFIsoDBAndrea_B= makeMuNu("muNuRelPFIsoDBAndrea_B",
      "lepton.pfIsolationR04().sumNeutralHadronEt +lepton.pfIsolationR04().sumPhotonEt",
      source,True),
    muNuRelPFIsoDB_B = makeMuNu("muNuRelPFIsoDB_B",
      "lepton.photonIso()+lepton.neutralHadronIso()",source,True),
    muNuRelPFIsoDBAndrea_C= makeMuNu("muNuRelPFIsoDBAndrea_C",
      "lepton.pfIsolationR04().sumPUPt",source,True),
    muNuRelPFIsoDB_C = makeMuNu("muNuRelPFIsoDB_C","lepton.puChargedHadronIso",source,True),

    dz = makeMuNu("dz",'abs(lepton.userFloat("dz"))',source,True),
    ht = makeMuNu("ht","ht",source,True),
 )
 return commonCollections

def addMuNuEventTreePtData(process,name,source = 'wCandsJets',sourceZ = 'diMuonsSorted'):
   process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
   eventTree = cms.EDAnalyzer('EventTreeMaker',
      makeCollections(source,sourceZ),
      makeJetIDInfo(source = 'wCandsJets'),
      makeJetCorrectionInfo(),
      makeJetSVInfo(source='wCandsJets'),
      coreCollections = cms.VInputTag( cms.InputTag(source) ),
      J1_eta_L23 = makeJetUserFloat('eta_L23')[0],
      J1_phi_L23 = makeJetUserFloat('phi_L23')[0],
      J2_eta_L23 = makeJetUserFloat('eta_L23')[1],
      J2_phi_L23 = makeJetUserFloat('phi_L23')[1],
      J3_eta_L23 = makeJetUserFloat('eta_L23')[2],
      J3_phi_L23 = makeJetUserFloat('phi_L23')[2],
      J4_eta_L23 = makeJetUserFloat('eta_L23')[3],
      J4_phi_L23 = makeJetUserFloat('phi_L23')[3]
   )
   setattr(process, name, eventTree)
   p = cms.Path(getattr(process,name))
   setattr(process, name+'Path', p)

def addMuNuEventTreePtMC(process,name,source = 'wCandsJets',sourceZ = 'diMuonsSorted',lhep="externalLHEProducer"):
   process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
   eventTree = cms.EDAnalyzer('EventTreeMaker',
      makeCollections(source,sourceZ),
      makeJetIDInfo(source = 'wCandsJets'),
      makeJetCorrectionInfo(),
      makeJetSVInfo(source='wCandsJets'),
      coreCollections = cms.VInputTag(
           cms.InputTag(source)
      ),
      topweight= cms.PSet(
           pluginType = cms.string("TopWeight"),
           src = cms.InputTag("genParticles")
      ),
      muNuLHEProduct = cms.PSet(
          pluginType = cms.string("LHEProductFiller"),
          src = cms.InputTag(lhep), 
          tag = cms.string("LHEProduct"),
      ),
      bHadronsPt = makeSimBHad("bhadrons","bHadronsPt","pt"),
      bHadronsEta = makeSimBHad("bhadrons","bHadronsEta","eta"),
      bHadronsPhi = makeSimBHad("bhadrons","bHadronsPhi","phi"),
  
      # Jet pT Reco + Gen from RecoTools/plugins/PATJetSmearer.h
      J1_pt_gen = makeJetUserFloat('pt_gen')[0],
      J1_pt_gen_two = makeJetUserFloat('pt_gen_two')[0],
      J1_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[0], 
      J1_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[0], 
      J1_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[0], 
      J2_pt_gen = makeJetUserFloat('pt_gen')[1],
      J2_pt_gen_two = makeJetUserFloat('pt_gen_two')[1],
      J2_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[1], 
      J2_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[1], 
      J2_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[1], 
      J3_pt_gen = makeJetUserFloat('pt_gen')[2],
      J3_pt_gen_two = makeJetUserFloat('pt_gen_two')[2],
      J3_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[2], 
      J3_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[2], 
      J3_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[2], 
      J4_pt_gen = makeJetUserFloat('pt_gen')[3],
      J4_pt_gen_two = makeJetUserFloat('pt_gen_two')[3],
      J4_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[3], 
      J4_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[3], 
      J4_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[3], 

      bCandidate1Pt = makeIVFBs("bCandidate1Pt","BC1PT"),
      bCandidate2Pt = makeIVFBs("bCandidate2Pt","BC2PT"),
      bCandidate1Eta = makeIVFBs("bCandidate1Eta","BC1ETA"),
      bCandidate2Eta = makeIVFBs("bCandidate2Eta","BC2ETA"),
      bCandidate1Phi = makeIVFBs("bCandidate1Phi","BC1PHI"),
      bCandidate2Phi = makeIVFBs("bCandidate2Phi","BC2PHI"),
      bCandidateDeltaR = makeIVFBs("bCandidateDeltaR","BCDeltaR"),
      bCandidateDeltaPhi = makeIVFBs("bCandidateDeltaPhi","BDeltaPHI"),

      nbHadrons = makeCollSize("bhadrons","nbHadrons"),
      genTs = makeCollSize("gentCands","genTs"),
      genTbars = makeCollSize("gentbarCands","genTbars"),
      genBs = makeCollSize("genbbCands","genBs"),
      genCs = makeCollSize("genccCands","genCs"),
      genDs = makeCollSize("genddCands","genDs"),
      genUs = makeCollSize("genuuCands","genUs"),
      genSs = makeCollSize("genSSCands","genSs"),
      genWs = makeCollSize("genWs","genWs"),
 
      EffWEIGHTCSVL = makeMuNu("EffWEIGHTCSVL","SFCSVL1",source,True),
      EffWEIGHTCSVL2 = makeMuNu("EffWEIGHTCSVL2","SFCSVL2",source,True),
      EffWEIGHTCSVM = makeMuNu("EffWEIGHTCSVM","SFCSVM1",source,True),
      EffWEIGHTCSVM2 = makeMuNu("EffWEIGHTCSVM2","SFCSVM2",source,True),
      EffWEIGHTCSVT = makeMuNu("EffWEIGHTCSVT","SFCSVT1",source,True),
      EffWEIGHTCSVT2 = makeMuNu("EffWEIGHTCSVT2","SFCSVT2",source,True),
      EffWEIGHTSSVHEM = makeMuNu("EffWEIGHTSSVHEM","SFSSVHE1",source,True),
      EffWEIGHTSSVHEM2 = makeMuNu("EffWEIGHTSSVHEM2","SFSSVHE2",source,True),
      weightEtaMuonIso = makeMuNu("weightEtaMuonIso","EffWEIGHTeta_IS",source,True),
      weightEtaMuonID = makeMuNu("weightEtaMuonID","EffWEIGHTeta_ID",source,True),
      weightEtaMuonTrig = makeMuNu("weightEtaMuonTrig","EffWEIGHTeta_TR",source,True),
   )
   setattr(process, name, eventTree)
   p = cms.Path(getattr(process,name))
   setattr(process, name+'Path', p)
