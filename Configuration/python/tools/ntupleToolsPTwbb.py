import FWCore.ParameterSet.Config as cms

from UWAnalysis.Configuration.tools.analysisToolsPT import TriggerPaths
jetRanks = [0,1,2,3,4,5,6,7]
jetNames = ['J1_','J2_','J3_','J4_','J5_','J6_','J7_','J8_']

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

def makeElectron(tagName,methodName,sourceElectrons='smearedElectrons',lo=True):
  PSet = cms.PSet(
        pluginType  = cms.string("PATElectronFiller"),
        src         = cms.InputTag(sourceElectrons),
        tag         = cms.string(tagName),
        method      = cms.string(methodName),
        leadingOnly = cms.untracked.bool(lo)
  )
  return PSet

def makeMuon(tagName,methodName,sourceMuons='selectedPatMuons',lo=True):
   PSet = cms.PSet(
         pluginType  = cms.string("PATMuonFiller"),
         src         = cms.InputTag(sourceMuons),
         tag         = cms.string(tagName),
         method      = cms.string(methodName),
         leadingOnly = cms.untracked.bool(lo)
   )
   return PSet

def makeMET(tagName,methodName,sourceMET='smearedMET',lo=True):
  PSet = cms.PSet(
         pluginType  = cms.string("PATMETFiller"),
         src         = cms.InputTag(sourceMET),
         tag         = cms.string(tagName),
         method      = cms.string(methodName),
         leadingOnly = cms.untracked.bool(lo)
  )
  return PSet

#jetRanks = [0,1,2,3,4,5,6,7]
#jetNames = ['J1_','J2_','J3_','J4_','J5_','J6_','J7_','J8_']
#
#def makeJetStringPar(strName,xn='',source='wCandsJets'):
# PSet_List = []
# for rank,name in zip(jetRanks,jetNames):
#  nameTag = name+strName+xn
#  PSet_List.append(cms.PSet(
#        pluginType = cms.string("PATMuonNuPairPtJetVarFiller"),
#        src        = cms.InputTag(source),
#        tag        = cms.string(nameTag),
#        method     = cms.string(strName+'()'),
#        rank       = cms.untracked.double(rank)
#  ))
# return PSet_List


def makeJetList(strName,methodName,xn='',sourceJets='smearedJets'):
   PSet_List = []
   for rank,name in zip(jetRanks,jetNames):
    nameTag = name+strName+xn
    PSet_List.append(cms.PSet(
          pluginType  = cms.string("PATJetFiller"),
          src         = cms.InputTag(sourceJets),
          tag         = cms.string(nameTag),
          method      = cms.string(methodName),
          leadingOnly = cms.untracked.bool(False),
          myrank      = cms.untracked.double(rank)
    ))
   return PSet_List

def makeJet(tagName,methodName,sourceJets='smearedJets',lo=True):
   PSet = cms.PSet(
         pluginType  = cms.string("PATJetFiller"),
         src         = cms.InputTag(sourceJets),
         tag         = cms.string(tagName),
         method      = cms.string(methodName),
         leadingOnly = cms.untracked.bool(lo)
   )
   return PSet

def makeEleNu(tagName,methodName,source='wCandsJets',lo=True):
  PSet = cms.PSet(
        pluginType  = cms.string("PATElectronNuPairFiller"),
        src         = cms.InputTag(source),
        tag         = cms.string(tagName),
        method      = cms.string(methodName),
        leadingOnly = cms.untracked.bool(lo)
  )
  return PSet

def makeMuNu(tagName,methodName,source='wCandsJets',lo=True):
  PSet = cms.PSet(
        pluginType  = cms.string("PATMuonNuPairFiller"),
        src         = cms.InputTag(source),
        tag         = cms.string(tagName),
        method      = cms.string(methodName),
        leadingOnly = cms.untracked.bool(lo)
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
        leadingOnly = cms.untracked.bool(True)
  )
  return PSet

def makeJetsLegit(sourceJets = 'smearedJets'):
  theJets = cms.PSet(
    Jet_legit_pt = makeJet('Jet_legit_pt','pt',sourceJets,False),
    Jet_legit_eta = makeJet('Jet_legit_eta','eta',sourceJets,False),
    Jet_legit_phi = makeJet('Jet_legit_phi','phi',sourceJets,False),
    Jet_legit_ID = makeJet('Jet_legit_ID','userFloat("idLoose")',sourceJets,False),
    Jet_legit_CSV = makeJet('Jet_legit_CSV','bDiscriminator("combinedSecondaryVertexBJetTags")',sourceJets,False),
  )
  return theJets

def makeMETsInd(sourceMETs = 'smearedMET'):
   theMETsInd = cms.PSet(
    met_pt_Ind = makeMET("met_pt_Ind","pt()",sourceMET=sourceMETs,lo=False),
    met_eta_Ind = makeMET("met_eta_Ind","eta()",sourceMET=sourceMETs,lo=False),
   )
   return theMETsInd


def makeMETsColl(sourceMETs = 'wCandsSelEle'):
   theMETsColl = cms.PSet(
    met_pt_coll = makeEleNu("met_pt_coll","met().pt",sourceMETs,True),
    met_eta_coll = makeEleNu("met_eta_coll","met().eta",sourceMETs,True),
   )
   return theMETsColl


def makeJetsInd(sourceJets = 'smearedJets'):
  theJetsInd = cms.PSet(
    J1_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[0],
    J2_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[1],
    J3_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[2],
    J4_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[3],
    J5_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[4],
    J6_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[5],
    J7_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[6],
    J8_pt_Ind = makeJetList('pt','pt()','_Ind',sourceJets)[7],
    J1_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[0],
    J2_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[1],
    J3_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[2],
    J4_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[3],
    J5_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[4],
    J6_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[5],
    J7_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[6],
    J8_eta_Ind = makeJetList('eta','eta()','_Ind',sourceJets)[7],
    J1_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[0],
    J2_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[1],
    J3_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[2],
    J4_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[3],
    J5_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[4],
    J6_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[5],
    J7_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[6],
    J8_phi_Ind = makeJetList('phi','phi()','_Ind',sourceJets)[7],
    J1_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[0],
    J2_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[1],
    J3_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[2],
    J4_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[3],
    J5_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[4],
    J6_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[5],
    J7_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[6],
    J8_idLoose_Ind = makeJetList('idLoose','userFloat("idLoose")','_Ind',sourceJets)[7],
    J1_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[0],
    J2_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[1],
    J3_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[2],
    J4_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[3],
    J5_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[4],
    J6_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[5],
    J7_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[6],
    J8_CSV_Ind = makeJetList('CSV','bDiscriminator("combinedSecondaryVertexBJetTags")','_Ind',sourceJets)[7],
  )
  return theJetsInd


def makeGenBWeight(tagName,methodName):
  PSet = cms.PSet(
        pluginType  = cms.string("GenBDWeightFiller"),
        src         = cms.InputTag('GBWeightProducer','GenBDWeight'),
        tag         = cms.string(tagName),
        method      = cms.string(methodName+"()"),
        leadingOnly = cms.untracked.bool(True)
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

  # from RecoTools/plugins/PATJetCSVreweight.h
   J1_CSVreweight = makeJetUserFloat('CSVreweight','',source)[0],
   J2_CSVreweight = makeJetUserFloat('CSVreweight','',source)[1],
   J3_CSVreweight = makeJetUserFloat('CSVreweight','',source)[2],
   J4_CSVreweight = makeJetUserFloat('CSVreweight','',source)[3],

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
    #J1_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[0], 
    #J2_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[1], 
    #J3_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[2], 
    #J4_fullDiscriminantPu = makeJetUserFloat('fullDiscriminant','Pu',source)[3], 
    #J1_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[0], 
    #J2_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[1], 
    #J3_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[2], 
    #J4_philv1DiscriminantPu = makeJetUserFloat('philv1Discriminant','Pu',source)[3], 
    #J1_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[0], 
    #J2_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[1], 
    #J3_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[2], 
    #J4_simpleDiscriminantPu = makeJetUserFloat('simpleDiscriminant','Pu',source)[3], 
    J1_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[0],
    J2_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[1],
    J3_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[2],
    J4_idBetaPu = makeJetUserFloat('idBeta','Pu',source)[3],
    #J1_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[0],
    #J2_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[1],
    #J3_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[2],
    #J4_idBetaClassicPu = makeJetUserFloat('idBetaClassic','Pu',source)[3],
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
    #J1_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[0],
    #J2_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[1],
    #J3_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[2],
    #J4_idClosestDzPu = makeJetUserFloat('idClosestDz','Pu',source)[3],
    #J1_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[0],
    #J2_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[1],
    #J3_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[2],
    #J4_idClosestDxyPu = makeJetUserFloat('idClosestDxy','Pu',source)[3],
   )
   return jetIDInfo

def makeJetCorrectionInfo(source = 'wCandsJets'):
    # Jets at Various Levels of Correction from RecoTools/plugins/PATJetOverloader.h
   jetCorrectionInfo = cms.PSet(
    J1_pt_L1  = makeJetUserFloat('pt_L1' ,'',source)[0],
    J1_eta_L1 = makeJetUserFloat('eta_L1','',source)[0],
    J1_phi_L1 = makeJetUserFloat('phi_L1','',source)[0],
    J1_pt_L2  = makeJetUserFloat('pt_L2' ,'',source)[0],
    J1_eta_L2 = makeJetUserFloat('eta_L2','',source)[0],
    J1_phi_L2 = makeJetUserFloat('phi_L2','',source)[0],
    J1_pt_L3  = makeJetUserFloat('pt_L3' ,'',source)[0],
    J1_eta_L3 = makeJetUserFloat('eta_L3','',source)[0],
    J1_phi_L3 = makeJetUserFloat('phi_L3','',source)[0],
    J2_pt_L1  = makeJetUserFloat('pt_L1' ,'',source)[1],
    J2_eta_L1 = makeJetUserFloat('eta_L1','',source)[1],
    J2_phi_L1 = makeJetUserFloat('phi_L1','',source)[1],
    J2_pt_L2  = makeJetUserFloat('pt_L2' ,'',source)[1],
    J2_eta_L2 = makeJetUserFloat('eta_L2','',source)[1],
    J2_phi_L2 = makeJetUserFloat('phi_L2','',source)[1],
    J2_pt_L3  = makeJetUserFloat('pt_L3' ,'',source)[1],
    J2_eta_L3 = makeJetUserFloat('eta_L3','',source)[1],
    J2_phi_L3 = makeJetUserFloat('phi_L3','',source)[1],
    J3_pt_L1  = makeJetUserFloat('pt_L1' ,'',source)[2],
    J3_eta_L1 = makeJetUserFloat('eta_L1','',source)[2],
    J3_phi_L1 = makeJetUserFloat('phi_L1','',source)[2],
    J3_pt_L2  = makeJetUserFloat('pt_L2' ,'',source)[2],
    J3_eta_L2 = makeJetUserFloat('eta_L2','',source)[2],
    J3_phi_L2 = makeJetUserFloat('phi_L2','',source)[2],
    J3_pt_L3  = makeJetUserFloat('pt_L3' ,'',source)[2],
    J3_eta_L3 = makeJetUserFloat('eta_L3','',source)[2],
    J3_phi_L3 = makeJetUserFloat('phi_L3','',source)[2],
    J4_pt_L1  = makeJetUserFloat('pt_L1' ,'',source)[3],
    J4_eta_L1 = makeJetUserFloat('eta_L1','',source)[3],
    J4_phi_L1 = makeJetUserFloat('phi_L1','',source)[3],
    J4_pt_L2  = makeJetUserFloat('pt_L2' ,'',source)[3],
    J4_eta_L2 = makeJetUserFloat('eta_L2','',source)[3],
    J4_phi_L2 = makeJetUserFloat('phi_L2','',source)[3],
    J4_pt_L3  = makeJetUserFloat('pt_L3' ,'',source)[3],
    J4_eta_L3 = makeJetUserFloat('eta_L3','',source)[3],
    J4_phi_L3 = makeJetUserFloat('phi_L3','',source)[3],
   )
   return jetCorrectionInfo

def makeGenJetInfo(source='wCandJets'):
    # Gen Jet information embedded from RecoTools/plugins/PATJetSmearer.h 
   jetGenInfo = cms.PSet(
    J1_pt_gen_NoNu = makeJetUserFloat('pt_gen_NoNu','',source)[0],
    J2_pt_gen_NoNu = makeJetUserFloat('pt_gen_NoNu','',source)[1],
    J3_pt_gen_NoNu = makeJetUserFloat('pt_gen_NoNu','',source)[2],
    J4_pt_gen_NoNu = makeJetUserFloat('pt_gen_NoNu','',source)[3],
    J1_eta_gen_NoNu = makeJetUserFloat('eta_gen_NoNu','',source)[0],
    J2_eta_gen_NoNu = makeJetUserFloat('eta_gen_NoNu','',source)[1],
    J3_eta_gen_NoNu = makeJetUserFloat('eta_gen_NoNu','',source)[2],
    J4_eta_gen_NoNu = makeJetUserFloat('eta_gen_NoNu','',source)[3],
    J1_phi_gen_NoNu = makeJetUserFloat('phi_gen_NoNu','',source)[0],
    J2_phi_gen_NoNu = makeJetUserFloat('phi_gen_NoNu','',source)[1],
    J3_phi_gen_NoNu = makeJetUserFloat('phi_gen_NoNu','',source)[2],
    J4_phi_gen_NoNu = makeJetUserFloat('phi_gen_NoNu','',source)[3],
    J1_pt_gen_Nu = makeJetUserFloat('pt_gen_Nu','',source)[0],
    J2_pt_gen_Nu = makeJetUserFloat('pt_gen_Nu','',source)[1],
    J3_pt_gen_Nu = makeJetUserFloat('pt_gen_Nu','',source)[2],
    J4_pt_gen_Nu = makeJetUserFloat('pt_gen_Nu','',source)[3],
    J1_eta_gen_Nu = makeJetUserFloat('eta_gen_Nu','',source)[0],
    J2_eta_gen_Nu = makeJetUserFloat('eta_gen_Nu','',source)[1],
    J3_eta_gen_Nu = makeJetUserFloat('eta_gen_Nu','',source)[2],
    J4_eta_gen_Nu = makeJetUserFloat('eta_gen_Nu','',source)[3],
    J1_phi_gen_Nu = makeJetUserFloat('phi_gen_Nu','',source)[0],
    J2_phi_gen_Nu = makeJetUserFloat('phi_gen_Nu','',source)[1],
    J3_phi_gen_Nu = makeJetUserFloat('phi_gen_Nu','',source)[2],
    J4_phi_gen_Nu = makeJetUserFloat('phi_gen_Nu','',source)[3],
   )
   return jetGenInfo

def makeLeptons(srcGoodMu='goodMuons',srcVetoMu='vetoMuons',srcQCDMu='qcdMuons',
                srcGoodEle='goodElectrons',srcVetoEle='vetoElectrons',srcQCDEle='qcdElectrons'):

    electronLooseID = \
        '((' \
        'abs(superCluster.eta) < 1.442 &&' \
        'abs(deltaEtaSuperClusterTrackAtVtx) < 0.007 &&' \
        'abs(deltaPhiSuperClusterTrackAtVtx) < 0.15 &&' \
        'sigmaIetaIeta < 0.01 &&' \
        'hadronicOverEm < 0.12' \
        ')||(' \
        'abs(superCluster.eta) > 1.566 && abs(superCluster.eta) < 2.5 &&' \
        'abs(deltaEtaSuperClusterTrackAtVtx) < 0.009 &&' \
        'abs(deltaPhiSuperClusterTrackAtVtx) < 0.10 &&' \
        'sigmaIetaIeta < 0.03 &&' \
        'hadronicOverEm < 0.10' \
        ')) &&' \
        'abs(dB) < 0.02 &&' \
        'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &&' \
        'passConversionVeto &&' \
        'gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 &&' \
        'triggerObjectMatches.size >= 0'

    electronTightID =  \
        '((' \
        'abs(superCluster.eta) < 1.442 &&' \
        'abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &&' \
        'abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 &&' \
        'sigmaIetaIeta < 0.01 &&' \
        'hadronicOverEm < 0.12' \
        ')||(' \
        'abs(superCluster.eta) > 1.566 && abs(superCluster.eta) < 2.5 &&' \
        'abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 &&' \
        'abs(deltaPhiSuperClusterTrackAtVtx) < 0.02 &&' \
        'sigmaIetaIeta < 0.03 &&' \
        'hadronicOverEm < 0.10' \
        ')) &&' \
        'abs(dB) < 0.02 &&' \
        'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &&' \
        'passConversionVeto &&' \
        'gsfTrack.trackerExpectedHitsInner.numberOfHits < 1'

    electronIso = '(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et'
    muonIso = "(pfIsolationR04().sumChargedHadronPt + max((pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt),0.0))/pt"

    leptons = cms.PSet(
      allMu_pt_vec = makeMuon("allMu_pt_vec","pt",sourceMuons="allMuons",lo=False),
      allMu_eta_vec = makeMuon("allMu_eta_vec","eta",sourceMuons="allMuons",lo=False),
      allMu_phi_vec = makeMuon("allMu_phi_vec","phi",sourceMuons="allMuons",lo=False),
      allMu_charge_vec = makeMuon("allMu_charge_vec","charge",sourceMuons="allMuons",lo=False),
      allMu_looseID_vec = makeMuon("allMu_looseID_vec","userFloat('looseID')",sourceMuons="allMuons",lo=False),
      allMu_Iso_vec = makeMuon("allMu_Iso_vec",muonIso,sourceMuons="allMuons",lo=False),
      allMu_chIso_vec = makeMuon("allMu_chIso_vec","pfIsolationR04().sumChargedHadronPt",sourceMuons="allMuons",lo=False),
      allMu_nhIso_vec = makeMuon("allMu_nhIso_vec","pfIsolationR04().sumNeutralHadronEt",sourceMuons="allMuons",lo=False),
      allMu_phIso_vec = makeMuon("allMu_phIso_vec","pfIsolationR04().sumPhotonEt",sourceMuons="allMuons",lo=False),
      allMu_puIso_vec = makeMuon("allMu_puIso_vec","pfIsolationR04().sumPUPt",sourceMuons="allMuons",lo=False),

      goodMu_pt_vec = makeMuon("goodMu_pt_vec","pt",sourceMuons=srcGoodMu,lo=False),
      goodMu_eta_vec = makeMuon("goodMu_eta_vec","eta",sourceMuons=srcGoodMu,lo=False),
      goodMu_phi_vec = makeMuon("goodMu_phi_vec","phi",sourceMuons=srcGoodMu,lo=False),
      goodMu_charge_vec = makeMuon("goodMu_charge_vec","charge",sourceMuons=srcGoodMu,lo=False),
 
      vetoMu_pt_vec = makeMuon("vetoMu_pt_vec","pt",sourceMuons=srcVetoMu,lo=False),
      vetoMu_eta_vec = makeMuon("vetoMu_eta_vec","eta",sourceMuons=srcVetoMu,lo=False),
      vetoMu_phi_vec = makeMuon("vetoMu_phi_vec","phi",sourceMuons=srcVetoMu,lo=False),
      vetoMu_charge_vec = makeMuon("vetoMu_charge_vec","charge",sourceMuons=srcVetoMu,lo=False),
 
      qcdMu_pt_vec = makeMuon("qcdMu_pt_vec","pt",sourceMuons=srcQCDMu,lo=False),
      qcdMu_eta_vec = makeMuon("qcdMu_eta_vec","eta",sourceMuons=srcQCDMu,lo=False),
      qcdMu_phi_vec = makeMuon("qcdMu_phi_vec","phi",sourceMuons=srcQCDMu,lo=False),
      qcdMu_charge_vec = makeMuon("qcdMu_charge_vec","charge",sourceMuons=srcQCDMu,lo=False),
 
      allEle_pt_vec = makeElectron("allEle_pt_vec","pt",sourceElectrons="allElectrons",lo=False),
      allEle_eta_vec = makeElectron("allEle_eta_vec","eta",sourceElectrons="allElectrons",lo=False),
      allEle_phi_vec = makeElectron("allEle_phi_vec","phi",sourceElectrons="allElectrons",lo=False),
      allEle_charge_vec = makeElectron("allEle_charge_vec","charge",sourceElectrons="allElectrons",lo=False),
      #allEle_looseID_vec = makeElectron("allEle_looseID_vec",electronLooseID,sourceElectrons="allElectrons",lo=False),
      #allEle_tightID_vec = makeElectron("allEle_tightID_vec",electronTightID,sourceElectrons="allElectrons",lo=False),
      allEle_iso_vec = makeElectron("allEle_iso_vec",electronIso,sourceElectrons="allElectrons",lo=False),
      allEle_scEta_vec = makeElectron("allEle_scEta_vec","superCluster.eta",sourceElectrons="allElectrons",lo=False),
      allEle_dEscTratVtx_vec = makeElectron("allEle_dEscTratVtx_vec","deltaEtaSuperClusterTrackAtVtx",sourceElectrons="allElectrons",lo=False),
      allEle_dPscTratVtx_vec = makeElectron("allEle_dPscTratVtx_vec","deltaPhiSuperClusterTrackAtVtx",sourceElectrons="allElectrons",lo=False),
      allEle_sIeIe_vec = makeElectron("allEle_sIeIe_vec","sigmaIetaIeta",sourceElectrons="allElectrons",lo=False),
      allEle_hoEM_vec = makeElectron("allEle_hoEM_vec","hadronicOverEm",sourceElectrons="allElectrons",lo=False),
      allEle_dB_vec = makeElectron("allEle_dB_vec","dB",sourceElectrons="allElectrons",lo=False),
      allEle_ecalE_vec = makeElectron("allEle_ecalE_vec","ecalEnergy",sourceElectrons="allElectrons",lo=False),
      allEle_escoP_vec = makeElectron("allEle_escoP_vec","eSuperClusterOverP",sourceElectrons="allElectrons",lo=False),
      allEle_pConvVto_vec = makeElectron("allEle_pConvVto_vec","passConversionVeto",sourceElectrons="allElectrons",lo=False),
      allEle_trkNrH_vec = makeElectron("allEle_trkNrH_vec","gsfTrack.trackerExpectedHitsInner.numberOfHits",sourceElectrons="allElectrons",lo=False),
      allEle_IoEmIoP_vec = makeElectron("allEle_IoEmIoP_vec","1./ecalEnergy - eSuperClusterOverP/ecalEnergy",sourceElectrons="allElectrons",lo=False),
      allEle_chIso_vec = makeElectron("allEle_chIso_vec","chargedHadronIso",sourceElectrons="allElectrons",lo=False),
      allEle_nhIso_vec = makeElectron("allEle_nhIso_vec","neutralHadronIso",sourceElectrons="allElectrons",lo=False),
      allEle_phIso_vec = makeElectron("allEle_phIso_vec","photonIso",sourceElectrons="allElectrons",lo=False),
      allEle_puIso_vec = makeElectron("allEle_puIso_vec","puChargedHadronIso",sourceElectrons="allElectrons",lo=False),

      goodEle_pt_vec = makeElectron("goodEle_pt_vec","pt",sourceElectrons=srcGoodEle,lo=False),
      goodEle_eta_vec = makeElectron("goodEle_eta_vec","eta",sourceElectrons=srcGoodEle,lo=False),
      goodEle_phi_vec = makeElectron("goodEle_phi_vec","phi",sourceElectrons=srcGoodEle,lo=False),
      goodEle_charge_vec = makeElectron("goodEle_charge_vec","charge",sourceElectrons=srcGoodEle,lo=False),
 
      vetoEle_pt_vec = makeElectron("vetoEle_pt_vec","pt",sourceElectrons=srcVetoEle,lo=False),
      vetoEle_eta_vec = makeElectron("vetoEle_eta_vec","eta",sourceElectrons=srcVetoEle,lo=False),
      vetoEle_phi_vec = makeElectron("vetoEle_phi_vec","phi",sourceElectrons=srcVetoEle,lo=False),
      vetoEle_charge_vec = makeElectron("vetoEle_charge_vec","charge",sourceElectrons=srcVetoEle,lo=False),
 
      qcdEle_pt_vec = makeElectron("qcdEle_pt_vec","pt",sourceElectrons=srcQCDEle,lo=False),
      qcdEle_eta_vec = makeElectron("qcdEle_eta_vec","eta",sourceElectrons=srcQCDEle,lo=False),
      qcdEle_phi_vec = makeElectron("qcdEle_phi_vec","phi",sourceElectrons=srcQCDEle,lo=False),
      qcdEle_charge_vec = makeElectron("qcdEle_charge_vec","charge",sourceElectrons=srcQCDEle,lo=False),
    )
    return leptons


def makeCollections(source = 'wCandsJets', sourceZ = 'diMuonsSorted',sourceE = 'weCandsJets', srcElectrons="smearedElectrons",srcMuons="smearedMuons"):
 commonCollections = cms.PSet(         
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
    J3_partonFlavour = makeJetStringPar('partonFlavour','',source)[2],
    J4_partonFlavour = makeJetStringPar('partonFlavour','',source)[3],
    J1_muonMultiplicity = makeJetStringPar('muonMultiplicity','',source)[0],
    J2_muonMultiplicity = makeJetStringPar('muonMultiplicity','',source)[1],
    J3_muonMultiplicity = makeJetStringPar('muonMultiplicity','',source)[2],
    J4_muonMultiplicity = makeJetStringPar('muonMultiplicity','',source)[3],
    J1_chargedMultiplicity = makeJetStringPar('chargedMultiplicity','',source)[0],
    J2_chargedMultiplicity = makeJetStringPar('chargedMultiplicity','',source)[1],
    J3_chargedMultiplicity = makeJetStringPar('chargedMultiplicity','',source)[2],
    J4_chargedMultiplicity = makeJetStringPar('chargedMultiplicity','',source)[3],
    J1_electronMultiplicity = makeJetStringPar('electronMultiplicity','',source)[0],
    J2_electronMultiplicity = makeJetStringPar('electronMultiplicity','',source)[1],
    J1_photonMultiplicity = makeJetStringPar('photonMultiplicity','',source)[0],
    J2_photonMultiplicity = makeJetStringPar('photonMultiplicity','',source)[1],
    J1_jetCharge = makeJetString('jetCharge','',source)[0],
    J2_jetCharge = makeJetString('jetCharge','',source)[1],
    J3_jetCharge = makeJetString('jetCharge','',source)[2],
    J4_jetCharge = makeJetString('jetCharge','',source)[3],
    J1_mass = makeJetString('mass','',source)[0],
    J2_mass = makeJetString('mass','',source)[1],
    #J1_photonEnergy = makeJetString('photonEnergy','',source)[0],
    #J2_photonEnergy = makeJetString('photonEnergy','',source)[1],
    #J1_photonEnergyFraction = makeJetString('photonEnergyFraction','',source)[0],
    #J2_photonEnergyFraction = makeJetString('photonEnergyFraction','',source)[1],
    #J1_electronEnergy = makeJetString('electronEnergy','',source)[0],
    #J2_electronEnergy = makeJetString('electronEnergy','',source)[1],
    #J1_chargedEmEnergyFraction = makeJetString('chargedEmEnergyFraction','',source)[0],
    #J2_chargedEmEnergyFraction = makeJetString('chargedEmEnergyFraction','',source)[1],
    #J1_chargedMuEnergy = makeJetString('chargedMuEnergy','',source)[0],
    #J2_chargedMuEnergy = makeJetString('chargedMuEnergy','',source)[1],
    #J1_muonEnergyFraction = makeJetString('muonEnergyFraction','',source)[0],
    #J2_muonEnergyFraction = makeJetString('muonEnergyFraction','',source)[1],
    #J1_chargedHadronMultiplicity = makeJetStringPar('chargedHadronMultiplicity','',source)[0],
    #J2_chargedHadronMultiplicity = makeJetStringPar('chargedHadronMultiplicity','',source)[1],
# Jet Kinematic Variables
    J1_pt = makeJetStringPar('pt','',source)[0], # formerly highestJetPt
    J2_pt = makeJetStringPar('pt','',source)[1],
    J3_pt = makeJetStringPar('pt','',source)[2],
    J4_pt = makeJetStringPar('pt','',source)[3],
    J5_pt = makeJetStringPar('pt','',source)[4], # formerly highestJetPt
    J6_pt = makeJetStringPar('pt','',source)[5],
    J7_pt = makeJetStringPar('pt','',source)[6],
    J8_pt = makeJetStringPar('pt','',source)[7],
    J1_phi = makeJetStringPar('phi','',source)[0],
    J2_phi = makeJetStringPar('phi','',source)[1],
    J3_phi = makeJetStringPar('phi','',source)[2],
    J4_phi = makeJetStringPar('phi','',source)[3],
    J5_phi = makeJetStringPar('phi','',source)[4],
    J6_phi = makeJetStringPar('phi','',source)[5],
    J7_phi = makeJetStringPar('phi','',source)[6],
    J8_phi = makeJetStringPar('phi','',source)[7],
    J1_eta = makeJetStringPar('eta','',source)[0],
    J2_eta = makeJetStringPar('eta','',source)[1],
    J3_eta = makeJetStringPar('eta','',source)[2],
    J4_eta = makeJetStringPar('eta','',source)[3],
    J5_eta = makeJetStringPar('eta','',source)[4],
    J6_eta = makeJetStringPar('eta','',source)[5],
    J7_eta = makeJetStringPar('eta','',source)[6],
    J8_eta = makeJetStringPar('eta','',source)[7],

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

# Few Electron Variables for MuEle control region
# BIG WARNING: This ID is super outdated. We need to go over the new EGamma ID and actualize it. WP80 is similar to the new "medium" one
    electron_pt_vec = makeElectron("electron_pt_vec","pt",sourceElectrons=srcElectrons,lo=False),
    electron_eta_vec = makeElectron("electron_eta_vec","eta",sourceElectrons=srcElectrons,lo=False),
    electron_phi_vec = makeElectron("electron_phi_vec","phi",sourceElectrons=srcElectrons,lo=False),

    ptEle = makeElectron("ele_pt","pt",sourceElectrons=srcElectrons),
    phiEle = makeElectron("ele_phi","phi",sourceElectrons=srcElectrons),
    etaEle = makeElectron("ele_eta","eta",sourceElectrons=srcElectrons),
    wp80Ele = makeElectron("ele_wp80","userFloat('wp80')",sourceElectrons=srcElectrons),
    chargeEle = makeElectron("ele_charge","charge",sourceElectrons=srcElectrons),
    isoEleDB = makeElectron("isoEleDB","(userIso(0)+max(userIso(1)+neutralHadronIso()-0.5*userIso(2),0.0))/pt",sourceElectrons=srcElectrons),

# Z Variables
    DiMuonMass = makeZColl("DiMuonMass","mass",sourceZ),
    DiMuonPt = makeZColl("DiMuonPt","pt",sourceZ),
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

    nrvtx = makeMuon("nrvtx",'userFloat("nvtx")',sourceMuons=srcMuons),
    nrPU  = makeMuon("nrPU",'userFloat("nrPU")',sourceMuons=srcMuons),
    preWeight = makeMuon("preWeight",'userFloat("preWeight")',sourceMuons=srcMuons),

    nrMu = makeCollSize('selectedPatMuons','nrMu'),
    nrEle = makeCollSize('selectedPatElectrons','nrEle'),
    nrMuLoose = makeCollSize('preselectedPatMuons','nrMuLoose'),
    nrEleLoose = makeCollSize('preselectedPatElectrons','nrEleLoose'),
    nrGoodMu = makeCollSize('goodMuons','nrGoodMu'),
    nrVetoMu = makeCollSize('vetoMuons','nrVetoMu'),
    nrQCDMu = makeCollSize('qcdMuons','nrQCDMu'),
    nrGoodEle = makeCollSize('goodElectrons','nrGoodEle'),
    nrVetoEle = makeCollSize('vetoElectrons','nrVetoEle'),
    nrQCDEle = makeCollSize('qcdElectrons','nrQCDEle'),
 

    #muon_pt = makeMuNu("muon_pt","lepton().pt",source),
    #muon_eta = makeMuNu("muon_eta","lepton.eta",source),
    #muon_phi = makeMuNu("muon_phi","lepton.phi",source),
    #muon_charge = makeMuNu("muon_charge","lepton.charge()",source),
    muon_pt_vec = makeMuNu("muon_pt_vec","lepton().pt",source,lo=False),
    muon_eta_vec  = makeMuNu("muon_eta_vec","lepton.eta",source,lo=False),
    muon_phi_vec  = makeMuNu("muon_phi_vec","lepton.phi",source,lo=False),
    muon_charge_vec  = makeMuNu("muon_charge_vec","lepton.charge()",source,lo=False),
    Wpt = makeMuNu("Wpt","corPt()",source),
    muon_pt_vec_2 = makeMuon("muon_pt_vec_2","pt",sourceMuons=srcMuons,lo=False),
    muon_eta_vec_2 = makeMuon("muon_eta_vec_2","eta",sourceMuons=srcMuons,lo=False),
    muon_phi_vec_2 = makeMuon("muon_phi_vec_2","phi",sourceMuons=srcMuons,lo=False),
    muon_charge_vec_2 = makeMuon("muon_charge_vec_2","charge",sourceMuons=srcMuons,lo=False),

    eleNuRelPFIsoDB_A = makeElectron("eleNuRelPFIsoDB_A",
      "(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et",
      srcElectrons,lo=True),
    eleNuRelPFIsoDB_A_vec = makeElectron("eleNuRelPFIsoDB_A_vec",
      "(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et",
      srcElectrons,lo=False),

    muNuRelPFIsoDB_A = makeMuNu("muNuRelPFIsoDB_A",
      "(lepton.pfIsolationR04().sumChargedHadronPt + max((lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5*lepton.pfIsolationR04().sumPUPt),0.0))/lepton.pt",
      source,lo=True),
    muNuRelPFIsoDB_A_vec = makeMuNu("muNuRelPFIsoDB_A_vec",
      "(lepton.pfIsolationR04().sumChargedHadronPt + max((lepton.pfIsolationR04().sumNeutralHadronEt + lepton.pfIsolationR04().sumPhotonEt - 0.5*lepton.pfIsolationR04().sumPUPt),0.0))/lepton.pt",
      source,lo=False),

    ipDXY = makeMuNu("ipDXY","lepton.userFloat('ipDXY')",source,True),
    met_pt = makeMuNu("met_pt","met().pt",source,True),
    met_phi = makeMuNu("met_phi","met().phi",source,True),
    mt_muon_vec = makeMuon("mt_muon_vec","userFloat('mt')",sourceMuons=srcMuons,lo=False),
    mt_electron_vec = makeElectron("mt_electron_vec","userFloat('mt')",sourceElectrons=srcElectrons,lo=False),
    
    metJJ = makeMuNu("metjj","metjj",source),
    leptonjj = makeMuNu("leptonjj","leptonjj",source),
    muNuDPhi = makeMuNu("muNuDPhi","dPhi",source),
    muNuRecoil = makeMuNu("muNuRecoil","recoil().pt()",source),

    dz = makeMuNu("dz",'abs(lepton.userFloat("dz"))',source,True),
    ht = makeMuNu("ht","ht",source,True),
 )
 return commonCollections

def addMuNuEventTreePtData(process,name,source = 'wCandsJets',sourceZ = 'diMuonsSorted'):
   process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
   eventTree = cms.EDAnalyzer('EventTreeMaker',
      makeCollections(source,sourceZ),
      makeJetIDInfo(source),
      #makeJetCorrectionInfo(),
      makeJetSVInfo(source),
      coreCollections = cms.VInputTag( cms.InputTag(source) ),
      bCandidatenBC = makeIVFBs("bCandidatenBC","nBC"),
      bCandidate1Pt = makeIVFBs("bCandidate1Pt","BC1PT"),
      bCandidate2Pt = makeIVFBs("bCandidate2Pt","BC2PT"),
      bCandidate1Eta = makeIVFBs("bCandidate1Eta","BC1ETA"),
      bCandidate2Eta = makeIVFBs("bCandidate2Eta","BC2ETA"),
      bCandidate1Phi = makeIVFBs("bCandidate1Phi","BC1PHI"),
      bCandidate2Phi = makeIVFBs("bCandidate2Phi","BC2PHI"),
      bCandidateDeltaR = makeIVFBs("bCandidateDeltaR","BCDeltaR"),
      bCandidateDeltaPhi = makeIVFBs("bCandidateDeltaPhi","BDeltaPHI"),
      bCandidateBC1nTracks = makeIVFBs("bCandidateBC1nTracks","BC1nTracks"),
      bCandidateBC2nTracks = makeIVFBs("bCandidateBC2nTracks","BC2nTracks"),
      bCandidateBCEnergySum = makeIVFBs("bCandidateBCEnergySum","BCEnergySum"),
      bCandidateBCGenEnergySum = makeIVFBs("bCandidateBCGenEnergySum","BCGenEnergySum"),
      bCandidatenJetAssociatedToBC1 = makeIVFBs("bCandidatenJetAssociatedToBC1","nJetAssociatedToBC1"),
      bCandidatenJetAssociatedToBC2 = makeIVFBs("bCandidatenJetAssociatedToBC2","nJetAssociatedToBC2"),
      bCandidatenJetAssociatedToBC = makeIVFBs("bCandidatenJetAssociatedToBC","nJetAssociatedToBC"),
      bCandidateJet1BC1PT = makeIVFBs("bCandidateJet1BC1PT","Jet1BC1PT"),
      bCandidateJet2BC1PT = makeIVFBs("bCandidateJet2BC1PT","Jet2BC1PT"),
      bCandidateJet3BC1PT = makeIVFBs("bCandidateJet3BC1PT","Jet3BC1PT"),
      bCandidateJet4BC1PT = makeIVFBs("bCandidateJet4BC1PT","Jet4BC1PT"),
      bCandidateJet1BC1ETA = makeIVFBs("bCandidateJet1BC1ETA","Jet1BC1ETA"),
      bCandidateJet2BC1ETA = makeIVFBs("bCandidateJet2BC1ETA","Jet2BC1ETA"),
      bCandidateJet3BC1ETA = makeIVFBs("bCandidateJet3BC1ETA","Jet3BC1ETA"),
      bCandidateJet4BC1ETA = makeIVFBs("bCandidateJet4BC1ETA","Jet4BC1ETA"),
      bCandidateJet1BC1PHI = makeIVFBs("bCandidateJet1BC1PHI","Jet1BC1PHI"),
      bCandidateJet2BC1PHI = makeIVFBs("bCandidateJet2BC1PHI","Jet2BC1PHI"),
      bCandidateJet3BC1PHI = makeIVFBs("bCandidateJet3BC1PHI","Jet3BC1PHI"),
      bCandidateJet4BC1PHI = makeIVFBs("bCandidateJet4BC1PHI","Jet4BC1PHI"),
      bCandidateJet1BC2PT = makeIVFBs("bCandidateJet1BC2PT","Jet1BC2PT"),
      bCandidateJet2BC2PT = makeIVFBs("bCandidateJet2BC2PT","Jet2BC2PT"),
      bCandidateJet3BC2PT = makeIVFBs("bCandidateJet3BC2PT","Jet3BC2PT"),
      bCandidateJet4BC2PT = makeIVFBs("bCandidateJet4BC2PT","Jet4BC2PT"),
      bCandidateJet1BC2ETA = makeIVFBs("bCandidateJet1BC2ETA","Jet1BC2ETA"),
      bCandidateJet2BC2ETA = makeIVFBs("bCandidateJet2BC2ETA","Jet2BC2ETA"),
      bCandidateJet3BC2ETA = makeIVFBs("bCandidateJet3BC2ETA","Jet3BC2ETA"),
      bCandidateJet4BC2ETA = makeIVFBs("bCandidateJet4BC2ETA","Jet4BC2ETA"),
      bCandidateJet1BC2PHI = makeIVFBs("bCandidateJet1BC2PHI","Jet1BC2PHI"),
      bCandidateJet2BC2PHI = makeIVFBs("bCandidateJet2BC2PHI","Jet2BC2PHI"),
      bCandidateJet3BC2PHI = makeIVFBs("bCandidateJet3BC2PHI","Jet3BC2PHI"),
      bCandidateJet4BC2PHI = makeIVFBs("bCandidateJet4BC2PHI","Jet4BC2PHI"),
      bCandidateBCDeltaR = makeIVFBs("bCandidateBCDeltaR","BCDeltaR"),
      bCandidateBDeltaPHI = makeIVFBs("bCandidateBDeltaPHI","BDeltaPHI"),
      bCandidateBC1MASS = makeIVFBs("bCandidateBC1MASS","BC1MASS"),
      bCandidateBC2MASS = makeIVFBs("bCandidateBC2MASS","BC2MASS"),
   )
   setattr(process, name, eventTree)
   p = cms.Path(getattr(process,name))
   setattr(process, name+'Path', p)

def addMuNuEventTreePtMC(process,name,source = 'MuNuSel',sourceZ = 'diMuonsSorted',lhep="externalLHEProducer",srcMuons="smearedMuonsall", srcElectrons="smearedElectronsall"):
   process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
   eventTree = cms.EDAnalyzer('EventTreeMaker',
      makeCollections(source = 'MuNuAll', sourceZ = 'diMuonsSorted',sourceE = 'EleNuAll', srcElectrons="smearedElectronsiall",srcMuons="smearedMuonsall"),
      makeJetIDInfo(source="MuNuAll"),
      #makeJetCorrectionInfo(),
      makeGenJetInfo(source="MuNuAll"),
      #makeLeptons(srcGoodMu='smearedMuonsall',srcVetoMu='smearedMuonsall',srcQCDMu='smearedMuonsall',srcGoodEle='smearedElectronsall',srcVetoEle='smearedElectronsall',srcQCDEle='smearedElectronsall'),
      #makeLeptons(srcGoodMu='smearedMuonsgood',srcVetoMu='smearedMuonsveto',srcQCDMu='smearedMuonsqcd',srcGoodEle='smearedElectronsgood',srcVetoEle='smearedElectronsveto',srcQCDEle='smearedElectronsqcd'),
      makeLeptons(srcGoodMu='goodMuons',srcVetoMu='vetoMuons',srcQCDMu='qcdMuons',srcGoodEle='goodElectrons',srcVetoEle='vetoElectrons',srcQCDEle='qcdElectrons'),
      makeJetsLegit(sourceJets='smearedJetsall'),
      makeJetsInd(sourceJets='smearedJetsall'),
      makeMETsColl(sourceMETs = 'EleNuAll'),
      makeMETsInd(sourceMETs = 'smearedMETall'),
      makeJetSVInfo(source='MuNuSel'),
      coreCollections = cms.VInputTag(
           cms.InputTag(source),
           cms.InputTag(srcElectrons),
           cms.InputTag("allMuons"),
           cms.InputTag("allElectrons"),
           cms.InputTag('metCorrected')
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
  
      ## Jet pT Reco + Gen from RecoTools/plugins/PATJetSmearer.h
      #J1_pt_gen = makeJetUserFloat('pt_gen')[0],
      #J1_pt_gen_two = makeJetUserFloat('pt_gen_two')[0],
      #J1_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[0], 
      #J1_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[0], 
      #J1_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[0], 
      #J2_pt_gen = makeJetUserFloat('pt_gen')[1],
      #J2_pt_gen_two = makeJetUserFloat('pt_gen_two')[1],
      #J2_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[1], 
      #J2_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[1], 
      #J2_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[1], 
      #J3_pt_gen = makeJetUserFloat('pt_gen')[2],
      #J3_pt_gen_two = makeJetUserFloat('pt_gen_two')[2],
      #J3_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[2], 
      #J3_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[2], 
      #J3_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[2], 
      #J4_pt_gen = makeJetUserFloat('pt_gen')[3],
      #J4_pt_gen_two = makeJetUserFloat('pt_gen_two')[3],
      #J4_pt_uncorr = makeJetUserFloat('pt_uncorr','',source)[3], 
      #J4_pt_smearedUp = makeJetUserFloat('pt_smearedUp','',source)[3], 
      #J4_pt_smearedDown = makeJetUserFloat('pt_smearedDown','',source)[3], 

      bCandidatenBC = makeIVFBs("bCandidatenBC","nBC"),
      bCandidate1Pt = makeIVFBs("bCandidate1Pt","BC1PT"),
      bCandidate2Pt = makeIVFBs("bCandidate2Pt","BC2PT"),
      bCandidate1Eta = makeIVFBs("bCandidate1Eta","BC1ETA"),
      bCandidate2Eta = makeIVFBs("bCandidate2Eta","BC2ETA"),
      bCandidate1Phi = makeIVFBs("bCandidate1Phi","BC1PHI"),
      bCandidate2Phi = makeIVFBs("bCandidate2Phi","BC2PHI"),
      bCandidateDeltaR = makeIVFBs("bCandidateDeltaR","BCDeltaR"),
      bCandidateDeltaPhi = makeIVFBs("bCandidateDeltaPhi","BDeltaPHI"),
      bCandidateBC1nTracks = makeIVFBs("bCandidateBC1nTracks","BC1nTracks"),
      bCandidateBC2nTracks = makeIVFBs("bCandidateBC2nTracks","BC2nTracks"),
      bCandidateBCEnergySum = makeIVFBs("bCandidateBCEnergySum","BCEnergySum"),
      bCandidateBCGenEnergySum = makeIVFBs("bCandidateBCGenEnergySum","BCGenEnergySum"),
      bCandidatenJetAssociatedToBC1 = makeIVFBs("bCandidatenJetAssociatedToBC1","nJetAssociatedToBC1"),
      bCandidatenJetAssociatedToBC2 = makeIVFBs("bCandidatenJetAssociatedToBC2","nJetAssociatedToBC2"),
      bCandidatenJetAssociatedToBC = makeIVFBs("bCandidatenJetAssociatedToBC","nJetAssociatedToBC"),
      bCandidateJet1BC1PT = makeIVFBs("bCandidateJet1BC1PT","Jet1BC1PT"),
      bCandidateJet2BC1PT = makeIVFBs("bCandidateJet2BC1PT","Jet2BC1PT"),
      bCandidateJet3BC1PT = makeIVFBs("bCandidateJet3BC1PT","Jet3BC1PT"),
      bCandidateJet4BC1PT = makeIVFBs("bCandidateJet4BC1PT","Jet4BC1PT"),
      bCandidateJet1BC1ETA = makeIVFBs("bCandidateJet1BC1ETA","Jet1BC1ETA"),
      bCandidateJet2BC1ETA = makeIVFBs("bCandidateJet2BC1ETA","Jet2BC1ETA"),
      bCandidateJet3BC1ETA = makeIVFBs("bCandidateJet3BC1ETA","Jet3BC1ETA"),
      bCandidateJet4BC1ETA = makeIVFBs("bCandidateJet4BC1ETA","Jet4BC1ETA"),
      bCandidateJet1BC1PHI = makeIVFBs("bCandidateJet1BC1PHI","Jet1BC1PHI"),
      bCandidateJet2BC1PHI = makeIVFBs("bCandidateJet2BC1PHI","Jet2BC1PHI"),
      bCandidateJet3BC1PHI = makeIVFBs("bCandidateJet3BC1PHI","Jet3BC1PHI"),
      bCandidateJet4BC1PHI = makeIVFBs("bCandidateJet4BC1PHI","Jet4BC1PHI"),
      bCandidateJet1BC2PT = makeIVFBs("bCandidateJet1BC2PT","Jet1BC2PT"),
      bCandidateJet2BC2PT = makeIVFBs("bCandidateJet2BC2PT","Jet2BC2PT"),
      bCandidateJet3BC2PT = makeIVFBs("bCandidateJet3BC2PT","Jet3BC2PT"),
      bCandidateJet4BC2PT = makeIVFBs("bCandidateJet4BC2PT","Jet4BC2PT"),
      bCandidateJet1BC2ETA = makeIVFBs("bCandidateJet1BC2ETA","Jet1BC2ETA"),
      bCandidateJet2BC2ETA = makeIVFBs("bCandidateJet2BC2ETA","Jet2BC2ETA"),
      bCandidateJet3BC2ETA = makeIVFBs("bCandidateJet3BC2ETA","Jet3BC2ETA"),
      bCandidateJet4BC2ETA = makeIVFBs("bCandidateJet4BC2ETA","Jet4BC2ETA"),
      bCandidateJet1BC2PHI = makeIVFBs("bCandidateJet1BC2PHI","Jet1BC2PHI"),
      bCandidateJet2BC2PHI = makeIVFBs("bCandidateJet2BC2PHI","Jet2BC2PHI"),
      bCandidateJet3BC2PHI = makeIVFBs("bCandidateJet3BC2PHI","Jet3BC2PHI"),
      bCandidateJet4BC2PHI = makeIVFBs("bCandidateJet4BC2PHI","Jet4BC2PHI"),
      bCandidateBCDeltaR = makeIVFBs("bCandidateBCDeltaR","BCDeltaR"),
      bCandidateBDeltaPHI = makeIVFBs("bCandidateBDeltaPHI","BDeltaPHI"),
      bCandidateBC1MASS = makeIVFBs("bCandidateBC1MASS","BC1MASS"),
      bCandidateBC2MASS = makeIVFBs("bCandidateBC2MASS","BC2MASS"),
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
      GBWeight = makeGenBWeight("GBWeight","theWeight"),
   )
   setattr(process, name, eventTree)
   p = cms.Path(getattr(process,name))
   setattr(process, name+'Path', p)
