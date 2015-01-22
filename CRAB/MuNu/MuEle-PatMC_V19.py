import FWCore.ParameterSet.Config as cms
import sys

sys.setrecursionlimit(10000)

process = cms.Process("ANALYSIS")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V27::All' # for global tag with re-reco data

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(
"file:///afs/hep.wisc.edu/cms/mcepeda/UpdateRecipe/CMSSW_5_3_22/src/FinalStateAnalysis/PatTools/test/test.root"
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-00277FF2-7B84-E211-9475-782BCB27B958.root", # "source"
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-002D70FC-7F84-E211-A7E6-782BCB6E0938.root" # "source"

#"root://cmsxrootd.hep.wisc.edu//store/user/tperry/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-00065798-2704-E211-B308-0025901D4C44.root" # "source"
#"root://cmsxrootd.hep.wisc.edu//store/user/tperry/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-0004D7A5-74E4-E211-9F79-002618943963.root"  #"externalLHEProducer"

# synch files
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-00277FF2-7B84-E211-9475-782BCB27B958.root", 
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-5ACBC7A4-7D84-E211-A0B6-782BCB27B958.root", 
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-AED5F4C3-6E84-E211-B83B-0002C90EEE6E.root", 
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-002D70FC-7F84-E211-A7E6-782BCB6E0938.root", 
#"root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/Fall2014PATTuples_V1/patTuple_cfg-5ADB775B-A284-E211-B0A2-90B11C18B19E.root", 
 ),
 inputCommands=cms.untracked.vstring(
  'keep *',
  'drop *_finalState*_*_*',
  'drop *_patFinalStateEvent*_*_*'
 )
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 10000000

process.load("PhysicsTools.PatAlgos.patSequences_cff")
from UWAnalysis.Configuration.tools.analysisToolsPT import *

defaultReconstructionPT(process,
 'HLT',
 ['HLT_IsoMu24_eta2p1_v','HLT_Mu40_eta2p1_v','HLT_Ele27_WP80'],
 itsMC=True,itsData=False)

createGeneratedParticles(process,
 'genWs',
 ["keep++ pdgId = {W+}","keep++ pdgId = {W-}"]
)
createGeneratedParticles(process,
 'genbbCands',
 ["keep abs(pdgId) = 5"]
)
createGeneratedParticles(process,
 'genbbCandsMinus',
 ["keep pdgId = -5"]
)
createGeneratedParticles(process,
'genccCands',
 ["keep abs(pdgId) = 4"]
)
createGeneratedParticles(process,
'gentCands',
 ["keep pdgId = 6"]
)
createGeneratedParticles(process,
'gentbarCands',
 ["keep pdgId = -6"]
)

process.load("UWAnalysis.Configuration.wMuNuAnalysisPT_cff")
process.eventSelection = cms.Path(process.selectionSequence) ##changing to multiples see below

# createSystematics in analysisToolsPT
process.selectionSequenceMuonUp    = createSystematics(process,process.selectionSequence,'MuonUp'  ,1.01, 1.0, 1.0, 0, 1.0)
process.selectionSequenceMuonDown  = createSystematics(process,process.selectionSequence,'MuonDown',0.99, 1.0, 1.0, 0, 1.0)
process.selectionSequenceJetUp     = createSystematics(process,process.selectionSequence,'JetUp'   ,1.00, 1.0, 1.0, 1, 1.0)
process.selectionSequenceJetDown   = createSystematics(process,process.selectionSequence,'JetDown' ,1.00, 1.0, 1.0,-1, 1.0)

process.eventSelectionMuonUp    =  cms.Path(process.selectionSequenceMuonUp)
process.eventSelectionMuonDown  =  cms.Path(process.selectionSequenceMuonDown)
process.eventSelectionJetUp     =  cms.Path(process.selectionSequenceJetUp)
process.eventSelectionJetDown   =  cms.Path(process.selectionSequenceJetDown)

from UWAnalysis.Configuration.tools.ntupleToolsPTwbbClean import *
LHESource = "externalLHEProducer"
addEventTreeMC(process,'muEleEventTree',
      srcGMu='smearedGoodMuons',
      srcVMu='smearedVetoMuons',
      srcQMu='smearedQCDMuons',
      srcAMu='smearedAllMuons',
      srcGEle='smearedGoodElectrons',
      srcVEle='smearedVetoElectrons',
      srcQEle='smearedQCDElectrons',
      srcAEle='smearedAllElectrons',
      srcGJet='smearedGoodJets',
      srcFJet='smearedFwdJets',
      srcAJet='smearedAllJets',
      srcCJet='smearedCleanJets',
      lhep=LHESource
      )

addEventSummary(process,True)

addEventTreeMC(process,'muEleEventTreeJetUp',
      srcGMu='smearedGoodMuonsJetUp',
      srcVMu='smearedVetoMuonsJetUp',
      srcQMu='smearedQCDMuonsJetUp',
      srcAMu='smearedAllMuonsJetUp',
      srcGEle='smearedGoodElectronsJetUp',
      srcVEle='smearedVetoElectronsJetUp',
      srcQEle='smearedQCDElectronsJetUp',
      srcAEle='smearedAllElectronsJetUp',
      srcGJet='smearedGoodJetsJetUp',
      srcFJet='smearedFwdJetsJetUp',
      srcAJet='smearedAllJetsJetUp',
      srcCJet='smearedCleanJetsJetUp',
      lhep=LHESource
      )

addEventTreeMC(process,'muEleEventTreeJetDown',
      srcGMu='smearedGoodMuonsJetDown',
      srcVMu='smearedVetoMuonsJetDown',
      srcQMu='smearedQCDMuonsJetDown',
      srcAMu='smearedAllMuonsJetDown',
      srcGEle='smearedGoodElectronsJetDown',
      srcVEle='smearedVetoElectronsJetDown',
      srcQEle='smearedQCDElectronsJetDown',
      srcAEle='smearedAllElectronsJetDown',
      srcGJet='smearedGoodJetsJetDown',
      srcFJet='smearedFwdJetsJetDown',
      srcAJet='smearedAllJetsJetDown',
      srcCJet='smearedCleanJetsJetDown',
      lhep=LHESource
      )

addEventTreeMC(process,'muEleEventTreeMuonUp',
      srcGMu='smearedGoodMuonsMuonUp',
      srcVMu='smearedVetoMuonsMuonUp',
      srcQMu='smearedQCDMuonsMuonUp',
      srcAMu='smearedAllMuonsMuonUp',
      srcGEle='smearedGoodElectronsMuonUp',
      srcVEle='smearedVetoElectronsMuonUp',
      srcQEle='smearedQCDElectronsMuonUp',
      srcAEle='smearedAllElectronsMuonUp',
      srcGJet='smearedGoodJetsMuonUp',
      srcFJet='smearedFwdJetsMuonUp',
      srcAJet='smearedAllJetsMuonUp',
      srcCJet='smearedCleanJetsMuonUp',
      lhep=LHESource
      )

addEventTreeMC(process,'muEleEventTreeMuonDown',
      srcGMu='smearedGoodMuonsMuonDown',
      srcVMu='smearedVetoMuonsMuonDown',
      srcQMu='smearedQCDMuonsMuonDown',
      srcAMu='smearedAllMuonsMuonDown',
      srcGEle='smearedGoodElectronsMuonDown',
      srcVEle='smearedVetoElectronsMuonDown',
      srcQEle='smearedQCDElectronsMuonDown',
      srcAEle='smearedAllElectronsMuonDown',
      srcGJet='smearedGoodJetsMuonDown',
      srcFJet='smearedFwdJetsMuonDown',
      srcAJet='smearedAllJetsMuonDown',
      srcCJet='smearedCleanJetsMuonDown',
      lhep=LHESource
      )

process.TFileService.fileName = cms.string('MariaTest.root') 

## makes EDM output of all collections
#process.out = cms.OutputModule("PoolOutputModule",
# fileName = cms.untracked.string('/scratch/tperry/analysis.root'),
# outputCommands = cms.untracked.vstring(
#                       'keep *'),
# )
#process.e = cms.EndPath(process.out)
#
#Pax = open("myConfig.py","w")
#Pax.write(process.dumpPython())
