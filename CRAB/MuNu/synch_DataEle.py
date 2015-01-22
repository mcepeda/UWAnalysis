import FWCore.ParameterSet.Config as cms
import sys

theFile = sys.argv[2]
nr = sys.argv[3]

sys.setrecursionlimit(1000)
process = cms.Process("ANALYSIS")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_53_V21_AN6::All' #new GT

process.maxEvents = cms.untracked.PSet(
 input = cms.untracked.int32(2500)
)

process.load("UWAnalysis/Configuration/2012FullJSON")
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/hep.wisc.edu/cms/tperry/Wbb_CMSSW_5_3_14_patch1/src/UWAnalysis/CRAB/MuNu/weights/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.json').getVLuminosityBlockRange()

process.source.fileNames= cms.untracked.vstring(
theFile
)
process.source.inputCommands=cms.untracked.vstring(
 'keep *', 'drop *_finalState*_*_*',
 'drop *_patFinalStateEvent*_*_*'
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("PhysicsTools.PatAlgos.patSequences_cff")
from UWAnalysis.Configuration.tools.analysisToolsPT import *

defaultReconstructionPT(process,'HLT',['HLT_IsoMu24_eta2p1_v','HLT_IsoMu24','HLT_Ele27_WP80'],
                        itsMC=False,itsData=True)

process.load("UWAnalysis.Configuration.wMuNuAnalysisPT_cff")
process.eventSelection = cms.Path(process.selectionSequence) 

from UWAnalysis.Configuration.tools.ntupleToolsPTwbbClean import *

addEventTreeData(process,'electronEventTree',
      'smearedGoodElectrons',
      'smearedQCDElectrons',
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
      srcCJet='smearedCleanJets'
   )

addEventSummary(process,True)
process.TFileService.fileName = cms.string('data_synch%s.root'%(nr))

## makes EDM output of all collections
#process.out = cms.OutputModule("PoolOutputModule",
# fileName = cms.untracked.string('/scratch/tperry/data.root'),
# outputCommands = cms.untracked.vstring(
#                       'keep *'),
# )
#process.e = cms.EndPath(process.out)
