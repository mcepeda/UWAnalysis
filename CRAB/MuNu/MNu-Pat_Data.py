import FWCore.ParameterSet.Config as cms
import sys
sys.setrecursionlimit(1000)
process = cms.Process("ANALYSIS")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_53_V21_AN6::All' #new GT

process.maxEvents = cms.untracked.PSet(
 input = cms.untracked.int32(10)
)

process.load("UWAnalysis/Configuration/2012FullJSON")
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/hep.wisc.edu/cms/tperry/Wbb_CMSSW_5_3_14_patch1/src/UWAnalysis/CRAB/MuNu/weights/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.json').getVLuminosityBlockRange()

process.source.fileNames= cms.untracked.vstring(
 #'root://cmsxrootd.hep.wisc.edu//store/user/mverzett/SingleMu/Run2012A-22Jan2013-v1/AOD/Spring2014PATTuples_v2/patTuple_cfg-589969E7-646F-E211-998B-1CC1DE1CEDB2.root'

# for synching
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-000A8B8E-2875-E211-BC1E-00259073E3D6.root', 
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-0012575F-4C75-E211-A333-00259073E4D6.root', 
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-001833BF-0575-E211-B2C8-E0CB4E19F9AF.root', 
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-001DD9EE-3975-E211-8FE1-90E6BA19A1FE.root', 
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-002723BE-2575-E211-A28F-90E6BA19A231.root', 
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-0033FBE6-0E75-E211-91C6-00259073E4BE.root', 

#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A4B3C43-C575-E211-8960-90E6BA0D09AD.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A4BE637-F777-E211-BA8E-485B39800BB5.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A542E64-4377-E211-A930-00259073E31C.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A5B339D-C974-E211-9242-002590747DE2.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A6226B7-CA78-E211-82F3-E0CB4E553643.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A68F191-7478-E211-A32A-00261834B51D.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A69AD6D-6477-E211-99CD-00259074AE44.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A6B0C87-2375-E211-B043-E0CB4EA0A926.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A6CF6B1-2675-E211-8CBD-90E6BA19A245.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A7A05AC-AB74-E211-AAEE-20CF3027A5C4.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A7ABD3B-2577-E211-9E3C-00259073E3D2.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A80B8A0-5775-E211-A3E6-E0CB4E29C505.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A82620F-4478-E211-9FF8-E0CB4E19F96E.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A870C87-2377-E211-8950-00259073E4C4.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A8FFDA5-B275-E211-B536-00261834B5D2.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A9005DE-4C77-E211-9C4F-00259073E340.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A9670CA-C875-E211-ACFC-20CF305616E1.root',
#'root://cmsxrootd.hep.wisc.edu//store/user/mcepeda/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples/patTuple_cfg-2A9B8F4B-EF77-E211-A0AA-485B39800B75.root',

#'root://cmsxrootd.hep.wisc.edu//store/user/tperry/SingleMu/Run2012C-22Jan2013-v1/AOD/Spring2014PATTuples_v3/patTuple_cfg-040AE0DE-8875-E211-AF34-001EC9D8D993.root',
'file:///afs/hep.wisc.edu/cms/tperry/FSAv2_CMSSW_5_3_14/src/FinalStateAnalysis/data_14ocbr14_2A4B3C43-C575-E211-8960-90E6BA0D09AD.root'
#'file:///afs/hep.wisc.edu/cms/tperry/FSAv2_CMSSW_5_3_14/src/FinalStateAnalysis/data2_2A4B3C43-C575-E211-8960-90E6BA0D09AD.root'
#'file:///afs/hep.wisc.edu/cms/tperry/FSAv2_CMSSW_5_3_14/src/FinalStateAnalysis/PatTools/test/data.root'
)
process.source.inputCommands=cms.untracked.vstring(
 'keep *', 'drop *_finalState*_*_*',
 'drop *_patFinalStateEvent*_*_*'
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load("PhysicsTools.PatAlgos.patSequences_cff")
from UWAnalysis.Configuration.tools.analysisToolsPT import *

defaultReconstructionPT(process,'HLT',
                        ['HLT_IsoMu24_eta2p1_v','HLT_Mu40_eta2p1_v','HLT_Ele27_WP80'],
                        itsMC=False,itsData=True)

process.load("UWAnalysis.Configuration.wMuNuAnalysisPT_cff")
process.eventSelection = cms.Path(process.selectionSequence) 

from UWAnalysis.Configuration.tools.ntupleToolsPTwbb import *

addMuNuEventTreePtData(process,'muNuEventTree')
addEventSummary(process,True)
process.TFileService.fileName = cms.string('data_T1Txy.root')

## makes EDM output of all collections
#process.out = cms.OutputModule("PoolOutputModule",
# fileName = cms.untracked.string('/scratch/tperry/data.root'),
# outputCommands = cms.untracked.vstring(
#                       'keep *'),
# )
#process.e = cms.EndPath(process.out)
