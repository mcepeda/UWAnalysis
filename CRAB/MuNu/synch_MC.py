import FWCore.ParameterSet.Config as cms
import sys

theFile = sys.argv[2]
nr = sys.argv[3]

sys.setrecursionlimit(10000)

process = cms.Process("ANALYSIS")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V27::All' # for global tag with re-reco data

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2500)
)

process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(
   theFile
 ),
 inputCommands=cms.untracked.vstring(
  'keep *',
  'drop *_finalState*_*_*',
  'drop *_patFinalStateEvent*_*_*'
 ),
#eventsToProcess = cms.untracked.VEventRange('1:76947168','1:75784819','1:63254787','1:80376621','1:66903586','1:63257967','1:22134922','1:80365288','1:5523470','1:24561900','1:76947901','1:5521023','1:45723622','1:15627512','1:66862177','1:80813693','1:75786530','1:75772927','1:14715666','1:39788037','1:63254235','1:63253537','1:24858743','1:75787408','1:45723620','1:57931470','1:80367443','1:88787628','1:76947521','1:3197986','1:76950837','1:1580935','1:57926203','1:12280574','1:24563655','1:6472380','1:5522241','1:80374198','1:30792789','1:80375891','1:75787820','1:66861167','1:65469963','1:11571789','1:87455940','1:44401745','1:6474255','1:75786451','1:75774206','1:6473729','1:80364915','1:24561998','1:66904402','1:66903614','1:22131089','1:70329889','1:79335097','1:76947159','1:45727079','1:24565066','1:24858091','1:70327196','1:66900560','1:75775179','1:22132393')
#eventsToProcess = cms.untracked.VEventRange('1:76947168','1:75784819','1:63254787','1:76947672','1:80376621','1:57926203','1:66903586','1:12280574','1:80375891','1:22134922','1:80365288','1:6472380','1:24561900','1:41519977','1:80374198','1:30792789','1:5523470','1:76947901','1:5521023','1:45723622','1:65469963','1:15627512','1:87455940','1:44401745','1:6471598','1:41520479','1:24563655','1:6474255','1:66862177','1:80813693','1:75786530','1:80365903','1:80364915','1:75772927','1:1580935','1:75774206','1:66904402','1:14715666','1:66903614','1:63254235','1:11571789','1:22131089','1:70329889','1:5522241','1:66903328','1:63253537','1:76947159','1:24858743','1:24565066','1:24858091','1:70327196','1:22132947','1:45723620','1:57931470','1:88787628','1:76947521','1:3197986','1:66861167','1:24561998','1:76950837','1:75775179','1:22132393')
#eventsToProcess = cms.untracked.VEventRange('1:65464584','1:65029282','1:65013712','1:65023142','1:65032541','1:65468610','1:65023339','1:65011185','1:65019452','1:65029644','1:52246427','1:60690779','1:52247853','1:60690573','1:65030360','1:52248063','1:65032514','1:60690757','1:52249945','1:52248269','1:68491059','1:65464430','1:68490345','1:75786530','1:52244957','1:52245809','1:65464976','1:55345112')
#eventsToProcess = cms.untracked.VEventRange('1:55344984')
)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 1

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

from UWAnalysis.Configuration.tools.ntupleToolsPTwbbClean import *

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
      lhep="externalLHEProducer"
      )

addEventSummary(process,True)

process.TFileService.fileName = cms.string('tt_synch%s.root'%(nr)) 
