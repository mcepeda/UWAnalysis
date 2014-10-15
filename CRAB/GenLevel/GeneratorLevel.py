import FWCore.ParameterSet.Config as cms
import sys
import os
from FWCore.ParameterSet.VarParsing import VarParsing

sys.setrecursionlimit(100)
process = cms.Process("ANALYSIS")

#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_R_42_V19::All'

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

options = VarParsing ('analysis')
options.outputFile = "2.root"
options.register(
    'CrossSection',
#    12234.4*3*0.05626,
        46.3*3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "CrossSection")

options.register(
    'selectNUP',
    -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "selectNUP")

options.parseArguments()


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#"file:/scratch/mcepeda/WJets_noMPI.root"
#"file:/scratch/mcepeda/BB_Gen.root"

#"/store/mc/Summer12_DR53X/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola/AODSIM/PU_S10_START53_V7A-v1/00000/001C0D8C-F40D-E211-9302-002590200814.root",
options.inputFiles
#"file:/nfs_scratch/mcepeda/testbb.root"
)
)
process.source.inputCommands=cms.untracked.vstring(
                        'keep *', 'drop *_finalState*_*_*',
                        'drop *_patFinalStateEvent*_*_*',

)

from UWAnalysis.Configuration.tools.analysisToolsPT import *
TriggerPaths="null"

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.JetProducers.ak5GenJets_cfi")

process.genParticlesNoNeutrinos = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *", # this is the default
        "keep status = 1 ",
        "drop abs(pdgId) = 12", "drop abs(pdgId) = 14", "drop abs(pdgId) = 16" 
    )
)

process.neutrinos = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *", # this is the default
        "keep abs(pdgId) = 12 & status =1 ", "keep abs(pdgId) = 14 & status =  1", "keep abs(pdgId) = 16 & status = 1"
    )
)



process.GenJetParameters.src=cms.InputTag("genParticlesNoNeutrinos")

process.ak5GenJetsNoNu = cms.EDProducer(
     "FastjetJetProducer",
     GenJetParameters,
     AnomalousCellParameters,
     jetAlgorithm = cms.string("AntiKt"),
     rParam       = cms.double(0.5)
     )

process.partonJetsSeq=cms.Sequence(process.genParticlesNoNeutrinos*process.ak5GenJetsNoNu*process.neutrinos)

process.LHEFilter = cms.EDFilter("GenFilterLHE",
        PartonMultiplicity=cms.untracked.int32(0),
        PartonMultiplicityMAX=cms.untracked.int32(20),
)

# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
      # Fix MADGRAPH if buggy (this PDF set will also appear on output, 
      # so only two more PDF sets can be added in PdfSetNames if not "")
      FixMADGRAPH = cms.untracked.int32(3),
      GenTag = cms.untracked.InputTag("genParticles"),
      PdfInfoTag = cms.untracked.InputTag("generator"),
      PdfSetNames = cms.untracked.vstring(
             "CT10.LHgrid",
#             "MSTW2008nlo68cl.LHgrid",
#             "cteq66.LHgrid",
             "MSTW2008nnlo68cl.LHgrid",
             "cteq6ll.LHpdf" ,
      )
)

process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                  quarkId = cms.uint32(5)
                                  )

process.btaggedGenJets= cms.EDProducer('GENJetBtagger',
                                src=cms.InputTag("ak5GenJetsNoNu"),
                                srcSIMB=cms.InputTag("bhadrons")
)

process.cleanGenJets= cms.EDProducer('GenJetCleaner',
                                src=cms.InputTag("btaggedGenJets"),
                                srcDressedLepton=cms.InputTag("dressedLeptons"),
                                onlyGoodLeptons=cms.bool(True)
)

process.dressedLeptons = cms.EDProducer("DressedGenLeptonsProducer",
                                src= cms.InputTag("genParticles"),
#                                genPDGID= cms.untracked.int32(0), # 0 --> muons and electrons,  13 --> muons  11--> electrons
#                                genSTATUS= cms.untracked.int32(1)
                                 minPT = cms.untracked.double(0),
                                 maxETA = cms.untracked.double(10),
)


process.cleanGenJetsAndrea= cms.EDProducer('GenJetCleaner',
                                src=cms.InputTag("btaggedGenJets"),
                                srcDressedLepton=cms.InputTag("dressedLeptons"),
                                onlyGoodLeptons=cms.bool(False),
                                excludeTaus=cms.untracked.bool(True)
)




process.Filter=cms.Sequence(process.partonJetsSeq*process.dressedLeptons*process.bhadrons*process.btaggedGenJets*process.cleanGenJets*process.cleanGenJetsAndrea*process.LHEFilter)#*process.pdfWeights)
process.eventSelectionMC = cms.Path(process.Filter) ##changing to multiples see below

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("genParticles")
)

#process.list=cms.Path(process.printTree)



process.crossSectionJet25Lepton30=cms.EDFilter("GenStudy",
      DressedLepton=cms.untracked.InputTag("dressedLeptons"),
      Jets=cms.untracked.InputTag("cleanGenJets"),
      SIMB=cms.untracked.InputTag("bhadrons"),
      Neutrinos=cms.untracked.InputTag("neutrinos"),
      filterByNUP=cms.untracked.int32(options.selectNUP),
      selLeptonPDGID=cms.untracked.int32(13),
      ptCutJetis=cms.untracked.double(25),
      ptCutLepton=cms.untracked.double(30),
      inputCrossSection=cms.untracked.double(options.CrossSection)
)

process.crossSectionJet25Lepton30Andrea=cms.EDFilter("GenStudyAndrea",
      DressedLepton=cms.untracked.InputTag("dressedLeptons"),
      Jets=cms.untracked.InputTag("cleanGenJets"),
      SIMB=cms.untracked.InputTag("bhadrons"),
      Neutrinos=cms.untracked.InputTag("neutrinos"),
      filterByNUP=cms.untracked.int32(options.selectNUP),
      selLeptonPDGID=cms.untracked.int32(13),
      ptCutJetis=cms.untracked.double(25),
      ptCutLepton=cms.untracked.double(30),
      inputCrossSection=cms.untracked.double(options.CrossSection)
)

#process.pAll=cms.Path(process.crossSectionJet25Lepton30)
process.eventSelectionMC = cms.Path(process.Filter*process.crossSectionJet25Lepton30Andrea)

from UWAnalysis.Configuration.tools.ntupleToolsPTGENwbb import *
addMCEventTree(process,'MC',src="genParticles",srcMuons="dressedLeptons",srcTaggedGenJets="btaggedGenJets")





process.out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring( 
         'keep *'),
     fileName = cms.untracked.string('test.root')
)
 
#process.outpath = cms.EndPath(process.out)
