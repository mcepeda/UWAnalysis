import FWCore.ParameterSet.Config as cms
from CommonTools.ParticleFlow.Isolation.tools_cfi import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *
import sys
from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from JetMETCorrections.Type1MET.pfMETCorrections_cff import *

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

def defaultReconstructionPT(process,triggerProcess = 'HLT',triggerPaths = ['HLT_Mu9','HLT_Mu11_PFTau15_v1','HLT_Mu11_PFTau15_v1','HLT_Mu11_PFTau15_v2','HLT_Mu15_v1','HLT_Mu15_v2'],itsMC=False,itsData=False):
  process.load("UWAnalysis.Configuration.startUpSequence_cff")
  process.load("Configuration.Geometry.GeometryIdeal_cff")
  process.load("Configuration.StandardSequences.MagneticField_cff")
  process.load("Configuration.StandardSequences.Services_cff")
  process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
  process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
  process.load("DQMServices.Core.DQM_cfg")
  process.load("DQMServices.Components.DQMEnvironment_cfi")

  #Make the TriggerPaths Global variable to be accesed by the ntuples
  global TriggerPaths
  TriggerPaths= triggerPaths
  process.analysisSequence = cms.Sequence()

  #Add trigger Matching
  muonTriggerMatchPT(process,triggerProcess) 
  electronTriggerMatchPT(process,triggerProcess)  
  tauTriggerMatchPT(process,triggerProcess)    
  goodVertexFilter(process)       
  
  # leptons
  #rochesterCorrector(process,muons="cleanPatMuons")
  #muScleCorrector(process,muons="cleanPatMuons",isMC=itsMC)
  #electronEnergyCorrector(process,'cleanPatElectrons')

  muonIDer(process,muons="cleanPatMuons")
  eleIsolater(process,electrons="cleanPatElectrons")
  leptonSFer(process,muons="IDedMuons",electrons="IsoedElectrons")
  #leptonSFer(process,muons="IDedMuons",electrons="cleanPatElectrons")

  # met
  metAndmT(process,met='systematicsMET',muons="muAvecSF",electrons="eleAvecSF",isMC=itsMC)
  wgtPUembedding(process,"metCorrected","primaryVertexFilter","addPileupInfo")
  #vertexEmbedding(process,"IDedMuons","primaryVertexFilter","addPileupInfo")

  # jets
  if itsMC:
   genLevel(process)
   resolutionSmearJets(process,jets='selectedPatJetsAK5chsPF')
  elif itsData:
   ReNameJetColl(process,inputJets='selectedPatJetsAK5chsPF')

  jetOverloading(process,"resolutionSmearedJets")
  jetPUEmbedding(process,"patOverloadedJets")

   # btag
  SVReconstruction(process,"patPUEmbeddedJets",isMC=itsMC,isData=itsData)
  #SVReconstruction(process,"patPUEmbeddedJets","rochCorMuons",isMC=itsMC,isData=itsData)

   # cleaning
  applyDefaultSelectionsPT(process,"patBRecoJets","mtMuons","mtEles")
  #applyDefaultSelectionsPT(process,"patCSVreweightedJets","mtMuons","mtEles")

  process.runAnalysisSequence = cms.Path(process.analysisSequence)
  #jetReading(process,jets='patBRecoJets')

  #mvaMet(process) #Build MVA MET
  IVF(process)
  if itsMC:
   genHadrons(process) 
   GenBWeight(process) 


def muonTriggerMatchPT(process,triggerProcess):
   process.triggeredPatMuons = cms.EDProducer("MuonTriggerMatcher",
    src = cms.InputTag("cleanPatMuons"),
    trigEvent = cms.InputTag("hltTriggerSummaryAOD"),
    filters = cms.VInputTag(
     cms.InputTag('hltSingleMuIsoL3IsoFiltered12','',triggerProcess),
     cms.InputTag('hltSingleMuIsoL3IsoFiltered15','',triggerProcess),
     cms.InputTag('hltSingleMuIsoL3IsoFiltered24','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu15IsoPFTau15','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu12IsoPFTau10','',triggerProcess),
     cms.InputTag('hltSingleMuIsoL3IsoFiltered17','',triggerProcess),
     cms.InputTag('hltSingleMuIsoL1s14L3IsoFiltered15eta2p1',"",triggerProcess),
     cms.InputTag('hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f24L3IsoFiltered','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu15IsoPFTau20','',triggerProcess),
     cms.InputTag('hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoFiltered10','',triggerProcess),
     cms.InputTag('hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu18LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu17LooseIsoPFTau20','',triggerProcess)
    ),
    pdgId = cms.int32(13)
   )  
   process.triggeredPatMuonSeq = cms.Sequence(process.triggeredPatMuons)
   process.triggeredPatMuonPath=cms.Path(process.triggeredPatMuonSeq)


def electronTriggerMatchPT(process,triggerProcess):
   process.triggeredPatElectronsL = cms.EDProducer("ElectronTriggerMatcher",
    src = cms.InputTag("cleanPatElectrons"),
    trigEvent = cms.InputTag("hltTriggerSummaryAOD"),
    filters = cms.VInputTag(
     cms.InputTag('hltEle17CaloIdLCaloIsoVLPixelMatchFilterDoubleEG125','',triggerProcess),
    ),
    pdgId = cms.int32(0)
   )
   process.triggeredPatElectrons = cms.EDProducer("ElectronTriggerMatcher",
    src = cms.InputTag("triggeredPatElectronsL"),
    trigEvent = cms.InputTag("hltTriggerSummaryAOD"),
    filters = cms.VInputTag(
     cms.InputTag('hltOverlapFilterIsoEle15IsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle15TightIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle18MediumIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle18TightIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle18IsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20MediumIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20WP90LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter','',triggerProcess),
     cms.InputTag('hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1IsoEG18OrEG20','',triggerProcess)
    ),
    pdgId = cms.int32(11)
   )   
   process.triggeredPatElectronsSeq = cms.Sequence(process.triggeredPatElectronsL*process.triggeredPatElectrons)
   process.triggeredPatElectronsPath=cms.Path(process.triggeredPatElectronsSeq)


def tauTriggerMatchPT(process,triggerProcess):
   strTrig=''
   for i in TriggerPaths:
    if i==TriggerPaths[0]:
      strTrig+='path(\"'+i+'\")'
    else:
      strTrig+='|| path(\"'+i+'\")'
   #Match With The triggers
   process.preTriggeredPatTaus = cms.EDProducer("TauTriggerMatcher",
    src = cms.InputTag("cleanPatTaus"),
    trigEvent = cms.InputTag("hltTriggerSummaryAOD"),
    filters = cms.VInputTag(
     cms.InputTag('hltFilterDoubleIsoPFTau20Trk5LeadTrack5IsolationL1HLTMatched','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu18LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20WP90LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu17LooseIsoPFTau20','',triggerProcess)
    ),
    pdgId = cms.int32(0)
   )
   process.triggeredPatTaus = cms.EDProducer("TauTriggerMatcher",
    src = cms.InputTag("preTriggeredPatTaus"),
    trigEvent = cms.InputTag("hltTriggerSummaryAOD"),
    filters = cms.VInputTag(
     cms.InputTag('hltOverlapFilterIsoMu15IsoPFTau15','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu15IsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu15MediumIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu15TightIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterMu15IsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu12IsoPFTau10','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu15IsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle15IsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle15TightIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle18MediumIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20MediumIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu18LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoEle20WP90LooseIsoPFTau20','',triggerProcess),
     cms.InputTag('hltOverlapFilterIsoMu17LooseIsoPFTau20','',triggerProcess)
    ),
    pdgId = cms.int32(15)
   )                                            
   process.triggeredPatTausSeq = cms.Sequence(process.preTriggeredPatTaus*process.triggeredPatTaus)
#   process.triggeredPatTausPath=cms.Path(process.triggeredPatTausSeq)


def goodVertexFilter(process):
  process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag('offlinePrimaryVertices'),
   cut = cms.string('ndof()>4&&position().rho()<2&&abs(z())<24'),
   filter = cms.bool(False)
  )
  process.goodVertexFilterSeq = cms.Sequence(process.primaryVertexFilter)
  process.goodVertexFilterPath=cms.Path(process.goodVertexFilterSeq)


def genLevel(process):
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
  process.genSeq = cms.Sequence(process.partonJetsSeq*process.dressedLeptons*process.bhadrons*process.btaggedGenJets*process.cleanGenJets*process.cleanGenJetsAndrea*process.LHEFilter)
  process.genPath = cms.Path(process.genSeq)
  return process.genPath


 


def reRunJets(process,isMC=False,isData=False,L1Only=False): 
  process.load("RecoJets.Configuration.RecoPFJets_cff")
  import PhysicsTools.PatAlgos.tools.jetTools as jettools
  process.load("PhysicsTools.PatAlgos.patSequences_cff")
  process.load("UWAnalysis.Configuration.tools.patJetPUId_cfi")
  process.pileupJetIdProducer.applyJec = cms.bool(True)

  process.simpleSecondaryVertex = cms.ESProducer("SimpleSecondaryVertexESProducer",
      use3d = cms.bool(True),
      unBoost = cms.bool(False),
      useSignificance = cms.bool(True),
      minTracks = cms.uint32(2)
  )

  # define the b-tag squences for offline reconstruction
  process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
  process.load("RecoBTau.JetTagComputer.combinedMVA_cff")
  process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
  process.load('RecoBTag/SecondaryVertex/bToCharmDecayVertexMerger_cfi')

  btag_options = {'doBTagging': True}
  btag_options['btagInfo'] = [
     'impactParameterTagInfos',
     'secondaryVertexTagInfos',
     'softMuonTagInfos'
  ]
  btag_options['btagdiscriminators'] = [
     'trackCountingHighEffBJetTags',
     'simpleSecondaryVertexHighEffBJetTags',
     'combinedSecondaryVertexMVABJetTags',
     'combinedSecondaryVertexBJetTags',
     'jetBProbabilityBJetTags', 
     'jetProbabilityBJetTags',
  ]

  if isMC:
   jec = ['L1FastJet', 'L2Relative', 'L3Absolute'] 
  elif isData:
   jec = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

  jettools.switchJetCollection(
        process,
        cms.InputTag('ak5PFJets'),
        doJTA=False,
        jetCorrLabel=('AK5PF', jec), 
        doType1MET=False,
        doJetID=True,
        genJetCollection=cms.InputTag("ak5GenJets"),
        outputModules=[],
        **btag_options
  )
  process.patJets.embedPFCandidates = False
  process.patJets.embedCaloTowers = False
  if isMC:
    process.patJets.embedGenJetMatch = True
  elif isData:
    process.patJets.embedGenJetMatch = False
  process.patJets.addAssociatedTracks = False
  if isMC:
    process.patJets.embedGenPartonMatch = True
  elif isData:
    process.patJets.embedGenPartonMatch = False
  process.patJets.tagInfoSources = cms.VInputTag(
   cms.InputTag("impactParameterTagInfos"),
   cms.InputTag("secondaryVertexTagInfos"),
   cms.InputTag("secondaryVertexNegativeTagInfos")
  ) 
  process.patJets.discriminatorSources = cms.VInputTag(
   cms.InputTag("jetBProbabilityBJetTags"),
   cms.InputTag("jetProbabilityBJetTags"),
   cms.InputTag("trackCountingHighPurBJetTags"),
   cms.InputTag("trackCountingHighEffBJetTags"),
   cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),
   cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
   cms.InputTag("combinedSecondaryVertexBJetTags"),
   cms.InputTag("combinedSecondaryVertexMVABJetTags"),
   cms.InputTag("simpleInclusiveSecondaryVertexHighEffBJetTags"),
   cms.InputTag("simpleInclusiveSecondaryVertexHighPurBJetTags")
  )
  process.patJets.trackAssociationSource = cms.InputTag("ak5JetTracksAssociatorAtVertex")
  process.patJets.addBTagInfo = cms.bool(True)
  process.patJets.addDiscriminators = cms.bool(True)
  process.patJets.addTagInfos = cms.bool(True)

  process.ak5JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
              tracks       = cms.InputTag("generalTracks"),
              jets         = cms.InputTag("ak5PFJets"),
              coneSize     = cms.double(0.5)
  )

  process.patJetCharge.src = cms.InputTag("ak5JetTracksAssociatorAtVertex")
  process.load("PhysicsTools.PatAlgos.patSequences_cff")
  process.NewSelectedPatJets = process.selectedPatJets.clone(src = cms.InputTag("patJetId"))
  process.L1OnlyJets = process.selectedPatJets.clone(src = cms.InputTag("patJetId"))
  print('')

#  process.load('UWAnalysis.Configuration.tools.RecoBTag_cff')
#  All the btagging bits should go to a btag.cff 

  process.btagging = cms.Sequence(
   (
    # impact parameters and IP-only algorithms
    impactParameterTagInfos *
    (trackCountingHighEffBJetTags +
     trackCountingHighPurBJetTags +
     jetProbabilityBJetTags +
     jetBProbabilityBJetTags +

     # SV tag infos depending on IP tag infos, and SV (+IP) based algos
     secondaryVertexTagInfos *
     (simpleSecondaryVertexHighEffBJetTags +
      simpleSecondaryVertexHighPurBJetTags +
      combinedSecondaryVertexBJetTags +
      combinedSecondaryVertexMVABJetTags
     )+
     secondaryVertexNegativeTagInfos*simpleSecondaryVertexNegativeHighEffBJetTags
     +
     ghostTrackVertexTagInfos *
     ghostTrackBJetTags
    )
   )
  )

  process.makePatJetsData = cms.Sequence(
   process.patJetCorrFactors+
   process.patJetCharge+
   process.patJets)

  if isMC:
   process.reRunJetSeq = cms.Sequence(
    process.inclusiveVertexing *
    process.inclusiveMergedVerticesFiltered *
    process.bToCharmDecayVertexMerged*
    process.ak5PFJets*
    process.ak5JetTracksAssociatorAtVertex*
    process.btagging*
    process.inclusiveSecondaryVertexFinderTagInfosFiltered*
    process.simpleInclusiveSecondaryVertexHighEffBJetTags *
    process.simpleInclusiveSecondaryVertexHighPurBJetTags *
    process.pileupJetIdProducer*
    process.makePatJets*
    process.patJetsPUID*
    process.patJetId*
    process.NewSelectedPatJets
   )
  elif isData:
   process.reRunJetSeq = cms.Sequence(
    process.inclusiveVertexing *
    process.inclusiveMergedVerticesFiltered *
    process.bToCharmDecayVertexMerged*
    process.ak5PFJets*
    process.ak5JetTracksAssociatorAtVertex*
    process.btagging*
    process.inclusiveSecondaryVertexFinderTagInfosFiltered*
    process.simpleInclusiveSecondaryVertexHighEffBJetTags *
    process.simpleInclusiveSecondaryVertexHighPurBJetTags *
    process.pileupJetIdProducer*
    process.makePatJetsData*
    process.patJetsPUID*
    process.patJetId*
    process.NewSelectedPatJets
   )
  process.reRunJetPath = cms.Path(process.reRunJetSeq)
  return process.reRunJetPath


#def resolutionSmearJets(process,jets='NewSelectedPatJets',genJets='ak5GenJets'):
def resolutionSmearJets(process,jets='NewSelectedPatJets',genJets='ak5GenJetsNoNu'):
  process.resolutionSmearedJets = cms.EDProducer("PATJetSmearer",
   src = cms.InputTag(jets),
   genJet = cms.InputTag(genJets)
  )
  process.resolutionSmearedJetsSeq = cms.Sequence(process.resolutionSmearedJets)
  process.resolutionSmearedJetsPath = cms.Path(process.resolutionSmearedJetsSeq)
  return process.resolutionSmearedJetsPath


def ReNameJetColl(process, inputJets="selectedPatJets"):
  process.load("PhysicsTools.PatAlgos.patSequences_cff")
  process.resolutionSmearedJets = process.selectedPatJets.clone(src = cms.InputTag(inputJets))
  process.reNameJetSeq = cms.Sequence(process.resolutionSmearedJets)
  process.reNameJetPath = cms.Path(process.reNameJetSeq)
  return process.reNameJetPath


def PFtoPAT(process,inputJets='ak5PFchsJets'):
  process.tomsPatJets = cms.EDProducer("PFtoPAT",
   src=cms.InputTag(inputJets)
  )
  process.tomsPatJetsSeq = cms.Sequence(process.tomsPatJets)
  process.tomsPatJetsPath= cms.Path(process.tomsPatJetsSeq)
  return process.tomsPatJetsPath


def metRenamer(process,met='systematicsMET'):
  process.metTypeOne = cms.EDProducer("PATMETRenamer",
   src = cms.InputTag(met),
  )
  process.metTypeOneSeq = cms.Sequence(process.metTypeOne)
  process.metTypeOnePath= cms.Path(process.metTypeOneSeq)
  return process.metTypeOnePath


def leptonSFer(process,muons="IDedMuons",electrons="cleanPatElectrons"):
  process.eleAvecSF = cms.EDProducer("PATeleSFer",
   src = cms.InputTag(electrons),
  )
  process.muAvecSF = cms.EDProducer("PATmuonSFer",
   src = cms.InputTag(muons),
  )
  
  process.leptonSFSeq = cms.Sequence(process.eleAvecSF + process.muAvecSF) 
  process.leptonSFPath= cms.Path(process.leptonSFSeq)
  return process.leptonSFPath

def metAndmT(process,met='systematicsMET',muons="IDedMuons",electrons="cleanPatElectrons",isMC=False):
  process.metCorrected = cms.EDProducer("PATMETCorrector",
   srcMET = cms.InputTag(met),
   isMC = cms.bool(isMC)
  )
  process.mtMuons = cms.EDProducer("PATmTmuonCalculator",
   srcMuon = cms.InputTag( muons ),
   srcMET = cms.InputTag( "metCorrected" )
  )
  process.mtEles = cms.EDProducer("PATmTeleCalculator",
   srcEle = cms.InputTag( electrons ),
   srcMET = cms.InputTag( "metCorrected" ) 
  )
  process.metCorrectedSeq = cms.Sequence(process.metCorrected)
  process.mtMuonsSeq = cms.Sequence(process.mtMuons)
  process.mtElesSeq = cms.Sequence (process.mtEles)
  process.metAndmTPath = cms.Path(process.metCorrected + process.mtMuonsSeq + process.mtElesSeq)
  return process.metAndmTPath


def metCorrector(process,met='systematicsMET',jets123='NewSelectedPatJets',isMC=False):
  process.selectedVerticesForMEtCorr = cms.EDFilter("VertexSelector",
   src = cms.InputTag('offlinePrimaryVertices'),
   cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
   filter = cms.bool(False)
  )

  process.metCorrected = cms.EDProducer("PATMETCorrector",
   srcMET = cms.InputTag(met),
   srcJ123 = cms.InputTag(jets123),
   srcVert = cms.InputTag("selectedVerticesForMEtCorr"),
   mc = cms.bool(isMC)
  )
  process.metCorrectedSeq = cms.Sequence(process.selectedVerticesForMEtCorr * process.metCorrected)
  process.metCorrectedPath= cms.Path(process.metCorrectedSeq)
  return process.metCorrectedPath


def officialMetCorrector(process,met='systematicsMET',jets='NewSelectedPatJets',isMC=False):
  process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
  #process.patPFMet.metSource = cms.InputTag(met)
  process.selectedPatJetsForMETtype1p2Corr.src = cms.InputTag(jets)
  process.selectedPatJetsForMETtype2Corr.src = cms.InputTag(jets)
  if isMC: process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L3Absolute")
  else: process.patPFJetMETtype1p2Corr.jetCorrLabel = cms.string("L2L3Residual")
  process.patType1CorrectedPFMet.src = cms.InputTag(met)
  process.patType1p2CorrectedPFMet.src = cms.InputTag(met)
  process.metCorrSeq = cms.Sequence(
  #  process.patPFMet
  # * process.pfCandsNotInJet
    process.selectedPatJetsForMETtype1p2Corr
   * process.selectedPatJetsForMETtype2Corr
   * process.patPFJetMETtype1p2Corr
   * process.patPFJetMETtype2Corr
   #* process.pfCandMETcorr
   * process.patType1CorrectedPFMet
   #* process.patType1p2CorrectedPFMet
  )
  process.metCorrPath= cms.Path(process.metCorrSeq)
  return process.metCorrPath


def updatedOfficialMetCorrector(process,isMC=False):
  process.load("JetMETCorrections.Type1MET.correctionTermsCaloMet_cff")
  process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
  if isMC: process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
  else: process.corrPfMetType1.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
  process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff")
  process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff")
  process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
  if isMC: process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
  else: process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
  process.load("JetMETCorrections.Type1MET.correctedMet_cff")
  process.metCorrPath = cms.Path(
  process.correctionTermsPfMetType1Type2 +
  process.correctionTermsPfMetType0RecoTrack +
  process.correctionTermsPfMetType0PFCandidate +
  process.correctionTermsPfMetShiftXY +
  process.correctionTermsCaloMet +
  process.caloMetT1 +
  process.caloMetT1T2 +
  process.pfMetT0rt +
  process.pfMetT0rtT1 +
  process.pfMetT0pc +
  process.pfMetT0pcT1 +
  process.pfMetT0rtTxy +
  process.pfMetT0rtT1Txy +
  process.pfMetT0pcTxy +
  process.pfMetT0pcT1Txy +
  process.pfMetT1 +
  process.pfMetT1Txy
  )
  return process.metCorrPath


def jetMCMatching(process,jets):
  process.patJetMCMatched = cms.EDProducer('MCBHadronJetProducer',
   src = cms.InputTag(jets)
  )
  process.jetMCMatchedSeq = cms.Sequence(process.patJetMCMatched)
  process.jetMCMatchedPath=cms.Path(process.jetMCMatchedSeq)
  return process.jetMCMatchedPath


def jetOverloading(process,jets):
  process.patOverloadedJets = cms.EDProducer('PATJetOverloader',
   src = cms.InputTag(jets),
   leptons = cms.InputTag("cleanPatMuons"),
   vertices=cms.InputTag("offlinePrimaryVertices")
  )
  process.jetOverloadingSeq = cms.Sequence(process.patOverloadedJets)
  process.jetOverloadingPath=cms.Path(process.jetOverloadingSeq)
  return process.jetOverloadingPath


def jetPUEmbedding(process,jets):
  process.patPUEmbeddedJets = cms.EDProducer('PATSimplePUID',
   src = cms.InputTag(jets),
   vertices=cms.InputTag("offlinePrimaryVertices")
  )
  process.jetPUEmbeddingSeq = cms.Sequence(process.patPUEmbeddedJets)
  process.jetPUEmbeddingPath=cms.Path(process.jetPUEmbeddingSeq)
  return process.jetPUEmbeddingPath


def rochesterCorrector(process,muons="cleanPatMuons",rochCor="RochCor2012"):
  process.rochCorMuons = cms.EDProducer("PATMuonRochesterCorrector",
   src = cms.InputTag( muons ),
   corr_type = cms.string( "p4_" + rochCor )
  )
  process.rochCorMuonsSeq = cms.Sequence(process.rochCorMuons)
  process.rochCorMuonsPath= cms.Path(process.rochCorMuonsSeq)
  return process.rochCorMuonsPath


def muScleCorrector(process,muons="cleanPatMuons",isMC=False):
  if isMC: 
    muScleID = "Summer12_DR53X_smearReReco"
  else:
    muScleID = "Data2012_53X_ReReco"
  process.muScleMuons = cms.EDProducer("MuScleFitPATMuonCorrector",
                        src = cms.InputTag(muons),
                        debug = cms.bool(False),
                        identifier = cms.string(muScleID),
                        applySmearing = cms.bool(isMC),
                        fakeSmearing = cms.bool(False)
                        )
  process.muScleMuonsSeq = cms.Sequence(process.muScleMuons)
  process.muScleMuonsPath= cms.Path(process.muScleMuonsSeq)
  return process.muScleMuonsPath


def electronEnergyCorrector(process,electrons = 'cleanPatElectrons'):
  process.energyCorrectedElectrons = cms.EDProducer("RegressionEnergyPatElectronProducer",
    nameEnergyErrorReg = cms.string('eneErrorRegForGsfEle'),
    rhoCollection = cms.InputTag("kt6PFJets","rho","RECO"),
    inputElectronsTag = cms.InputTag(electrons),
    nameEnergyReg = cms.string('eneRegForGsfEle'),
    produceValueMaps = cms.bool(False),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    inputCollectionType = cms.uint32(1),
    regressionInputFile = cms.string('EgammaAnalysis/ElectronTools/data/eleEnergyRegWeights_WithSubClusters_VApr15.root'),
    energyRegressionType = cms.uint32(2),
    debug = cms.untracked.bool(False),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    useRecHitCollections = cms.bool(False)
    #useRecHitCollections = cms.bool(True)
  )
  process.energyCorrectedElectronsSeq = cms.Sequence(process.energyCorrectedElectrons)
  process.energyCorrectedElectronsPath= cms.Path(process.energyCorrectedElectronsSeq)
  return process.energyCorrectedElectronsPath


def muonIDer(process,muons="cleanPatMuons"):
  process.IDedMuons = cms.EDProducer("PATMuonIDer",
   src = cms.InputTag( muons ),
  )
  process.IDedMuonSeq = cms.Sequence(process.IDedMuons)
  process.IDedMuonPath= cms.Path(process.IDedMuonSeq)
  return process.IDedMuonPath


def eleIsolater(process,electrons="cleanPatElectrons"):
  process.IsoedElectrons = cms.EDProducer("PATeleIsoEmbedder",
   src = cms.InputTag( electrons ),
  )
  process.IsoedElectronsSeq = cms.Sequence(process.IsoedElectrons)
  process.IsoedElectronsPath= cms.Path(process.IsoedElectronsSeq)
  return process.IsoedElectronsPath


#def vertexEmbedding(process,muons="IDedMuons",vertices="primaryVertexFilter",pu="addPileupInfo"):
#  process.vertexEmbeddedMuons = cms.EDProducer("PATMuonVertexWeighter",
#   srcMuon = cms.InputTag( muons ),
#   srcVert = cms.InputTag( vertices ),
#   srcPU   = cms.InputTag( pu )
#  )
#  process.vertexEmbeddedMuonSeq = cms.Sequence(process.vertexEmbeddedMuons)
#  process.vertexEmbeddedMuonPath= cms.Path(process.vertexEmbeddedMuonSeq)
#  return process.vertexEmbeddedMuonPath

def wgtPUembedding(process,met="",vertices="primaryVertexFilter",pu="addPileupInfo"):
  process.wgtPUembedMET = cms.EDProducer("PATMETVertexWeighter",
   srcMET = cms.InputTag( met ),
   srcVert = cms.InputTag( vertices ),
   srcPU   = cms.InputTag( pu )
  )
  process.wgtPUembedMETSeq = cms.Sequence(process.wgtPUembedMET)
  process.wgtPUembedMETPath= cms.Path(process.wgtPUembedMETSeq)
  return process.wgtPUembedMETPath


def mtMaker(process,muons="IDedMuons", electrons="cleanPatElectrons", met="metCorrected"):
  process.mtMuons = cms.EDProducer("PATmTmuonCalculator",
   srcMuon = cms.InputTag( muons ),
   srcMET = cms.InputTag( met )
  )
  process.mtEles = cms.EDProducer("PATmTeleCalculator",
   srcEle = cms.InputTag( electrons ),
   srcMET = cms.InputTag( met ) 
  )
  process.mtMuonsSeq = cms.Sequence(process.mtMuons)
  process.mtElesSeq = cms.Sequence (process.mtEles)
  process.mtMakerPath= cms.Path(process.mtMuonsSeq + process.mtElesSeq)
  return process.mtMakerPath


#def SVReconstruction(process,jets,muons,isMC=False,isData=False): 
def SVReconstruction(process,jets,isMC=False,isData=False): 
  process.SVFoundJets=cms.EDProducer("bTagSFer",
    src = cms.InputTag(jets),
    #leptons = cms.InputTag(muons),
    vertices=cms.InputTag("offlinePrimaryVertices") 
  )
  process.patSSVJets=cms.EDProducer("PATSSVJetEmbedder", 
   src = cms.InputTag("SVFoundJets") 
  ) 
  process.patSSVJets2=cms.EDProducer("PATCSVVertex",
   src = cms.InputTag("patSSVJets")
  )
  process.patBpmRecoJets = cms.EDProducer('PATJetBpmReco', 
   src = cms.InputTag("patSSVJets2"), 
   vertices=cms.InputTag("offlinePrimaryVertices") 
  )
  process.patBRecoJets = cms.EDProducer('PATJetBReco', 
   src = cms.InputTag("patBpmRecoJets"), 
   vertices=cms.InputTag("offlinePrimaryVertices") 
  ) 
  #process.patCSVreweightedJets = cms.EDProducer('PATJetCSVreweight',
  # src = cms.InputTag("patBRecoJets")
  #)
  
  process.BReconstruction = cms.Sequence(
   process.SVFoundJets *
   process.patSSVJets *
   process.patSSVJets2 *
   process.patBpmRecoJets *
   process.patBRecoJets 
   #process.patCSVreweightedJets
  ) 
  process.createBRecoJets=cms.Path(process.BReconstruction) 
  return process.createBRecoJets 


def applyDefaultSelectionsPT(process,jets,muons,electrons):
  process.selectedPatTaus = cms.EDFilter("PATTauSelector",
   src = cms.InputTag("cleanPatTaus"),
   cut = cms.string('pt>15&&tauID("byLooseIsolationMVA")&&tauID("againstElectronLoose")&&tauID("againstMuonLoose")'),
   filter = cms.bool(False)
   )

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
  #electronIso = '(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et'
  #electronIso = 'userFloat("Iso04")'
  electronIso = 'userFloat("Iso03")'
  electronLooseIso = electronIso+' < 0.15'
  electronTightIso = electronIso+' < 0.10 '
  electronAntiIso  = electronIso+' > 0.15 '

  muonLooseID = ' userFloat("looseID")>0 '
  muonTightID = ' userInt("tightID")==1 '
  muonIso = '(pfIsolationR04().sumChargedHadronPt + max((pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt),0.0))/pt'
  muonLooseIso = muonIso+' < 0.20 '
  muonTightIso = muonIso+' < 0.12 '
  muonAntiIso  = muonIso+' > 0.20 '

  process.allMuons = cms.EDFilter("PATMuonSelector",
   src = cms.InputTag( muons ),
   cut = cms.string(""),
   filter = cms.bool(False)
  )
  process.allElectrons = cms.EDFilter("PATElectronSelector",
   src = cms.InputTag( electrons ),
   cut = cms.string(""),
   filter = cms.bool(False)
  )
  process.goodMuons = cms.EDFilter("PATMuonSelector",
   src = cms.InputTag( muons ),
   cut = cms.string('pt>30 && abs(eta)<2.1 && %s && %s'%(muonTightIso,muonTightID)),
   filter = cms.bool(False),
  )
  process.goodElectrons = cms.EDFilter("PATElectronSelector",
   src = cms.InputTag( electrons ),
   cut = cms.string('pt>30 && abs(eta)<2.1 && %s && %s'%(electronTightIso,electronTightID)),
   #cut = cms.string('pt>30 && ( abs(eta)<1.4442 || ( abs(eta)>1.566 && abs(eta)<2.1 ) ) && %s && %s'%(electronTightIso,electronTightID)),
   filter = cms.bool(False),
  )
  process.vetoMuons = cms.EDFilter("PATMuonSelector",
   src = cms.InputTag( muons ),
   cut = cms.string('pt>10 && abs(eta)<2.4 && %s && %s'%(muonLooseIso,muonLooseID)),
   filter = cms.bool(False),
  )
  process.vetoElectrons = cms.EDFilter("PATElectronSelector",
   src = cms.InputTag( electrons ),
   cut = cms.string('pt>10 && abs(eta)<2.4 && %s && %s'%(electronLooseIso,electronLooseID)),
   filter = cms.bool(False),
  )
  process.qcdMuons = cms.EDFilter("PATMuonSelector",
   src = cms.InputTag( muons ),
   cut = cms.string('pt>30 && abs(eta)<2.1 && %s && %s'%(muonAntiIso,muonTightID)),
   filter = cms.bool(False),
  )
  process.qcdElectrons = cms.EDFilter("PATElectronSelector",
   src = cms.InputTag( electrons ),
   cut = cms.string('pt>30 && abs(eta)<2.1 && %s && %s'%(electronAntiIso,electronTightID)),
   #cut = cms.string('pt>30 && ( abs(eta)<1.4442 || ( abs(eta)>1.566 && abs(eta)<2.1 ) ) && %s && %s'%(electronAntiIso,electronTightID)),
   filter = cms.bool(False),
  )

  process.allJets = cms.EDProducer("PATJetCleaner",
   src = cms.InputTag(jets),
   preselection = cms.string(''),
   checkOverlaps = cms.PSet( 
       muons = cms.PSet(
        src = cms.InputTag(muons),
        algorithm = cms.string("byDeltaR"),
        preselection = cms.string(
          'pt > 30 && abs(eta) < 2.1 && %s && %s'%(muonTightIso,muonTightID)
        ),
        deltaR = cms.double(0.5),
        checkRecoComponents = cms.bool(False),
        pairCut = cms.string(""),
        requireNoOverlaps = cms.bool(False),
       ),
       electrons = cms.PSet(
        src = cms.InputTag( electrons ),
        algorithm = cms.string("byDeltaR"),
        preselection = cms.string(
          'pt > 30 && abs(eta) < 2.1 && %s && %s'%(electronTightIso,electronTightID)
        ),
        deltaR = cms.double(0.5),
        checkRecoComponents = cms.bool(False),
        pairCut = cms.string(""),
        requireNoOverlaps = cms.bool(False),
       ),
   ),
   finalCut = cms.string('')
   )

  process.cleanJets = cms.EDProducer("PATJetCleaner",
   src = cms.InputTag(jets),
   preselection = cms.string('abs(eta)<5.0&&pt>20&&userFloat("idLoose")>0'),
   checkOverlaps = cms.PSet( 
       muons = cms.PSet(
        src = cms.InputTag(muons),
        algorithm = cms.string("byDeltaR"),
        preselection = cms.string(
          'pt > 30 && abs(eta) < 2.1 && %s && %s'%(muonTightIso,muonTightID)
        ),
        deltaR = cms.double(0.5),
        checkRecoComponents = cms.bool(False),
        pairCut = cms.string(""),
        requireNoOverlaps = cms.bool(True),
       ),
       electrons = cms.PSet(
        src = cms.InputTag( electrons ),
        algorithm = cms.string("byDeltaR"),
        preselection = cms.string(
          #'pt > 30 && abs(eta) < 2.1 && %s && %s'%(electronTightIso,electronTightID)
          'pt > 30 && ( abs(eta)<1.4442 || ( abs(eta)>1.566 && abs(eta)<2.1 ) ) && %s && %s'%(electronTightIso,electronTightID)
        ),
        deltaR = cms.double(0.5),
        checkRecoComponents = cms.bool(False),
        pairCut = cms.string(""),
        requireNoOverlaps = cms.bool(True),
       ),
   ),
   finalCut = cms.string('')
   )	

  process.goodJets = process.cleanJets.clone(
   preselection = cms.string("abs(eta)<2.4 && pt>25 && userFloat('idLoose')>0")
  )

  process.fwdJets = process.cleanJets.clone(
   preselection = cms.string("abs(eta)>=2.4 && abs(eta)<5.0 && pt>25 && userFloat('idLoose')>0")
  )

  process.selectedObjectsForSyst = cms.Sequence(
     process.allMuons 
   + process.goodMuons
   + process.vetoMuons 
   + process.qcdMuons 
   + process.allElectrons
   + process.goodElectrons 
   + process.vetoElectrons 
   + process.qcdElectrons 
   + process.allJets
   + process.cleanJets
   + process.goodJets
   + process.fwdJets
) 
  #process.selectedObjectsForSyst = cms.Sequence(process.preselectedPatMuons+process.selectedPatMuons+process.preselectedPatElectrons+process.selectedPatElectrons+process.cleanPatJets + process.goodMuons + process.goodElectrons + process.vetoMuons + process.vetoElectrons + process.qcdMuons + process.qcdElectrons + process.allMuons + process.allElectrons) 
  process.selectedObjPath = cms.Path(process.selectedObjectsForSyst)
#  process.preselectedPatElectrons = cms.EDFilter("PATElectronSelector",
#   src = cms.InputTag( electrons ),
#   cut = cms.string('pt>10 && abs(eta)<2.4 &&'
#     '(('
#      'abs(superCluster.eta) < 1.442 &'
#      'abs(deltaEtaSuperClusterTrackAtVtx) < 0.007 &'
#      'abs(deltaPhiSuperClusterTrackAtVtx) < 0.15 &'
#      'sigmaIetaIeta < 0.01 &'
#      'hadronicOverEm < 0.12'
#     ')|('
#      'abs(superCluster.eta) > 1.566 & abs(superCluster.eta) < 2.5 &'
#      'abs(deltaEtaSuperClusterTrackAtVtx) < 0.009 &'
#      'abs(deltaPhiSuperClusterTrackAtVtx) < 0.10 &'
#      'sigmaIetaIeta < 0.03 &'
#      'hadronicOverEm < 0.10'
#     ')) &'
#     'abs(dB) < 0.02 &'
#     'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &'
#     '(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et < 0.15 &'
#     'passConversionVeto &'
#     'gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 &'
#     'triggerObjectMatches.size >= 0'
#   ), # 
#   filter = cms.bool(False)
#   )
#
#  process.preselectedPatMuons = cms.EDFilter("PATMuonSelector",
#   src = cms.InputTag(muons),
#   cut = cms.string('pt>10 && abs(eta)<2.4 && userFloat("looseID")>0 && (pfIsolationR04().sumChargedHadronPt + max((pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt),0.0))/pt<0.2'), # LooseID, Loose Iso for Muon Veto
#   filter = cms.bool(False)
#   )
#
#  process.selectedPatElectrons = cms.EDFilter("PATElectronSelector",
#   src = cms.InputTag( electrons ),
#   cut = cms.string(
#           'pt > 30 & abs(eta) < 2.4 &'
#           '(('
#           'abs(superCluster.eta) < 1.442 &'
#           'abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &'
#           'abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 &'
#           'sigmaIetaIeta < 0.01 &'
#           'hadronicOverEm < 0.12'
#          ')|('
#           'abs(superCluster.eta) > 1.566 & abs(superCluster.eta) < 2.5 &'
#           'abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 &'
#           'abs(deltaPhiSuperClusterTrackAtVtx) < 0.02 &'
#           'sigmaIetaIeta < 0.03 &'
#           'hadronicOverEm < 0.10'
#          ')) &'
#          'abs(dB) < 0.02 &'
#          'abs(1./ecalEnergy - eSuperClusterOverP/ecalEnergy) < 0.05 &'
#          '(chargedHadronIso + max((neutralHadronIso + photonIso - 0.5*puChargedHadronIso),0.0))/et < 0.10 &'
#          'passConversionVeto &'
#          'gsfTrack.trackerExpectedHitsInner.numberOfHits < 1'
#   ),
#   #cut = cms.string('pt>10 && userFloat("wp95")>0&&(userIso(0)+max(photonIso+neutralHadronIso()-0.5*userIso(2),0.0))/pt()<0.15'),
#   filter = cms.bool(False)
#   )
#  process.selectedPatMuons = cms.EDFilter("PATMuonSelector",
#   src = cms.InputTag(muons),
#   cut = cms.string('pt>10 && userInt("tightID")==1'), # Tight ID, No Iso for GOOD, QCD
#   filter = cms.bool(False)
#   )

  return process.selectedObjPath


def jetReading(process,jets):
  process.patOverloadedJets = cms.EDProducer('PATJetReader',
   src = cms.InputTag(jets),
  )
  process.jetReadingSeq = cms.Sequence(process.patOverloadedJets)
  process.jetReadingPath=cms.Path(process.jetReadingSeq)
  return process.jetReadingPath


def mvaMet(process):
  process.load("RecoMET.METProducers.mvaPFMET_cff")
  process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
  
  process.metEleTaus = cms.EDFilter("PATTauSelector",
                                           src = cms.InputTag("cleanPatTaus"),
                                           cut = cms.string('abs(userFloat("dz"))<0.2&&pt>20&&tauID("byLooseIsolationMVA")&&tauID("againstElectronMedium")&&tauID("againstElectronMVA")&&tauID("againstMuonLoose")&&abs(eta())<2.3'),
                                           filter = cms.bool(False)
   )
  process.metMuTaus = cms.EDFilter("PATTauSelector",
                                           src = cms.InputTag("cleanPatTaus"),
                                           cut = cms.string('abs(userFloat("dz"))<0.2&&pt>20&&tauID("byLooseIsolationMVA")&&tauID("againstElectronLoose")&&tauID("againstMuonTight")&&abs(eta())<2.3'),
                                           filter = cms.bool(False)
   )
  process.metElectrons = cms.EDFilter("PATElectronSelector",
                                           src = cms.InputTag("cleanPatElectrons"),
                                           cut = cms.string('abs(userFloat("dz"))<0.2&&pt>20&&userFloat("wp95V")>0&&(userIso(0)+max(photonIso+neutralHadronIso()-0.5*userIso(2),0.0))/pt()<0.1&&!(userFloat("hasConversion")>0)&&userInt("missingHits")==0&&abs(userFloat("ipDXY"))<0.045&&abs(eta())<2.1'),
                                           filter = cms.bool(False)
   )
  process.metMuons = cms.EDFilter("PATMuonSelector",
                                           src = cms.InputTag("cleanPatMuons"),
                                           cut = cms.string('abs(userFloat("dz"))<0.2&&pt>17&&userInt("tightID")>0&&(userIso(0)+max(photonIso+neutralHadronIso()-0.5*userIso(2),0.0))/pt()<0.1&&abs(eta())<2.1&&abs(userFloat("ipDXY"))<0.045'),
                                           filter = cms.bool(False)
   )
  process.pfMEtMVA2.inputFileNames.U = cms.FileInPath('pharris/MVAMet/data/gbrmet_52.root')
  process.pfMEtMVA2.inputFileNames.DPhi = cms.FileInPath('pharris/MVAMet/data/gbrmetphi_52.root')
  process.pfMEtMVA2.inputFileNames.CovU1 = cms.FileInPath('pharris/MVAMet/data/gbrmetu1cov_52.root')
  process.pfMEtMVA2.inputFileNames.CovU2 = cms.FileInPath('pharris/MVAMet/data/gbrmetu2cov_52.root')
  
  process.mvaMetMuTau = process.pfMEtMVA2.clone()
  process.mvaMetMuTau.srcLeptons = cms.VInputTag('metMuons', 'metMuTaus')
  process.patMVAMetMuTau = process.patMETs.clone(
   metSource = cms.InputTag('mvaMetMuTau'),
   addMuonCorrections = cms.bool(False),
   addGenMET = cms.bool(False)
  )
  process.mvaMetEleTau = process.pfMEtMVA2.clone()
  process.mvaMetEleTau.srcLeptons = cms.VInputTag('metElectrons', 'metEleTaus')
  process.patMVAMetEleTau = process.patMETs.clone(
   metSource = cms.InputTag('mvaMetEleTau'),
   addMuonCorrections = cms.bool(False),
   addGenMET = cms.bool(False)
  )
  process.mvaMetMuTauSequence = cms.Sequence(process.metMuTaus*process.metMuons*process.mvaMetMuTau*process.patMVAMetMuTau)
  process.mvaMetEleTauSequence = cms.Sequence(process.metEleTaus*process.metElectrons*process.mvaMetEleTau*process.patMVAMetEleTau)
  process.analysisSequence = cms.Sequence(process.analysisSequence*process.mvaMetMuTauSequence*process.mvaMetEleTauSequence)


def genHadrons(process):
  process.bhadrons = cms.EDProducer('MCBHadronProducer',
                                    quarkId = cms.uint32(5)
                                    )
  process.cbarCands=cms.EDProducer("GenParticlePruner",
                                   src = cms.InputTag("genParticles"),
                                   select = cms.vstring("keep pdgId = -4",
                                                        )
                                   )
  process.cCands=cms.EDProducer("GenParticlePruner",
                                src = cms.InputTag("genParticles"),
                                select = cms.vstring("keep pdgId= 4 ",
                                                     )
                                )
  process.createSimCollections=cms.Path(process.bhadrons*process.cbarCands*process.cCands)


def GenBWeight(process):
  process.GBWeightProducer=cms.EDProducer("GenBWeightProducer",
	pprop=cms.FileInPath('UWAnalysis/RecoTools/data/bHProp.txt'),
  	p511=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_511_decaytable.txt'),
  	p521=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_521_decaytable.txt'),
  	p531=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_531_decaytable.txt'),
  	p541=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_541_decaytable.txt'),
  	p411=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_411_decaytable.txt'),
  	p421=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_421_decaytable.txt'),
  	p431=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_431_decaytable.txt'),
  	p441=cms.FileInPath('UWAnalysis/RecoTools/data/fHadron_441_decaytable.txt')

  )
  process.GetBDWeights=cms.Path(process.GBWeightProducer)

def IVF(process):
  process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
  process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
  process.load("RecoBTag.SecondaryVertex.bToCharmDecayVertexMerger_cfi")

  process.IVF = process.inclusiveVertexFinder.clone();
  process.IVF.clusterScale = cms.double(1.)
  process.IVF.clusterMinAngleCosine = cms.double(0.5)
  process.IVF.primaryVertices = cms.InputTag("offlinePrimaryVertices")
  process.IVF.seedMin3DIPSignificance = cms.double(1.2)
  process.IVF.vertexMinAngleCosine = cms.double(0.95)

  process.IMV1 = process.vertexMerger.clone()
  process.IMV1.secondaryVertices = cms.InputTag("IVF")

  process.TVA=process.trackVertexArbitrator.clone()
  process.TVA.secondaryVertices = cms.InputTag("IMV1")
  process.TVA.primaryVertices = cms.InputTag("offlinePrimaryVertices")

  process.IMV2 = process.vertexMerger.clone()
  process.IMV2.secondaryVertices = cms.InputTag("TVA")
  process.IMV2.maxFraction = 0.2
  process.IMV2.minSignificance = cms.double(10.)

  process.IMVF = process.bVertexFilter.clone()
  process.IMVF.vertexFilter.multiplicityMin = 2
  process.IMVF.secondaryVertices = cms.InputTag("IMV2")

  process.bCandidates=process.bToCharmDecayVertexMerged.clone()
  process.bCandidates.secondaryVertices = cms.InputTag("IMVF")

  process.LCProducer=cms.EDProducer("BCProducer")
  process.LCProducer.src= cms.InputTag("bCandidates")
  process.LCProducer.primaryVertices= cms.InputTag("offlinePrimaryVertices")
  process.LCProducer.jetSource= cms.InputTag("ak5PFJets")
  process.LCProducer.jetGenSource= cms.InputTag("ak5GenJets")
  process.dump = cms.EDAnalyzer("EventContentAnalyzer")
  process.LCProducer.isMC = cms.bool(False)
  process.createBCandidateCollection=cms.Path(process.IVF*process.IMV1*process.TVA*process.IMV2*process.IMVF*process.bCandidates*process.LCProducer)


  



#def BTAGGING(process):
#
#  process.load("RecoJets.Configuration.RecoPFJets_cff")
#  import PhysicsTools.PatAlgos.tools.jetTools as jettools
#  process.load("PhysicsTools.PatAlgos.patSequences_cff")
#
#  process.simpleSecondaryVertex = cms.ESProducer("SimpleSecondaryVertexESProducer",
#      use3d = cms.bool(True),
#      unBoost = cms.bool(False),
#      useSignificance = cms.bool(True),
#      minTracks = cms.uint32(2)
#  )
#
#
#  # define the b-tag squences for offline reconstruction
#  process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
#  process.load("RecoBTau.JetTagComputer.combinedMVA_cff")
#
#  process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
#  process.load('RecoBTag/SecondaryVertex/bToCharmDecayVertexMerger_cfi')
#    
#  process.ak5JetTracksAssociatorAtVertex= cms.EDProducer("JetTracksAssociatorAtVertex",
#      tracks = cms.InputTag("generalTracks"),
#      jets = cms.InputTag("ak5PFJets"),
#      coneSize = cms.double(0.5)
#  )
#
#  process.load('RecoBTag/ImpactParameter/impactParameter_cff')
#  process.load('RecoBTag/SecondaryVertex/secondaryVertex_cff')
#  
#  process.btagging = cms.Sequence(
#   (
#   # impact parameters and IP-only algorithms
#    impactParameterTagInfos *
#    (trackCountingHighEffBJetTags +
#     trackCountingHighPurBJetTags +
#     jetProbabilityBJetTags +
#     jetBProbabilityBJetTags +
#     # SV tag infos depending on IP tag infos, and SV (+IP) based algos
#     secondaryVertexTagInfos *
#     (simpleSecondaryVertexHighEffBJetTags +
#      simpleSecondaryVertexHighPurBJetTags +
#      combinedSecondaryVertexBJetTags +
#      combinedSecondaryVertexMVABJetTags
#     )+
#     secondaryVertexNegativeTagInfos*
#     simpleSecondaryVertexNegativeHighEffBJetTags+
#     ghostTrackVertexTagInfos *
#     ghostTrackBJetTags
#    )
#   ))
#  
#  process.patJetsForBTagging = cms.EDProducer("PATJetProducer",
#    addJetCharge = cms.bool(False),
#    addGenJetMatch = cms.bool(False),
#    embedPFCandidates = cms.bool(False),
#    embedGenJetMatch = cms.bool(False),
#    addAssociatedTracks = cms.bool(False),
#    addGenPartonMatch = cms.bool(False),
#    genPartonMatch = cms.InputTag(""),
#    addTagInfos = cms.bool(True),
#    addPartonJetMatch = cms.bool(False),
#    embedGenPartonMatch = cms.bool(False),
#    jetSource = cms.InputTag("ak5PFJets"),
#    addEfficiencies = cms.bool(False),
#    trackAssociationSource = cms.InputTag("ak5JetTracksAssociatorAtVertex"),
#    tagInfoSources = cms.VInputTag(
#      #cms.InputTag("impactParameterTagInfos"), 
#      cms.InputTag("secondaryVertexTagInfos"), 
#      cms.InputTag("secondaryVertexNegativeTagInfos"),
#      #cms.InputTag("softMuonTagInfos")
#    ),
#    discriminatorSources = cms.VInputTag(
#      cms.InputTag("jetBProbabilityBJetTags"), 
#      cms.InputTag("jetProbabilityBJetTags"), 
#      cms.InputTag("trackCountingHighPurBJetTags"), 
#      cms.InputTag("trackCountingHighEffBJetTags"),
#      cms.InputTag("simpleSecondaryVertexHighEffBJetTags"), 
#      cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
#      cms.InputTag("combinedSecondaryVertexBJetTags"),
#      cms.InputTag("combinedSecondaryVertexMVABJetTags"),
#      cms.InputTag("simpleInclusiveSecondaryVertexHighEffBJetTags"),
#      cms.InputTag("simpleInclusiveSecondaryVertexHighPurBJetTags"),
#      cms.InputTag("doubleSecondaryVertexHighEffBJetTags"),
#      #cms.InputTag("softMuonBJetTags"),
#      # cms.InputTag("softMuonByPtBJetTags"), cms.InputTag("softMuonByIP3dBJetTags")
#    ),
#    addBTagInfo = cms.bool(True),
#    embedCaloTowers = cms.bool(False),
#    addResolutions = cms.bool(False),
#    getJetMCFlavour = cms.bool(False),
#    addDiscriminators = cms.bool(True),
#    jetChargeSource = cms.InputTag("patJetCharge"),
#    addJetCorrFactors = cms.bool(False),
#    jetIDMap = cms.InputTag("ak5JetID"),
#    addJetID = cms.bool(False)
#  )
#  
#  process.reDOBTAGGING = cms.Path(
#     process.inclusiveVertexing *
#     process.inclusiveMergedVerticesFiltered *
#     process.bToCharmDecayVertexMerged * 
#     process.btagging * 
#     process.inclusiveSecondaryVertexFinderTagInfosFiltered *
#     process.simpleInclusiveSecondaryVertexHighEffBJetTags *
#     process.simpleInclusiveSecondaryVertexHighPurBJetTags *
#     process.doubleSecondaryVertexHighEffBJetTags * 
#     process.patJetsForBTagging
#  )
#
#def createGeneratedTaus(process,decayMode,fiducialCuts):
#  process.generatedTaus = cms.EDFilter("TauGenJetDecayModeSelector",
#                                       src = cms.InputTag("tauGenJets"),
#                                       select = cms.vstring(decayMode),
#                                       filter = cms.bool(False)
#                                       )
#  process.generatedTausInAcceptance = cms.EDFilter("GenJetSelector",
#                                           src = cms.InputTag("generatedTaus"),
#                                           cut = cms.string(fiducialCuts),
#                                           filter = cms.bool(False)
#                                           )
#  process.analysisSequence*= process.generatedTaus
#  process.analysisSequence*= process.generatedTausInAcceptance
#
#def cloneAndReplaceInputTag(process,sequence,oldValue,newValue,postfix):
#  #First Clone the sequence
#  p = cloneProcessingSnippet(process, sequence, postfix)
#  massSearchReplaceAnyInputTag(p,oldValue,newValue )
#  modules = listModules(p)
#  #Change the labels of the counters
#  for mod in modules:
#    if(hasattr(mod,'label')):
#      if mod.label().find('Counter') !=-1 :
#        if(hasattr(mod,'name')):
#          newValue = mod.name.value()+postfix
#          mod.name=cms.string(newValue)
#  return p
#
#def createRecoilSystematics(process,sequence,postfix,metScale,metResolution):
#  #First Clone the sequence
#  p = cloneProcessingSnippet(process, sequence, postfix)
#  modules = listModules(p)
#  #Change the labels of the counters and the smearign modules
#  for mod in modules:
#    if(hasattr(mod,'label')):
#      if mod.label().find('Counter') !=-1 :
#        if(hasattr(mod,'name')):
#          newValue = mod.name.value()+postfix
#          mod.name=cms.string(newValue)
#    if(hasattr(mod,'metCalibration')):
#          mod.metCalibration.shiftScale = cms.untracked.double(metScale)
#          mod.metCalibration.shiftRes = cms.untracked.double(metResolution)
#  return cms.Path(p)
#
#def addTagAndProbePlotter(process,type,name,src,ref,selections,methods,triggers,triggersProbe):
#  process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
#  muonAnalysis = cms.EDAnalyzer(type,
#                                src=cms.InputTag(src),
#                                vertices=cms.InputTag("primaryVertexFilter"),
#                                ref=cms.InputTag(ref),
#                                patTrigger = cms.InputTag("patTrigger"),
#                                id= cms.vstring(selections),
#                                methods= cms.vstring(methods),
#                                triggers = cms.vstring(triggers),
#                                triggersProbe = cms.vstring(triggersProbe),
#  )
#  setattr(process,'tagAndProbe'+name,muonAnalysis)
#  p = cms.EndPath(getattr(process,'tagAndProbe'+name))
#  setattr(process,'tagAndProbePath'+name,p)


def createGeneratedParticles(process,name,commands):
  refObjects = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop * "
    )
   )
  refObjects.select.extend(commands)
  setattr(process,name,refObjects)
  process.analysisSequence*= getattr(process,name)


def createSystematics(process,sequence,postfix,muScale,eScale,tauScale,jetScale,unclusteredScale,electronresb = 0.0, electronrese = 0.0):
  #First Clone the sequence
  p = cloneProcessingSnippet(process, sequence, postfix)
  modules = listModules(p)
  #Change the labels of the counters and the smearign modules
  for mod in modules:
    if(hasattr(mod,'label')):
      if mod.label().find('Counter') !=-1 :
        if(hasattr(mod,'name')):
          newValue = mod.name.value()+postfix
          mod.name=cms.string(newValue)
      if mod.label().find('smearedVetoMuons') !=-1 :
          mod.energyScale = cms.double(muScale)
      if mod.label().find('smearedQCDMuons') !=-1 :
          mod.energyScale = cms.double(muScale)
      if mod.label().find('smearedGoodMuons') !=-1 :
          mod.energyScale = cms.double(muScale)
      if mod.label().find('smearedTaus') !=-1 :
          mod.energyScale = cms.double(tauScale)
      if mod.label().find('smearedVetoElectrons') !=-1 :
          mod.energyScale = cms.double(eScale)
          mod.deltaPtB = cms.double(electronresb)
          mod.deltaPtE = cms.double(electronrese)
      if mod.label().find('smearedQCDElectrons') !=-1 :
          mod.energyScale = cms.double(eScale)
          mod.deltaPtB = cms.double(electronresb)
          mod.deltaPtE = cms.double(electronrese)
      if mod.label().find('smearedGoodElectrons') !=-1 :
          mod.energyScale = cms.double(eScale)
          mod.deltaPtB = cms.double(electronresb)
          mod.deltaPtE = cms.double(electronrese)
      if mod.label().find('smearedGoodJets') !=-1 :
          mod.energyScaleDB = cms.double(jetScale) ##changed from int32 to double
      if mod.label().find('smearedFwdJets') !=-1 :
          mod.energyScaleDB = cms.double(jetScale) ##changed from int32 to double
      if mod.label().find('smearedMET') !=-1 :
          mod.unclusteredScale= cms.double(unclusteredScale)
  return cms.Sequence(p)


def addEventSummary(process,onSkim = False,name = 'summary',path = 'eventSelection'):
   process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
   summary = cms.EDAnalyzer('EventSummary',
                            src =getAllEventCounters(process,getattr(process,path),onSkim)
   )
   setattr(process,name,summary)
   if onSkim:
        process.EDMtoMEConverter = cms.EDAnalyzer("EDMtoMEConverter",
           Name = cms.untracked.string('EDMtoMEConverter'),
           Verbosity = cms.untracked.int32(1), # 0 provides no output # 1 provides basic output
           Frequency = cms.untracked.int32(50),
           convertOnEndLumi = cms.untracked.bool(True),
           convertOnEndRun = cms.untracked.bool(True),
	   runInputTag = cms.InputTag('MEtoEDMConverter', 'MEtoEDMConverterRun'),
	   lumiInputTag = cms.InputTag('MEtoEDMConverter', 'MEtoEDMConverterLumi')
        )
        eventSummaryPath=cms.EndPath(process.EDMtoMEConverter+getattr(process,name))
        setattr(process,name+"Path",eventSummaryPath)
   else:
        eventSummaryPath=cms.EndPath(getattr(process,name))
        setattr(process,name+"Path",eventSummaryPath)


def getAllEventCounters(process,path,onSkim = False):
        stringList = []
        if onSkim:
          stringList.append('processedEvents')

        modules = listModules(path)
    
        for mod in modules:
            if(hasattr(mod,'label')):
                if mod.label().find('Counter') !=-1 :
                    stringList.append(mod.name.value())
        print 'List Of Filters'
        print stringList

        return cms.untracked.vstring(stringList)

