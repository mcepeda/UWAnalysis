import FWCore.ParameterSet.Config as cms

#Import tool that creates the cut sequence
from UWAnalysis.Configuration.tools.CutSequenceProducer import *

analysisConfigurator = CutSequenceProducer(initialCounter = 'initialEvents',
                                  pyModuleName = __name__,
                                  pyNameSpace = locals())

analysisConfigurator.addSmearedAll('allMuons','goodMuons','vetoMuons','qcdMuons','allElectrons','goodElectrons','vetoElectrons','qcdElectrons','allJets','cleanJets','goodJets','fwdJets',mpost='')

#analysisConfigurator.addSmearing('patOverloadedTaus','goodMuons','goodElectrons','cleanPatJets','metCorrected',mpost="good")
#analysisConfigurator.addSmearing('patOverloadedTaus','allMuons','allElectrons','cleanPatJets','metCorrected',mpost="all")
#analysisConfigurator.addSmearing('patOverloadedTaus','vetoMuons','vetoElectrons','cleanPatJets','metCorrected',mpost="veto")
#analysisConfigurator.addSmearing('patOverloadedTaus','qcdMuons','qcdElectrons','cleanPatJets','metCorrected',mpost="qcd")
#analysisConfigurator.addSmearing('patOverloadedTaus','selectedPatMuons','selectedPatElectrons','cleanPatJets','metCorrected')

#DiMuonPreSel='leg1.isGlobalMuon && leg1.isTrackerMuon && leg2.isGlobalMuon && leg2.isTrackerMuon && leg1.pt()>20 && leg2.pt()>20 && leg1.charge * leg2.charge<0'
#DiMuonPreSel2='(leg1.isolationR03().sumPt+leg1.isolationR03().emEt+leg1.isolationR03().hadEt)/leg1.pt()<0.15 '
#DiMuonPreSel3='(leg2.isolationR03().sumPt+leg2.isolationR03().emEt+leg2.isolationR03().hadEt)/leg2.pt()<0.15 '
#
##Make Di Muons to VETO DY
##analysisConfigurator.addDiCandidateModule('diMuons','PATMuPairProducer', 'smearedMuons','smearedMuons','smearedMET','smearedJets',0,9999,text = '',leadingObjectsOnly = True,dR = 0.15,recoMode = "")
#analysisConfigurator.addDiCandidateModule('diMuons','PATMuPairProducer', 'smearedGoodMuons','smearedMuonsall','smearedMETall','smearedJetsall',0,9999,text = '',leadingObjectsOnly = True,dR = 0.15,recoMode = "")
##analysisConfigurator.addDiCandidateModule('diMuons','PATMuPairProducer', 'smearedMuonsgood','smearedMuonsgood','smearedMETgood','smearedJetsgood',0,9999,text = '',leadingObjectsOnly = True,dR = 0.15,recoMode = "")
#analysisConfigurator.addSelector('diMuonsSel','PATMuPairSelector',DiMuonPreSel+'&&'+DiMuonPreSel2+'&&'+DiMuonPreSel3,'diMuonSel',0,999)
#analysisConfigurator.addSorter('diMuonsSorted','PATMuPairSorter')
#
##Make Muons+MET cleanPatJets
##analysisConfigurator.addCandidateMETModule('MuNuSel','PATMuonNuPairProducer','smearedMuonsgood','smearedMETgood','smearedJetsgood',0,9999,'MuNuSel',genParticles="genDaughters")
##analysisConfigurator.addCandidateMETModule('MuNuVeto','PATMuonNuPairProducer','smearedMuonsveto','smearedMETveto','smearedJetsveto',0,9999,'MuNuVeto',genParticles="genDaughters")
#analysisConfigurator.addCandidateMETModule('MuNuAll','PATMuonNuPairProducer','smearedMuonsall','smearedMETall','smearedJetsall',0,9999,'MuNuAll',genParticles="genDaughters")
#
##analysisConfigurator.addSelector('wCandsKIN','PATMuonNuPairSelector','lepton.pt()>25 && abs(lepton.eta())<2.1' ,'wCandsKIN',0)
##analysisConfigurator.addSelector('wCandsJets','PATMuonNuPairSelector','','wCandsSel',0)
#
#
##elenuKIN='lepton.pt()>25 && abs(lepton.eta())<2.5&&(abs(lepton.eta)<1.4442||abs(lepton.eta)>1.5666)'
##elenuSel=''
##
##analysisConfigurator.addCandidateMETModule('wCandsEle','PATElectronNuPairProducer','smearedElectronsgood','smearedMETgood','smearedJetsgood',0,9999,'AtLeastOneWCandidateELE',genParticles="genDaughters")
##analysisConfigurator.addSelector('wCandsKINEle','PATElectronNuPairSelector',elenuKIN,'wCandsKINEle',0)
##analysisConfigurator.addSelector('wCandsSelEle','PATElectronNuPairSelector',elenuSel,'wCandsSelEle',0)
#
#
##Make Electrons+MET
##analysisConfigurator.addCandidateMETModule('EleNuSel','PATElectronNuPairProducer','smearedMuonsgood','smearedMETgood','smearedJetsgood',0,9999,'EleNuSel',genParticles="genDaughters")
##analysisConfigurator.addCandidateMETModule('EleNuVeto','PATElectronNuPairProducer','smearedMuonsveto','smearedMETveto','smearedJetsveto',0,9999,'EleNuVeto',genParticles="genDaughters")
#analysisConfigurator.addCandidateMETModule('EleNuAll','PATElectronNuPairProducer','smearedMuonsall','smearedMETall','smearedJetsall',0,9999,'EleNuAll',genParticles="genDaughters")
#

selectionSequence = analysisConfigurator.returnSequence()
