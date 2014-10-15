import FWCore.ParameterSet.Config as cms

def addMCEventTree(process,name,src="genParticles",srcMuons="dressedLeptons",srcTaggedGenJets="cleanGenJets"):
   process.TFileService = cms.Service("TFileService", fileName = cms.string("analysis.root") )
   eventTree = cms.EDAnalyzer('EventTreeMaker',
                              coreCollections = cms.VInputTag(
                                   cms.InputTag(src)
                              ),
			     # gen Level identification
	                genID = cms.PSet(
                	        pluginType = cms.string("GenID"),
                        	src        = cms.InputTag("genParticles"),
        	        ),

                genID2 = cms.PSet(
                        pluginType = cms.string("GenID2"),
       		        verbose    = cms.untracked.bool(False),
                        saveCs     = cms.untracked.bool(True),
                        saveBs     = cms.untracked.bool(True),
       		        ptMin      = cms.untracked.double(0),
                        etaMax     = cms.untracked.double(10),
                         src        = cms.InputTag("genParticles"),

                ),

		 	     simBHadrons = cms.PSet(
                                    pluginType = cms.string("SimBHadronsFiller"),
                                    src        = cms.InputTag("bhadrons"),
                                    tag        = cms.string("bHadronsID"),
				    method     = cms.string("pdgId()"),
				    leadingOnly=cms.untracked.bool(False)	
                              ),

                             simBHadronsPT = cms.PSet(
                                    pluginType = cms.string("SimBHadronsFiller"),
                                    src        = cms.InputTag("bhadrons"),
                                    tag        = cms.string("bHadronsPT"),
                                    method     = cms.string("pt()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             simBHadronsETA = cms.PSet(
                                    pluginType = cms.string("SimBHadronsFiller"),
                                    src        = cms.InputTag("bhadrons"),
                                    tag        = cms.string("bHadronsETA"),
                                    method     = cms.string("eta()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             simBHadronsPHI = cms.PSet(
                                    pluginType = cms.string("SimBHadronsFiller"),
                                    src        = cms.InputTag("bhadrons"),
                                    tag        = cms.string("bHadronsPHI"),
                                    method     = cms.string("phi()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),



                             simBHadronsSIZE = cms.PSet(
                                    pluginType = cms.string("CollectionSizeFiller"),
                                    src        = cms.InputTag("bhadrons"),
                                    tag        = cms.string("bHadronsSIZE"),
                              ),

                           jetPT = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPt"),
                                   method     = cms.string("pt()"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),

                           jetPHI = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPhi"),
                                   method     = cms.string("phi()"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),

                           jetETA = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetEta"),
                                   method     = cms.string("eta()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),
                           jetPTDRHadron = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPtDRHadron"),
                                   method     = cms.string("userFloat('minDRBHadron')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),

                           jetPTLeadB = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPtLeadB"),
                                   method     = cms.string("userFloat('leadJetBHadronPt')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),
                           jetPTSecondB = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPtSecondB"),
                                   method     = cms.string("userFloat('secondJetBHadronPt')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),
                           jetETALeadB = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetEtaLeadB"),
                                   method     = cms.string("userFloat('leadJetBHadronEta')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),
                           jetETASecondB = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetEtaSecondB"),
                                   method     = cms.string("userFloat('secondJetBHadronEta')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),

                           jetPhiLeadB = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPhiLeadB"),
                                   method     = cms.string("userFloat('leadJetBHadronPhi')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),
                           jetPhiSecondB = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetPhiSecondB"),
                                   method     = cms.string("userFloat('secondJetBHadronPhi')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),


                           jetAngleBs = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetAngleBs"),
                                   method     = cms.string("userFloat('DR2BHadrons')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),

                           jetMatched = cms.PSet(
                                   pluginType = cms.string("PATJetFiller"),
                                   src        = cms.InputTag(srcTaggedGenJets),
                                   tag        = cms.string("genJetBsMatches"),
                                   method     = cms.string("userFloat('matchedHadron')"),
                                   leadingOnly=cms.untracked.bool(False)
                             ),

                             dressedLeptons = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonPDGID"),
                                    method     = cms.string("pdgId()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),


                             dressedLeptonsType = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonType"),
                                    method     = cms.string("Type()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             dressedLeptonsMOM = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonMotherPDGID"),
                                    method     = cms.string("MotherPDGID()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             dressedLeptonsOriginalPt = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonOriginalPT"),
                                    method     = cms.string("originalPt()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             dressedLeptonsPT = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonPT"),
                                    method     = cms.string("pt()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             dressedLeptonsETA = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonETA"),
                                    method     = cms.string("eta()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             dressedLeptonsPHI = cms.PSet(
                                    pluginType = cms.string("DressedLeptonFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonPHI"),
                                    method     = cms.string("phi()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             NeutrinosPDGID = cms.PSet(
                                    pluginType = cms.string("PATGenParticleFiller"),
                                    src        = cms.InputTag("neutrinos"),
                                    tag        = cms.string("neutrinoPDGID"),
                                    method     = cms.string("pdgId()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ), 

                             NeutrinosPT = cms.PSet(
                                    pluginType = cms.string("PATGenParticleFiller"),
                                    src        = cms.InputTag("neutrinos"),
                                    tag        = cms.string("neutrinoPT"),
                                    method     = cms.string("pt()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             NeutrinosETA = cms.PSet(
                                    pluginType = cms.string("PATGenParticleFiller"),
                                    src        = cms.InputTag("neutrinos"),
                                    tag        = cms.string("neutrinoETA"),
                                    method     = cms.string("eta()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),

                             NeutrinosPHI = cms.PSet(
                                    pluginType = cms.string("PATGenParticleFiller"),
                                    src        = cms.InputTag("neutrinos"),
                                    tag        = cms.string("neutrinoPHI"),
                                    method     = cms.string("phi()"),
                                    leadingOnly=cms.untracked.bool(False)
                              ),



                             dressedLeptonsSIZE = cms.PSet(
                                    pluginType = cms.string("CollectionSizeFiller"),
                                    src        = cms.InputTag(srcMuons),
                                    tag        = cms.string("dressedLeptonSIZE"),
                              ),



   )

   setattr(process, name, eventTree)
   p = cms.Path(getattr(process,name))
   setattr(process, name+'Path', p)



