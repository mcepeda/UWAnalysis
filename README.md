UWAnalysis
==========
UW Analysis for 53x (8Tev)


download this using the recipe below:

```bash

scram pro -n Wbb_slc6_472_CMSSW_5_3_14_patch1 CMSSW CMSSW_5_3_14_patch1
cd Wbb_slc6_472_CMSSW_5_3_14_patch1/src
cmsenv

eval `ssh-agent -s`
ssh-add

git cms-addpkg DataFormats/PatCandidates     
git cms-addpkg HiggsAnalysis/CombinedLimit

git cms-addpkg PhysicsTools/PatAlgos
git cms-addpkg DataFormats/StdDictionaries
git cms-addpkg PhysicsTools/PatUtils 
git cms-addpkg CommonTools/ParticleFlow
git cms-addpkg FWCore/GuiBrowsers
git cms-addpkg RecoParticleFlow/PFProducer
git cms-addpkg CommonTools/RecoUtils
git cms-addpkg DataFormats/HLTReco
git cms-addpkg JetMETCorrections/Type1MET

git cms-addpkg RecoJets/JetProducers
git cms-cvs-history import METPU_5_3_X_v4 RecoJets/JetProducers 

git cms-addpkg PhysicsTools/Utilities
git cms-addpkg RecoBTag/SecondaryVertex
git cms-addpkg RecoVertex/AdaptiveVertexFinder
git cms-addpkg RecoEgamma/EgammaTools

git cms-addpkg EgammaAnalysis/ElectronTools
pushd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
popd

git cms-addpkg DataFormats/METReco  
git cms-addpkg RecoTauTag/RecoTau
git cms-addpkg RecoTauTag/Configuration
git cms-addpkg CondFormats/EgammaObjects

git clone https://github.com/cms-analysis/AnalysisDataFormats-TauAnalysis.git AnalysisDataFormats/TauAnalysis
pushd $CMSSW_BASE/src/AnalysisDataFormats/TauAnalysis
git checkout AnalysisDataFormats-TauAnalysis-bMinimalSVfit-08-03-11
popd

git clone https://github.com/cms-analysis/TauAnalysis-CandidateTools.git TauAnalysis/CandidateTools
pushd $CMSSW_BASE/src/TauAnalysis/CandidateTools
git checkout TauAnalysis-CandidateTools-bMinimalSVfit_2012May13
popd

git cms-addpkg RecoLuminosity/LumiDB
git cms-addpkg PhysicsTools/PythonAnalysis
git cms-addpkg FWCore/PythonUtilities
git-cms-cvs-history import V05-00-16 DataFormats/JetReco
 
git cms-addpkg DataFormats/MuonReco

git clone https://github.com/sdevissc/ZSV
git clone https://github.com/tmrhombus/UWAnalysis.git

```
