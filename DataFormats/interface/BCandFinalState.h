#ifndef BCandFinalState_H
#define BCandFinalState_H
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"



class BCandFinalState{
public:

  BCandFinalState(){}

  //set info
  void nBC(int val)                 {nBC_=val;}
  void BC1PT(float val)             {BC1PT_=val;}
  void BC2PT(float val)             {BC2PT_=val;}
  void BC1ETA(float val)            {BC1ETA_=val;}
  void BC2ETA(float val)            {BC2ETA_=val;}
  void BC1PHI(float val)            {BC1PHI_=val;}
  void BC2PHI(float val)            {BC2PHI_=val;}
  void BC1nTracks(int val)          {BC1nTracks_=val;}
  void BC2nTracks(int val)          {BC2nTracks_=val;}
  void BCEnergySum(float val)       {BCEnergySum_=val;}
  void BCGenEnergySum(float val)       {BCGenEnergySum_=val;}
  void nJetAssociatedToBC1(int val) {nJetAssociatedToBC1_=val;}
  void nJetAssociatedToBC2(int val) {nJetAssociatedToBC2_=val;}
  void nJetAssociatedToBC(int val) {nJetAssociatedToBC_=val;}
  void Jet1BC1PT(float val)         {Jet1BC1PT_=val;}
  void Jet2BC1PT(float val)         {Jet2BC1PT_=val;}
  void Jet3BC1PT(float val)         {Jet3BC1PT_=val;}
  void Jet4BC1PT(float val)         {Jet4BC1PT_=val;}
  void Jet1BC1ETA(float val)        {Jet1BC1ETA_=val;}
  void Jet2BC1ETA(float val)        {Jet2BC1ETA_=val;}
  void Jet3BC1ETA(float val)        {Jet3BC1ETA_=val;}
  void Jet4BC1ETA(float val)        {Jet4BC1ETA_=val;}
  void Jet1BC1PHI(float val)        {Jet1BC1PHI_=val;}
  void Jet2BC1PHI(float val)        {Jet1BC1PHI_=val;}
  void Jet3BC1PHI(float val)        {Jet1BC1PHI_=val;}	
  void Jet4BC1PHI(float val)        {Jet1BC1PHI_=val;}
  void Jet1BC2PT(float val)         {Jet1BC2PT_=val;}
  void Jet2BC2PT(float val)         {Jet2BC2PT_=val;}
  void Jet3BC2PT(float val)         {Jet3BC2PT_=val;}
  void Jet4BC2PT(float val)         {Jet4BC2PT_=val;}
  void Jet1BC2ETA(float val)        {Jet1BC2ETA_=val;}
  void Jet2BC2ETA(float val)        {Jet2BC2ETA_=val;}
  void Jet3BC2ETA(float val)        {Jet3BC2ETA_=val;}
  void Jet4BC2ETA(float val)        {Jet4BC2ETA_=val;}
  void Jet1BC2PHI(float val)        {Jet1BC2PHI_=val;}
  void Jet2BC2PHI(float val)        {Jet1BC2PHI_=val;}
  void Jet3BC2PHI(float val)        {Jet1BC2PHI_=val;}
  void Jet4BC2PHI(float val)        {Jet1BC2PHI_=val;}
  void BCDeltaR(float val)          {BCDeltaR_=val;}
  void BDeltaPHI(float val)         {BDeltaPHI_=val;}
  void BC1MASS(float val)          {BC1MASS_=val;}
  void BC2MASS(float val)         {BC2MASS_=val;}


  //get info
   int nBC()		     const {return nBC_;}	
   float BC1PT()             const {return BC1PT_;}
   float BC2PT()             const {return BC2PT_;}
   float BC1ETA()            const {return BC1ETA_;}
   float BC2ETA()            const {return BC2ETA_;}
   float BC1PHI()            const {return BC1PHI_;}
   float BC2PHI()            const {return BC2PHI_;}
   int BC1nTracks()          const {return BC1nTracks_;}
   int BC2nTracks()          const {return BC2nTracks_;}
   float BCEnergySum()       const {return BCEnergySum_;}
   float BCGenEnergySum()       const {return BCGenEnergySum_;}
   int nJetAssociatedToBC1() const {return nJetAssociatedToBC1_;}
   int nJetAssociatedToBC2() const {return nJetAssociatedToBC2_;}
   int nJetAssociatedToBC() const {return nJetAssociatedToBC_;}
   float Jet1BC1PT()         const {return Jet1BC1PT_;}
   float Jet2BC1PT()         const {return Jet2BC1PT_;}
   float Jet3BC1PT()         const {return Jet3BC1PT_;}
   float Jet4BC1PT()         const {return Jet4BC1PT_;}
   float Jet1BC1ETA()        const {return Jet1BC1ETA_;}
   float Jet2BC1ETA()        const {return Jet2BC1ETA_;}
   float Jet3BC1ETA()        const {return Jet3BC1ETA_;}
   float Jet4BC1ETA()        const {return Jet4BC1ETA_;}
   float Jet1BC1PHI()        const {return Jet1BC1PHI_;}
   float Jet2BC1PHI()        const {return Jet2BC1PHI_;}
   float Jet3BC1PHI()        const {return Jet3BC1PHI_;}
   float Jet4BC1PHI()        const {return Jet4BC1PHI_;}
   float Jet1BC2PT()         const {return Jet1BC2PT_;}
   float Jet2BC2PT()         const {return Jet2BC2PT_;}
   float Jet3BC2PT()         const {return Jet3BC2PT_;}
   float Jet4BC2PT()         const {return Jet4BC2PT_;}
   float Jet1BC2ETA()        const {return Jet1BC2ETA_;}
   float Jet2BC2ETA()        const {return Jet2BC2ETA_;}
   float Jet3BC2ETA()        const {return Jet3BC2ETA_;}
   float Jet4BC2ETA()        const {return Jet4BC2ETA_;}
   float Jet1BC2PHI()        const {return Jet1BC2PHI_;}
   float Jet2BC2PHI()        const {return Jet2BC2PHI_;}
   float Jet3BC2PHI()        const {return Jet3BC2PHI_;}
   float Jet4BC2PHI()        const {return Jet4BC2PHI_;}
   float BCDeltaR()          const {return BCDeltaR_;}
   float BDeltaPHI()         const {return BDeltaPHI_;}
   float BC1MASS()         const {return BC1MASS_;}
   float BC2MASS()         const {return BC2MASS_;}

private:
   int nBC_;
   float BC1PT_;
   float BC2PT_;
   float BC1ETA_;
   float BC2ETA_;
   float BC1PHI_;
   float BC2PHI_;
   int BC1nTracks_;
   int BC2nTracks_;	
   float BCEnergySum_;
   float BCGenEnergySum_;
   int nJetAssociatedToBC1_;
   int nJetAssociatedToBC2_;
   int nJetAssociatedToBC_;
   float Jet1BC1PT_;
   float Jet2BC1PT_;
   float Jet3BC1PT_;
   float Jet4BC1PT_;
   float Jet1BC1ETA_;
   float Jet2BC1ETA_;
   float Jet3BC1ETA_;
   float Jet4BC1ETA_;
   float Jet1BC1PHI_;
   float Jet2BC1PHI_;
   float Jet3BC1PHI_;
   float Jet4BC1PHI_;
   float Jet1BC2PT_;
   float Jet2BC2PT_;
   float Jet3BC2PT_;
   float Jet4BC2PT_;
   float Jet1BC2ETA_;
   float Jet2BC2ETA_;
   float Jet3BC2ETA_;
   float Jet4BC2ETA_;
   float Jet1BC2PHI_;
   float Jet2BC2PHI_;
   float Jet3BC2PHI_;
   float Jet4BC2PHI_;
   float BCDeltaR_;
   float BDeltaPHI_;
   float BC1MASS_;
   float BC2MASS_;


};

typedef std::vector<BCandFinalState> BCandFinalStateCollection;
typedef edm::Ref<BCandFinalStateCollection> BCandFinalStateRef;
typedef edm::RefProd<BCandFinalStateCollection> BCandFinalStateRefProd;
typedef edm::RefVector<BCandFinalStateCollection> BCandFinalStateRefVector;

#endif
