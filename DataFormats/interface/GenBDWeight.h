#ifndef GenBDWeight_H
#define GenBDWeight_H
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>
#include "DataFormats/Math/interface/Vector3D.h"



class GenBDWeight{
public:

  GenBDWeight(){}

  //set info
  void theWeight(float val)                 {theWeight_=val;}

   float theWeight()	     const {return theWeight_;}	

private:
   float theWeight_;

};

typedef std::vector<GenBDWeight> GenBDWeightCollection;
typedef edm::Ref<GenBDWeightCollection> GenBDWeightRef;
typedef edm::RefProd<GenBDWeightCollection> GenBDWeightRefProd;
typedef edm::RefVector<GenBDWeightCollection> GenBDWeightRefVector;

#endif
