#ifndef _myiptools_h_
#define _myiptools_h_

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//get impact parameter for a track
template<class T>
std::pair<bool,Measurement1D> getImpactParameter(const T &trkRef, reco::VertexRef &vtx, const edm::EventSetup &iSetup, bool is3d=true)
{
  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
  reco::TransientTrack tt = trackBuilder->build(trkRef.get());
  if(is3d) return IPTools::absoluteImpactParameter3D(tt, *vtx);
  else     return IPTools::absoluteTransverseImpactParameter(tt, *vtx);
}

#endif
