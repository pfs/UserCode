#ifndef _object_tools_h_
#define _object_tools_h_

#include <vector>
#include <algorithm>

#include "TLorentzVector.h"

class Particle {
  
  public:
    Particle(TLorentzVector p4, int charge, int id, double puppi = 1)
      : p4_(p4), charge_(charge), id_(id), puppi_(puppi) {}
    
    double px()     { return p4_.Px();  }
    double py()     { return p4_.Py();  }
    double pz()     { return p4_.Pz();  }
    double e()      { return p4_.E();   }
    double pt()     { return p4_.Pt();  }
    double eta()    { return p4_.Eta(); }
    double phi()    { return p4_.Phi(); }
    double energy() { return p4_.E(); }
    double m()      { return p4_.M(); }
    double mass()   { return p4_.M(); }
    TLorentzVector p4()       { return p4_; }
    TLorentzVector momentum() { return p4_; }
    int charge()    { return charge_; }
    int id()        { return id_; }
    double puppi()  { return puppi_; }
  
  private:
    TLorentzVector p4_;
    int charge_;
    int id_;
    double puppi_;
};

/**
   @short summarizes the information on a jet needed for the charmed meson analysis
 */
typedef std::pair<TLorentzVector,int> IdTrack;

class Jet {

  public:
    Jet(TLorentzVector p4, int flavor, int idx)
      : p4_(p4), flavor_(flavor), idx_(idx), overlap_(0) {}
    Jet(TLorentzVector p4, float csv, int idx)
      : p4_(p4), csv_(csv), idx_(idx) {}
    ~Jet() {}
    
    double pt()     { return p4_.Pt();  }
    TLorentzVector p4()       { return p4_; }
    TLorentzVector momentum() { return p4_; }
    std::vector<Particle> particles() { return particles_; }
    int flavor()  { return flavor_; }
    int overlap() { return overlap_; }
    
    void addParticle(Particle p) { particles_.push_back(p); }
    void setFlavor(int flavor)   { flavor_ = flavor; }
    void setOverlap(int overlap) { overlap_ = overlap; }
    
    void addTrack(TLorentzVector p4, int pfid) { trks_.push_back( IdTrack(p4,pfid) ); }
    TLorentzVector &getVec() { return p4_; }
    float &getCSV() { return csv_; }
    int &getJetIndex() { return idx_; }
    std::vector<IdTrack> &getTracks() { return trks_; }
    void sortTracksByPt() { sort(trks_.begin(),trks_.end(), sortIdTracksByPt); }
    
    static bool sortJetsByPt(Jet i, Jet j)  { return i.getVec().Pt() > j.getVec().Pt(); }
    static bool sortJetsByCSV(Jet i, Jet j) { return i.getCSV() > j.getCSV(); }
  
  private:
    static bool sortIdTracksByPt(IdTrack i, IdTrack j)  { return i.first.Pt() > j.first.Pt(); }
    
    TLorentzVector p4_;
    std::vector<Particle> particles_;
    std::vector<IdTrack> trks_;
    float csv_;
    int flavor_;
    int idx_;
    int overlap_;
};
#endif
