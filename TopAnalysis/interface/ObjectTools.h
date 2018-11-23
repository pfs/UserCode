#ifndef _object_tools_h_
#define _object_tools_h_

#include <vector>
#include <algorithm>

#include "TLorentzVector.h"

class Particle : public TLorentzVector {
  
  public:
    Particle(TLorentzVector p4, int charge, int id, int qualityFlags, int origRef, double puppi = 1,double unc=0)
      : TLorentzVector(p4), charge_(charge), id_(id), qualityFlags_(qualityFlags), origRef_(origRef), puppi_(puppi), unc_(unc) {}
   
    Particle( const Particle &p) 
      : TLorentzVector(p.px(),p.py(),p.pz(),p.e()), charge_(p.charge_), id_(p.id_), qualityFlags_(p.qualityFlags_), origRef_(p.origRef_), puppi_(p.puppi_), unc_(p.unc_) {}

    double px() const    { return TLorentzVector::Px();  }
    double py() const    { return TLorentzVector::Py();  }
    double pz() const    { return TLorentzVector::Pz();  }
    double e()  const    { return TLorentzVector::E();   }
    double pt()     { return TLorentzVector::Pt();  }
    double eta()    { return TLorentzVector::Eta(); }
    double phi()    { return TLorentzVector::Phi(); }
    double energy() { return TLorentzVector::E(); }
    double m()      { return TLorentzVector::M(); }
    double mass()   { return TLorentzVector::M(); }
    TLorentzVector p4()       { return TLorentzVector(TLorentzVector::Px(),TLorentzVector::Py(),TLorentzVector::Pz(),TLorentzVector::E()); }
    TLorentzVector momentum() { return TLorentzVector(TLorentzVector::Px(),TLorentzVector::Py(),TLorentzVector::Pz(),TLorentzVector::E()); }
    int charge()    { return charge_; }
    int id()        { return id_; }
    int qualityFlags() { return qualityFlags_; }
    bool hasQualityFlag(int bit) { return ((qualityFlags_>>bit)&0x1); }
    int originalReference() { return origRef_; }
    void setOriginalReference(int origRef) { origRef_=origRef; }
    double puppi()  { return puppi_; }
    double scaleUnc() { return unc_; }

  private:
    int charge_, id_, qualityFlags_,origRef_;
    double puppi_;
    double unc_;
};

/**
   @short summarizes the information on a jet needed for the charmed meson analysis
 */
typedef std::pair<TLorentzVector,int> IdTrack;

class Jet : public TLorentzVector {

  public:
    Jet(TLorentzVector p4, int flavor, int idx)
      : TLorentzVector(p4), flavor_(flavor), idx_(idx), overlap_(0),scaleUnc_(0) {}
    Jet(TLorentzVector p4, float csv, int idx)
      : TLorentzVector(p4), csv_(csv), idx_(idx),scaleUnc_(0) {}
    Jet(const Jet &j) 
      : TLorentzVector(j.Px(),j.Py(),j.Pz(),j.E()), particles_(j.particles_), trks_(j.trks_), csv_(j.csv_), pumva_(j.pumva_), flavor_(j.flavor_), idx_(j.idx_), overlap_(j.overlap_), scaleUnc_(j.scaleUnc_) {}
    ~Jet() {}

    double pt()     { return TLorentzVector::Pt();  }
    double eta()    { return TLorentzVector::Eta();  }
    TLorentzVector p4()       { return TLorentzVector(TLorentzVector::Px(),TLorentzVector::Py(),TLorentzVector::Pz(),TLorentzVector::E()); }
    TLorentzVector momentum() { return TLorentzVector(TLorentzVector::Px(),TLorentzVector::Py(),TLorentzVector::Pz(),TLorentzVector::E()); }
    std::vector<Particle> &particles() { return particles_; }
    int flavor()  { return flavor_; }
    int overlap() { return overlap_; }
    
    void addParticle(Particle p) { particles_.push_back(p); }
    void setFlavor(int flavor)   { flavor_ = flavor; }
    void setOverlap(int overlap) { overlap_ = overlap; }
    
    void addTrack(TLorentzVector p4, int pfid) { trks_.push_back( IdTrack(p4,pfid) ); }
    //TLorentzVector &getVec() { return p4_; }
    float &getCSV() { return csv_; }
    void setCSV(float csv) { csv_=csv; }
    float &getPUMVA() { return pumva_; }
    void setPUMVA(float pumva) { pumva_=pumva; }
    void setDeepCSV(float deepcsv) { deepcsv_=deepcsv; }
    int &getJetIndex() { return idx_; }
    std::vector<IdTrack> &getTracks() { return trks_; }
    void sortTracksByPt() { sort(trks_.begin(),trks_.end(), sortIdTracksByPt); }
    void setScaleUnc(float scaleUnc) { scaleUnc_=scaleUnc; }
    float getScaleUnc() { return scaleUnc_; }   
    static bool sortJetsByPt(Jet i, Jet j)  { return i.Pt() > j.Pt(); }
    static bool sortJetsByCSV(Jet i, Jet j) { return i.getCSV() > j.getCSV(); }
  
  private:
    static bool sortIdTracksByPt(IdTrack i, IdTrack j)  { return i.first.Pt() > j.first.Pt(); }
    
    //TLorentzVector p4_;
    std::vector<Particle> particles_;
    std::vector<IdTrack> trks_;
    float csv_,deepcsv_;
    float pumva_;    
    int flavor_;
    int idx_;
    int overlap_;
    float scaleUnc_;
};

#endif
