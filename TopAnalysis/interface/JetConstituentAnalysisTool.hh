#ifndef _JetConstituentAnalysisTool_hh_
#define _JetConstituentAnalysisTool_hh_
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <map>

class MiniEvent_t;

using namespace std;
class JetConstituentAnalysisTool
{
public:
  
  const MiniEvent_t        * event_ptr_;
  const TLorentzVector     * jet_ptr_;
  unsigned char index_;
  map<TString, TH1*> * plots_ptr_;
  map<TString, TH2*> * plots2D_ptr_;

  float weight_;
  void AssignHistograms()  const;
  void AnalyseAllJets()    const;
  void AnalyseLightJets()  const;
  void AnalyseBJets()      const;
  void Do();

};

#endif
