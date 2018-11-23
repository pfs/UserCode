#ifndef _TemplatedFitTools_h_
#define _TemplatedFitTools_h_

#include <vector>

#include "TH1F.h"
#include "TObjArray.h"
#include "TString.h"

struct TemplatedFitResult_t
{
  float nExp,nExpUnc,nObs,nObsUnc,sf,sfUnc;
  int minuitStatus;
};

class TemplatedFitTools
{
 public:
  TemplatedFitTools();
  TemplatedFitResult_t fit(TObjArray &fracTempl,TH1F *data,Int_t idxOfInterest=0,TString saveResultIn="");
  ~TemplatedFitTools();
};

#endif
