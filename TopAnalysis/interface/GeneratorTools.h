#ifndef _GENERATORTOOLS_H
#define _GENERATORTOOLS_H

#include <vector>
#include <string>
#include "TFile.h"
#include "TString.h"

typedef std::pair<TString,int> WeightSysts_t;
std::vector< WeightSysts_t > getWeightSysts(TFile *);

#endif
