/*#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ReadTree.h"
#include "TopLJets2015/TopAnalysis/interface/TemplatedFitTools.h"
#include "TopLJets2015/TopAnalysis/interface/myfunction.h"*/

#ifdef __CINT__

#pragma link off all class; 
#pragma link off all function; 
#pragma link off all global; 
#pragma link off all struct; 
#pragma link off all typedef;
#pragma link off all enum;

#pragma link C++ function ReadTree+;
#pragma link C++ function myfunction+;
#pragma link C++ function RunOverSamples;

//#pragma link C++ defined_in "src/TopAnalysisLinkDef.h";
#pragma link C++ function attachToMiniEventTree;
#pragma link C++ function createMiniEventTree;
#pragma link C++ enum FlavourSplitting;
#pragma link C++ class TemplatedFitTools;
#pragma link C++ struct TemplatedFitResult_t;

#endif

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 8
// End:

