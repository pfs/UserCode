#ifndef _VBFVectorBoson_h_
#define _VBFVectorBoson_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TRandom3.h"
#include "TMath.h"
using namespace std;

//Vector boson will be either Z or photon at the moment

struct Category{
  float MM,A,VBF,HighPt,HighPtVBF,V1J,HighPtOfflineVBF,HighMJJ,LowMJJ,HighMJJLP,LowMJJLP;
  Category(){ reset(); }
  Category(std::vector<bool> &cat){
    reset();
    set(cat);
  };  
  void reset(){
    std::vector<bool> cat(9,false);
    set(cat);
  };
  void set(std::vector<bool> &cat){
    MM               =(float) cat[0];
    A                = (float)cat[1];
    VBF              = (float)cat[2];
    HighPt           = (float)cat[3];
    HighPtVBF        = (float)cat[4];
    V1J              = (float)cat[5];
    HighPtOfflineVBF = (float)cat[6];
    HighMJJ        = (float)cat[7];
    LowMJJ         = (float)cat[8];
    HighMJJLP       = (float)cat[9];
    LowMJJLP        = (float)cat[10];
  };
  std::vector<TString> getChannelTags() {
    std::vector<TString> chTags;
    TString chTag("");

    if(MM>0) chTag="MM";
    if(A>0)  chTag="A";
    if(chTag=="") return chTags;

    if(VBF>0)               chTags.push_back("VBF"+chTag);
    if(HighPt>0)            chTags.push_back("HighPt"+chTag);
    if(HighPtVBF>0)         chTags.push_back("HighPtVBF"+chTag);
    if(V1J>0)               chTags.push_back("V1J"+chTag);
    if(HighPtOfflineVBF>0)  chTags.push_back("HighPtOfflineVBF"+chTag);
    if(HighMJJ>0)           chTags.push_back("HighMJJ"+chTag);
    if(LowMJJ>0)            chTags.push_back("LowMJJ"+chTag);
    if(HighMJJLP>0)         chTags.push_back("HighMJJLP"+chTag);
    if(LowMJJLP>0)          chTags.push_back("LowMJJLP"+chTag);
    return chTags;
  }
};

class VBFVectorBoson{
public:
	VBFVectorBoson(TString filename_,
                       TString outname_,
                       Int_t anFlag_,
                       TH1F *normH_, 
                       TH1F *genPU_,
                       TString era_,
                       Bool_t debug_=false, Bool_t CR_ =false, Bool_t QCDTemp_ =true, Bool_t SRfake_ = false, Bool_t skimtree_=false, bool blind =true):
  filename(filename_),outname(outname_),anFlag(anFlag_), normH(0), genPU(0), era(era_), debug(debug_), CR(CR_), QCDTemp(QCDTemp_), SRfake(SRfake_), skimtree(skimtree_), doBlindAnalysis(blind)
	{
          if(normH_) normH = (TH1F*)normH_->Clone("normH_c");
	  if(genPU_) genPU = (TH1F*)genPU_->Clone("genPu_c");
	  fMVATree = NULL;
	  newTree = NULL;
	  init();
	  setXsecs();
          rnd.SetSeed(123456789);
	};
	~VBFVectorBoson(){}
	void init(){
		this->readTree();
		cout << "...producing " << outname << " from " << nentries << " events" << endl;
		this->prepareOutput();
		this->bookHistograms();
		this->loadCorrections();
		if(skimtree){
			this->addMVAvars();
		}
		selector = new SelectionTool(filename, debug, triggerList,SelectionTool::VBF);
		std::cout << "init done" << std::endl;
	}

	void saveHistos();
	void readTree();
	void prepareOutput();
	void bookHistograms();
	void setGammaZPtWeights();
	void loadCorrections();
	void addMVAvars();
	void fill(MiniEvent_t ev, TLorentzVector boson, std::vector<Jet> jets, std::vector<double> cplotwgts, TString c, std::map<TString, int> mults, std::vector<Particle> fakePhotonCR ={}, std::vector<Particle> tightPhotonCR={});
	void RunVBFVectorBoson();
	void initVariables(std::vector<Jet>);
	

private:
	TString filename, outname, baseName;
	Int_t anFlag, nentries;
	TH1F * normH, * genPU;
	TString era;
	TH1* triggerList;
	Bool_t debug, CR, QCDTemp, SRfake, skimtree, isQCDEMEnriched;
	TFile * f /*inFile*/, *fMVATree, *fOut;
	TTree * t /*inTree*/, *newTree /*MVA*/;
  	HistTool * ht;
	MiniEvent_t ev;
        float sihih,chiso,r9,hoe, ystar,relbpt,dphibjj;
	double mindrl;

	TRandom3 rnd;
	
  	//LUMINOSITY+PILEUP
	LumiTools * lumi;
  
  	//LEPTON EFFICIENCIES
	EfficiencyScaleFactorsWrapper * gammaEffWR;
  
        //L1 prefire efficiencies
        L1PrefireEfficiencyWrapper *l1PrefireWR;

  	//JEC/JER
  	//JECTools * jec;
	
        //Photon/Z pt weights
  	std::map<TString,TGraph *> photonPtWgts;
  	std::map<TString,std::pair<double,double> > photonPtWgtCtr;

	//Fake rate
	FakeRateTool * fr;

	//Variables to be added to the MVA Tree
	float centraleta, forwardeta, jjetas, centjy, ncentjj, dphivj0, dphivj1, dphivj2, dphivj3;
	float evtWeight, mjj, detajj , dphijj ,jjpt;
	float isotropy, circularity,sphericity,	aplanarity, C, D;
	float scalarht,balance, mht, training, lead_qg;
        float leadj_c1_05,subleadj_c1_05;
        float leadj_gawidth,subleadj_gawidth, subleadj_c2_02, subleadj_pt;
        float vbfmva,vbffisher;
        bool doBlindAnalysis;
        
	std::vector<Particle> relaxedTightPhotons, photons, tmpPhotons; 
	/////////////////////////////////////
	// Categorie for VBF:              //
	//   MM:A:VBF:HighPt:HighPtVBF:V1J // 
	/////////////////////////////////////
	Category category;

	SelectionTool * selector;

	/////////////////////////////////////////
	// To select events for training       //
	/////////////////////////////////////////

	int useForTraining(){
	  double myRnd = rnd.Rndm();
	  if (myRnd >= 0.5) return 1;
	  return 0;
	}

	/////////////////////////////////////////
	// Quick and ugly cross section getter //
	// To be updated such as to use the    //
	// json file                           //
	/////////////////////////////////////////
	std::map<TString, float> xsecRefs;

	void setXsecs(){
	  xsecRefs["MC13TeV_TTJets"        ] = 832;
	  xsecRefs["MC13TeV_ZZ"            ] = 0.5644; 
	  xsecRefs["MC13TeV_WZ"            ] = 47.13;  
	  xsecRefs["MC13TeV_WW"            ] = 12.178; 
	  xsecRefs["MC13TeV_SingleTbar_tW" ] = 35.85;  
	  xsecRefs["MC13TeV_SingleT_tW"    ] = 35.85;  
	  xsecRefs["MC13TeV_DY50toInf"     ] = 5765.4; 
	  xsecRefs["MC13TeV_QCDEM_15to20"  ] = 2302200;
	  xsecRefs["MC13TeV_QCDEM_20to30"  ] = 5352960;
	  xsecRefs["MC13TeV_QCDEM_30to50"  ] = 9928000;
	  xsecRefs["MC13TeV_QCDEM_50to80"  ] = 2890800;
	  xsecRefs["MC13TeV_QCDEM_80to120" ] = 350000; 
	  xsecRefs["MC13TeV_QCDEM_120to170"] = 62964;  
	  xsecRefs["MC13TeV_QCDEM_170to300"] = 18810;  
	  xsecRefs["MC13TeV_QCDEM_300toInf"] = 1350;  
	  xsecRefs["MC13TeV_GJets_HT40to100" ] = 20790 ;
	  xsecRefs["MC13TeV_GJets_HT100to200"] = 9238;  
	  xsecRefs["MC13TeV_GJets_HT200to400"] = 2305;  
	  xsecRefs["MC13TeV_GJets_HT400to600"] = 274.4; 
	  xsecRefs["MC13TeV_GJets_HT600toInf"] = 93.46; 
	  xsecRefs["MC13TeV_EWKZJJ"        ] = 4.32;
	  xsecRefs["MC13TeV_AJJ_EWK_LO"    ] = 32.49;
	  xsecRefs["MC13TeV_AJJ_EWK_INT_LO"] = 8.3;
	}
	
	float getXsec(){
	  for (auto const& x : xsecRefs){
	      if (filename.Contains(x.first))
		return x.second;
	  }
	  return 1;
	}
	
};
#endif
