#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace std;

float LEADINGJETCUT=50.;
bool APPLYEECLEANING=true;
bool APPLYPHOTONTRIGSAFE=true;
int PUIDX=1;

//
void VBFVectorBoson::RunVBFVectorBoson()
{
  bool is2018(filename.Contains("2018"));
  bool isJetHT(filename.Contains("JetHT"));
  float highMJJcut(1000.);
  if (isJetHT) highMJJcut = 50.;
  float minMJJ(150);
  if(isJetHT) minMJJ = 50;
  float minBosonHighPt(200.);
  TString vbfPhotonTrigger = "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v";
  TString highPtPhotonTrigger("HLT_Photon200_v");
  SelectionTool::QualityFlags offlinePhoton(SelectionTool::TIGHT);
  if(is2018) {
    cout << "[VBFVectorBoson::RunVBFVectorBoson] this is 2018, adapting" << endl;
    vbfPhotonTrigger="HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v";
    highPtPhotonTrigger="HLT_Photon165_R9Id90_HE10_IsoM_v"; //HLT_Photon300_NoHE_v
    offlinePhoton=SelectionTool::LOOSE;
    minBosonHighPt=165.;
  }
  std::vector<TString> photonTriggerPaths={vbfPhotonTrigger,highPtPhotonTrigger};
  if(isJetHT){
    cout <<"[VBFVectorBoson::RunVBFVectorBoson] this is JetHT data, adapting the trigger" <<endl;
    photonTriggerPaths.clear();
    photonTriggerPaths = {"HLT_PFJet40_v","HLT_PFJet60_v","HLT_PFJet80_v","HLT_PFJet140_v",
                          "HLT_PFJet200_v","HLT_PFJet260_v","HLT_PFJet320_v","HLT_PFJet400_v",
			  "HLT_PFJet450_v","HLT_PFJet500_v","HLT_PFJet550_v","HLT_PFJetFwd40_v",
			  "HLT_PFJetFwd60_v","HLT_PFJetFwd80_v","HLT_PFJetFwd140_v","HLT_PFJetFwd200_v",
			  "HLT_PFJetFwd260_v","HLT_PFJetFwd320_v","HLT_PFJetFwd400_v","HLT_PFJetFwd450_v",
			  "HLT_PFJetFwd500_v"};

  }
  selector->setPhotonSelection(photonTriggerPaths,offlinePhoton);


  //TMVA configuration
  std::map<TString,TMVA::Reader *> readers;
  std::map<TString,TGraph *> mvaCDFinv;
  TString method("BDT_VBF0HighMJJ");
  TString weightFile("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/HighMJJ_BDT_VBF0HighMJJ.weights.xml");
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("D",             &vbfVars_.D);
  readers[method]->AddVariable("C",             &vbfVars_.C);
  readers[method]->AddVariable("circularity",   &vbfVars_.circularity);
  readers[method]->AddVariable("balance",       &vbfVars_.balance);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("j_pt[1]",       &vbfVars_.subleadj_pt);
  readers[method]->AddVariable("j_c2_00[0]",    &vbfVars_.leadj_c2_02);
  readers[method]->AddVariable("j_c2_00[1]",    &vbfVars_.subleadj_c2_02);
  readers[method]->AddVariable("j_gawidth[0]",  &vbfVars_.leadj_gawidth);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->BookMVA(method,weightFile);
  
  method="BDT_VBF0LowMJJ";
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/LowMJJ_BDT_VBF0LowMJJ.weights.xml";
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("D",             &vbfVars_.D);
  readers[method]->AddVariable("C",             &vbfVars_.C);
  readers[method]->AddVariable("balance",       &vbfVars_.balance);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("dphivj1",       &vbfVars_.dphivj1);
  readers[method]->AddVariable("j_pt[1]",       &vbfVars_.subleadj_pt);
  readers[method]->AddVariable("detajj",        &vbfVars_.detajj);
  readers[method]->AddVariable("j_c2_00[0]",    &vbfVars_.leadj_c2_02);
  readers[method]->AddVariable("j_c2_00[1]",    &vbfVars_.subleadj_c2_02);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("j_qg[0]",       &vbfVars_.leadj_qg);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->BookMVA(method,weightFile);

  //read the transformations based on CDF^{-1}
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/inverse_cdfs.root";
  gSystem->ExpandPathName(weightFile);
  TFile *fcdf=TFile::Open(weightFile);
  for(std::map<TString,TMVA::Reader *>::iterator it=readers.begin(); it!=readers.end(); it++) {
    TString key=it->first;
    mvaCDFinv[key]=(TGraph *)fcdf->Get(key+"_cdfinv");
  }


  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  float xsec = getXsec();
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(debug) cout << "Number of event: "<<iev<<endl;
      if(iev%10000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      std::vector<double>plotwgts(1,1.0);
      ht->fill("qscale",ev.g_qscale,plotwgts);
      
      //assign randomly a run period
      TString period = lumi->assignRunPeriod();
      
      //////////////////
      // CORRECTIONS //
      ////////////////      
      //jec->smearJetEnergies(ev);
             
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector->flagFinalState(ev,{},{},CR, QCDTemp, SRfake);
      photons                            = selector->getSelPhotons(); 
      relaxedTightPhotons                = selector->getRelaxedTightPhotons();
      tmpPhotons                         = selector->getQCDTmpPhotons();
      std::vector<Particle> &leptons     = selector->getSelLeptons(); 
      std::vector<Jet>       alljets     = selector->getGoodJets(ev,30.,4.7,leptons,photons);
      int nTotalJets = alljets.size();
      std::vector<Jet> jets;
      std::vector<Particle> fakeACR;
      std::vector<Particle> tightACR;
      std::map<TString, int> mults;
      mults["loosefake"]   = 0;
      mults["tightfake"]   = 0;
      mults["looseprompt"] = 0;
      mults["tightprompt"] = 0;
      std::vector<Particle> QCDTemplate;

      //Pileup jet id
      for(auto j : alljets) {
        int idx=j.getJetIndex();
        int jid=ev.j_id[idx];
        //cleanup for ecal noise
        if(APPLYEECLEANING && fabs(j.Eta())>2.7 && fabs(j.Eta())<3 && ev.j_emf[idx]>0.55) continue;
	//bool passPu((jid>>2)&0x1);
	bool passPu((jid>>PUIDX)&0x1);
	if(CR){
	  if(jets.size() == 0 && passPu)
	    continue;
	} else if(!passPu) continue;
      	jets.push_back(j);
      }
     //Fake and tight photons in CR
      int nLPGamma(0);
      for(auto a : photons) {
        int idx=a.originalReference();
	if (selector->isFakePhoton(ev,idx)){
	  fakeACR.push_back(a);
	  if (ev.gamma_isPromptFinalState[idx])
	    mults["looseprompt"]++;
	  else
	    mults["loosefake"]++;
	} else if (a.hasQualityFlag(SelectionTool::TIGHT)){
	  tightACR.push_back(a);
	  if (ev.gamma_isPromptFinalState[idx])
	    mults["tightprompt"]++;
	  else
	    mults["tightfake"]++;
	}
        if(a.Pt()>50 && fabs(a.Eta()) > 2.25 && fabs(a.Eta()) < 3.0)
          nLPGamma++;
      }

      //Category selection
      if(chTag!="A" && chTag!="MM") continue;

      category.reset();
      std::vector<bool> cat(8,0);
      if(chTag == "MM") cat[0] = true;
      if(chTag == "A") cat[1] = true;

      //L1-prefire jet candidates
      int nLPJets(0);
      for(auto j : jets) {
	if(j.Pt()>100 && fabs(j.Eta()) > 2.25 && fabs(j.Eta()) < 3.0) nLPJets++;
      }

      //define the boson
      TLorentzVector boson(0,0,0,0);
      float bosonScaleUnc(0.);
      sihih = 0, chiso = 0 ,r9 = 0, hoe = 0, eveto=0, pixelseed=0;
      if(chTag=="A") {
        boson += photons[0];
        bosonScaleUnc = photons[0].scaleUnc()/photons[0].Pt();
        int pidx=photons[0].originalReference();
        sihih         = ev.gamma_sieie[pidx];
        chiso         = ev.gamma_chargedHadronIso[pidx];
        r9            = ev.gamma_r9[pidx];
        hoe           = ev.gamma_hoe[pidx];
        eveto         = ev.gamma_idFlags[pidx] & 0x1;
        pixelseed     = (ev.gamma_idFlags[pidx]>>1) & 0x1;
      }else {
        boson   += leptons[0];
        boson   += leptons[1];
        bosonScaleUnc= TMath::Sqrt( pow(leptons[0].Pt()*leptons[1].scaleUnc(),2)+
                                    pow(leptons[1].Pt()*leptons[0].scaleUnc(),2) )/boson.Pt();
      }

      //leptons,jets and boson
      mindrl=9999.;
      for(auto &l: leptons) mindrl=min(l.DeltaR(boson),mindrl);
      mindrj=9999.;
      for(auto &j: jets) mindrj=min(j.DeltaR(boson),mindrj);
      
      //vbf-dedicated
      vbfVars_.fillDiscriminatorVariables(boson,jets,ev);

      //final categories
      bool passJetMult(jets.size()>=2 && jets[0].Pt()>LEADINGJETCUT);
      bool passMJJ(passJetMult && vbfVars_.mjj>highMJJcut);
      bool passJets(passJetMult && vbfVars_.mjj>minMJJ);
      bool passVBFJetsTrigger(passJets && vbfVars_.detajj>3.0);
      //L1-prefiltering check
      bool passLP(nLPJets > 0 || nLPGamma > 0);
      bool isHighPt(false),isVBF(false),isHighPtAndVBF(false),isHighPtAndOfflineVBF(false),
        isBosonPlusOneJet(false),isHighMJJ(false),isLowMJJ(false), isHighMJJLP(false),isLowMJJLP(false);
      if(chTag=="A") {        
        bool isPhotonTrigSafe(r9>0.9 && hoe<0.01 && eveto>0 && pixelseed<1);
        if(!APPLYPHOTONTRIGSAFE) isPhotonTrigSafe=true;
        isVBF    = (selector->hasTriggerBit(vbfPhotonTrigger, ev.triggerBits) 
                    && photons[0].Pt()>75 
                    && fabs(photons[0].Eta())<1.442
                    && isPhotonTrigSafe
                    && passVBFJetsTrigger);
        isHighPt = ( selector->hasTriggerBit(highPtPhotonTrigger, ev.triggerBits) 
                     && isPhotonTrigSafe
                     && photons[0].Pt()>minBosonHighPt);
        isHighPtAndOfflineVBF = (isHighPt && fabs(photons[0].Eta())<1.442 && passVBFJetsTrigger);
        isHighPtAndVBF = (isHighPt && isVBF);
        isBosonPlusOneJet=(isHighPt && nTotalJets==1);
	// A very simple categorization based on MJJ and boson Pt
	isHighMJJ = (isVBF && (photons[0].Pt() < minBosonHighPt) && passMJJ);
	isLowMJJ  = (passJets && isHighPt);
	//L1 Prefiltering check
	isHighMJJLP = (isHighMJJ & !passLP);
	isLowMJJLP = (isLowMJJ & !passLP);
        //veto prompt photons on the QCDEM enriched sample
        if( vetoPromptPhotons && ev.gamma_isPromptFinalState[ photons[0].originalReference() ] ) {
          isVBF                 = false;
          isHighPt              = false;
          isHighPtAndVBF        = false;
          isHighPtAndOfflineVBF = false;
          isBosonPlusOneJet     = false;
	  isHighMJJ             = false;
	  isLowMJJ              = false;
	  isHighMJJLP           = false;
	  isLowMJJLP            = false;
        }
      } else {
        isVBF    = boson.Pt()>75 && fabs(boson.Rapidity())<1.442 && passVBFJetsTrigger;
        isHighPt = boson.Pt()>minBosonHighPt;
        isHighPtAndVBF = (isHighPt && isVBF);
        isBosonPlusOneJet=(isHighPt && nTotalJets==1);
	// A very simple categorization based on MJJ and boson Pt
	isHighMJJ = (isVBF && boson.Pt() < minBosonHighPt && passMJJ);
	isLowMJJ  = (passJets && isHighPt);
	//L1 Prefiltering check
	isHighMJJLP = (isHighMJJ & !passLP);
	isLowMJJLP = (isLowMJJ & !passLP);
      }

      if(isVBF)                 cat[2]=true;
      if(isHighPt)	        cat[3]=true;
      if(isHighPtAndVBF)        cat[4]=true;
      if(isBosonPlusOneJet)     cat[5]=true;
      if(isHighPtAndOfflineVBF) cat[6]=true;
      if(isHighMJJ)             cat[7]=true;
      if(isLowMJJ)              cat[8]=true;
      if(isHighMJJLP)           cat[9]=true;
      if(isLowMJJLP)            cat[10]=true;
      category.set(cat);
      std::vector<TString> chTags( category.getChannelTags() );
      TString baseCategory(chTag);
      if(isHighPt)   baseCategory="HighPt"+chTag;
      else if(isVBF) baseCategory="VBF"+chTag;

      //evaluate discriminator MVA
      vbfmva = -1000;
      if(passJets) {
        TString key(isHighMJJ ?"BDT_VBF0HighMJJ":"BDT_VBF0LowMJJ");
        vbfmva = readers[key]->EvaluateMVA(key);
        if(mvaCDFinv[key]) vbfmva=max(0.,mvaCDFinv[key]->Eval(vbfmva));
        if(doBlindAnalysis && ev.isData && vbfmva>0.8) vbfmva=-1000;
      }
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<float>puWgts(3,1.0);
      EffCorrection_t trigSF(1.0,0.),selSF(1.0,0.),l1prefireProb(1.0,0.);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
            
        // pu weight
        ht->fill("puwgtctr",0,plotwgts);
        puWgts=lumi->pileupWeight(ev.g_pu,period);
        std::vector<double>puPlotWgts(1,puWgts[0]);
        ht->fill("puwgtctr",1,puPlotWgts);

        //l1 prefire probability
        EffCorrection_t l1prefireProb=l1PrefireWR->getJetBasedCorrection(jets);
        wgt *= l1prefireProb.first;

        // photon trigger*selection weights        
        if(chTag=="A")
          {
            trigSF = gammaEffWR->getTriggerCorrection({},photons,{}, period);
            selSF  = gammaEffWR->getOfflineCorrection(photons[0], period);
          }
        else
          {
            trigSF = gammaEffWR->getTriggerCorrection(leptons,{},{}, period);
            selSF  = gammaEffWR->getOfflineCorrection(leptons[0], period);
            EffCorrection_t sel2SF=gammaEffWR->getOfflineCorrection(leptons[1], period);
            selSF.first *= sel2SF.first;
            selSF.second = TMath::Sqrt( pow(selSF.second,2)+pow(sel2SF.second,2) );
          }
        wgt *= puWgts[0]*trigSF.first*selSF.first;
        
       
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //control histograms
      for( auto c : chTags) {        
        std::vector<double> cplotwgts(plotwgts);

        //photon pT weighting
        if(chTag=="A") {
          float photonPtWgt(1.0);
          if(photonPtWgts.find(c)!=photonPtWgts.end()) {
            photonPtWgt=photonPtWgts[c]->Eval(boson.Pt());
            if(photonPtWgt>0) photonPtWgt = 1./photonPtWgt;
            else              photonPtWgt = 1.0;
          }
          photonPtWgtCtr[c].first  += 1.0;
          photonPtWgtCtr[c].second += photonPtWgt;
          cplotwgts[0]*=photonPtWgt;

	  if(ev.isData && SRfake && (isHighMJJ || isLowMJJ))  {
	    cout << "The Fake Rate will be applied! " <<endl;
	    TString Cat = "LowMJJ";
	    if(isHighMJJ) Cat = "HighMJJ";
	    cplotwgts[0]*=fr->getWeight(Cat, vbfVars_.mjj, photons[0].Eta());
	  }
        }

        //fill plots
	evtWeight = cplotwgts[0]*xsec;
	training = useForTraining(); 
	fill( ev,  boson,  jets,  cplotwgts, c, mults, fakeACR, tightACR);
      }
        
      //experimental systs cycle: better not to do anything else after this...
      //final category selection is repeated ad nauseam with varied objects/weights and mva is re-evaluated several times
      if(ev.isData) continue;
      std::vector<std::pair<float,float> > mvaWithWeights;
      selector->setDebug(false);
      vbf::DiscriminatorInputs origVbfVars(vbfVars_);
      for(size_t is=0; is<expSysts_.size(); is++){
        
        //reset to the original values
        vbfVars_.assignValuesFrom(origVbfVars);

        //uncertainty
        TString sname=expSysts_[is];
        bool isUpVar(sname.Contains("up"));
        
        //base values and kinematics
        TString icat(baseCategory);
        float imva=vbfmva;
        float iwgt=(ev.g_nw>0 ? ev.g_w[0] : 1.0);
        TLorentzVector iBoson(boson);
        std::vector<Jet> ijets(jets);
        bool reSelect(false);        

        if(sname=="puup")        iwgt *= puWgts[1]*trigSF.first*selSF.first*l1prefireProb.first;
        if(sname=="pudn")        iwgt *= puWgts[2]*trigSF.first*selSF.first*l1prefireProb.first;
        if(sname=="trigup")      iwgt *= puWgts[0]*(trigSF.first+trigSF.second)*selSF.first*l1prefireProb.first;
        if(sname=="trigdn")      iwgt *= puWgts[0]*(trigSF.first-trigSF.second)*selSF.first*l1prefireProb.first;
        if(sname=="selup")       iwgt *= puWgts[2]*trigSF.first*(selSF.first+selSF.second)*l1prefireProb.first;
        if(sname=="seldn")       iwgt *= puWgts[2]*trigSF.first*(selSF.first+selSF.second)*l1prefireProb.first;
        if(sname=="l1prefireup") iwgt *= puWgts[1]*trigSF.first*selSF.first*(l1prefireProb.first+l1prefireProb.second);
        if(sname=="l1prefiredn") iwgt *= puWgts[1]*trigSF.first*selSF.first*(l1prefireProb.first+l1prefireProb.second);
        if(sname.Contains("aes") && chTag=="A")  {
          reSelect=true;
          iBoson *= (1+(isUpVar?1:-1)*bosonScaleUnc); 
        }
        if(sname.Contains("mes") && chTag=="MM") {
          //technically we should re-select the leptons but given we're looking to high pT Z's
          //assume effect is negligible and all that counts is the Z energy scale?
          reSelect=true;
          iBoson *= (1+(isUpVar?1:-1)*bosonScaleUnc); 
        }
        if(sname.Contains("JEC") || sname.Contains("JER") )  {
          reSelect=true;
          int jecIdx=-1;
          if(sname.Contains("AbsoluteStat"))     jecIdx=0;
          if(sname.Contains("AbsoluteScale"))    jecIdx=1; 
          if(sname.Contains("AbsoluteMPFBias"))  jecIdx=2; 
          if(sname.Contains("Fragmentation"))    jecIdx=3; 
          if(sname.Contains("SinglePionECAL"))   jecIdx=4; 
          if(sname.Contains("SinglePionHCAL"))   jecIdx=5; 
          if(sname.Contains("FlavorPureGluon"))  jecIdx=6; 
          if(sname.Contains("FlavorPureQuark"))  jecIdx=7; 
          if(sname.Contains("FlavorPureCharm"))  jecIdx=8; 
          if(sname.Contains("FlavorPureBottom")) jecIdx=9; 
          if(sname.Contains("TimePtEta"))        jecIdx=10; 
          if(sname.Contains("RelativeJEREC1"))   jecIdx=11; 
          if(sname.Contains("RelativeJEREC2"))   jecIdx=12; 
          if(sname.Contains("RelativeJERHF"))    jecIdx=13; 
          if(sname.Contains("RelativePtBB"))     jecIdx=14; 
          if(sname.Contains("RelativePtEC1"))    jecIdx=15; 
          if(sname.Contains("RelativePtEC2"))    jecIdx=16; 
          if(sname.Contains("RelativePtHF"))     jecIdx=17; 
          if(sname.Contains("RelativeBal"))      jecIdx=18; 
          if(sname.Contains("RelativeFSR"))      jecIdx=19; 
          if(sname.Contains("RelativeStatFSR"))  jecIdx=20; 
          if(sname.Contains("RelativeStatEC"))   jecIdx=21; 
          if(sname.Contains("RelativeStatHF"))   jecIdx=22; 
          if(sname.Contains("PileUpDataMC"))     jecIdx=23; 
          if(sname.Contains("PileUpPtRef"))      jecIdx=24; 
          if(sname.Contains("PileUpPtBB"))       jecIdx=25; 
          if(sname.Contains("PileUpPtEC1"))      jecIdx=26; 
          if(sname.Contains("PileUpPtEC2"))      jecIdx=27; 
          if(sname.Contains("PileUpPtHF"))       jecIdx=28;
          
          //re-scale and re-select jets
          std::vector<Jet> newJets = selector->getGoodJets(ev,20.,4.7,leptons,photons);
          ijets.clear();
          for(auto j : alljets) {

            int idx=j.getJetIndex();
            int jflav(abs(ev.j_flav[idx]));

            //shift jet energy
            float scaleVar(1.0);
            if(jecIdx<0) {
              scaleVar=isUpVar ? ev.j_jerUp[idx] : ev.j_jerDn[idx];
            } 
            else {
              bool flavorMatches(true);
              if(jecIdx==6 && jflav!=21) flavorMatches=false; //FlavorPureGluon
              if(jecIdx==7 && jflav>=4)  flavorMatches=false; //FlavorPureQuark
              if(jecIdx==8 && jflav!=4)  flavorMatches=false; //FlavorPureCharm
              if(jecIdx==9 && jflav!=5)  flavorMatches=false; //FlavorPureGluon
              if(flavorMatches)
                scaleVar=isUpVar ? ev.j_jecUp[jecIdx][idx] : ev.j_jecDn[jecIdx][idx];
            }
            j*=scaleVar;
            if(j.Pt()<30) continue;
            
            int jid=ev.j_id[idx];
            if(APPLYEECLEANING && fabs(j.Eta())>2.7 && fabs(j.Eta())<3 && ev.j_emf[idx]>0.55) continue;
            bool passPu((jid>>PUIDX)&0x1);
            if(!passPu) continue;
            ijets.push_back(j);
          }
        }
        
        //re-select if needed
        if(reSelect) {
          
          if (ijets.size()<2 || ijets[0].Pt()<LEADINGJETCUT) continue;

          vbfVars_.fillDiscriminatorVariables(iBoson,ijets,ev);
          if(vbfVars_.mjj<minMJJ) continue;
          
          //final event category
          bool passVBFJetsTrigger(vbfVars_.detajj>3.0 && vbfVars_.mjj>highMJJcut);
          bool isVBF(false),isHighPt(false);
          if(chTag=="A") {
            isVBF    = (selector->hasTriggerBit(vbfPhotonTrigger, ev.triggerBits) 
                        && iBoson.Pt()>75 
                        && fabs(iBoson.Eta())<1.442
                        && passVBFJetsTrigger);
            isHighPt = (selector->hasTriggerBit(highPtPhotonTrigger, ev.triggerBits) 
                        && iBoson.Pt()>minBosonHighPt);
          }
          else {
            isVBF    = (iBoson.Pt()>75 
                        && fabs(iBoson.Rapidity())<1.442 
                        && passVBFJetsTrigger);
            isHighPt = (iBoson.Pt()>minBosonHighPt);
          }
          
          //set the new tag
          if(isHighPt)   icat="HighPt"+chTag;
          else if(isVBF) icat="VBF"+chTag;
          else continue;
          
          //re-evaluate MVA         
          TString key(isHighMJJ ?"BDT_VBF0HighMJJ":"BDT_VBF0LowMJJ");
          imva = readers[key]->EvaluateMVA(key);
          if(mvaCDFinv[key]) imva=max(0.,mvaCDFinv[key]->Eval(imva));
        }

        //fill with new values/weights
        std::vector<double> eweights(1,iwgt);
        ht->fill2D("vbfmva_exp",imva,is,eweights,icat);
      }
      selector->setDebug(debug);
    }

  //close input file
  f->Close();

  //compute the scale factor needed to keep the normalization
  //due to photon pT weighting
  for(auto &wit : photonPtWgtCtr) {
    if(wit.second.second<=0) wit.second.first=1.0;
    else                     wit.second.first /= wit.second.second;
  }
  
  //save histos to file  
  saveHistos();
  if(skimtree){
    fMVATree->cd();
    newTree->Write();
    fMVATree->Close();
  }
}

void VBFVectorBoson::saveHistos(){
  fOut->cd();
  for (auto& it : ht->getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    for(auto &wit : photonPtWgtCtr){
      if(!it.first.Contains(wit.first)) continue;
      it.second->Scale(wit.second.first);
      break;
    }
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht->get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    for(auto &wit : photonPtWgtCtr){
      if(!it.first.Contains(wit.first)) continue;
      it.second->Scale(wit.second.first);
      break;
    }
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();		
}

void VBFVectorBoson::readTree(){
  f = TFile::Open(filename);
  triggerList=(TH1 *)f->Get("analysis/triggerList");
  t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  nentries = t->GetEntriesFast();
  if (debug) nentries = 10000; //restrict number of entries for testing
  //nentries = 10000;
  t->GetEntry(0);
  vetoPromptPhotons = filename.Contains("_QCDEM_") || filename.Contains("_TTJets");
  weightSysts_=getWeightSysts(f);
}

void VBFVectorBoson::prepareOutput(){
  baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  
  if(skimtree){
    fMVATree=TFile::Open(dirName+"/MVATree_"+baseName,"RECREATE");
    newTree = t->CloneTree(0);
  }
  fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
}

void VBFVectorBoson::bookHistograms(){
  ht = new HistTool(0);
  ht->addHist("puwgtctr",      new TH1F("puwgtctr",         ";Weight sums;Events",2,0,2));  
  ht->addHist("qscale",        new TH1F("qscale",           ";Q^{2} scale;Events",100,0,2000));  
  ht->addHist("nvtx",          new TH1F("nvtx",             ";Vertex multiplicity;Events",100,-0.5,99.5));  
  ht->addHist("vpt", 	       new TH1F("vectorbosonPt",    ";Boson p_{T}[GeV];Events",25,50,550));  
  ht->addHist("vy", 	       new TH1F("vectorbosony",     ";Boson rapidity;Events",25,-3,3));  
  ht->addHist("mindrl",        new TH1F("mindrl",           ";min #Delta R(boson,lepton);Events",25,0,6));  
  ht->addHist("mindrj",        new TH1F("mindrj",           ";min #Delta R(boson,jet);Events",25,0,6));  
  ht->addHist("sihih", 	       new TH1F("sihih",            ";#sigma(i#eta,i#eta);Events",50,0,0.1));  
  ht->addHist("hoe", 	       new TH1F("hoe",              ";h/e;Events",25,0,0.05));  
  ht->addHist("eveto", 	       new TH1F("eveto",            ";e veto;Events",2,0,2));
  ht->addHist("pixelseed",     new TH1F("pixelseed",        ";Pixel seed veto;Events",2,0,2) );
  ht->addHist("r9", 	       new TH1F("r9",               ";r9;Events",25,0.5,1.0));  
  ht->addHist("chiso", 	       new TH1F("chiso",            ";Charged isolation [GeV];Events",50,0,10));  
  ht->addHist("vystar",        new TH1F("vectorbosonystar", ";y-(1/2)(y_{j1}+y_{j2});Events",25,-5,5));  
  ht->addHist("njets",         new TH1F("njets",            ";Jet multiplicity;Events",10,-0.5,9.5));  
  ht->addHist("mjj", 	       new TH1F("mjj",              ";Dijet invariant mass [GeV];Events",40,0,4000));  
  ht->addHist("leadpt",        new TH1F("leadpt",           ";Leading jet p_{T} [GeV];Events",25,0,500));  
  ht->addHist("subleadpt",     new TH1F("subleadpt"   ,     ";Sub-leading jet p_{T} [GeV];Events",25,0,500));  
  ht->addHist("drj1b",         new TH1F("drj1b",            ";#DeltaR(j_{1},boson);Events",25,0,8));  
  ht->addHist("drj2b",         new TH1F("drj2b"   ,         ";#DeltaR(j_{2},boson);Events",25,0,8));  
  ht->addHist("leadpumva",     new TH1F("leadpumva",        ";Pileup MVA;Events",25,-1,1));  
  ht->addHist("subleadpumva",  new TH1F("subleadpumva"   ,  ";Pileup MVA;Events",25,-1,1));  
  ht->addHist("centraleta",    new TH1F("centraleta",       ";Most central jet |#eta|;Events",25,0,5));  
  ht->addHist("forwardeta",    new TH1F("forwardeta",       ";Most forward jet |#eta|;Events",25,0,5));  
  ht->addHist("dijetpt",       new TH1F("dijetpt",          ";Dijet p_{T} [GeV];Events",20,0,1000));  
  ht->addHist("detajj",        new TH1F("detajj" ,          ";#Delta#eta(J,J);Events",20,0,8));  
  ht->addHist("dphijj",        new TH1F("dphijj" ,          ";#Delta#phi(J,J) [rad];Events",20,-3.15,3.15));  
  ht->addHist("ht",            new TH1F("ht",               ";H_{T} [GeV];Events",20,0,4000));  
  ht->addHist("mht",           new TH1F("mht",              ";Missing H_{T} [GeV];Events",20,0,500));  
  ht->addHist("balance",       new TH1F("balance",          ";System p_{T} balance [GeV];Events",20,0,300));  
  ht->addHist("relbpt",        new TH1F("relbpt",           ";#Sigma p_{T}(j)/Boson p_{T};Events",25,0,2));  
  ht->addHist("dphibjj",       new TH1F("dphibjj",          ";#Delta#phi(JJ,boson) [rad];Events",20,-3.15,3.15));  
  ht->addHist("sphericity",    new TH1F("sphericity",       ";Sphericity;Events",20,0,1.0));  
  ht->addHist("aplanarity",    new TH1F("aplanarity",       ";Aplanarity;Events",20,0,1.0));  
  ht->addHist("C",             new TH1F("C",                ";C;Events",20,0,1.0));  
  ht->addHist("D",             new TH1F("D",                ";D;Events",20,0,1.0));  
  ht->addHist("isotropy",      new TH1F("isotropy",         ";Isotropy;Events",20,0,1.0));  
  ht->addHist("circularity",   new TH1F("circularity",      ";Circularity;;Events",20,0,1.0));
  //Photons in CR
  ht->addHist("allsihih",      new TH1F("allsihih",         ";All #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht->addHist("relaxedTightsihih",new TH1F("relaxedTightsihih",      ";Relaxed tight #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht->addHist("fakesihih",     new TH1F("fakesihih",        ";Fake #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht->addHist("tightsihih",    new TH1F("tightsihih",       ";Tight #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht->addHist("fakechiso",     new TH1F("fakechiso",        ";Fake ch. isolation [GeV];Photons",50,0,10));  
  ht->addHist("tightchiso",    new TH1F("tightchiso",       ";Tight ch. isolation [GeV];Photons",50,0,10));
  ht->addHist("fakeneutiso",   new TH1F("fakeneutiso",      ";Fake neut. isolation [GeV];Photons",50,0,10));  
  ht->addHist("tightneutiso",  new TH1F("tightneutiso",     ";Tight neut. isolation [GeV];Photons",50,0,10));
  ht->addHist("fakeaiso",      new TH1F("fakeaiso",         ";Fake #gamma isolation [GeV];Photons",50,0,10));  
  ht->addHist("tightaiso",     new TH1F("tightaiso",        ";Tight #gamma isolation [GeV];Photons",50,0,10));  
  ht->addHist("nloose",        new TH1F("nloose",           ";Number of loose #gamma; Events",20,-0.5,19.5));
  ht->addHist("ntight",        new TH1F("ntight",           ";Number of tight #gamma; Events",20,-0.5,19.5));
  ht->addHist("nloosefake",    new TH1F("nloosefake",       ";Number of fake loose #gamma; Events",20,-0.5,19.5));
  ht->addHist("ntightfake",    new TH1F("ntightfake",       ";Number of fake tight #gamma; Events",20,-0.5,19.5));
  ht->addHist("nlooseprompt",  new TH1F("nlooseprompt",     ";Number of prompt loose #gamma; Events",20,-0.5,19.5));
  ht->addHist("ntightprompt",  new TH1F("ntightprompt",     ";Number of prompt tight #gamma; Events",20,-0.5,19.5));
  //2D's for Mjj-binned FR
  double bins[]={0,500,1000,2000,4000};
  
  ht->addHist("relaxedTightMjjEB",  new TH2F("relaxedTightMjjEB",";Relaxed tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins)); //80,0,4000
  ht->addHist("tightMjjEB",         new TH2F("tightMjjEB",";Tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht->addHist("looseMjjEB",         new TH2F("looseMjjEB",";Loose #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht->addHist("allMjjEB",           new TH2F("allMjjEB",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht->addHist("tmpQCDMjjEB",        new TH2F("tmpQCDMjjEB",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));

  ht->addHist("relaxedTightMjjEE",  new TH2F("relaxedTightMjjEE",";Relaxed tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins)); //80,0,4000
  ht->addHist("tightMjjEE",         new TH2F("tightMjjEE",";Tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht->addHist("looseMjjEE",         new TH2F("looseMjjEE",";Loose #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht->addHist("allMjjEE",           new TH2F("allMjjEE",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht->addHist("tmpQCDMjjEE",        new TH2F("tmpQCDMjjEE",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  // Study of jet variables
  ht->addHist("etaphi",           new TH2F("etaphi",          ";Most central jet |#eta|V|#phi|;Events",30,-4.7,4.7,25,-TMath::Pi(),TMath::Pi()));
  ht->addHist("jet_raw_pt", 	  new TH1F("jet_raw_pt",      ";raw PT of jets;Jets",50,0,200));
  ht->addHist("jet_raw_empt", 	  new TH1F("jet_raw_empt",    ";raw e.m. PT of jets;Jets",50,0,200));
  ht->addHist("jet_emf", 	  new TH1F("jet_emf",         ";Electromagnetic fraction;Jets",50,0,1));
  ht->addHist("jet_qg", 	  new TH1F("jet_qg",          ";q/g discriminator;Jets",50,0,1));
  ht->addHist("jet_c2_00", 	  new TH1F("jet_c2_00",       ";Energy correlation C_{2}^{(0)};Jets",50,0,1));  
  ht->addHist("jet_c2_02", 	  new TH1F("jet_c2_02",       ";Energy correlation C_{2}^{(2)};Jets",50,0,1));  
  ht->addHist("jet_c2_05", 	  new TH1F("jet_c2_05",       ";Energy correlation C_{2}^{(5)};Jets",50,0,1));  
  ht->addHist("jet_zg", 	  new TH1F("jet_zg",          ";Jet z_{g};Jets",50,0,1));  
  ht->addHist("jet_gaptd", 	  new TH1F("jet_gaptd",       ";Jet p_{T}^{D};Jets",50,0,1));  
  ht->addHist("jet_gawidth",      new TH1F("jet_gawidth",     ";Jet width;Jets",50,0,1));
  //additional variables from https://link.springer.com/content/pdf/10.1140/epjc/s10052-017-5315-6.pdf
  ht->addHist("jjetas", 	  new TH1F("jjetas",          ";#eta_{j1}#eta_{j2};Events",50,-25,15));  
  ht->addHist("centjy",		  new TH1F("centjy",          ";Central jet rapidity;Jets",25,0,3));  
  ht->addHist("ncentj", 	  new TH1F("ncentjj",          ";Number of central jets;Events",10,-0.5,9.5));  
  ht->addHist("dphivj0", 	  new TH1F("dphivj0",          ";#Delta#phi(V,j0);Jets",20,0,TMath::Pi()));  
  ht->addHist("dphivj1", 	  new TH1F("dphivj1",          ";#Delta#phi(V,j1);Jets",20,0,TMath::Pi()));  
  ht->addHist("dphivj2", 	  new TH1F("dphivj2",          ";#Delta#phi(V,j2);Jets",20,0,TMath::Pi()));  
  ht->addHist("dphivj3", 	  new TH1F("dphivj3",          ";#Delta#phi(V,j3);Jets",20,0,TMath::Pi()));
  //final analyses distributions
  ht->addHist("evcount",         new TH1F("evcount",        ";Pass;Events",3,0,3));  
  ht->addHist("vbfmva",          new TH1F("vbfmva",         ";VBF MVA;Events",50,0,1));  

  TString expSystNames[]={"puup","pudn","trigup","trigdn","selup","seldn","l1prefireup","l1prefiredn",
                          "aesup","aesdn",
                          "mesup","mesdn",
                          "JERup","JERdn",
                          "AbsoluteStatJECup","AbsoluteScaleJECup","AbsoluteMPFBiasJECup","FragmentationJECup","SinglePionECALJECup","SinglePionHCALJECup","FlavorPureGluonJECup","FlavorPureQuarkJECup","FlavorPureCharmJECup","FlavorPureBottomJECup","TimePtEtaJECup","RelativeJEREC1JECup","RelativeJEREC2JECup","RelativeJERHFJECup","RelativePtBBJECup","RelativePtEC1JECup","RelativePtEC2JECup","RelativePtHFJECup","RelativeBalJECup","RelativeFSRJECup","RelativeStatFSRJECup","RelativeStatECJECup","RelativeStatHFJECup","PileUpDataMCJECup","PileUpPtRefJECup","PileUpPtBBJECup","PileUpPtEC1JECup","PileUpPtEC2JECup","PileUpPtHFJECup",
                          "AbsoluteStatJECdn","AbsoluteScaleJECdn","AbsoluteMPFBiasJECdn","FragmentationJECdn","SinglePionECALJECdn","SinglePionHCALJECdn","FlavorPureGluonJECdn","FlavorPureQuarkJECdn","FlavorPureCharmJECdn","FlavorPureBottomJECdn","TimePtEtaJECdn","RelativeJEREC1JECdn","RelativeJEREC2JECdn","RelativeJERHFJECdn","RelativePtBBJECdn","RelativePtEC1JECdn","RelativePtEC2JECdn","RelativePtHFJECdn","RelativeBalJECdn","RelativeFSRJECdn","RelativeStatFSRJECdn","RelativeStatECJECdn","RelativeStatHFJECdn","PileUpDataMCJECdn","PileUpPtRefJECdn","PileUpPtBBJECdn","PileUpPtEC1JECdn","PileUpPtEC2JECdn","PileUpPtHFJECdn"};
  
  size_t nexpSysts=sizeof(expSystNames)/sizeof(TString);
  expSysts_=std::vector<TString>(expSystNames,expSystNames+nexpSysts);  
  ht->addHist("vbfmva_exp",      new TH2F("vbfmva_exp",     ";VBF MVA;Systs;Events",20,0,1,nexpSysts,0,nexpSysts));
  for(size_t is=0; is<nexpSysts; is++)
    ht->get2dPlots()["vbfmva_exp"]->GetYaxis()->SetBinLabel(is+1,expSystNames[is]);
 
  size_t nthSysts(weightSysts_.size());
  if(nthSysts>0){
    ht->addHist("vbfmva_th",       new TH2F("vbfmva_th",      ";VBF MVA;Systs;Events",20,0,1,nthSysts,0,nthSysts));  
    for(size_t is=0; is<nthSysts; is++)
      ht->get2dPlots()["vbfmva_th"]->GetYaxis()->SetBinLabel(is+1,weightSysts_[is].first);
  }
}

void VBFVectorBoson::setGammaZPtWeights(){
  TString wgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson/raw/plots/ratio_plotter.root");
  gSystem->ExpandPathName(wgtUrl);
  TFile *wgtF=TFile::Open(wgtUrl);
  if(wgtF) {
    cout << "Reading photon/Z pT weights" << endl;
    TString pfix(baseName.Contains("Data13TeV_") ? "" : "_mc_MC");
    photonPtWgts["VBFA"]      = new TGraph((TH1* )wgtF->Get("VBFA_vectorbosonPt_ratio/VBFA_vectorbosonPt"+pfix));
    photonPtWgts["HighPtA"]   = new TGraph((TH1* )wgtF->Get("HighPtA_vectorbosonPt_ratio/HighPtA_vectorbosonPt"+pfix));
    photonPtWgtCtr["VBFA"]    = std::pair<double,double>(0.0,0.0);
    photonPtWgtCtr["HighPtA"] = std::pair<double,double>(0.0,0.0);
    wgtF->Close();
  } else {
    cout << "Requested to reweight photon spectrum but could not find " << wgtUrl << endl
         << "Proceeding without" << endl;
  }
}
void VBFVectorBoson::loadCorrections(){
  fr = new FakeRateTool(era, "fakeRatios.root");
  lumi = new LumiTools(era,genPU);

  std::map<TString,TString> cfgMap;
  cfgMap["g_id"]="Tight";
  cfgMap["m_id"]="TightID";
  cfgMap["m_iso"]="TightRelIso";
  cfgMap["m_id4iso"]="TightIDandIPCut";
  gammaEffWR = new EfficiencyScaleFactorsWrapper(filename.Contains("Data13TeV"),era,cfgMap);
  //jec = new JECTools(era);
  l1PrefireWR = new L1PrefireEfficiencyWrapper(filename.Contains("Data13TeV"),era);
  if(anFlag>0) this->setGammaZPtWeights();
}
void VBFVectorBoson::addMVAvars(){
  newTree->Branch("centralEta",       &vbfVars_.centraleta);
  newTree->Branch("forwardeta",       &vbfVars_.forwardeta);
  newTree->Branch("leadj_pt",         &vbfVars_.leadj_pt);
  newTree->Branch("subleadj_pt",      &vbfVars_.subleadj_pt);
  newTree->Branch("mjj",              &vbfVars_.mjj);
  newTree->Branch("detajj",           &vbfVars_.detajj);
  newTree->Branch("jjpt",             &vbfVars_.jjpt);
  newTree->Branch("dphijj",           &vbfVars_.dphijj);
  newTree->Branch("ystar",            &vbfVars_.ystar);
  newTree->Branch("relbpt",           &vbfVars_.relbpt);
  newTree->Branch("dphibjj",          &vbfVars_.dphibjj);
  newTree->Branch("balance",          &vbfVars_.balance);
  newTree->Branch("leadj_gawidth",    &vbfVars_.leadj_gawidth);
  newTree->Branch("subleadj_gawidth", &vbfVars_.subleadj_gawidth);
  newTree->Branch("subleadj_c2_02",   &vbfVars_.subleadj_c2_02);
  newTree->Branch("jjetas",           &vbfVars_.jjetas);
  newTree->Branch("centjy",           &vbfVars_.centjy);
  newTree->Branch("ncentjj",          &vbfVars_.ncentj);
  newTree->Branch("dphivj0",          &vbfVars_.dphivj0);
  newTree->Branch("dphivj1",          &vbfVars_.dphivj1);
  newTree->Branch("dphivj2",          &vbfVars_.dphivcentj[0]);
  newTree->Branch("dphivj3",          &vbfVars_.dphivcentj[0]);
  newTree->Branch("mht",              &vbfVars_.mht);
  newTree->Branch("ht",               &vbfVars_.scalarht);
  newTree->Branch("isotropy",         &vbfVars_.isotropy);
  newTree->Branch("circularity",      &vbfVars_.circularity);
  newTree->Branch("sphericity",       &vbfVars_.sphericity);
  newTree->Branch("aplanarity",       &vbfVars_.aplanarity);
  newTree->Branch("C",                &vbfVars_.C);
  newTree->Branch("D",                &vbfVars_.D);
  newTree->Branch("evtWeight",        &evtWeight);
  newTree->Branch("training",         &training);
  newTree->Branch("category",         &category, "MM:A:VBF:HighPt:HighPtVBF:V1J:HighPtOfflineVBF:HighMJJ:LowMJJ:HighMJJLP:LowMJJLP");
}


void VBFVectorBoson::fill(MiniEvent_t ev, TLorentzVector boson, std::vector<Jet> jets, std::vector<double> cplotwgts, TString c, std::map<TString, int> mults, std::vector<Particle> fakeACR, std::vector<Particle> tightACR){
  ht->fill("nvtx",   ev.nvtx,          cplotwgts,c);        

  //boson histos
  ht->fill("vpt",    boson.Pt(),       cplotwgts,c);
  ht->fill("vy",     boson.Rapidity(), cplotwgts,c);   
  ht->fill("sihih",  sihih,            cplotwgts,c);   
  ht->fill("r9",     r9,               cplotwgts,c);   
  ht->fill("hoe",    hoe,              cplotwgts,c);   
  ht->fill("eveto",      eveto,               cplotwgts,c);   
  ht->fill("pixelseed",    pixelseed,              cplotwgts,c);   
  ht->fill("chiso",  chiso,            cplotwgts,c);   
  ht->fill("mindrl", mindrl,           cplotwgts,c);   
  ht->fill("mindrj", mindrj,           cplotwgts,c);   

  for(auto a : photons) {
    int idx = a.originalReference();
    ht->fill("allsihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("allMjjEB"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
    else
      ht->fill2D("allMjjEE"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
  }
  for(auto a : relaxedTightPhotons) {
    int idx = a.originalReference();
    ht->fill("relaxedTightsihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("relaxedTightMjjEB"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
    else
      ht->fill2D("relaxedTightMjjEE"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
  }
  for(auto a : tmpPhotons) {
    int idx = a.originalReference();
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("tmpQCDMjjEB"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
    else
      ht->fill2D("tmpQCDMjjEE"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
  }
  //bosons in CR and fakes
  for(auto a : fakeACR) {
    int idx = a.originalReference();
    ht->fill("fakesihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    ht->fill("fakechiso",   ev.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
    ht->fill("fakeneutiso", ev.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
    ht->fill("fakeaiso",    ev.gamma_photonIso[idx]         ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("looseMjjEB"  ,  ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
    else
      ht->fill2D("looseMjjEE"  ,  ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
  }
  for(auto a : tightACR) { 
    int idx = a.originalReference();
    ht->fill("tightsihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    ht->fill("tightchiso",   ev.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
    ht->fill("tightneutiso", ev.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
    ht->fill("tightaiso",    ev.gamma_photonIso[idx]         ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("tightMjjEB"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
    else
      ht->fill2D("tightMjjEE"  ,   ev.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
  }
  ht->fill("nloose",        fakeACR.size(),       cplotwgts,c);
  ht->fill("ntight",        tightACR.size(),      cplotwgts,c);
  ht->fill("nloosefake",    mults["loosefake"],   cplotwgts,c);
  ht->fill("ntightfake",    mults["tightfake"],   cplotwgts,c);
  ht->fill("nlooseprompt",  mults["looseprompt"], cplotwgts,c);
  ht->fill("ntightprompt",  mults["tightprompt"], cplotwgts,c);

  //jet histos
  for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
    TString jtype(ij==0?"lead":"sublead");
    ht->fill(jtype+"pt",       jets[ij].Pt(),        cplotwgts,c);          
    ht->fill(jtype+"pumva",    jets[ij].getPUMVA(),  cplotwgts,c);
    ht->fill(Form("drj%db",(int)ij+1),   jets[ij].DeltaR(boson),  cplotwgts,c);
    ht->fill("jet_c2_00", 	ev.j_c2_00[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_c2_02", 	ev.j_c2_02[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_c2_05",	ev.j_c2_05[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_zg", 		ev.j_zg[jets[ij].getJetIndex()]	          ,  cplotwgts,c);
    ht->fill("jet_gaptd", 	ev.j_gaptd[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_gawidth",     ev.j_gawidth[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    float j_rawpt   = jets[ij].pt()/ev.j_rawsf[jets[ij].getJetIndex()];
    float j_emf     = ev.j_emf[jets[ij].getJetIndex()];
    float j_rawempt = j_rawpt*j_emf;
    ht->fill("jet_emf"     , j_emf        ,cplotwgts,c);
    ht->fill("jet_raw_pt" ,j_rawpt        , cplotwgts,c );
    ht->fill("jet_raw_empt" ,j_rawempt        , cplotwgts,c );
    ht->fill("jet_qg"   ,ev.j_qg[jets[ij].getJetIndex()]    , cplotwgts,c);
    ht->fill2D("etaphi",  jets[ij].Eta(),jets[ij].Phi() ,   cplotwgts,c);    
  }
 
  
  if(jets.size() >= 2){
    ht->fill("jjetas",  vbfVars_.jjetas,   cplotwgts,c);
    ht->fill("dphivj0", vbfVars_.dphivj0 ,  cplotwgts,c);
    ht->fill("dphivj1", vbfVars_.dphivj1 ,  cplotwgts,c);
  }

  //central jet activity
  if(vbfVars_.ncentj>0){
    ht->fill("ncentj", vbfVars_.ncentj, cplotwgts, c);
    ht->fill("dphivj2", vbfVars_.dphivcentj[0] ,  cplotwgts,c);
    if(vbfVars_.ncentj>1) 
      ht->fill("dphivj3", vbfVars_.dphivcentj[1] ,  cplotwgts,c);    
  }

  ht->fill("njets",        jets.size(), cplotwgts,c);
  ht->fill("ht",           vbfVars_.scalarht,    cplotwgts,c);
  ht->fill("mht",          vbfVars_.mht,         cplotwgts,c);
  ht->fill("centraleta",   vbfVars_.centraleta,  cplotwgts,c);
  ht->fill("forwardeta",   vbfVars_.forwardeta,  cplotwgts,c);
  ht->fill("dijetpt",      vbfVars_.jjpt,        cplotwgts,c);
  ht->fill("detajj",       vbfVars_.detajj,      cplotwgts,c);
  ht->fill("dphijj",       vbfVars_.dphijj,      cplotwgts,c);
  ht->fill("relbpt",       vbfVars_.relbpt,      cplotwgts,c);
  ht->fill("dphibjj",      vbfVars_.dphibjj,     cplotwgts,c);
  ht->fill("mjj", 	   vbfVars_.mjj,         cplotwgts,c);
	
  //visible system histos
  ht->fill("vystar",       vbfVars_.ystar,              cplotwgts,c);        
  ht->fill("balance",      vbfVars_.balance,            cplotwgts,c);
  ht->fill("isotropy",     vbfVars_.isotropy,     cplotwgts,c);
  ht->fill("circularity",  vbfVars_.circularity,  cplotwgts,c);
  ht->fill("sphericity",   vbfVars_.sphericity, cplotwgts,c);
  ht->fill("aplanarity",   vbfVars_.aplanarity, cplotwgts,c);
  ht->fill("C",            vbfVars_.C,          cplotwgts,c);
  ht->fill("D",            vbfVars_.D,          cplotwgts,c);

  //final analysis histograms
  ht->fill("evcount",  0, cplotwgts, c);
  if(vbfmva>=0) {
    ht->fill("vbfmva", vbfmva, cplotwgts,c);
    if(vbfmva>0.9)
      ht->fill("evcount",  1, cplotwgts, c);
    if(vbfmva>0.9 && vbfVars_.mjj>1000.)
      ht->fill("evcount",  2, cplotwgts, c);
  }
  if(!ev.isData){
    //replicas for theory systs
    for(size_t is=0; is<weightSysts_.size(); is++){
      std::vector<double> sweights(1,cplotwgts[0]);
      size_t idx=weightSysts_[is].second;
      sweights[0] *= (ev.g_w[idx]/ev.g_w[0])*(normH->GetBinContent(idx+1)/normH->GetBinContent(1));
      ht->fill2D("vbfmva_th",vbfmva,is,sweights,c);
    }
  }

  if(skimtree) newTree->Fill();
}
