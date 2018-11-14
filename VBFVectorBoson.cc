
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

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace std;

//TODOS
// pedro: jet distribution before pu id in boson+1jet
//

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
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  reader->AddVariable("ht",&scalarht);
  reader->AddVariable("forwardeta",&forwardeta);
  //  reader->AddVariable("j_c2_02[1]",&subleadj_c2_02);
  reader->AddVariable("balance",&balance);
  //  reader->AddVariable("aplanarity",&aplanarity);
  reader->AddVariable("dphivj0",&dphivj0);
  reader->AddVariable("mjj", &mjj);
  reader->AddVariable("j_qg[0]", &lead_qg) ;
  //  reader->AddVariable("D",&D);
  //  reader->AddVariable("j_pt[1]",&subleadj_pt);
  //  reader->AddVariable("j_gawidth[0]",&leadj_gawidth);
  //  reader->AddVariable("j_c1_05[0]",&leadj_c1_05);






  TString weightFiles[]={//"${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/BDT.weights.xml",
                         //"${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/Fisher.weights.xml",
			 //"${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/FisherHighMJJ.weights.xml",
                         //"${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/FisherLowMJJ.weights.xml"};
    "${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/BDTHighMJJ.weights.xml",
    "${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/BDTLowMJJ.weights.xml"};
  TString mvaMethod[]={//"BDT",
		       //"Fisher",
		       // "FisherHighMJJ",
		       // "FisherLowMJJ"};
		       "BDT_VBF0HighMJJ",
		       "BDT_VBF0LowMJJ"};

  for(size_t i=0; i<2; i++){
    gSystem->ExpandPathName(weightFiles[i]);
    reader->BookMVA(mvaMethod[i], weightFiles[i] );
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
      std::vector<Jet>      &alljets     = selector->getJets();
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
	bool passLoosePu((jid>>2)&0x1);
	if(CR){
	  if(jets.size() == 0 && passLoosePu)
	    continue;
	} else if(!passLoosePu) continue;
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

      //jet related variables and selection
      initVariables(jets);

      scalarht = 0.;
      TLorentzVector mhtP4(0,0,0,0);
      mht = 0;
      int nLPJets(0);
      for(auto j : jets) {
        scalarht += j.Pt();
	//jet_pt_raw=j.pt()/j.rawsf();
	//jet_emf=j.emf();
	jet_phi= j.Phi();
        mhtP4 += j;
	if(j.Pt()>100 && fabs(j.Eta()) > 2.25 && fabs(j.Eta()) < 3.0)
	  nLPJets++;
      }
      mht = mhtP4.Pt();
      bool passJetMult(jets.size()>=2);
      bool passMJJ(passJetMult && mjj>highMJJcut);
      bool passJets(passJetMult && mjj>minMJJ);
      bool passVBFJetsTrigger(passJets && detajj>3.0);
      //L1-prefiltering check
      bool passLP(nLPJets > 0 || nLPGamma > 0);

      //categorize the event according to the boson kinematics
      //for the photon refine also the category according to the trigger  bit
      TLorentzVector boson(0,0,0,0);     
      bool isHighPt(false),isVBF(false),isHighPtAndVBF(false),isHighPtAndOfflineVBF(false),isBosonPlusOneJet(false),isHighMJJ(false),isLowMJJ(false), isHighMJJLP(false),isLowMJJLP(false);
      sihih = 0, chiso = 0 ,r9 = 0, hoe = 0;
      if(chTag=="A") {        
        boson += photons[0];
        sihih = ev.gamma_sieie[photons[0].originalReference()];
        chiso = ev.gamma_chargedHadronIso[photons[0].originalReference()];
        r9    = ev.gamma_r9[photons[0].originalReference()];
        hoe   = ev.gamma_hoe[photons[0].originalReference()];        
        isVBF    = (selector->hasTriggerBit(vbfPhotonTrigger, ev.triggerBits) 
                    && photons[0].Pt()>75 
                    && fabs(photons[0].Eta())<1.442
                    && passVBFJetsTrigger);
        isHighPt = ( selector->hasTriggerBit(highPtPhotonTrigger, ev.triggerBits) 
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
        if( isQCDEMEnriched && ev.gamma_isPromptFinalState[ photons[0].originalReference() ] ) {
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
        boson   += leptons[0];
        boson   += leptons[1];
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

      //if(!isVBF && !isHighPt && !isBosonPlusOneJet) continue;      

      //leptons and boson
      double mindrl(9999.);
      for(auto &l: leptons) mindrl=min(l.DeltaR(boson),mindrl);

      if(isVBF)                 cat[2]=true;
      if(isHighPt)	        cat[3]=true;
      if(isHighPtAndVBF)        cat[4]=true;
      if(isBosonPlusOneJet)     cat[5]=true;
      if(isHighPtAndOfflineVBF) cat[6]=true;
      if(isHighMJJ)             cat[7]=true;
      if(isLowMJJ)              cat[8]=true;
      if(isHighMJJLP)           cat[9]=true;
      if(isLowMJJLP)              cat[10]=true;
      category.set(cat);
      std::vector<TString> chTags( category.getChannelTags() );
      
      //leptons and boson
      mindrl = 9999.;
      for(auto &l: leptons) mindrl=min(l.DeltaR(boson),mindrl);
      //system variables and event shapes
      
      if(passJets) {
        ystar=boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
        balance=(boson+jets[0]+jets[1]).Pt();
        relbpt=(jets[0].Pt()+jets[1].Pt())/boson.Pt();
        dphibjj=boson.DeltaPhi( jets[0]+jets[1] );
      }
      else if(jets.size()>0){
        balance=(boson+jets[0]).Pt();
        relbpt=jets[0].Pt()/boson.Pt();
        dphibjj=boson.DeltaPhi(jets[0]);
      }
        
      std::vector<math::XYZVector> inputVectors;
      inputVectors.push_back( math::XYZVector(boson.Px(),boson.Py(),boson.Pz()) );
      for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
        inputVectors.push_back( math::XYZVector(jets[ij].Px(),jets[ij].Py(),jets[ij].Pz()) );
      }
      EventShapeVariables esv(inputVectors);
      isotropy    = esv.isotropy();
      circularity = esv.circularity();
      sphericity  = esv.sphericity(1.);
      aplanarity  = esv.aplanarity(1.);
      C           = esv.C(1.);
      D           = esv.D(1.);
      dphivj0 = (jets.size() >= 2) ? fabs(jets[0].DeltaPhi(boson)) : -99;
      
      vbffisher = passJets ? reader->EvaluateMVA(mvaMethod[0]) : -99;
      if(isHighMJJ) 
	vbfmva = passJets ? reader->EvaluateMVA(mvaMethod[0]) : -99;
      else
	vbfmva = passJets ? reader->EvaluateMVA(mvaMethod[1]) : -99;
      if(doBlindAnalysis && ev.isData && vbfmva>0.1)    vbfmva=-1000;
      if(doBlindAnalysis && ev.isData && vbffisher>0.1) vbffisher=-1000;

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
            
        // pu weight
        ht->fill("puwgtctr",0,plotwgts);
        double puWgt(lumi->pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht->fill("puwgtctr",1,puPlotWgts);

        //l1 prefire probability
        float l1prefireProb=l1PrefireWR->getJetBasedCorrection(jets).first;
        wgt *= l1prefireProb;

        // photon trigger*selection weights
        float trigSF(1.0), selSF(1.0);
        if(chTag=="A")
          {
            trigSF *= gammaEffWR->getTriggerCorrection({},photons,{}, period).first;
            selSF  *= gammaEffWR->getOfflineCorrection(photons[0], period).first;
          }
        else
          {
            trigSF *=gammaEffWR->getTriggerCorrection(leptons,{},{}, period).first;
            selSF  *=gammaEffWR->getOfflineCorrection(leptons[0], period).first;
            selSF  *=gammaEffWR->getOfflineCorrection(leptons[1], period).first;
          }
        wgt *= puWgt*trigSF*selSF;
        
       
        // generaor level weights
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
	    cplotwgts[0]*=fr->getWeight(Cat, mjj, photons[0].Eta());
	  }
        }


	//What is the final weight? 0 or 1 in the array?
	evtWeight = cplotwgts[0]*xsec;
	training = useForTraining(); 
	fill( ev,  boson,  jets,  cplotwgts, c, mults, fakeACR, tightACR);
      }
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
  isQCDEMEnriched = filename.Contains("_QCDEM_");
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
  ht->addHist("sihih", 	       new TH1F("sihih",            ";#sigma(i#eta,i#eta);Events",50,0,0.1));  
  ht->addHist("hoe", 	       new TH1F("hoe",              ";h/e;Events",25,0,0.1));  
  ht->addHist("r9", 	       new TH1F("r9",               ";r9;Events",25,0,1.0));  
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
  int Xbin=30, Ybin=25;
  double X1=0.0 , X2=5.0,y1=-3.0 , y2=3.0;
  ht->addHist("etaphi",           new TH2F("etaphi",       ";Most central jet |#eta|V|#phi|;Events",Xbin,X1,X2,Ybin,y1,y2));
  ht->addHist("jet_raw_pt", 	  new TH1F("jet_raw_pt",          ";raw PT of jets;Jets",100,-1,1));
  ht->addHist("jet_emf", 	  new TH1F("jet_emf",          ";EM effect of jets;Jets",100,0,1));
  ht->addHist("jet_qg", 	  new TH1F("jet_qg",          ";qg of jets;Jets",100,-1,1));
  // ht->addHist("jet_pumva", 	  new TH1F("jet_pumva",          ";pileup mva of jets;Jets",100,-1,1));
  ht->addHist("jet_c2_00", 	  new TH1F("jet_c2_00",          ";Jet shape var. c2_00;Jets",100,-1,1));  
  ht->addHist("jet_c2_02", 	  new TH1F("jet_c2_02",          ";Jet shape var. c2_02;Jets",100,-1,1));  
  ht->addHist("jet_c2_05", 	  new TH1F("jet_c2_05",          ";Jet shape var. c2_05;Jets",100,-1,1));  
  ht->addHist("jet_zg", 	  new TH1F("jet_zg",          ";Jet shape var. zg;Jets",100,-1,1));  
  ht->addHist("jet_gaptd", 	  new TH1F("jet_gaptd",          ";Jet shape var. gaptd;Jets",100,-1,1));  
  ht->addHist("jet_gawidth",      new TH1F("jet_gawidth",          ";Jet shape var. gawidth;Jets",100,-1,1));
  //additional variables from https://link.springer.com/content/pdf/10.1140/epjc/s10052-017-5315-6.pdf
  ht->addHist("jjetas", 	  new TH1F("jjetas",          ";#eta_{j1}#eta_{j2};Events",200,-25,25));  
  ht->addHist("centjy",		  new TH1F("centjy",          ";Central jet rapidity;Jets",25,0,3));  
  ht->addHist("ncentj", 	  new TH1F("ncentjj",          ";Number of central jets;Events",10,-0.5,9.5));  
  ht->addHist("dphivj0", 	  new TH1F("dphivj0",          ";#Delta#phi(V,j0);Jets",20,0,4));  
  ht->addHist("dphivj1", 	  new TH1F("dphivj1",          ";#Delta#phi(V,j1);Jets",20,0,4));  
  ht->addHist("dphivj2", 	  new TH1F("dphivj2",          ";#Delta#phi(V,j2);Jets",20,0,4));  
  ht->addHist("dphivj3", 	  new TH1F("dphivj3",          ";#Delta#phi(V,j3);Jets",20,0,4));
  //final analyses distributions
  ht->addHist("evcount",         new TH1F("evcount",        ";Pass;Events",1,0,1));  
  ht->addHist("vbfmva",          new TH1F("vbfmva",         ";VBF MVA;Events",20,-1,1));  
  ht->addHist("vbffisher",       new TH1F("vbffisher",      ";VBF Fisher;Events",40,-2,3));  
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
  gammaEffWR = new EfficiencyScaleFactorsWrapper(filename.Contains("Data13TeV"),era);
  //jec = new JECTools(era);
  l1PrefireWR = new L1PrefireEfficiencyWrapper(filename.Contains("Data13TeV"),era);
  if(anFlag>0) this->setGammaZPtWeights();
}
void VBFVectorBoson::addMVAvars(){
  newTree->Branch("centralEta", &centraleta);
  newTree->Branch("etaphi",&etaphi);
  newTree->Branch("subleadj_pt", &subleadj_pt);
  newTree->Branch("mjj", &mjj);
  newTree->Branch("detajj", &detajj);
  newTree->Branch("jjpt", &jjpt);
  newTree->Branch("dphijj", &dphijj);
  newTree->Branch("ystar", &ystar);
  newTree->Branch("relbpt", &relbpt);
  newTree->Branch("dphibjj", &dphibjj);
  newTree->Branch("balance", &balance);
  newTree->Branch("forwardeta", &forwardeta);
  newTree->Branch("leadj_gawidth", &leadj_gawidth);
  newTree->Branch("subleadj_gawidth", &subleadj_gawidth);
  //newTree->Branch("leadj_c1_05", &leadj_c1_05);
  // newTree->Branch("subleadj_c1_05", &subleadj_c1_05);
  newTree->Branch("subleadj_c2_02", &subleadj_c2_02);
  newTree->Branch("jjetas", &jjetas);
  newTree->Branch("centjy", &centjy);
  newTree->Branch("ncentjj", &ncentjj);
  newTree->Branch("dphivj0", &dphivj0);
  newTree->Branch("dphivj1", &dphivj1);
  newTree->Branch("dphivj2", &dphivj2);
  newTree->Branch("dphivj3", &dphivj3);
  newTree->Branch("evtWeight", &evtWeight);
  newTree->Branch("mht", &mht);
  newTree->Branch("balance", &balance);
  newTree->Branch("ht", &scalarht);
  newTree->Branch("isotropy", &isotropy);
  newTree->Branch("circularity",&circularity);
  newTree->Branch("sphericity",&sphericity);
  newTree->Branch("aplanarity",&aplanarity);
  newTree->Branch("C",&C);
  newTree->Branch("D",&D);
  newTree->Branch("training",&training);
  newTree->Branch("category", &category, "MM:A:VBF:HighPt:HighPtVBF:V1J:HighPtOfflineVBF:HighMJJ:LowMJJ:HighMJJLP:LowMJJLP");
}

void VBFVectorBoson::initVariables(std::vector<Jet> jets){
  mjj    = (jets.size()>=2 ?  (jets[0]+jets[1]).M() : 0.);
  detajj = (jets.size()>=2 ? fabs(jets[0].Eta()-jets[1].Eta()) : -99.);
  dphijj = (jets.size()>=2 ? jets[0].DeltaPhi(jets[1]) : -99.);
  jjpt   = (jets.size()>=2 ? (jets[0]+jets[1]).Pt() : 0.);
  leadj_gawidth    = (jets.size()>=1 ? ev.j_gawidth[jets[0].getJetIndex()] : -99);
  //leadj_c2_05      = (jets.size()>1 ? ev.j_c2_05[jets[0].getJetIndex()] : -99);
  subleadj_gawidth = (jets.size()>=2 ? ev.j_gawidth[jets[1].getJetIndex()] : -99);
  //subleadj_c1_05   = (jets.size()>2 ? ev.j_c1_05[jets[1].getJetIndex()] : -99);
  subleadj_c2_02   = (jets.size()>=2 ? ev.j_c2_02[jets[1].getJetIndex()] : -99);
  subleadj_pt      = (jets.size()>=2 ? ev.j_pt[jets[1].getJetIndex()] : -99);
  lead_qg          = (jets.size()>=1 ? ev.j_qg[jets[0].getJetIndex()] : -99);
  ystar       = -99;
  balance     = -99;
  relbpt      = -99;
  dphibjj     = -99;
  isotropy    = -99;
  circularity = -99;
  sphericity  = -99;
  aplanarity  = -99;
  C           = -99;
  D           = -99;
  jjetas      = -99;
  dphivj0     = -99; 
  dphivj1     = -99;
}

void VBFVectorBoson::fill(MiniEvent_t ev, TLorentzVector boson, std::vector<Jet> jets, std::vector<double> cplotwgts, TString c, std::map<TString, int> mults, std::vector<Particle> fakeACR, std::vector<Particle> tightACR){
  ht->fill("nvtx",   ev.nvtx,          cplotwgts,c);        

  //boson histos
  ht->fill("vpt",    boson.Pt(),       cplotwgts,c);
  ht->fill("vy",     boson.Rapidity(), cplotwgts,c);   
  ht->fill("sihih",  sihih,            cplotwgts,c);   
  ht->fill("r9",     r9,               cplotwgts,c);   
  ht->fill("hoe",    hoe,              cplotwgts,c);   
  ht->fill("chiso",  chiso,            cplotwgts,c);   
  ht->fill("mindrl", mindrl,           cplotwgts,c);   

  for(auto a : photons) {
    int idx = a.originalReference();
    ht->fill("allsihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("allMjjEB"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c); 
    else
      ht->fill2D("allMjjEE"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c); 
  }
  for(auto a : relaxedTightPhotons) {
    int idx = a.originalReference();
    ht->fill("relaxedTightsihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("relaxedTightMjjEB"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c); 
    else
      ht->fill2D("relaxedTightMjjEE"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c); 
  }
  for(auto a : tmpPhotons) {
    int idx = a.originalReference();
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("tmpQCDMjjEB"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c); 
    else
      ht->fill2D("tmpQCDMjjEE"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c); 
  }
  //bosons in CR and fakes
  for(auto a : fakeACR) {
    int idx = a.originalReference();
    ht->fill("fakesihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    ht->fill("fakechiso",   ev.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
    ht->fill("fakeneutiso", ev.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
    ht->fill("fakeaiso",    ev.gamma_photonIso[idx]         ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("looseMjjEB"  ,  ev.gamma_sieie[idx], mjj        ,cplotwgts,c);
    else
      ht->fill2D("looseMjjEE"  ,  ev.gamma_sieie[idx], mjj        ,cplotwgts,c);
  }
  for(auto a : tightACR) { 
    int idx = a.originalReference();
    ht->fill("tightsihih",   ev.gamma_sieie[idx]             ,cplotwgts,c);
    ht->fill("tightchiso",   ev.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
    ht->fill("tightneutiso", ev.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
    ht->fill("tightaiso",    ev.gamma_photonIso[idx]         ,cplotwgts,c);
    if (fabs(ev.gamma_eta[idx]) < 1.442)
      ht->fill2D("tightMjjEB"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c);
    else
      ht->fill2D("tightMjjEE"  ,   ev.gamma_sieie[idx], mjj        ,cplotwgts,c);
  }
  ht->fill("nloose",        fakeACR.size(),       cplotwgts,c);
  ht->fill("ntight",        tightACR.size(),      cplotwgts,c);
  ht->fill("nloosefake",    mults["loosefake"],   cplotwgts,c);
  ht->fill("ntightfake",    mults["tightfake"],   cplotwgts,c);
  ht->fill("nlooseprompt",  mults["looseprompt"], cplotwgts,c);
  ht->fill("ntightprompt",  mults["tightprompt"], cplotwgts,c);

  //jet histos
  centraleta = 9999; 
  forwardeta = -9999;
  for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
    TString jtype(ij==0?"lead":"sublead");
    ht->fill(jtype+"pt",       jets[ij].Pt(),        cplotwgts,c);          
    ht->fill(jtype+"pumva",    jets[ij].getPUMVA(),  cplotwgts,c);
    ht->fill(Form("drj%db",(int)ij+1),   jets[ij].DeltaR(boson),  cplotwgts,c);
    ht->fill("jet_c2_00", 	ev.j_c2_00[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_c2_02", 	ev.j_c2_02[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_c2_05",	ev.j_c2_05[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_zg", 		ev.j_zg[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_gaptd", 	ev.j_gaptd[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_gawidth",     ev.j_gawidth[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht->fill("jet_emf"     , ev.j_emf[jets[ij].getJetIndex()]               ,cplotwgts,c);
    ht->fill("jet_raw_pt" ,jets[ij].pt()/ev.j_rawsf[jets[ij].getJetIndex()]        , cplotwgts,c );
    ht->fill("jet_qg"   ,ev.j_qg[jets[ij].getJetIndex()]    , cplotwgts,c);
    centraleta=min(centraleta,float(fabs(jets[ij].Eta())));
    forwardeta=max(forwardeta,float(fabs(jets[ij].Eta())));
    //int nnnn = (int) ij;
    ht->fill2D("etaphi",    centraleta,jet_phi ,   cplotwgts,c);
    
  }
 
  
  if(jets.size() >= 2){
    jjetas = jets[0].Eta()*jets[1].Eta();
    //dphivj0 = fabs(jets[0].DeltaPhi(boson));
    dphivj1 = fabs(jets[1].DeltaPhi(boson));
    ht->fill("jjetas",  jjetas,   cplotwgts,c);
    ht->fill("dphivj0", dphivj0 ,  cplotwgts,c);
    ht->fill("dphivj1", dphivj1 ,  cplotwgts,c);
  }
  dphivj2 = 9999; dphivj3 = 9999;
  centjy = 9999; ncentjj = 0;  
  if(jets.size() > 2){
    dphivj2 = fabs(jets[2].DeltaPhi(boson));
    ht->fill("dphivj2", dphivj2 ,  cplotwgts,c);
    for(unsigned int iJet = 2; iJet < jets.size(); iJet++){	
      float dy = fabs(jets[0].Rapidity() - jets[1].Rapidity())/2;
      float sumy = (jets[0].Rapidity() + jets[1].Rapidity())/2;
      if(fabs(jets[iJet].Rapidity() - sumy) < dy){
        centjy =  jets[iJet].Rapidity();
        ht->fill("centjy",centjy,  cplotwgts,c);
        ncentjj++;
      }
    }
    ht->fill("ncentj", ncentjj, cplotwgts, c);
  }
  if(jets.size() > 3){
    dphivj3 =  fabs(jets[3].DeltaPhi(boson));
    ht->fill("dphivj3", dphivj3 , cplotwgts,c);
  }
  ht->fill("njets",        jets.size(), cplotwgts,c);
  ht->fill("ht",           scalarht,    cplotwgts,c);
  ht->fill("mht",          mht,         cplotwgts,c);
  ht->fill("centraleta",   centraleta,  cplotwgts,c);
  ht->fill("forwardeta",   forwardeta,  cplotwgts,c);
  ht->fill("dijetpt",      jjpt,        cplotwgts,c);
  ht->fill("detajj",       detajj,      cplotwgts,c);
  ht->fill("dphijj",       dphijj,      cplotwgts,c);
  ht->fill("relbpt",       relbpt,      cplotwgts,c);
  ht->fill("dphibjj",      dphibjj,     cplotwgts,c);
  ht->fill("mjj", 	   mjj,         cplotwgts,c);
	
  //visible system histos
  ht->fill("vystar",       ystar,              cplotwgts,c);        
  ht->fill("balance",      balance,            cplotwgts,c);
  ht->fill("isotropy",     isotropy,     cplotwgts,c);
  ht->fill("circularity",  circularity,  cplotwgts,c);
  ht->fill("sphericity",   sphericity, cplotwgts,c);
  ht->fill("aplanarity",   aplanarity, cplotwgts,c);
  ht->fill("C",            C,          cplotwgts,c);
  ht->fill("D",            D,          cplotwgts,c);

  //final analysis histograms
  ht->fill("evcount",    0, cplotwgts, c);
  if(!(doBlindAnalysis && vbfmva<-99))
    ht->fill("vbfmva",       vbfmva, cplotwgts,c);
  if(!(doBlindAnalysis && vbffisher<-99))
    ht->fill("vbffisher",       vbffisher, cplotwgts,c);

  if(skimtree) newTree->Fill();
}
