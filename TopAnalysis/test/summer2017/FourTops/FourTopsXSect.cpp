#include <stdio.h>
#include <TFile.h>
#include <string.h>

void FourTopsXSect()
{
    string dir("/afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTops/");

    string sigFileName("MC13TeV_TTTT.root");
    std::vector<string> bgFileName =
    {
        "MC13TeV_TTJets.root",
        "MC13TeV_SingleTbar_tW.root",
        "MC13TeV_SingleT_tW.root",
        "MC13TeV_SingleTbar_t.root",
        "MC13TeV_SingleT_t.root",
        "MC13TeV_DY50toInf_mlm.root",
        "MC13TeV_DY10to50.root",
        "MC13TeV_ZZTo2L2Nu.root",
        "MC13TeV_ZZTo2L2Q.root)",
        "MC13TeV_WWToLNuQQ.root",
        "MC13TeV_WWTo2L2Nu.root",
        "MC13TeV_WZTo3LNu.root",
        "MC13TeV_TTWToQQ.root",
        "MC13TeV_TTZToQQ.root",
        "MC13TeV_TTZToLLNuNu.root",
        "MC13TeV_tZq.root",
        "MC13TeV_ZZ4L.root",
        "MC13TeV_GGHZZ4L.root",
        "MC13TeV_VHtoNonbb.root",
        "MC13TeV_WWW.root",
        "MC13TeV_WWZ.root",
        "MC13TeV_ZZZ.root",
        "MC13TeV_TTWW.root",
        "MC13TeV_TTTW.root",
        "MC13TeV_TTWZ.root",
        "MC13TeV_TTZZ.root",
        "MC13TeV_TTWH.root",
        "MC13TeV_TTZH.root",
        "MC13TeV_TTHH.root",
        "MC13TeV_TTTJ.root"
    };

    std::vector<string> dataFileName = 
    {
        "Data13TeV_DoubleEG_2016B.root",
        "Data13TeV_DoubleEG_2016C.root",
        "Data13TeV_DoubleEG_2016D.root",
        "Data13TeV_DoubleEG_2016E.root",
        "Data13TeV_DoubleEG_2016F.root",
        "Data13TeV_DoubleEG_2016G.root",
        "Data13TeV_DoubleEG_2016Hv2.root",
        "Data13TeV_DoubleEG_2016Hv3.root",
        "Data13TeV_DoubleMuon_2016B.root",
        "Data13TeV_DoubleMuon_2016C.root",
        "Data13TeV_DoubleMuon_2016D.root",
        "Data13TeV_DoubleMuon_2016E.root",
        "Data13TeV_DoubleMuon_2016F.root",
        "Data13TeV_DoubleMuon_2016G.root",
        "Data13TeV_DoubleMuon_2016Hv2.root",
        "Data13TeV_DoubleMuon_2016Hv3.root",
        "Data13TeV_MuonEG_2016B.root",
        "Data13TeV_MuonEG_2016C.root",
        "Data13TeV_MuonEG_2016D.root",
        "Data13TeV_MuonEG_2016E.root",
        "Data13TeV_MuonEG_2016F.root",
        "Data13TeV_MuonEG_2016G.root",
        "Data13TeV_MuonEG_2016Hv2.root",
        "Data13TeV_MuonEG_2016Hv3.root",
        "Data13TeV_SingleElectron_2016B.root"
        "Data13TeV_SingleElectron_2016C.root",
        "Data13TeV_SingleElectron_2016D.root",
        "Data13TeV_SingleElectron_2016E.root",
        "Data13TeV_SingleElectron_2016F.root",
        "Data13TeV_SingleElectron_2016G.root",
        "Data13TeV_SingleElectron_2016Hv2.root",
        "Data13TeV_SingleElectron_2016Hv3.root",
        "Data13TeV_SingleMuon_2016B.root",
        "Data13TeV_SingleMuon_2016C.root",
        "Data13TeV_SingleMuon_2016D.root",
        "Data13TeV_SingleMuon_2016E.root",
        "Data13TeV_SingleMuon_2016F.root",
        "Data13TeV_SingleMuon_2016G.root",
        "Data13TeV_SingleMuon_2016Hv2.root",
        "Data13TeV_SingleMuon_2016Hv3.root"
    };

    TFile *sigFile = new TFile((dir+sigFileName).c_str());
    std::vector<TFile*> bgFilePtr();
    std::vector<TFile*> dataFilePtr();

    for (int i=0;i<bgFileName.size();i++) 
        bgFilePtr.push_back(new TFile((dir+bgFileName[i]).c_str()));

    for (int i=0;i<dataFileName.size();i++)
        dataFilePtr.push_back(new TFile((dir+dataFileName[i]).c_str()));
    
    std::vector<string> histograms = {"bdt","mlp_ann"};
    double factor[histograms.size()];

    TH1F *histCache;

    for (int i=0;i<histograms.size();i++)
    {
        histCache = sigFile->Get(histograms[i].c_str());

        int numBins = histCache->GetNbinsX();
        
        double y_sig [numBins];
        double y_bg  [numBins];
        double y_data[numBins];

        for (int j=0;j<numBins;j++)
        {
            y_sig[j] = histCache->GetBinContent(j+1);
            y_bg[j] = 0.;
            y_data[j] = 0.;
        }

        for (int k=0;k<bgFilePtr.size();k++) 
        {
            histCache = bgFilePtr[k]->Get(histograms[i].c_str());
            for (int j=0;j<numBins;j++) y_bg[j] += histCache->GetBinContent(j+1);
        }

        for (int k=0;k<dataFilePtr.size();k++)
        {
            histCache = dataFilePtr[k]->Get(histograms[i].c_str());
            for (int j=0;j<numBins;j++) y_data[j] += histCache->GetBinContent(j+1);
        }

        printf("%s",histograms[i]);
        printf("\n---------------------------------------\n");
        printf("SIGNAL\n");
        for (int i=0;i<numBins;i++) printf("%lf ",y_sig[i]);
        printf("\nBACKGROUND\n");
        for (int i=0;i<numBins;i++) printf("%lf ",y_bg[i]);
        printf("\nDATA\n");
        for (int i=0;i<numBins;i++) printf("%lf ",y_data[i]);
    }
}