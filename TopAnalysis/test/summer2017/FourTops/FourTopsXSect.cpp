#include <stdio.h>
#include <TFile.h>
#include <string.h>

void FourTopsXSect()
{
    string dir("/afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTops/");

    TFile *plotterFile = new TFile("plots/plotter.root");
    string sigProcess = "t#bar{t}t#bar{t} x 100";
    std::vector<string> bgProcess =
    {
        "t#bar{t}",
        "Single top",
        "W",
        "DY",
        "Multiboson",
        "t#bar{t}+V",
        "Rares"
    };

    std::vector<string> histograms = {"bdt","mlp_ann"};
    double factor[histograms.size()];

    TH1F *histCache;

    for (int i=0;i<histograms.size();i++)
    {
        //printf("Now working on %s\n",histograms[i].c_str());
        histCache = (TH1F*) plotterFile->Get((histograms[i]+"/"+histograms[i]+"_"+sigProcess).c_str());

        int numBins = histCache->GetNbinsX();
        
        double y_sig [numBins];
        double y_bg  [numBins];
        double y_data[numBins];
        double y_data_nobg[numBins];

        for (int j=0;j<numBins;j++)
        {
            y_sig[j] = histCache->GetBinContent(j+1);
            y_bg[j] = 0.;
            y_data[j] = 0.;
        }
        //printf("y_sig filled\n");

        for (int k=0;k<bgProcess.size();k++) 
        {
            histCache = (TH1F*) plotterFile->Get((histograms[i]+"/"+histograms[i]+"_"+bgProcess[k]).c_str());
            //printf("histograms[i] obtained\n");
            for (int j=0;j<numBins;j++) y_bg[j] += histCache->GetBinContent(j+1);
        }
        //printf("y_bg filled\n");

        histCache = (TH1F*) plotterFile->Get((histograms[i]+"/"+histograms[i]).c_str());
        //printf("histograms[i] obtained\n");
        for (int j=0;j<numBins;j++) y_data[j] += histCache->GetBinContent(j+1);
        //printf("y_data filled\n");

        printf("%s",histograms[i].c_str());
        printf("\n---------------------------------------\n");
        printf("SIGNAL\n");
        for (int i=0;i<numBins;i++) printf("%lf ",y_sig[i]);
        printf("\nBACKGROUND\n");
        for (int i=0;i<numBins;i++) printf("%lf ",y_bg[i]);
        printf("\nDATA\n");
        for (int i=0;i<numBins;i++) printf("%lf ",y_data[i]);

        for (int j=0;j<numBins;j++) 
        {
            y_data_nobg[j] = y_data[j]-y_bg[j];
            if (y_data_nobg[j] < 0) y_data_nobg[j] = 0;
        }

        double sig_sq = 0.;
        double sig_data = 0.;

        for (int j=0;j<numBins;j++) 
        {
            sig_sq += y_sig[j]*y_sig[j];
            sig_data += y_sig[j]*y_data[j];
        }

        factor[i] = sig_data/sig_sq;

    }

    printf("\n==================RESULTS=====================\n");
    for (int i=0;i<histograms.size();i++) printf("%s\t%lf\n",histograms[i].c_str(),factor[i]);
}