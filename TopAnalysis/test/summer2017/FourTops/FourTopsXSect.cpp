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
    double factor_stat_err_pos[histograms.size()];
    double factor_stat_err_neg[histograms.size()];

    TH1F *histCache;

    double SDGoal = 1.; // 68% CL
    SDGoal = 1.645; // 90% CL
    SDGoal = 1.96; // 95% CL
    SDGoal = 2.576; // 99% CL

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

        double sum_sig = 0.;
        double sum_data = 0.;
        double sum_data_log_sig = 0.;

        for (int j=0;j<numBins;j++)
        {
            if (y_sig[j] <= 0) continue;
            sum_sig += y_sig[j];
            sum_data += y_data_nobg[j];
            sum_data_log_sig += y_data_nobg[j]*TMath::Log(y_sig[j]);
        }

        factor[i] = sum_data / sum_sig;
        double lowest_F = -TMath::Log(factor[i])*sum_data + factor[i]*sum_sig - sum_data_log_sig;

        double x1 = factor[i]+0.1;
        double error = 0.;

	printf("sum_sig: %lf\n",sum_sig);	
	printf("sum_data: %lf\n",sum_data);	
	printf("sum_data_log_sig: %lf\n",sum_data_log_sig);	
	printf("factor[%d]: %lf\n",i,factor[i]);
	printf("lowest_F: %lf\n",lowest_F);

        do
        {
            error = -TMath::Log(x1)*sum_data + x1*sum_sig - sum_data_log_sig - lowest_F - SDGoal/2.;
            if (TMath::Abs(error) < 1e-6) break;
            x1 = factor[i] + (x1 - factor[i])/(error+SDGoal/2.)*(SDGoal/2.);
            printf("\nfactor_stat_err_pos: %lf\terror: %lf",x1,error);
	    getchar();
        } while (TMath::Abs(error) >= 1e-6);

        factor_stat_err_pos[i] = x1;

        x1 = 0.0001;
        double F_at_zero = -TMath::Log(x1)*sum_data + x1*sum_sig - sum_data_log_sig;
        x1 = 0.0001+(factor[i]-0.0001)/(lowest_F-F_at_zero)*(-F_at_zero);

        do
        {
            error = -TMath::Log(x1)*sum_data + x1*sum_sig - sum_data_log_sig - lowest_F - SDGoal/2;
            if (TMath::Abs(error) < 1e-6) break;
            x1 = 0.0001 + (x1-0.0001)/(error-F_at_zero+lowest_F+SDGoal/2.)*(-F_at_zero+lowest_F+SDGoal/2.);
            printf("\nfactor_stat_err_pos: %lf\terror: %lf",x1,error);
	    getchar(); 
        } while (TMath::Abs(error) >= 1e-6);

        factor_stat_err_neg[i] = x1;
    }

    printf("\n==================RESULTS=====================\n");
    for (int i=0;i<histograms.size();i++) printf("%s\t%lf\t+ %lf\t- %lf\n",histograms[i].c_str(),factor[i],factor_stat_err_pos[i],factor_stat_err_neg[i]);
}
