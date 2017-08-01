#include <stdio.h>
#include <stdlib.h>

void plotXSect()
{
    TFile *fileIndex[2];
    fileIndex[0] = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    fileIndex[1] = new TFile("datacard_4tops/fourtops_bdt_Asymptotic.root");
    //TFile *file = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    TTree *tree;

    TCanvas *canvas = new TCanvas("canvas","Upper limits for four tops",800,600);
    canvas->Draw();
    canvas->cd();
    TH1F *histFrame = gPad->DrawFrame(0,0,80,4);
    histFrame->SetXTitle("95% CL Limit on #mu = #sigma_{obs} / #sigma_{SM}");
    histFrame->GetYaxis()->SetTickSize(0);
    histFrame->GetYaxis()->SetLabelSize(0);

    TBox *box1sig, *box2sig;
    TLine *lineexp, *lineobs;
    double limitVal;
    double limit2sigdown, limit1sigdown, limit0, limit1sigup, limit2sigup;
    double limitobs;
    double *limitAddr[] = {&limit2sigdown, &limit1sigdown, &limit0, &limit1sigup, &limit2sigup, &limitobs};

    for (int s = 0; s < 2; s++)
    {
        tree = (TTree*) fileIndex[s]->Get("limit");
        tree->SetBranchAddress("limit",&limitVal);
        for (int i=0;i<tree->GetEntries();i++)
        {
            tree->GetEntry(i);
            *limitAddr[i] = limitVal;
            //printf("limitVal = %lf\n",limitVal);
        }

        printf("limit2sigdown = %lf\n",limit2sigdown);
        printf("limit2sigup   = %lf\n",limit2sigup);
        printf("limit1sigdown = %lf\n",limit1sigdown);
        printf("limit1sigup   = %lf\n",limit1sigup);
        printf("limit0        = %lf\n",limit0);
        printf("limitobs      = %lf\n",limitobs);

        canvas->cd();
        box2sig = new TBox(limit2sigdown, s, limit2sigup, s+1);
        box2sig->SetFillColor(17);
        box2sig->Draw();
        box1sig = new TBox(limit1sigdown, s, limit1sigup, s+1);
        box1sig->SetFillColor(15);
        box1sig->Draw();

        //lineexp = new TLine(limit0, s, limit0, s+1);
        //lineexp->Draw();
        //lineexp->SetLineWidth((Width_t)2);
        lineobs = new TLine(limitobs, s, limitobs, s+1);
        lineobs->SetLineStyle(7);
        lineobs->SetLineWidth((Width_t)2);
        lineobs->Draw();

        TMarker *marker = new TMarker(limit0, s + 0.5, 21);
        marker->Draw();
    }

    TLine *lineSM = new TLine(1, 0, 1, 4);
    lineSM->SetLineColor(kRed);
    lineSM->SetLineWidth((Width_t)5);
    lineSM->Draw();

    TLatex *tex = new TLatex();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetNDC();
    tex->DrawLatex(0.15,0.9,"#bf{CMS} #it{Preliminary}");
    tex->DrawLatex(0.75,0.97,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
}
