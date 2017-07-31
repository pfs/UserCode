#include <stdio.h>
#include <stdlib.h>

void plotXSect()
{
    TFile *file = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    TTree *tree = (TTree*) file->Get("limit");

    TCanvas canvas("canvas","Upper limits for four tops",800,600);
    canvas.cd();
    TH1F *histFrame = gPad->DrawFrame(0,0,4,50);
    histFrame->SetXTitle("95% CL Limit on #mu = #sigma_{obs} / #sigma_{SM}");
    histFrame->GetYaxis()->SetTickSize(0);
    histFrame->GetYaxis()->SetLabelSize(0);

    TBox *box1sig, *box2sig;
    TLine *lineexp, *lineobs;
    double limitVal;
    double limit2sigdown, limit1sigdown, limit0, limit1sigup, limit2sigup;
    double limitobs;
    double *limitAddr[] = {&limit2sigdown, &limit1sigdown, &limit0, &limit1sigup, &limit2sigup, &limitobs};

    tree->SetBranchAddress("limit",&limitVal);
    for (int i=0;i<tree->GetEvents();i++)
    {
        tree->GetEvent(i);
        *limitAddr[0] = limitVal;
    }

    box2sig = new TBox(limit2sigdown, 0, limit2sigup, 1);
    box2sig->SetFillColor(17);
    box2sig->Draw();
    box1sig = new TBox(limit1sigdown, 0, limit1sigup, 1);
    box1sig->SetFillColor(15);
    box1sig->Draw();

    lineexp = new TLine(limit0,0,limit0,1);
    lineexp.Draw();
    lineobs = new TLine(limitobs,0,limitobs,1);
    lineobs.SetLineStyle(7);
    lineobs.Draw();

    TLine lineSM(1,0,1,4);
    lineSM.SetLineColor(kRed);
    lineSM.SetLineWidth(5);
    lineSM.Draw();
}