#include <stdio.h>
#include <stdlib.h>

void plotXSect()
{
    const int numberOfDatasets = 2;
    TFile *fileIndex[numberOfDatasets];
    fileIndex[0] = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    fileIndex[1] = new TFile("datacard_4tops/fourtops_bdt_Asymptotic.root");
    //TFile *file = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    TTree *tree;
    string fileName[] = {"#bf{Neural Network}", "#bf{BDT}"};

    TCanvas *canvas = new TCanvas("canvas","Upper limits for four tops",800,600);
    canvas->Draw();
    canvas->cd();
    TH1F *histFrame = gPad->DrawFrame(0,-0.15,80,4);
    histFrame->SetXTitle("95% CL Limit on #mu = #sigma_{obs} / #sigma_{SM}");
    histFrame->GetYaxis()->SetTickSize(0);
    histFrame->GetYaxis()->SetLabelSize(0);
    
    TLatex *tex = new TLatex();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    //tex->SetNDC();

    TBox *box1sig, *box2sig;
    TLine *lineexp, *lineobs;
    double limitVal, limitErrVal;
    double limit2sigdown, limit1sigdown, limit0, limit1sigup, limit2sigup;
    double limitObs[numberOfDatasets], limitObsError[numberOfDatasets], limitObsCache;
    double *limitAddr[] = {&limit2sigdown, &limit1sigdown, &limit0, &limit1sigup, &limit2sigup, &limitObsCache};

    double limitYPointCoord[numberOfDatasets], limitYPointCoordErr[numberOfDatasets];

    for (int s = 0; s < 2; s++)
    {
        tree = (TTree*) fileIndex[s]->Get("limit");
        tree->SetBranchAddress("limit", &limitVal);
        tree->SetBranchAddress("limitErr", &limitErrVal);
        for (int i=0;i<tree->GetEntries();i++)
        {
            tree->GetEntry(i);
            *limitAddr[i] = limitVal;
            //printf("limitVal = %lf\n",limitVal);
        }

        limitObs[s] = limitObsCache;
        tree->GetEntry(5);
        limitObsError[s] = limitErrVal;
        limitYPointCoord[s] = s + 0.5;
        limitYPointCoordErr[s] = 0.5;

        printf("limit2sigdown = %lf\n",limit2sigdown);
        printf("limit2sigup   = %lf\n",limit2sigup);
        printf("limit1sigdown = %lf\n",limit1sigdown);
        printf("limit1sigup   = %lf\n",limit1sigup);
        printf("limit0        = %lf\n",limit0);
        printf("limitObs      = %lf\n",limitObs[s]);

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
        //lineobs = new TLine(limitObs, s, limitObs, s+1);
        //lineobs->SetLineWidth((Width_t)2);
        //lineobs->Draw();

        //TMarker *marker = new TMarker(limitObs, s + 0.5, 21);
        //marker->Draw();

        tex->DrawLatex(50, s + 0.5, fileName[s].c_str());
    }

    TLine *lineSM = new TLine(1, -0.15, 1, 4);
    lineSM->SetLineColor(kRed);
    lineSM->SetLineWidth((Width_t)5);
    lineSM->Draw();

    TGraphErrors *obsPoints = new TGraphErrors(numberOfDatasets, limitObs, limitYPointCoord, limitObsError, limitYPointCoordErr);
    obsPoints->SetLineColor(1);
    obsPoints->SetLineWidth(2);
    obsPoints->SetMarkerColor(1);
    obsPoints->SetMarkerStyle(21);
    obsPoints->SetMarkerSize(1.6);
    obsPoints->SetMarkerColor(1);
    obsPoints->Draw("SAME P");

    tex->DrawLatex(0,4.1,"#bf{CMS} #it{Preliminary}");
    tex->DrawLatex(64,4.1,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
    tex->DrawLatex(1.5,3.8,"#scale[0.75]{#color[2]{#bf{SM Here}}}");

    auto hist1sigDummy = TH1F("hist1sigDummy","expected #pm 1 #sigma",1,0,1);
    auto hist2sigDummy = TH1F("hist2sigDummy","expected #pm 2 #sigma",1,0,1);
    hist1sigDummy->SetFillColor(15);
    hist1sigDummy->SetLineWidth(0);
    hist2sigDummy->SetFillColor(17);
    hist2sigDummy->SetLineWidth(0);

    TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetHeader("limits on #sigma_{tttt}");
    legend->AddEntry(obsPoints,"observed","lep");
    legend->AddEntry(hist1sigDummy,"expected #pm 1 #sigma","f");
    legend->AddEntry(hist2sigDummy,"expected #pm 2 #sigma","f");
    legend->Draw();
}
