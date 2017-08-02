#include <stdio.h>
#include <stdlib.h>

void plotXSect()
{
    const int numberOfDatasets = 5;
    TFile *fileIndex[numberOfDatasets];
    fileIndex[0] = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    fileIndex[1] = new TFile("datacard_4tops/fourtops_bdt_Asymptotic.root");
    //TFile *file = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    TTree *tree;
    string fileName[] = {"#splitline{#bf{Four lepton NN}}{#scale[0.7]{This analysis}}", "#splitline{#bf{Four lepton BDT}}{#scale[0.7]{This analysis}}", "#splitline{#bf{Dilepton}}{#scale[0.7]{TOP-16-016}}",  "#splitline{#bf{Single lepton}}{#scale[0.7]{TOP-16-016}}", "#splitline{#bf{Combined}}{#scale[0.7]{TOP-16-016}}"};

    TCanvas *canvas = new TCanvas("canvas","Upper limits for four tops",800,600);
    canvas->Draw();
    canvas->cd();
    TH1F *histFrame = gPad->DrawFrame(0,-0.15,100,7);
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

    limitObs[2] = 16.1;
    limitObs[3] = 14.9;
    limitObs[4] = 10.2;
    limitObsError[2] = 0.;
    limitObsError[3] = 0.;
    limitObsError[4] = 0.;

    for (int y = 2; y < 5; y++)
    {
        limitYPointCoord[y] = y + 0.5;
        limitYPointCoordErr[y] = 0.5;
    }

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

        tex->DrawLatex(70, s + 0.5, fileName[s].c_str());
    }

    for (int s = 2; s < 5; s++)
    {
        switch (s)
        {
            case 2: // Dilepton
                limit2sigdown = 9.4;
                limit1sigdown = 13.9;
                limit1sigup = 38.5;
                limit2sigup = 64.0;
                break;
            case 3: // Single lepton
                limit2sigdown = 5.9;
                limit1sigdown = 8.3;
                limit1sigup = 20.5;
                limit2sigup = 32.8;
                break;
            case 4: // Combined
                limit2sigdown = 5.0;
                limit1sigdown = 7.0;
                limit1sigup = 17.5;
                limit2sigup = 28.3;
                break;
        }

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

        tex->DrawLatex(70, s + 0.5, fileName[s].c_str());
    }

    TLine *lineSM = new TLine(1, -0.15, 1, 7);
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

    tex->DrawLatex(0,7.1,"#bf{CMS} #it{Work in progress}");
    tex->DrawLatex(80,7.1,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
    tex->DrawLatex(2,6.6,"#scale[0.75]{#color[2]{#bf{SM Here}}}");

    auto hist1sigDummy = new TH1F("hist1sigDummy","expected #pm 1 #sigma",1,0,1);
    auto hist2sigDummy = new TH1F("hist2sigDummy","expected #pm 2 #sigma",1,0,1);
    hist1sigDummy->SetFillColor(15);
    hist1sigDummy->SetLineColor(0);
    hist2sigDummy->SetFillColor(17);
    hist2sigDummy->SetLineColor(0);

    TLegend *legend = new TLegend(0.65,0.7,0.85,0.85);
    legend->SetLineColor(0);
    legend->SetHeader("limits on #sigma_{tttt}");
    legend->AddEntry(obsPoints,"observed","lep");
    legend->AddEntry(hist1sigDummy,"expected #pm 1 #sigma","f");
    legend->AddEntry(hist2sigDummy,"expected #pm 2 #sigma","f");
    legend->Draw();
}
