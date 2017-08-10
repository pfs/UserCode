#include <stdio.h>
#include <stdlib.h>

void plotXSect()
{
    const int numberOfDatasets = 6;
    TFile *fileIndex[numberOfDatasets];
    fileIndex[0] = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    fileIndex[1] = new TFile("datacard_4tops/fourtops_bdt_Asymptotic.root");
    //TFile *file = new TFile("datacard_4tops/fourtops_nn_Asymptotic.root");
    TTree *tree;
    
    // Dataset name here
    string fileName[] = {"#splitline{#bf{Four lepton NN}}{#scale[0.7]{This analysis}}", "#splitline{#bf{Four lepton BDT}}{#scale[0.7]{This analysis}}", "#splitline{#bf{SS Dilepton}}{#scale[0.7]{EPJC 76 (2016) 439}}", "#splitline{#bf{OS Dilepton}}{#scale[0.7]{TOP-16-016}}",  "#splitline{#bf{Single lepton}}{#scale[0.7]{TOP-16-016}}", "#splitline{#bf{Combined}}{#scale[0.4]{TOP-16-016 + EPJC 76 (2016) 439}}"};

    TCanvas *canvas = new TCanvas("canvas","Upper limits for four tops",800,600);
    canvas->Draw();
    canvas->cd();
    TH1F *histFrame = gPad->DrawFrame(0,-0.15,100,8);
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
    double limit2sigdown, limit1sigdown, limit0Cache, limit1sigup, limit2sigup;
    double limitObs[numberOfDatasets], limitObsError[numberOfDatasets], limitObsCache;
    double limit0[numberOfDatasets];
    double *limitAddr[] = {&limit2sigdown, &limit1sigdown, &limit0Cache, &limit1sigup, &limit2sigup, &limitObsCache};

    double limitYPointCoord[numberOfDatasets], limitYPointCoordErr[numberOfDatasets];

    // Hardcoding values for plots
    limitObs[2] = 12.9; // SS Dilepton
    limitObs[3] = 14.5; // OS Dilepton
    limitObs[4] = 17.2; // Single lepton
    limitObs[5] = 7.4; // Combined
    limitObsError[2] = 0.;
    limitObsError[3] = 0.;
    limitObsError[4] = 0.;
    limitObsError[5] = 0.;
    limit0[2] = 11.0;
    limit0[3] = 24.7;
    limit0[4] = 16.4;
    limit0[5] = 7.7;

    for (int y = 2; y < 6; y++)
    {
        limitYPointCoord[y] = 5 - y + 0.5;
        limitYPointCoordErr[y] = 0.5;
    }

    // Read data from file to plot
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

        limit0[s] = limit0Cache;
        limitObs[s] = limitObsCache;
        tree->GetEntry(5);
        limitObsError[s] = limitErrVal;
        limitYPointCoord[s] = (5.5-s) + 0.5;
        limitYPointCoordErr[s] = 0.5;

        /*
        printf("limit2sigdown = %lf\n",limit2sigdown);
        printf("limit2sigup   = %lf\n",limit2sigup);
        printf("limit1sigdown = %lf\n",limit1sigdown);
        printf("limit1sigup   = %lf\n",limit1sigup);
        printf("limit0        = %lf\n",limit0);
        printf("limitObs      = %lf\n",limitObs[s]);
        */

        // Draw 2-sig box first ...
        canvas->cd();
        box2sig = new TBox(limit2sigdown, (5.5-s), limit2sigup, (5.5-s)+1);
        box2sig->SetFillColor(17);
        box2sig->Draw();
        // ... and 1-sig box on top of it.
        box1sig = new TBox(limit1sigdown, (5.5-s), limit1sigup, (5.5-s)+1);
        box1sig->SetFillColor(15);
        box1sig->Draw();

        // Dataset name
        tex->DrawLatex(70, (5.5-s) + 0.35, fileName[s].c_str());
    }

    // Hardcode the values for plotting. You can also use arrays to do the job, but I'm too lazy here.
    for (int s = 2; s < 6; s++)
    {
        switch (s)
        {
            case 2: // SS Dilepton
                limit2sigdown = 7.2;
                limit1sigdown = 7.2;
                limit1sigup = 17.2;
                limit2sigup = 17.2;
                break;
            case 3: // OS Dilepton
                limit2sigdown = 10.7;
                limit1sigdown = 15.5;
                limit1sigup = 41.4;
                limit2sigup = 68.5;
                break;
            case 4: // Single lepton
                limit2sigdown = 7.6;
                limit1sigdown = 10.7;
                limit1sigup = 26.2;
                limit2sigup = 42.8;
                break;
            case 5: // Combined
                limit2sigdown = 3.7;
                limit1sigdown = 5.1;
                limit1sigup = 11.8;
                limit2sigup = 17.4;
                break;
        }

        // Again, 2-sig box
        canvas->cd();
        box2sig = new TBox(limit2sigdown, 5-s, limit2sigup, (5-s)+1);
        box2sig->SetFillColor(17);
        box2sig->Draw();
        // 1-sig box
        box1sig = new TBox(limit1sigdown, 5-s, limit1sigup, (5-s)+1);
        box1sig->SetFillColor(15);
        box1sig->Draw();

        // Dataset name
        tex->DrawLatex(70, (5-s) + 0.35, fileName[s].c_str());
    }

    // Separate between our analysis and TOP-16-016
    TLine *lineSep = new TLine(0, 4.5, 100, 4.5);
    lineSep->Draw();

    // Line indicating where SM is
    TLine *lineSM = new TLine(1, -0.15, 1, 8);
    lineSM->SetLineColor(kRed);
    lineSM->SetLineWidth((Width_t)5);
    lineSM->Draw();

    // Put in graphs
    TGraphErrors *obsPoints = new TGraphErrors(numberOfDatasets, limitObs, limitYPointCoord, limitObsError, limitYPointCoordErr);
    obsPoints->SetLineColor(1);
    obsPoints->SetLineWidth(2);
    obsPoints->SetMarkerColor(1);
    obsPoints->SetMarkerStyle(21);
    obsPoints->SetMarkerSize(1.6);
    obsPoints->SetMarkerColor(1);
    obsPoints->Draw("SAME P");
    TGraphErrors *expPoints = new TGraphErrors(numberOfDatasets, limit0, limitYPointCoord, limitObsError, limitYPointCoordErr);
    expPoints->SetLineColor(1);
    expPoints->SetLineWidth(1);
    expPoints->SetLineStyle(7);
    expPoints->SetMarkerColor(1);
    expPoints->SetMarkerStyle(0);
    expPoints->SetMarkerSize(1.6);
    expPoints->SetMarkerColor(1);
    expPoints->Draw("SAME P");

    // A couple of texts
    tex->DrawLatex(0,8.1,"#bf{CMS} #it{Work in progress}");
    tex->DrawLatex(70,8.1,"#scale[0.8]{35.9 fb^{-1} (13 TeV)}");
    tex->DrawLatex(70,4.15,"#scale[0.8]{2.6 fb^{-1} (13 TeV)}");
    tex->DrawLatex(2,7.6,"#scale[0.75]{#color[2]{#bf{SM Here}}}");

    // Dummy TH1F objects for TLegend
    auto hist1sigDummy = new TH1F("hist1sigDummy","expected #pm 1 #sigma",1,0,1);
    auto hist2sigDummy = new TH1F("hist2sigDummy","expected #pm 2 #sigma",1,0,1);
    hist1sigDummy->SetFillColor(15);
    hist1sigDummy->SetLineColor(0);
    hist2sigDummy->SetFillColor(17);
    hist2sigDummy->SetLineColor(0);

    // Legend
    TLegend *legend = new TLegend(0.65,0.75,0.85,0.88);
    legend->SetLineColor(0);
    legend->SetHeader("limits on #sigma_{tttt}");
    legend->AddEntry(obsPoints,"observed","lep");
    legend->AddEntry(expPoints,"expected","l");
    legend->AddEntry(hist1sigDummy,"expected #pm 1 #sigma","f");
    legend->AddEntry(hist2sigDummy,"expected #pm 2 #sigma","f");
    legend->Draw();
}
