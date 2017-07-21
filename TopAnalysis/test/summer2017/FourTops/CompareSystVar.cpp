#include <stdio.h>
#include <TFile.h>
#include <string.h>

void CompareSystVar()
{
    TFile *FourTops_nosys = new TFile("MC13TeV_TTTT.root");
    TFile *FourTops_fsrdn = new TFile("MC13TeV_TTTT_fsrdn.root");
    TFile *FourTops_fsrup = new TFile("MC13TeV_TTTT_fsrup.root");
    TFile *FourTops_isrdn = new TFile("MC13TeV_TTTT_isrdn.root");
    TFile *FourTops_isrup = new TFile("MC13TeV_TTTT_isrup.root");

    TFile *FourTops_addr[] = {FourTops_nosys, FourTops_fsrdn, FourTops_isrdn, FourTops_fsrup, FourTops_isrup};

    std::string fileName[] = {"FourTops_nosys", "FourTops_fsrdn", "FourTops_isrdn", "FourTops_fsrup", "FourTops_isrup"};
    std::string histograms[] = {"nbjets", "njets", "nleptons", "nvtx"};
    int markerStyle[] = {24, 20, 21, 22, 23};
    int markerColour[] = {1, 2, 3, 4, 5};

    TFile *out = new TFile("Compare_TTTT.root","RECREATE");
    TCanvas *canvas;
    TH1F *hist;
    TLegend *legend;

    for (short i = 0; i < 4; i++)
    {
        canvas = new TCanvas(histograms[i].c_str(),histograms[i].c_str(),800,600);
        canvas->cd();
        legend = new TLegend(0.5, 0.7, 0.9, 0.9);
        for (short j = 0; j < 5; j++)
        {
            hist = (TH1F*) FourTops_addr[j]->Get(histograms[i].c_str());
            hist->SetMarkerStyle(markerStyle[j]);
            hist->SetLineColor(markerColour[j]);
            hist->SetMarkerColor(markerColour[j]);
            hist->Draw("SAME");
            legend->AddEntry(hist, fileName[j], "lep");
        }
        legend->Draw();
        out->cd();
        canvas->Write();
    }
}
