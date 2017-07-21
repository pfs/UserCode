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

    std::string histograms[] = {"nbjets", "njets", "nleptons", "nvtx"};

    TFile *out = new TFile("Compare_TTTT.root","RECREATE");
    TCanvas *canvas;
    TH1F *hist;

    for (short i = 0; i < 4; i++)
    {
        canvas = new TCanvas(histograms[i].c_str(),histograms[i].c_str(),800,600);
        canvas->cd();
        for (short j = 0; j < 5; j++)
        {
            hist = (TH1F*) FourTops_addr[j]->Get(histograms[i].c_str());
            hist->Draw("SAME");
        }
        out->cd();
        canvas->Write();
    }
}
