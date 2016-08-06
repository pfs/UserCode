#include "CMS_lumi.C"
#include "tdrstyle.C"

TF1 *gnull = new TF1("gnull","gaus");
TF1 *galt  = new TF1("galt", "gaus");

Double_t ginter_fn(Double_t *x, Double_t *par) {
    return TMath::Abs(gnull->EvalPar(x,par) - galt->EvalPar(x,par));
}

RooStats::HypoTestResult *readLepFile(TDirectory *toyDir,  double rValue) {
    TString prefix = TString::Format("HypoTestResult_r%g_",rValue);
    RooStats::HypoTestResult *ret = 0;
    TIter next(toyDir->GetListOfKeys()); 
    TKey *k;
    while ((k = (TKey *) next()) != 0) {
        if (TString(k->GetName()).Index(prefix) != 0) continue;
        RooStats::HypoTestResult *toy = (RooStats::HypoTestResult *)(toyDir->Get(k->GetName()));
        if (toy == 0) continue;
        if (ret == 0) {
            ret = new RooStats::HypoTestResult(*toy);
        } else {
            ret->Append(toy);
        }
    }
    return ret;
}

void hypoTestResultTreeTopWid(TString fOutName,
                              double mass, 
                              double rValue=1.0,
                              const char *poiName="r", 
                              double numToys=1000,
                              const char *lfs="",
                              TString wid="1p0w", 
                              const char *dist="mlb", 
                              bool unblind = false) 
{
    if(wid == "1p0w") return;
    
    if (gROOT->GetListOfFiles()->GetSize() == 0) {
        std::cerr << "ERROR: you have to open at least one root file" << std::endl;
    }

    TFile *fOut = new TFile(fOutName, "RECREATE");
    TTree *tree = new TTree("q","Test statistics");

    float q, 
          mh = mass,
          r = rValue, 
          weight; 
    int   type;

    tree->Branch("q", &q, "q/F");
    tree->Branch("mh", &mh, "mh/F");
    tree->Branch("weight", &weight, "weight/F");
    tree->Branch("type", &type, "type/I");
    tree->Branch("r", &r, "r/F");

    TString prefix1 = TString::Format("HypoTestResult_mh%g_%s%g_",mass,poiName,rValue);
    TString prefix2 = TString::Format("HypoTestResult_%s%g_",poiName,rValue);
    long int nS = 0, 
             nB = 0;

    TCanvas *c = new TCanvas("","");

    for (int i = 0, n = gROOT->GetListOfFiles()->GetSize()-1;  i < n; ++i) {
        TDirectory *toyDir = ((TFile*) gROOT->GetListOfFiles()->At(i))->GetDirectory("toys");

        if (toyDir == 0) {
            std::cerr << "Error in file " 
                      << gROOT->GetListOfFiles()->At(i)->GetName() 
                      << ": directory /toys not found" << std::endl;
            continue;
        }

        TIter next(toyDir->GetListOfKeys()); TKey *k;
        while ((k = (TKey *) next()) != 0) {
            if (TString(k->GetName()).Index(prefix1) != 0 
                    && TString(k->GetName()).Index(prefix2) != 0) continue;

            RooStats::HypoTestResult *toy = dynamic_cast<RooStats::HypoTestResult *>(toyDir->Get(k->GetName()));

            if (toy == 0) continue;

            std::cout << " - " << k->GetName() << std::endl;
            RooStats::SamplingDistribution * bDistribution = toy->GetNullDistribution(), 
                                           * sDistribution = toy->GetAltDistribution();
            const std::vector<Double_t> & bdist   = bDistribution->GetSamplingDistribution();
            const std::vector<Double_t> & bweight = bDistribution->GetSampleWeights();
            for (int j = 0, nj = bdist.size(); j < nj; ++j) {
                q = bdist[j]; 
                weight = bweight[j]; 
                type = -1;
                tree->Fill(); 
                nB++;
            }

            const std::vector<Double_t> & sdist   = sDistribution->GetSamplingDistribution();
            const std::vector<Double_t> & sweight = sDistribution->GetSampleWeights();
            for (int j = 0, nj = sdist.size(); j < nj; ++j) {
                q = sdist[j]; 
                weight = sweight[j]; 
                type = 1;
                tree->Fill(); 
                nS++;
            }

            Double_t data = toy->GetTestStatisticData();
            weight = 1.0; 
            q = data; 
            type = 0;
            tree->Fill();
        }
    }

    tree->Write();

    /*
     * Storing plots
     */

    c->cd();
    setTDRStyle();

    std::cout << "Creating plots" << std::endl;
    TH1D *hnull = new TH1D("hnull","Null Hypothesis",100,-20,20);
    TH1D *halt  = new TH1D("halt" ,"Alternate Hypothesis",100,-20,20);
    TH1D *hnullstat = new TH1D("hnullstat","Null Hypothesis",5000,-1000,1000);
    TH1D *haltstat  = new TH1D("haltstat" ,"Alternate Hypothesis",5000,-1000,1000);

    tree->Draw("-2*q>>hnull","type>0","goff");
    tree->Draw("-2*q>>halt" ,"type<0","goff");

    hnull->SetLineColor(kBlue);
    hnull->SetStats(false);
    hnull->GetXaxis()->SetTitle("-2*ln(L_{alt}#/L_{null})");
    hnull->GetXaxis()->SetTitleSize(0.04);
    hnull->GetXaxis()->SetTitleOffset(1.2);
    hnull->GetYaxis()->SetTitle("Pseudoexperiments");
    halt->SetLineColor(kOrange);
    halt->SetStats(false);
    hnull->Draw();
    halt->Draw("SAME");

    c->BuildLegend(0.7,0.8,0.88,0.88);
    TString plotName = TString(lfs)+"_"+TString(wid)+"_"+TString(dist);

    c->SetTitle(plotName+" Toys");

    CMS_lumi(c,4,0); 
    c->SaveAs(plotName+".pdf");
    c->SaveAs(plotName+".png");

    /*
     * Outputting LaTeX table with useful statistics
     */
    
    tree->Draw("-2*q>>hnullstat","type>0","goff");
    tree->Draw("-2*q>>haltstat" ,"type<0","goff");

    //
    // Get x position where toy histograms approximately intersect
    // 
    gnull = new TF1("gnull","gaus",-1000,1000);
    galt  = new TF1("galt", "gaus",-1000,1000);
    hnullstat->Fit(gnull); 
     haltstat->Fit(galt); 

    TF1 *ginter = new TF1("",ginter_fn);
    Double_t xint = ginter->GetMinimumX();
    Int_t intbin  = gnull->GetXaxis()->FindBin(xint);

    std::cout << " - intersection x   is " << xint   << std::endl;
    std::cout << " - intersection bin is " << intbin << std::endl;

    //
    // get separation by integrating fits
    //
    std::cout << " - getting separation... " << std::endl;

    Double_t galtNorm = galt->Integral(-1000,1000);
    Double_t gnullNorm=gnull->Integral(-1000,1000);
    
    if(galtNorm == 0 or gnullNorm == 0) {
        std::cout<< " \n\nWARNING : INTEGRAL OF ONE FIT IS 0 \n\n " << std::endl;
    }

    Double_t separation = 0;
    if(hnullstat->GetMean() <= haltstat->GetMean()) {
       separation = galt->Integral(-1000,xint)/galtNorm
                 + gnull->Integral(xint,1000)/gnullNorm;
    } else {
       separation = gnull->Integral(-1000,xint)/gnullNorm 
                   + galt->Integral(xint,1000)/galtNorm;
    }

    separation = 1 - separation;

    //
    // get CLs from file (highly specific to our analysis)
    // 
    Double_t clsObs,clsbObs,clbObs,
             clsObsErr,clsbObsErr,clbObsErr;
    if(unblind) {
        std::cout << " - getting CLs... " << std::endl;
        TDirectory *toyDir = ((TFile*) gROOT->GetListOfFiles()->At(0))->GetDirectory("toys");
        if (toyDir == 0) {
            std::cerr << "Error in file " << gROOT->GetListOfFiles()->At(0)->GetName() 
                      << ": directory /toys not found" << std::endl;
        }
        std::cout << " - reading lep file... " << std::endl;
        RooStats::HypoTestResult *res = readLepFile(toyDir,rValue);
        std::cout << " - read lep file... " << std::endl;

        clsObs     = res->CLs(); 
        clsObsErr  = res->CLsError();
        clsbObs    = res->CLsplusb();
        clsbObsErr = res->CLsplusbError();
        clbObs     = res->CLb();
        clbObsErr  = res->CLbError();
        std::cout << " - got CLs... " << std::endl;
    }

    //
    // get fluctuations past median
    // 
    std::cout << " - getting density exceeding medians... " << std::endl;
    Double_t nullExceededDensity=0;
    Double_t  altExceededDensity=0;
    Double_t nullFitMedian = gnull->GetParameter("Mean"); 
    Double_t  altFitMedian =  galt->GetParameter("Mean"); 
    if(nullFitMedian < altFitMedian) {
      nullExceededDensity=gnull->Integral(altFitMedian,1000)/gnullNorm;
       altExceededDensity= galt->Integral(-1000,nullFitMedian)/galtNorm;
    } else {
      nullExceededDensity=gnull->Integral(-1000,altFitMedian)/gnullNorm;
       altExceededDensity= galt->Integral(nullFitMedian,1000)/galtNorm;
    }

    Double_t CLsExp = altExceededDensity/0.5;

    //
    // get quantiles for 1,2,3 sigma
    //
    std::cout << " - opening file... " << std::endl; 
    const int nquants = 7;
    Double_t quantpos[nquants] = { 0.0015, 0.023, 0.16, 0.50, 0.841, 0.977, 0.9985 };
    Double_t outq_nul[nquants];
    Double_t outq_alt[nquants];

    hnullstat->GetQuantiles(nquants,outq_nul,quantpos);
     haltstat->GetQuantiles(nquants,outq_alt,quantpos);

    //
    // store the information in a nice text file
    //
    std::ofstream ofs(TString(TString("stats__")+wid+
                                TString("_")+lfs+
                                TString("_")+dist+
                                TString(".txt")).Data(), 
                      std::ofstream::out);
    ofs << TMath::ErfInverse(separation)          << "$\\sigma$" << " # separation \n"
        << TMath::ErfInverse(nullExceededDensity) << "$\\sigma$" << " # null exceeded density \n"
        << TMath::ErfInverse(altExceededDensity)  << "$\\sigma$" << " # alt exceeded density\n"
        << xint                << " #  xint\n"
        << separation          << " #  sep\n"
        << altExceededDensity  << " #  aed\n"
        << nullExceededDensity << " #  ned\n"
        << CLsExp << " # cls expected \n";

    if(unblind) {
        ofs << clsObs  << " \\pm " << clsObsErr  << " # cls observed \n"
            << clbObs  << " \\pm " << clbObsErr  << " # clb observed \n"
            << clsbObs << " \\pm " << clsbObsErr << " # clsb observed \n";
    }

    ofs << "nulquant;"; 
        for(int i = 0; i < nquants; i++) {
            ofs << outq_nul[i] << ";";
        }
    ofs << "\naltquant;";
        for(int i = 0; i < nquants; i++) {
            ofs << outq_alt[i] << ";";
        }
    ofs.close();

    /* 
     * Cleanup
     */
    fOut->Close();
    std::cout << "Saved test statistics distributions for " 
              << nS << " signal toys and " 
              << nB << " background toys to " 
              << fOutName << "." << std::endl;
}
