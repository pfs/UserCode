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
        //if (TString(k->GetName()).Index(prefix) != 0) continue;
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
                              bool unblind = false,
                              TString prepost="") 
{
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

    TCanvas *c = new TCanvas("","",500,500);

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
     * Outputting LaTeX table with useful statistics
     */
    double maxToyValue = -2*tree->GetMinimum("q");
    double minToyValue = -2*tree->GetMaximum("q");
           maxToyValue = maxToyValue + 0.1*TMath::Abs(maxToyValue-minToyValue);
           minToyValue = minToyValue - 0.1*TMath::Abs(maxToyValue-minToyValue);

    std::cout<<"Max is: "<<maxToyValue<<" and min is "<<minToyValue<<std::endl;

    TH1D *hnullstat = new TH1D("hnullstat","Null Hypothesis",
                                    maxToyValue-minToyValue,
                                    minToyValue,maxToyValue);
    TH1D *haltstat  = new TH1D("haltstat" ,"Alternate Hypothesis",
                                    maxToyValue-minToyValue,
                                    minToyValue,maxToyValue);
    
    tree->Draw("-2*q>>hnullstat","type>0","goff");
    tree->Draw("-2*q>>haltstat" ,"type<0","goff");


    //
    // Get x position where toy histograms approximately intersect
    // 
    gnull = new TF1("gnull","gaus",minToyValue,maxToyValue);
    galt  = new TF1("galt", "gaus",minToyValue,maxToyValue);
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

    Double_t galtNorm = galt->Integral(minToyValue,maxToyValue);
    Double_t gnullNorm=gnull->Integral(minToyValue,maxToyValue);
    
    if(galtNorm == 0 or gnullNorm == 0) {
        std::cout<< " \n\nWARNING : INTEGRAL OF ONE FIT IS 0 \n\n " << std::endl;
    }

    Double_t separation = 0;
    if(hnullstat->GetMean() <= haltstat->GetMean()) {
       separation = galt->Integral(minToyValue,xint)/galtNorm
                 + gnull->Integral(xint,maxToyValue)/gnullNorm;
    } else {
       separation = gnull->Integral(minToyValue,xint)/gnullNorm 
                   + galt->Integral(xint,maxToyValue)/galtNorm;
    }

    separation = 1 - separation;

    //
    // get CLs from file (highly specific to our analysis)
    // 
    Double_t clsObs,clsbObs,clbObs,
             clsObsErr,clsbObsErr,clbObsErr,qobs;
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
        qobs       = res->GetTestStatisticData();
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
      nullExceededDensity=gnull->Integral(altFitMedian,maxToyValue)/gnullNorm;
       altExceededDensity= galt->Integral(minToyValue,nullFitMedian)/galtNorm;
    } else {
      nullExceededDensity=gnull->Integral(minToyValue,altFitMedian)/gnullNorm;
       altExceededDensity= galt->Integral(nullFitMedian,maxToyValue)/galtNorm;
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
    std::ofstream ofs(TString(prepost+TString("stats__")+wid+
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
            << clsbObs << " \\pm " << clsbObsErr << " # clsb observed \n"
            << "qobs;" << qobs << "\n";

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
     * Storing plots
     */

    c->cd();
    setTDRStyle();

    std::cout << "Creating plots" << std::endl;

    TH1D *hnull = new TH1D("hnull","Null Hypothesis",
                            TMath::CeilNint((maxToyValue-minToyValue)/4),
                            minToyValue,maxToyValue);
    TH1D *halt  = new TH1D("halt" ,"Alternate Hypothesis",
                            TMath::CeilNint((maxToyValue-minToyValue)/4),
                            minToyValue,maxToyValue);

    tree->Draw("-2*q>>hnull","type>0","goff");
    tree->Draw("-2*q>>halt" ,"type<0","goff");


    hnull->SetLineColor(kBlue-7);
    hnull->SetFillColor(kBlue-7);
    hnull->SetStats(false);
    hnull->GetXaxis()->SetTitle("-2 ln [ L(alt) #/ L(SM) ]");
    hnull->GetXaxis()->SetTitleSize(0.04);
    hnull->GetXaxis()->SetTitleOffset(1.2);
    hnull->GetYaxis()->SetTitle("Pseudoexperiments");
    hnull->GetYaxis()->SetTitle("Pseudoexperiments");
    hnull->GetYaxis()->SetTitleSize(0.04);
    hnull->GetYaxis()->SetTitleOffset(1.3);
    hnull->GetYaxis()->SetLabelSize(0.03);
    hnull->GetXaxis()->SetLabelSize(0.03);
    hnull->SetMaximum(hnull->GetMaximum()*1.3);
    halt->SetLineColor(kOrange);
    halt->SetFillColor(kOrange);
    halt->SetStats(false);
    hnull->Draw();
    halt->Draw("SAME");

    TLegend *leg = new TLegend(0.60,0.70,0.88,0.88);
    leg->AddEntry(hnull,TString("#Gamma_{SM} Hypothesis"),"F");
    leg->AddEntry(halt,TString(wid).ReplaceAll("p",".").ReplaceAll("w","")+TString("#times#Gamma_{SM} Hypothesis"),"F");
    leg->SetFillColor(kNone);
    leg->SetLineColor(kNone);
    leg->SetShadowColor(kNone);
    leg->Draw();

    TString plotName = TString(lfs)+"_"+TString(wid)+"_"+TString(dist);
    c->SetTitle(plotName+" Toys");

    relPosX=0.225;
    CMS_lumi(c,4,0); 
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SaveAs(prepost+plotName+".pdf");
    c->SaveAs(prepost+plotName+".png");

    /* 
     * Cleanup
     */
    fOut->Close();
    std::cout << "Saved test statistics distributions for " 
              << nS << " signal toys and " 
              << nB << " background toys to " 
              << fOutName << "." << std::endl;
}
