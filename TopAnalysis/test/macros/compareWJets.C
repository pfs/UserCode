void compareWJets()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString baseDir("~/work/LJets2015-arcrev/");
  
  TCanvas *c1 = new TCanvas("c1","c1",500,500);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.02);
  c1->SetBottomMargin(0.1);
  TString dirs[]={"muplus","muminus","eplus","eminus"};
  TString plots[]={"metpt_1j0t", "minmlb_1j1t",
		   "metpt_2j0t", "minmlb_2j1t", "minmlb_2j2t",
		   "metpt_3j0t", "minmlb_3j1t", "minmlb_3j2t",
		   "metpt_4j0t", "minmlb_4j1t", "minmlb_4j2t"};
  for(size_t ip=0; ip<sizeof(plots)/sizeof(TString); ip++)
    {
      c1->Clear();

      TH1 *aMCatNLO=0, *MG=0;
      for(size_t id=0; id<sizeof(dirs)/sizeof(TString); id++)
	{
	  TFile *fIn=TFile::Open(baseDir+"/analysis_"+dirs[id]+"/MC13TeV_WJets.root");
	  if(aMCatNLO==0)
	    {
	      aMCatNLO=(TH1 *)fIn->Get(plots[ip])->Clone("amcatnlo");
	      aMCatNLO->SetDirectory(0);
	    }
	  else
	    aMCatNLO->Add( (TH1 *) fIn->Get(plots[ip]) );
	  fIn->Close();
       
	  fIn=TFile::Open(baseDir+"/analysis_"+dirs[id]+"/MC13TeV_WJets_madgraph.root");
	  if(MG==0)
	    {
	      MG=(TH1 *)fIn->Get(plots[ip])->Clone("mg");
	      MG->SetDirectory(0);
	    }
	  else
	    MG->Add( (TH1 *) fIn->Get(plots[ip]) );
	  fIn->Close();
	}
      
      MG->SetTitle("Madgraph MLM");
      MG->SetLineColor(kRed);
      MG->Draw("hist");
      MG->GetYaxis()->SetRangeUser(0,MG->GetMaximum()*1.8);
      aMCatNLO->SetTitle("aMC@NLO");
      aMCatNLO->SetLineColor(1);
      aMCatNLO->SetMarkerStyle(20);
      aMCatNLO->SetMarkerSize(1.0);
      aMCatNLO->Draw("e1same");

      TLegend *leg = new TLegend(0.15,0.8,0.4,0.95);
      leg->SetHeader("#bf{CMS} #it{simulation} #sqrt{s}=13 TeV");
      leg->SetFillStyle ( 0);
      leg->SetFillColor ( 0);
      leg->SetBorderSize( 0);
      leg->SetLineWidth(1);
      leg->SetTextSize(0.035);
      leg->Draw();
      leg->AddEntry(aMCatNLO,aMCatNLO->GetTitle(),"p");
      leg->AddEntry(MG,MG->GetTitle(),"l");

      Int_t nbins=MG->GetNbinsX();
      TLatex *txt=new TLatex();
      txt->SetTextFont(42);
      txt->SetTextSize(0.035);
      txt->SetNDC();
      txt->DrawLatex(0.6,0.91,Form("NLO/LO=%3.2f",aMCatNLO->Integral(0,nbins+1)/MG->Integral(0,nbins+1)));
      Double_t chi2;
      Int_t ndf,  igood;
      TH1 *normAMCatNLO=(TH1 *)aMCatNLO->Clone("normamcnlo");
      normAMCatNLO->Scale(1./aMCatNLO->Integral(0,nbins+1));
      TH1 *normMG=(TH1 *)MG->Clone("normmg");
      normMG->Scale(1./MG->Integral(0,nbins+1));
      normAMCatNLO->Chi2TestX(normMG,chi2,ndf,igood," WW");
      txt->DrawLatex(0.6,0.86,Form("#chi^{2}/ndf (shape)=%3.2f/%d",chi2,ndf));
      
      c1->Modified();
      c1->Update();
      c1->SaveAs("wjets_"+plots[ip]+".png");
      c1->SaveAs("wjets_"+plots[ip]+".pdf");
    }
}
