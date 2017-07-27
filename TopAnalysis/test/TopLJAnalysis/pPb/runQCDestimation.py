import ROOT
import sys
import pickle

from runDataFit import lumi

QCDNORM={}
PROCNORM={('ttbar','e'):59./200000.,
          ('ttbar','mu'):59./200000.,
          ('wjets','e'):1970./1.99183e+06,
          ('wjets','mu'):1970./3.93920e+06,
          ('dy','mu'):224./951411.,
          ('dy','e'):224./1.00000e+06,
          }

"""
"""
def getHistos(args,var):

    #get templates from file
    histos={'qcd':{}}
    for arg in args:
        proc,url=arg.split(':')
        fIn=ROOT.TFile.Open(url)
        histos[proc]={}
        for ch in ['e','mu']:
            for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
                key=(ch,cat)
                try:
                    histos[proc][key]=fIn.Get(var+'_%s%s'%key).Clone(var+'_%s%s_'%key+proc)
                    histos[proc][key].SetDirectory(0)
                except:
                    pass
                if proc!='data' or cat!='1l4j2q': continue
                qcdcat=cat.replace('l','f')
                histos['qcd'][key]=fIn.Get('%s_%s%s'%(var,ch,qcdcat)).Clone('%s_%s%s_%s'%(var,ch,qcdcat,proc))
                histos['qcd'][key].SetDirectory(0)
        fIn.Close()

    #fix zeros
    for key in histos:
        for cat in histos[key]:
            histos[key][cat].Rebin()
            for xbin in xrange(1,histos[key][cat].GetNbinsX()+1):
                val=histos[key][cat].GetBinContent(xbin)
                histos[key][cat].SetBinContent(xbin,max(val,1e-3))

    return histos


"""
"""
def fitQCD(args,varName='mtw',c=None):
    
    #fit the QCD yields
    histos=getHistos(args,varName)
    for key in histos['data']:

        var=ROOT.RooRealVar(varName,
                            'M_{T} [GeV]' if varName=='mtw' else 'MET [GeV]',
                            0,200)
        data=ROOT.RooDataHist('data','data',ROOT.RooArgList(var),histos['data'][key])


        histos['ttbar'][key].Scale(lumi[0]*PROCNORM[('ttbar',key[0])])
        nonqcdH=histos['ttbar'][key].Clone('nonqcd')
        try:
            histos['wjets'][key].Scale(lumi[0]*PROCNORM[('wjets',key[0])])
            nonqcdH.Add(histos['wjets'][key])
        except:
            pass
        try:
            histos['dy'][key].Scale(lumi[0]*PROCNORM[('dy',key[0])])
            nonqcdH.Add(histos['dy'][key])
        except:
            pass

        nonQCD=ROOT.RooDataHist('nonQCD','nonQCD',ROOT.RooArgList(var),nonqcdH)
        nonQCDpdf=ROOT.RooHistPdf('nonQCDpdf','nonQCDpdf',ROOT.RooArgSet(var),nonQCD,0)
        nnonQCDVal=nonqcdH.Integral()
        nnonQCD=ROOT.RooRealVar('nnonQCD','nnonQCD',nnonQCDVal,nnonQCDVal*0.7,nnonQCDVal*1.3)

        qcdKey=(key[0],'1l4j2q')
        qcd=ROOT.RooDataHist('qcd','qcd',ROOT.RooArgList(var),histos['qcd'][qcdKey])
        qcdpdf=ROOT.RooHistPdf('qcdpdf','qcdpdf',ROOT.RooArgSet(var),qcd,0)
        nqcd=ROOT.RooRealVar('nqcd','nqcd',0,1e6)

        model=ROOT.RooAddPdf('model','model',ROOT.RooArgList(nonQCDpdf,qcdpdf),ROOT.RooArgList(nnonQCD,nqcd) )

        model.fitTo(data,ROOT.RooFit.Save(),ROOT.RooFit.Minos(1))

        #
        qcdEst,qcdEstUnc=nqcd.getVal(),nqcd.getError()
        if not key in QCDNORM:
            QCDNORM[key]=(qcdEst,qcdEstUnc)
        else:
            diff=abs(qcdEst-QCDNORM[key][0])
            QCDNORM[key]=(QCDNORM[key][0],
                          ROOT.TMath.Sqrt(diff**2+QCDNORM[key][1]**2))

        #show fit results
        c.Clear()
        frame=var.frame()
        data.plotOn(frame,ROOT.RooFit.Name('data'))
        color=ROOT.kGray
        model.plotOn(frame,
                     ROOT.RooFit.Components('non*'),
                     ROOT.RooFit.Name('nonqcd'),
                     ROOT.RooFit.ProjWData(data),
                     ROOT.RooFit.LineColor(color),
                     ROOT.RooFit.FillColor(color),
                     ROOT.RooFit.FillStyle(1001),
                     ROOT.RooFit.DrawOption('f'),
                     ROOT.RooFit.LineWidth(1),
                     ROOT.RooFit.MoveToBack()
                     )
        model.plotOn(frame,ROOT.RooFit.Name('qcd'),ROOT.RooFit.ProjWData(data),ROOT.RooFit.MoveToBack())
        frame.Draw()

        leg=ROOT.TLegend(0.6,0.7,0.9,0.6)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.AddEntry(data,'Data','ep')
        leg.AddEntry('qcd','QCD','l')
        leg.AddEntry('nonqcd','t#bar{t},W,DY','f')
        leg.Draw()
    
        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.035)
        label.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
        label.DrawLatex(0.6,0.86,'174 nb^{-1} (#sqrt{s_{NN}}=8.16 TeV)')
        label.DrawLatex(0.6,0.82,'#it{%s%s}'%key)
        label.DrawLatex(0.6,0.78,'N_{QCD}=%3.0f#pm%3.0f'%(qcdEst,qcdEstUnc))

        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs(varName+'fit_%s%s.'%key+ext)
        #raw_input()

"""
"""
def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)
    
    fitQCD(sys.argv[1:],'met',c)
    fitQCD(sys.argv[1:],'mtw',c)

    #saturate the uncertainty at 50% at least but not more than 100%
    for key in QCDNORM:
        val,valUnc=QCDNORM[key]
        valUnc=min(val,max(0.5*val,valUnc))
        QCDNORM[key]=(val,valUnc)

    #dump results to a pickle file
    with open('qcdnorm.pck','wb') as fOut:
        pickle.dump(QCDNORM,fOut)

if __name__ == "__main__":
    main()
