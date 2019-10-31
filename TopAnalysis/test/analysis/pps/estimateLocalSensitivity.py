import ROOT
import os
import sys
import optparse
import numpy as np

CSILIST=np.arange(0.04,0.05,0.001)
PTLIST=np.arange(20,60,5)
NVTXLIST=np.arange(10,40,5)
PFSUMPZLIST=np.arange(4000,12000,1000)


def estimateLocalSensitivity(opt):

    """ steers the estimation of the local sensitivity for a given distribution """

    data=ROOT.TChain('data')
    for f in [os.path.join(opt.input,x) for x in os.listdir(opt.input) if 'Data13TeV' in x and 'DoubleMu' in x]:
        if 'MuonEG' in f : continue
        data.AddFile(f)

    cuts='xangle==%d && mixType==1 && cat==169 && l1pt>30 && l2pt>20 && bosonpt>40 && nvtx<20 && PFPzSumHF<12000'%opt.xangle
    finalCut=cuts+' && csi1>0.04 && csi2>0.04'
    data.Draw('mmiss >> (50,0,2500)',finalCut,'goff')
    h=data.GetHistogram()    
    h.SetDirectory(0)

    #local  variation graphs will be approximated by pol1
    csiEvol=[]   
    ptllEvol=[]
    vtxEvol=[]
    hfEvol=[]
    for xbin in range(50):
        csiEvol.append(ROOT.TGraph())
        ptllEvol.append(ROOT.TGraph())
        vtxEvol.append(ROOT.TGraph())
        hfEvol.append(ROOT.TGraph())
    gfunc=ROOT.TF1('grad','[0]*x+[1]',-100,100)

    #scan in csi
    print 'Scanning csi in',CSILIST
    for csi in CSILIST:
        finalCut='%s && csi1>%f && csi2>%f'%(cuts,csi,csi)
        data.Draw('mmiss >> (50,0,2500)',finalCut,'goff')
        hvar=data.GetHistogram()
        for xbin in range(h.GetNbinsX()):
            nomCts=h.GetBinContent(xbin+1)
            if nomCts==0: continue
            rel_diff=100.*(hvar.GetBinContent(xbin+1)/nomCts)
            csiEvol[xbin].SetPoint(csiEvol[xbin].GetN(),100*(csi/0.5-1),rel_diff)

    #scan in ptll
    print 'Scanning ptll in',PTLIST
    cuts='xangle==%d && mixType==1 && cat==169 && l1pt>30 && l2pt>20 && csi1>0.04 && csi2>0.04 && nvtx<20 && PFPzSumHF<12000'%opt.xangle
    for pt in PTLIST:
        finalCut='%s && bosonpt>%f'%(cuts,pt)
        data.Draw('mmiss >> (50,0,2500)',finalCut,'goff')
        hvar=data.GetHistogram()
        for xbin in range(h.GetNbinsX()):
            nomCts=h.GetBinContent(xbin+1)
            if nomCts==0: continue
            rel_diff=100.*(hvar.GetBinContent(xbin+1)/nomCts)
            ptllEvol[xbin].SetPoint(ptllEvol[xbin].GetN(),100*(pt/40.-1),rel_diff)

    #scan in vertices
    print 'Scanning nvtx in',NVTXLIST
    cuts='xangle==%d && mixType==1 && cat==169 && l1pt>30 && l2pt>20 && csi1>0.04 && csi2>0.04 && bosonpt>40 && PFPzSumHF<12000'%opt.xangle
    for n in NVTXLIST:
        finalCut='%s && nvtx<%f'%(cuts,n)
        data.Draw('mmiss >> (50,0,2500)',finalCut,'goff')
        hvar=data.GetHistogram()
        for xbin in range(h.GetNbinsX()):
            nomCts=h.GetBinContent(xbin+1)
            if nomCts==0: continue
            rel_diff=100.*(hvar.GetBinContent(xbin+1)/nomCts)
            vtxEvol[xbin].SetPoint(vtxEvol[xbin].GetN(),100*(n/20.-1),rel_diff)

    #scan in HF
    print 'Scanning PFSumPz in',PFSUMPZLIST
    cuts='xangle==%d && mixType==1 && cat==169 && l1pt>30 && l2pt>20 && csi1>0.04 && csi2>0.04 && bosonpt>40 && nvtx<20'%opt.xangle
    for hf in PFSUMPZLIST:
        finalCut='%s && PFPzSumHF<%f'%(cuts,hf)
        data.Draw('mmiss >> (50,0,2500)',finalCut,'goff')
        hvar=data.GetHistogram()
        for xbin in range(h.GetNbinsX()):
            nomCts=h.GetBinContent(xbin+1)
            if nomCts==0: continue
            rel_diff=100.*(hvar.GetBinContent(xbin+1)/nomCts)
            hfEvol[xbin].SetPoint(hfEvol[xbin].GetN(),100*(hf/12000.-1),rel_diff)

    #local sensitivities
    csils=h.Clone('csils')
    csils.Reset('ICE')
    csils.SetFillStyle(0)
    csils.SetLineWidth(2)
    csils.SetLineColor(ROOT.kAzure-2)    
    ptls=csils.Clone('ptls')
    ptls.SetLineColor(ROOT.kRed+1)
    vtxls=csils.Clone('vtxls')
    vtxls.SetLineColor(1)
    hfls=csils.Clone('hfls')
    hfls.SetLineColor(ROOT.kGreen+3)
    hfls.SetLineStyle(9)
    for xbin in range(h.GetNbinsX()):
        if csiEvol[xbin].GetN()>2 :
            csiEvol[xbin].Sort()
            csiEvol[xbin].Fit(gfunc,'MQ+')
            grad=gfunc.GetParameter(0)
            csils.SetBinContent(xbin+1,grad)
        if ptllEvol[xbin].GetN()>2:
            ptllEvol[xbin].Sort()
            ptllEvol[xbin].Fit(gfunc,'MQ+')
            grad=gfunc.GetParameter(0)
            ptls.SetBinContent(xbin+1,grad)
        if vtxEvol[xbin].GetN()>2:
            vtxEvol[xbin].Sort()
            vtxEvol[xbin].Fit(gfunc,'MQ+')
            grad=gfunc.GetParameter(0)
            vtxls.SetBinContent(xbin+1,grad)

        if hfEvol[xbin].GetN()>2:
            hfEvol[xbin].Sort()
            hfEvol[xbin].Fit(gfunc,'MQ+')
            grad=gfunc.GetParameter(0)
            hfls.SetBinContent(xbin+1,grad)

    #display results
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.04)
    c.SetGridy()
    h.Scale(100./h.Integral())
    h.SetFillStyle(1001)
    h.SetFillColor(ROOT.kGray)
    h.SetLineColor(ROOT.kGray)
    h.SetMarkerColor(ROOT.kGray)
    h.SetTitle('PDF')
    h.GetYaxis().SetTitle('d(N/N_{0}) / d(x/x_{0}) [%]')
    h.GetXaxis().SetTitle('Missing mass [GeV]')
    h.Draw('hist')
    h.GetYaxis().SetRangeUser(-10,10)
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetTitleSize(0.045)    

    csils.Draw('histsame')
    csils.SetTitle('#xi')
    ptls.Draw('histsame')
    ptls.SetTitle('p_{T}(ll)')
    vtxls.Draw('histsame')
    vtxls.SetTitle('N_{vtx}')
    hfls.Draw('histsame')
    hfls.SetTitle('#Sigma_{i#inHF} |p_{Z}|')

    leg=c.BuildLegend(0.2,0.4,0.4,0.2)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.SetFillStyle(0)
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.96,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.17,0.9,'[%d#murad]'%opt.xangle)
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,'13 TeV')  

    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']:
        c.SaveAs('%s/localsens_%d.%s'%(opt.output,opt.xangle,ext))

def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='input directory [%default]',  
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/',
                      type='string')
    parser.add_option('-o', '--out',          
                      dest='output',       
                      help='output directory [%default]',  
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/plots',
                      type='string')
    parser.add_option('--xangle',
                      dest='xangle',
                      default=120,
                      type=int,
                      help='crossing angle [%default]')
    (opt, args) = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    estimateLocalSensitivity(opt)


if __name__ == "__main__":
    sys.exit(main())
