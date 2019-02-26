import ROOT
import os,sys

url=sys.argv[1]
year=sys.argv[2]

HISTOLIST={'2017':[('offlinephoton_apt',    'photon200_apt'),
                   #('offlinephotonvbf_apt', 'photon75_vbf2017_apt'),
                   #('offlinephotonvbf_mjj', 'photon75_vbf2017_mjj')
                   ],
           '2016':[('offlinephoton_apt',    'photon175_apt'),
                   ('offlinephotonvbf_apt', 'photon75_vbf2016_apt'),
                   ('offlinephotonvbf_mjj', 'photon75_vbf2016_mjj')]
           }

def getData2MC(data,mc):

    gr=data.Clone('%s_SF'%data.GetName())
    gr.Set(0)
    x,y=ROOT.Double(0),ROOT.Double(0)
    for i in range(data.GetN()):
        data.GetPoint(i,x,y)
        
        xbin=mc.GetXaxis().FindBin(float(x))
        mcVal=mc.GetBinContent(xbin)

        if mcVal<0.05 : continue   

        ip=gr.GetN()
        gr.SetPoint(ip,float(x),float(y)/mcVal)
        gr.SetPointError(ip,0,0,data.GetErrorYlow(i),data.GetErrorYhigh(i))
    return gr

def getEffGraph(name,total_pass,total,isData):

    if not isData:
        gr=total_pass.Clone(name+'mc')
        gr.Divide(total)
        gr.SetFillStyle(0)
        gr.SetLineColor(ROOT.kBlue)
        gr.SetLineWidth(2)
        gr.SetTitle('MC')
    else:
        gr=ROOT.TGraphAsymmErrors()
        gr.BayesDivide(total_pass,total)
        gr.SetName(name)
        gr.SetTitle('Data')
        gr.SetFillStyle(0)
        gr.SetMarkerStyle(20)
        turnOn=ROOT.TF1('turnon',
                        '[0]+[1]/(1.+TMath::Exp([2]*(x-[3])))',
                        total.GetXaxis().GetXmin(),
                        total.GetXaxis().GetXmax())
        #turnOn.SetParLimits(0,0.,0.8)
        #turnOn.SetParLimits(1,0.5,1.2)
        #turnOn.SetParLimits(2,-2,2)
        # turnOn.SetParLimits(3,total.GetXaxis().GetXmin(),total.GetXaxis().GetXmax())
        gr.Fit(turnOn,'MR+')

    return gr


def getEffPlots(plist,isData=True):

    histos={}
    for a,b in HISTOLIST[year]:
        histos[a]=None
        histos[b]=None

    for p in plist:
        fIn=ROOT.TFile(p)
        for h in histos:

            if isData and year=='2017':
                if 'vbf' in h and not '2017F_SingleMuon' in p : continue
                
            h1=fIn.Get(h)
            try:
                h1.Integral()
                if histos[h] is None:
                    histos[h]=h1.Clone()
                    histos[h].SetDirectory(0)
                else:
                    histos[h].Add(h1)
            except:
                pass
        fIn.Close()
    
    print histos
    effGrs=[ getEffGraph(b,histos[a],histos[b],isData) for a,b in HISTOLIST[year] ]

    return effGrs

url=sys.argv[1]
effGr=getEffPlots(plist=[os.path.join(url,f) for f in os.listdir(url) if 'SingleMuon' in f])
mcEffGr=getEffPlots(plist=[os.path.join(url,'MC13TeV_'+year+'_EWKAJJ.root')],isData=False)


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

c=ROOT.TCanvas('c','c',500,500)
c.SetBottomMargin(0)
c.SetTopMargin(0)
c.SetLeftMargin(0)
c.SetRightMargin(0)

c.cd()
p1=ROOT.TPad('p1','p1',0,0.5,1,1.0)
p1.Draw()
p1.SetRightMargin(0.03)
p1.SetLeftMargin(0.12)
p1.SetTopMargin(0.1)
p1.SetBottomMargin(0.01)
p1.SetGridy()

c.cd()
p2=ROOT.TPad('p2','p2',0,0,1,0.5)
p2.SetRightMargin(0.03)
p2.SetLeftMargin(0.12)
p2.SetTopMargin(0.01)
p2.SetBottomMargin(0.18)
p2.SetGridy()
p2.Draw()


for i in xrange(len(effGr)):

    c.cd()
    p1.cd()
    p1.Clear()
    mcEffGr[i].Draw('hist')
    effGr[i].Draw('p')    
    mcEffGr[i].GetXaxis().SetTitleSize(0.0)
    mcEffGr[i].GetXaxis().SetLabelSize(0.0)
    mcEffGr[i].GetYaxis().SetTitleSize(0.08)
    mcEffGr[i].GetYaxis().SetLabelSize(0.08)
    mcEffGr[i].GetYaxis().SetTitleOffset(0.7)
    mcEffGr[i].GetYaxis().SetRangeUser(0.05,1.)

    leg1=p1.BuildLegend(0.15,0.88,0.6,0.66)
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.08)

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.08)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{preliminary}')
    p1.RedrawAxis()


    c.cd()
    p2.cd()
    frame=mcEffGr[i].Clone()
    frame.Reset('ICE')
    frame.Draw()
    frame.GetYaxis().SetNdivisions(5)
    frame.GetYaxis().SetRangeUser(0.72,1.08)
    frame.GetXaxis().SetTitleSize(0.08)
    frame.GetXaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleSize(0.08)
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().SetTitle('Data/MC')
    
    sfgr=getData2MC(effGr[i],mcEffGr[i])
    sfgr.Draw('p')
    p2.RedrawAxis()

    c.cd()
    c.Modified()
    c.Update()
    raw_input()
    for ext in ['png','pdf']:
        c.SaveAs('%s_%s.%s'%(effGr[i].GetName(),year,ext))
