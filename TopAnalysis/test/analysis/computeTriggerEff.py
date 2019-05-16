import ROOT
import os,sys

url=sys.argv[1]
year=sys.argv[2]

HISTOLIST=[('hptoff_apt','hpttrig_hptoff_apt'),
           ('lptoff_apt','lpttrig_lptoff_apt'),
           #('lpthmjjoff_apt','lpthmjjtrig_lpthmjjoff_apt'),
           ('lpthmjjoff_mjj','lpthmjjtrig_lpthmjjoff_mjj'),
           ]


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

    return gr

def scaleGr(gr,sf):
    x,y=ROOT.Double(0),ROOT.Double(0)
    for ip in range(gr.GetN()):
        gr.GetPoint(ip,x,y)
        gr.SetPoint(ip,float(x),float(y*sf))
        gr.SetPointEYhigh(ip,gr.GetErrorYhigh(ip)*sf)
        gr.SetPointEYlow(ip,gr.GetErrorYlow(ip)*sf)
        

def getEffPlots(url,year,isData=True):

    fIn=ROOT.TFile.Open(url)

    histos={}
    for a,b in HISTOLIST:
        histos[a]=None
        histos[b]=None

    for h in histos:
        tag=''
        if not isData: tag='_EWK #gammajj'
        try:
            histos[h]=fIn.Get('{0}/{0}{1}'.format(h,tag))
            histos[h].SetDirectory(0)
        except Exception as e:
            print e
            pass
        
    fIn.Close()    

    effGrs=[ getEffGraph(b,histos[b],histos[a],isData) for a,b in HISTOLIST ]
#    if isData:
#        for gr in effGrs:
#            grname=gr.GetName()
#            sf=1.
#            if year=='2016' and not 'hpt' in grname:
#                sf=2567./35879. if not 'hmjj' in grname else 28412./35879.
#            if year=='2017' and not 'hpt' in grname:
#                sf=1327./41367. if not 'hmjj' in grname else 7661./41367.
#            scaleGr(gr,sf)
#
    return effGrs

url=sys.argv[1]
effGr=getEffPlots(url,year)
mcEffGr=getEffPlots(url,year,isData=False)


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

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
    mcEffGr[i].GetYaxis().SetRangeUser(0.01,1.)

    leg1=ROOT.TLegend(0.15,0.88,0.6,0.66)
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.08)
    leg1.AddEntry(effGr[i],'Data','p')
    leg1.AddEntry(mcEffGr[i],'EWK #gammajj','l')
    leg1.Draw()

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
    frame.GetYaxis().SetRangeUser(0.62,1.28)
    frame.GetXaxis().SetTitleSize(0.08)
    frame.GetXaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleSize(0.08)
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().SetTitle('Data/MC')
    
    sfgr=getData2MC(effGr[i],mcEffGr[i])
    minX=150
    if 'vbf' in effGr[i].GetName(): minX=50
    sff=ROOT.TF1('sff','0.5*[0]*(1.+TMath::Erf((x-[1])/(TMath::Sqrt(2.)*[2])))',minX,frame.GetYaxis().GetXmax())    
    sff.SetParLimits(0,0.1,1)
    if 'hptoff_apt' in effGr[i].GetName():
        sff.SetParLimits(1,150,210)
        sff.SetParLimits(2,1,50)
    elif 'lptoff_apt' in effGr[i].GetName():
        sff.SetParLimits(1,70,80)
        sff.SetParLimits(2,1,50)
    else:
        sff.SetParLimits(1,200,1500)
        sff.SetParLimits(2,100,500)
        sfgr.Fit(sff,'M','')
        sfVal=sff.GetParameter(0)
        scaleGr(effGr[i],1./sfVal)
        sfgr=getData2MC(effGr[i],mcEffGr[i])

    sfgr.Fit(sff,'M+')


    sfgr.Draw('p')

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.07)
    txt.SetTextAlign(12)
    for ip in range(3):
        txt.DrawLatex(0.15,0.9-0.07*ip,'p_{%d}=%3.3f'%(ip+1,sff.GetParameter(ip)))
    p2.RedrawAxis()

    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s_%s.%s'%(effGr[i].GetName(),year,ext))
