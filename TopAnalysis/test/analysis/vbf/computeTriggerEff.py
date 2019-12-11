import ROOT
import os,sys
from array import array

def getData2MC(data,mc):

    """computes the ratio of two efficiency graphs"""

    gr=data.Clone('%s_SF'%data.GetName())
    gr.Set(0)
    x,y=ROOT.Double(0),ROOT.Double(0)
    xmc,ymc=ROOT.Double(0),ROOT.Double(0)
    for i in range(data.GetN()):

        data.GetPoint(i,x,y)

        #find closest (sometimes points are missing)
        closej=-1
        closedx=9999999.
        for j in range(mc.GetN()):
            mc.GetPoint(j,xmc,ymc)
            dx=abs(float(xmc)-float(x))
            if dx>closedx: continue
            closedx=dx
            closej=j
        if closedx>10: continue

        j=closej
        mc.GetPoint(j,xmc,ymc)
        if float(ymc)<0.01 : continue

        ip=gr.GetN()
        sf=float(y)/float(ymc)
        gr.SetPoint(ip,float(x),sf)
        sflo=(float(y)-data.GetErrorYlow(i))/(float(ymc)+mc.GetErrorYhigh(i))
        sfhi=(float(y)+data.GetErrorYhigh(i))/(float(ymc)-mc.GetErrorYlow(i))
        gr.SetPointError(ip,0,0,abs(sflo-sf),abs(sfhi-sf))        

    return gr


def fillTrigHisto(var,hdef,tag,probe,data,effOpt):

    """projects the tree to build a specific histogram before and after the probe cut"""

    histos=[]
    for hname,cut in [ ('tagdata',    tag),
                       ('probedata',  '%s && %s'%(tag,probe))]:
        data.Draw('{0} >> {1}'.format(var,hdef.GetName()),'wgt*({0})'.format(cut),'goff')
        histos.append( hdef.Clone(hname) )
        histos[-1].SetDirectory(0)
        histos[-1].Sumw2()
        hdef.Reset('ICE')

    #compute efficiency
    effgr=ROOT.TGraphAsymmErrors()
    effgr.SetFillStyle(0)
    effgr.Divide(histos[1],histos[0],effOpt)
    histos.append(effgr)

    return histos

def show2DSFs(h2d,extraTxt,outdir):

    """show the final 2D"""

    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.12)
    #c.SetLogx()

    h2d.Draw('colz') # text E')
    h2d.GetZaxis().SetTitleOffset(-0.3)
    #h2d.GetXaxis().GetMoreLogLabels()
    h2d.GetZaxis().SetRangeUser(0,1)

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{Preliminary}')
    txt.SetTextSize(0.04)
    txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    txt.DrawLatex(0.95,0.97,extraTxt)

    for ext in ['png','pdf']:
        c.SaveAs('%s/%s.%s'%(outdir,h2d.GetName(),ext))

def showTriggerEfficiency(histos,hname,extraTxt,xtitle,outdir):

    """plots a standard canvas with the trigger efficiencies and scale factors"""

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
    p1.cd()

    frame=histos['mc'][0].Clone('frame')
    frame.Reset('ICE')
    frame.Draw()
    frame.GetXaxis().SetTitleSize(0.0)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.GetYaxis().SetTitleSize(0.08)
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().SetTitle('Efficiency')
    frame.GetYaxis().SetRangeUser(0.04,1.)

    mg=ROOT.TMultiGraph()
    histos['mc'][-1].SetMarkerStyle(24)
    histos['mc'][-1].SetMarkerColor(ROOT.kAzure+1)
    histos['mc'][-1].SetLineColor(ROOT.kAzure+1)
    mg.Add(histos['mc'][-1],'p')

    histos['data'][-1].SetMarkerStyle(20)
    for i in range(histos['data'][-1].GetN()):
        histos['data'][-1].SetPointEXhigh(i,0)
        histos['data'][-1].SetPointEXlow(i,0)

    mg.Add(histos['data'][-1],'p')
    mg.Draw('p')

    leg=ROOT.TLegend(0.55,0.12,0.95,0.17)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.08)
    leg.SetNColumns(2)
    leg.AddEntry(histos['data'][-1], 'Data',         'ep')
    leg.AddEntry(histos['mc'][-1],   'EWK #gammajj', 'ep')
    leg.Draw()

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.08)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{Preliminary}')
    txt.SetTextSize(0.07)
    txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    txt.DrawLatex(0.95,0.95,extraTxt)
    p1.RedrawAxis()

    c.cd()
    p2=ROOT.TPad('p2','p2',0,0,1,0.5)
    p2.SetRightMargin(0.03)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.01)
    p2.SetBottomMargin(0.18)
    p2.SetGridy()
    p2.Draw()
    p2.cd()

    rframe=histos['mc'][0].Clone('rframe')
    rframe.Reset('ICE')
    rframe.Draw()
    rframe.GetXaxis().SetTitleSize(0.08)
    rframe.GetXaxis().SetLabelSize(0.08)
    rframe.GetYaxis().SetTitleSize(0.08)
    rframe.GetYaxis().SetLabelSize(0.08)
    rframe.GetYaxis().SetTitleOffset(0.7)
    rframe.GetYaxis().SetTitle('Data / MC')
    rframe.GetXaxis().SetTitle(xtitle)
    rframe.GetYaxis().SetRangeUser(0.0,1.04)

    sfgr=getData2MC(histos['data'][-1],histos['mc'][-1])
    sfgr.Draw('p')
    
    """  
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

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.07)
    txt.SetTextAlign(12)
    for ip in range(3):
        txt.DrawLatex(0.15,0.9-0.07*ip,'p_{%d}=%3.3f'%(ip+1,sff.GetParameter(ip)))
    p2.RedrawAxis()
    """

    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/%s.%s'%(outdir,hname,ext))


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPaintTextFormat('4.2f')
ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)


#build the chains with data fo trigger measuremnets
url=sys.argv[1]
outdir='%s/trigeff'%url
os.system('mkdir -p %s'%outdir)
year=2016
data=ROOT.TChain('tree')
mc=ROOT.TChain('tree')
for f in os.listdir(url):
    if 'Data' in f:
        data.AddFile(os.path.join(url,f))
        if '2017' in f: year=2017
        if '2018' in f: year=2018
    if 'EWKAJJ.root' in f:
        mc.AddFile(os.path.join(url,f))

binDef={'lowvpt':{'pt' :   [70,72.5,75,77.5,80,85,90,100,120,150,200],
                  'mjj':   [475,525,575,625,675,775,1000,2000],
                  'eta':   [0,1.442],
                  'tag':   'passLowPtCtrlTrig && lowPtHighMJJCtrTrigActive && mjj>500 && detajj>3.5 && fabs(veta)<1.5 && r9>0.9 && hoe<0.1',
                  'probe': 'passLowPtHighMJJTrig'
                  },
        'highvpt':{'pt' :   [170,172.5,175,177.5,180,185,190,195,200,205,210,220,230,240,250,300,400,500],
                   'mjj':   [475,525,575,625,775,1000,2000],
                   'eta':   [0,1.442,2.4],
                   'tag':   'passHighPtCtrlTrig',
                   'probe': 'passHighPtTrig'                   
                   },
        }
if year==2017 :
    binDef['highvpt']['pt']=[195,197.5,200,202.5,205,210,220,230,240,250,300,400,500]

varTitles={'vpt':'Photon p_{T} [GeV]',
           'mjj':'Dijet invariant mass [GeV]'}

finalSFs=[]

for cat in ['lowvpt']: # binDef:

    hbase={'vpt':ROOT.TH1F('hpt',  'hpt',  len(binDef[cat]['pt'])-1,  array('d',binDef[cat]['pt'])),
           'mjj':ROOT.TH1F('hmjj', 'hmjj', len(binDef[cat]['mjj'])-1, array('d',binDef[cat]['mjj']))}

    #loop over bins in eta (first bin is inclusive)
    for i in range(0,len(binDef[cat]['eta'])):

        if i==0:
            etamin=0
            etamax=binDef[cat]['eta'][-1]
        else:
            etamin=binDef[cat]['eta'][i-1]
            etamax=binDef[cat]['eta'][i]
    
        #get efficiency as function of photon pT in different mjj ranges
        hptvsmjj=ROOT.TH2F('hptvmjj_%d_%s'%(i,cat),  
                           ';%s;%s;Data / MC'%(varTitles['vpt'],varTitles['mjj']),  
                           len(binDef[cat]['pt'])-1,  array('d',binDef[cat]['pt']),
                           len(binDef[cat]['mjj'])-1, array('d',binDef[cat]['mjj']))    
        for im in range(1,len(binDef[cat]['mjj'])):        
            tagcut  = binDef[cat]['tag']+' && fabs(veta)>=%f && fabs(veta)<%f'%(etamin,etamax)
            tagcut += ' && mjj>=%f && mjj<%f'%(binDef[cat]['mjj'][im-1],binDef[cat]['mjj'][im])
            probecut = tagcut+' && ' + binDef[cat]['probe']
            histos={'data' : fillTrigHisto('vpt',hbase['vpt'],tagcut,probecut,data,'cp'), 
                    'mc'   : fillTrigHisto('vpt',
                                           hbase['vpt'],
                                           tagcut.replace('&& lowPtHighMJJCtrTrigActive',''),
                                           probecut.replace('&& lowPtHighMJJCtrTrigActive',''),
                                           mc,  
                                           'central')}
            sfgr=getData2MC(histos['data'][-1],histos['mc'][-1])

            #copy to the appropriate row in the final scale factor histogram
            xval,sfval=ROOT.Double(0),ROOT.Double(0)
            for ip in range(sfgr.GetN()):
                sfgr.GetPoint(ip,xval,sfval)
                sfvalunc=sfgr.GetErrorY(ip)
                xbin=hptvsmjj.GetXaxis().FindBin(float(xval))
                hptvsmjj.SetBinContent(xbin,im,float(sfval))
                hptvsmjj.SetBinError(xbin,im,sfvalunc)
                   
        show2DSFs(h2d=hptvsmjj,
                  extraTxt='%3.1f<|#eta|<%3.1f (%d)'%(etamin,etamax,year),
                  outdir=outdir)
        finalSFs.append(hptvsmjj)
        finalSFs[-1].SetDirectory(0)

        #plot vs pt or mjj inclusively
        for key in hbase:

            tagcut   = binDef[cat]['tag']+' && fabs(veta)>=%f && fabs(veta)<%f'%(etamin,etamax)

            #add low photon pt requirement
            if key=='mjj':
                lowvpt=75 if 'lowvpt' in cat else 200
                tagcut += ' && vpt>%f'%lowvpt

            probecut = tagcut+' && ' + binDef[cat]['probe']

            #compute trigger efficiency
            histos={'data' : fillTrigHisto(key,hbase[key],tagcut,probecut,data,'cp'), 
                    'mc'   : fillTrigHisto(key,
                                           hbase[key],
                                           tagcut.replace('&& lowPtHighMJJCtrTrigActive',''),
                                           probecut.replace('&& lowPtHighMJJCtrTrigActive',''),
                                           mc,  
                                           'central')}
            showTriggerEfficiency(histos=histos,
                                  hname='%s_%d_%s'%(key,i,cat),
                                  extraTxt='%3.1f<|#eta|<%3.1f (%d)'%(etamin,etamax,year),
                                  xtitle=varTitles[key],
                                  outdir=outdir)

fname='%s/trigger_sf.root'%outdir
fOut=ROOT.TFile(fname,'RECREATE')
for h in finalSFs:
    h.SetDirectory(fOut)
    h.Write()
fOut.Close()
print 'Scale factors for year=',year,'available in',fname
