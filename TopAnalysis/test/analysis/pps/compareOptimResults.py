import ROOT
import sys
import os
import pickle
import re

def showResults(plotColl,title,name):
    
    """shows a summary plot with the optimization results"""

    nPlots=len(plotColl)
    wid=1./nPlots

    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()
    leg=ROOT.TLegend(0.05,0.05,0.5,0.01)
    leg.SetNColumns(max(ROOT.TMath.FloorNint(nPlots/2),1))
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)
    leg.SetBorderSize(0)

    drawOpt='hbar'
    for i in range(nPlots):
        h,_=plotColl[i]
        h.SetBarWidth(wid)
        h.SetBarOffset(i*wid)
        h.GetYaxis().SetRangeUser(0,5 if 'sig' in name else 5)
        h.GetYaxis().SetTitleOffset(1)
        h.GetXaxis().SetTitleOffset(1)
        h.Draw(drawOpt)
        drawOpt='hbarsame'
        leg.AddEntry(h,h.GetTitle(),'f')

    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
    tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignBottom)
    tex.DrawLatex(0.95,0.96,title)

    c.Modified()
    c.Update()
    c.RedrawAxis()

    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(name,ext))


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

resultsDir=sys.argv[1]
resultsPerMass={}
colors={120:ROOT.kGray,130:ROOT.kOrange,140:ROOT.kGreen+1,150:ROOT.kAzure+3}

#search for optimization results pickle files in subdirectories
pckList=[]
for mres in os.listdir(resultsDir):
    if not 'stat_m' in mres:
        continue
    for url in os.listdir(os.path.join(resultsDir,mres)):
        if not 'optim_results' in url:
            continue
        pckList.append( os.path.join(resultsDir,mres,url) )

#loop over pickle files found
for url in pckList:
    print url
    regex=re.compile('_a(\d+)') 
    xangle=int(regex.findall(url)[0])
    with open(url,'r') as f:

        res=pickle.load(f)

        hr95=pickle.load(f)
        hr95.SetTitle('%d #murad'%xangle)
        hr95.SetLineColor(colors[xangle])
        hr95.SetFillColor(colors[xangle])

        hsig=pickle.load(f)
        hsig.SetTitle('%d #murad'%xangle)
        hsig.SetLineColor(colors[xangle])
        hsig.SetFillColor(colors[xangle])

        bestByR95=sorted(res,key=lambda x: x[-2],reverse=True)[0]
        bestBySig=sorted(res,key=lambda x: x[-1],reverse=True)[0]

        mass=bestByR95[1]
        if not mass in resultsPerMass:
            resultsPerMass[mass]={'r95':[],'sig':[],'xangle':[]}
        resultsPerMass[mass]['xangle'].append(xangle)
        resultsPerMass[mass]['r95'].append( (hr95.Clone('optim_r95_m%d_a%d'%(mass,xangle)), bestByR95[0:-2]) )
        resultsPerMass[mass]['sig'].append( (hsig.Clone('optim_sig_m%d_a%d'%(mass,xangle)), bestBySig[0:-2]) )
        resultsPerMass[mass]['r95'][-1][0].SetDirectory(0)
        resultsPerMass[mass]['sig'][-1][0].SetDirectory(0)

#save optimization results by significance
finalAnalysisCuts=[]
for key in ['r95','sig']:

    for mass in resultsPerMass:
        showResults(resultsPerMass[mass][key],
                    'm=%d GeV'%mass,
                    'optim_%s_m%d'%(key,mass))

        if key=='sig':
            nres=len(resultsPerMass[mass][key])
            for i in range(nres):
                _,res=resultsPerMass[mass][key][i]
                xangle=resultsPerMass[mass]['xangle'][i]
                finalAnalysisCuts.append( [mass,xangle]+res[2:6] )

print finalAnalysisCuts
with open(resultsDir+'/final_analysis_cuts.pck','w') as cache:
    pickle.dump(finalAnalysisCuts,cache,pickle.HIGHEST_PROTOCOL) 

