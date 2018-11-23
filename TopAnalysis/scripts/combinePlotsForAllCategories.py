from plotter import *
import sys

"""
"""
def doPlot(plotName,ch):

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    #open the files
    inFiles=[]
    baseDir='~/work/LJets2015-arcrev'
    if ch=='all' or ch=='plus' or ch=='mu':  inFiles.append(ROOT.TFile.Open('%s/analysis_muplus/plots/final_plotter.root'%baseDir))
    if ch=='all' or ch=='plus' or ch=='e':   inFiles.append(ROOT.TFile.Open('%s/analysis_eplus/plots/final_plotter.root'%baseDir))
    if ch=='all' or ch=='minus' or ch=='e':  inFiles.append(ROOT.TFile.Open('%s/analysis_eminus/plots/final_plotter.root'%baseDir))
    if ch=='all' or ch=='minus' or ch=='mu': inFiles.append(ROOT.TFile.Open('%s/analysis_muminus/plots/final_plotter.root'%baseDir))

    #baseDir='./'
    #inFiles.append(ROOT.TFile.Open('%s/analysis_mu/plots/plotter.root'%baseDir))
    #inFiles.append(ROOT.TFile.Open('%s/analysis_e/plots/plotter.root'%baseDir))

    plotsPerProc={}
    if plotName=='nbtags':

        #the template
        nbjetsnbtags=ROOT.TH1F('nbjetsnbtags',';Event category;Events',11,0,11)
        xbin=0
        for ijet in xrange(1,5):
            for itag in xrange(0,3):
                if itag>ijet: continue
                xbin+=1
                nbjetsnbtags.GetXaxis().SetBinLabel(xbin,'#splitline{%dj}{%dt}'%(ijet,itag))

        #sum up
        xbin=0
        for ijet in xrange(1,5):
            for itag in xrange(0,3):
                if itag>ijet: continue
                xbin+=1

                for f in inFiles:
                    pName='nbtags_%dj'%ijet
                    pdir=f.GetDirectory(pName)
                    for key in pdir.GetListOfKeys():
                        keyName=key.GetName()
                        if 'Graph' in keyName : continue

                        h=pdir.Get(keyName)
                        title=keyName.replace(pName+'_','')
                        if title==keyName : title='Data'
                        if not title in plotsPerProc:
                            plotsPerProc[title]=nbjetsnbtags.Clone(title)
                            plotsPerProc[title].SetTitle(title)
                            plotsPerProc[title].SetDirectory(0)
                            plotsPerProc[title].SetFillColor(h.GetFillColor())                        
                            plotsPerProc[title].SetLineColor(h.GetLineColor())                        
                            plotsPerProc[title].SetMarkerColor(h.GetMarkerColor())                        
                        newVal=h.GetBinContent(itag+1)+plotsPerProc[title].GetBinContent(xbin)
                        newErr=ROOT.TMath.Sqrt(h.GetBinError(itag+1)**2+plotsPerProc[title].GetBinError(xbin)**2)
                        plotsPerProc[title].SetBinContent(xbin,newVal)
                        plotsPerProc[title].SetBinError(xbin,newErr)
    else:

        for f in inFiles:
            pName=plotName
            pdir=f.GetDirectory(pName)
            for key in pdir.GetListOfKeys():
                keyName=key.GetName()
                if 'Graph' in keyName : continue

                h=pdir.Get(keyName)
                title=keyName.replace(pName+'_','')
                if title==keyName : title='Data'
                if not title in plotsPerProc:
                    plotsPerProc[title]=h.Clone(title)
                    plotsPerProc[title].Reset('ICE')
                    plotsPerProc[title].SetTitle(title)
                    plotsPerProc[title].SetDirectory(0)
                    plotsPerProc[title].SetFillColor(h.GetFillColor())                        
                    plotsPerProc[title].SetLineColor(h.GetLineColor())                        
                    plotsPerProc[title].SetMarkerColor(h.GetMarkerColor())   
                plotsPerProc[title].Add(h)
                    
    #show
    plot=Plot('%s%s'%(ch,plotName))    
    plot.savelog=True
    plot.wideCanvas=True if plotName=='nbtags' else False
    plot.ratiorange=(0.76,1.24)
    plot.plotformats=['root','pdf','png']
    for key in  plotsPerProc:
        isData=True if 'Data' in plotsPerProc[key].GetTitle() else False
        color=1 if isData else plotsPerProc[key].GetFillColor()
        plot.add(plotsPerProc[key],
                 plotsPerProc[key].GetTitle(),
                 color,
                 isData)
    plot.finalize()
    plot.show(outDir="plots/",lumi=2300)
    #plot.show(outDir="plots/",lumi=2.3,noStack=True,saveTeX=True)
    #raw_input()
                     
def main():

    os.system('mkdir -p plots')

    plots=sys.argv[1].split(',')
    ch='all'
    if len(sys.argv)>2: ch=sys.argv[2]
    for p in plots : doPlot(p,ch)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
