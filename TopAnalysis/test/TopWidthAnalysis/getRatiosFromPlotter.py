import optparse
import os,sys
import json
import ROOT
import math
import pickle

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default="./plotter.root",       type='string')
    parser.add_option('-d', '--divideFile',  dest='inDiv' ,      help='comparison input plotter',       default="./syst_plotter.root",  type='string')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi to print out',              default=35.9,                   type=float)
    parser.add_option(      '--saveLog',     dest='saveLog' ,    help='save log versions of the plots', default=False,             action='store_true')
    parser.add_option('--wids',  dest='wids'  , help='widths to compare to nominal in div', default="8.0",                   type='string')
    parser.add_option('--obs',   dest='obs'   , help='observables to process', default="incmlb", type='string')
    parser.add_option('--dists', dest='dists' , help='non-observable distributions to process', default="tmass",             type='string')
    parser.add_option('--sigs',  dest='sigs'  , help='signal processes to plot ratios for',     default="t#bar{t}",          type='string')
    parser.add_option('--ptCh',  dest='ptchs' , help='pt categories to consider',               default="highpt,lowpt",      type='string')
    parser.add_option('--ch',    dest='chs'   , help='final states to consider',                default="EE,EM,MM",          type='string')
    parser.add_option('--bcat',  dest='bcats' , help='b categories to consider',                default="1b,2b",             type='string')
    (opt, args) = parser.parse_args()

    wgtFin=ROOT.TFile(opt.inDir)
    divFin=ROOT.TFile(opt.inDiv)


    # get lists for distribution collection
    widList =opt.wids.split(',')
    ptChList=opt.ptchs.split(',')
    chList  =opt.chs.split(',')
    bcatList=opt.bcats.split(',')
    distList=opt.dists.split(',')
    obsList =opt.obs.split(',')
    sigList =opt.sigs.split(',')

    # useful settings
    divWid="1.0"

    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    outDir='./genvalidation'
    os.system('mkdir -p %s' % outDir)

    # collect a list of all distributions to plot
    plots = [("%s_%sw/%s_%sw_%s"%(dist,wid,dist,wid,sig),
        "%s_%sw/%s_%sw_%s widthx8"%(dist,divWid,dist,divWid,sig))
                for dist in distList
                for wid in widList
                for sig in sigList]
    plots+= [("%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s"%(ptC,ch,bc,obs,wid,ptC,ch,bc,obs,wid,sig),
        "%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s widthx8"%(ptC,ch,bc,obs,divWid,ptC,ch,bc,obs,divWid,sig))
                for ptC in ptChList
                for ch in chList
                for bc in bcatList
                for obs in obsList
                for wid in widList
                for sig in sigList]

    import tdrStyle
    import CMS_lumi

    tdrStyle.setTDRStyle()
    CMS_lumi.lumiTextSize=0.6
    CMS_lumi.extraText="Simulation Preliminary"
    CMS_lumi.extraOverCmsTextSize=0.5
    CMS_lumi.relPosX=0.18


    for p in plots :
        canvas=ROOT.TCanvas()
        if opt.saveLog : canvas.SetLogy()
        canvas.cd()

        wgtHist=wgtFin.Get(p[0]).Clone()
        divHist=divFin.Get(p[1]).Clone()

        divHist.Divide(wgtHist)

        divHist.GetYaxis().SetTitle("Events (4 #times#Gamma_{SM}) / Events (Reweighted)")
        divHist.GetYaxis().SetRangeUser(0.7,1.3)
        divHist.SetTitle("")
        divHist.SetMarkerColor(ROOT.kBlack)

        divHist.Draw()

        if "tmass" in p[0] :
            fit=ROOT.TF1("tm","pol1",150,200)
            divHist.Fit(fit)
            fit.SetLineColor(ROOT.kRed)

            fit.Draw("SAME")


        ROOT.gStyle.SetOptFit(0)

        CMS_lumi.CMS_lumi(canvas,4,0)
        canvas.SaveAs("%s/%s.pdf"%(outDir,p[0].split('/')[1].replace('.','w').replace('#','').replace('{','').replace('}','')))
        canvas.SaveAs("%s/%s.png"%(outDir,p[0].split('/')[1].replace('.','w').replace('#','').replace('{','').replace('}','')))

    for wid,obs,sig in [(a,b,c) for a in widList for b in obsList for c in sigList]:
        miniCanvas=ROOT.TCanvas()
        miniCanvas.SetCanvasSize(300*len(chList), 200*len(ptChList)*len(bcatList))
        miniCanvas.Clear()

        miniCanvas.cd()
        miniCanvas.Divide(len(chList),len(bcatList)*len(ptChList))

        info=ROOT.TLatex()
        info.SetTextSize(0.02)

        divHists={}
        fits={}
        for (iptC,ich,ibc) in [(iptC,ich,ibc)
                for iptC in range(len(ptChList))
                for ich  in range(len(chList))
                for ibc  in range(len(bcatList))] :

            miniCanvas.cd(len(chList)*len(bcatList)*iptC + len(chList)*ibc + (ich + 1))
            ROOT.gPad.SetGrid(1,1)
            ROOT.gPad.SetBottomMargin(ROOT.gPad.GetBottomMargin()*1.5)
            ROOT.gStyle.SetOptFit(0)
            ptC = ptChList[iptC]
            ch  = chList[ich]
            bc  = bcatList[ibc]

            wgtHist=wgtFin.Get("%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s"%(ptC,ch,bc,obs,wid,ptC,ch,bc,obs,wid,sig)).Clone()
            divHists["%i %i %i"%(iptC,ich,ibc)]=divFin.Get("%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s widthx8"%(ptC,ch,bc,obs,divWid,ptC,ch,bc,obs,divWid,sig)).Clone()

            divHist=divHists["%i %i %i"%(iptC,ich,ibc)]
            divHist.Divide(wgtHist)

            divHist.GetYaxis().SetTitle("Events / Reweighted Events")
            divHist.GetYaxis().SetTitleSize(0.065)
            divHist.GetYaxis().SetRangeUser(0.7,1.3)
            divHist.GetXaxis().SetTitleSize(0.065)
            divHist.SetTitle("")
            divHist.SetMarkerColor(ROOT.kBlack)

            divHist.Draw()


            fits["%i %i %i"%(iptC,ich,ibc)]=ROOT.TF1("%i %i %i"%(iptC,ich,ibc),"pol1",0,300)
            fits["%i %i %i"%(iptC,ich,ibc)].SetLineColor(ROOT.kRed)
            divHist.Fit(fits["%i %i %i"%(iptC,ich,ibc)],"NQ")

            fits["%i %i %i"%(iptC,ich,ibc)].Draw("SAME")


            info.SetTextSize(0.06)
            info.DrawLatexNDC(0.20,0.95,"%s %s %s"%(ptC,ch,bc))

            #CMS_lumi.CMS_lumi(miniCanvas,4,0)

        miniCanvas.SaveAs("%s/RatioGrid_%s_%s_%s.pdf"%(outDir,obs,sig,wid.replace('#','').replace('{','').replace('}','')))
        miniCanvas.SaveAs("%s/RatioGrid_%s_%s_%s.png"%(outDir,obs,sig,wid.replace('#','').replace('{','').replace('}','')))

    # use full dataset for comparison
    for wid,sig in [(a,c) for a in widList for c in sigList]:
        miniCanvas=ROOT.TCanvas()
        miniCanvas.SetCanvasSize(300*2, 200* int(1+ROOT.TMath.Ceil(len(obsList)/2)))
        miniCanvas.Clear()

        miniCanvas.cd()
        miniCanvas.Divide(2,int(ROOT.TMath.Ceil(len(obsList)/2))+1)

        info=ROOT.TLatex()
        info.SetTextSize(0.02)

        obsHists={}
        obsFits={}
        xFits={}

        i=1
        for obs in obsList :
            obsHists[obs]={}
            obsHists[obs]["1x"]=None
            obsHists[obs]["4x"]=None
            obsHists[obs]["Nom"]=None
            for (iptC,ich,ibc) in [(iptC,ich,ibc)
                    for iptC in range(len(ptChList))
                    for ich  in range(len(chList))
                    for ibc  in range(len(bcatList))] :

                ptC = ptChList[iptC]
                ch  = chList[ich]
                bc  = bcatList[ibc]

                wgtHist=wgtFin.Get("%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s"%(ptC,ch,bc,obs,wid,ptC,ch,bc,obs,wid,sig)).Clone()
                xHist=wgtFin.Get("%s%s%s_%s_1.0w/%s%s%s_%s_1.0w_%s"%(ptC,ch,bc,obs,ptC,ch,bc,obs,sig)).Clone()
                divHist=divFin.Get("%s%s%s_%s_%sw/%s%s%s_%s_%sw_%s widthx8"%(ptC,ch,bc,obs,divWid,ptC,ch,bc,obs,divWid,sig)).Clone()

                wgtHist.SetDirectory(0)
                xHist.SetDirectory(0)
                divHist.SetDirectory(0)

                if obsHists[obs]["1x"] is None :
                    obsHists[obs]["1x"]=wgtHist
                else :
                    obsHists[obs]["1x"].Add(wgtHist)

                if obsHists[obs]["4x"] is None :
                    obsHists[obs]["4x"]=divHist
                else :
                    obsHists[obs]["4x"].Add(divHist)

                if obsHists[obs]["Nom"] is None :
                    obsHists[obs]["Nom"]=xHist
                else :
                    obsHists[obs]["Nom"].Add(xHist)

            miniCanvas.cd(i)
            ROOT.gPad.SetGrid(1,1)
            ROOT.gPad.SetBottomMargin(ROOT.gPad.GetBottomMargin()*1.5)
            ROOT.gStyle.SetOptFit(0)
            totObsHist=obsHists[obs]["1x"]
            totObsHist.SetDirectory(0)
            totObsHist.Divide(obsHists[obs]["4x"])
            totObsHist.GetYaxis().SetTitle("1#timesSM Events / 8.0#times Events")
            totObsHist.GetYaxis().SetTitleSize(0.065)
            totObsHist.GetYaxis().SetRangeUser(0.7,1.3)
            totObsHist.GetXaxis().SetTitleSize(0.065)
            totObsHist.SetTitle("")
            totObsHist.SetMarkerColor(ROOT.kBlack)
            totObsHist.Draw()

            xObsHist=obsHists[obs]["Nom"]
            xObsHist.SetDirectory(0)
            xObsHist.Divide(obsHists[obs]["4x"])
            xObsHist.SetMarkerStyle(25)
            xObsHist.SetMarkerColor(ROOT.kRed)
            xObsHist.Draw("SAME")

            print '\n',obs
            obsFits[obs]=ROOT.TF1(obs,"pol1",0,300)
            obsFits[obs].SetLineColor(ROOT.kBlack)
            totObsHist.Fit(obsFits[obs])
            obsFits[obs].Draw("SAME")

            xFits[obs]=ROOT.TF1("1x"+obs,"pol1",0,300)
            xFits[obs].SetLineColor(ROOT.kRed)
            xObsHist.Fit(xFits[obs])
            xFits[obs].Draw("SAME")

            info.SetTextSize(0.06)
            info.DrawLatexNDC(0.20,0.95,"%s"%(obs))
            i+=1


        miniCanvas.SaveAs("%s/DistGrid_%s_%s.pdf"%(outDir,sig,wid.replace('#','').replace('{','').replace('}','')))
        miniCanvas.SaveAs("%s/DistGrid_%s_%s.png"%(outDir,sig,wid.replace('#','').replace('{','').replace('}','')))

    print '-'*50
    print 'Plots and summary ROOT file can be found in %s' % outDir
    print '-'*50


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

