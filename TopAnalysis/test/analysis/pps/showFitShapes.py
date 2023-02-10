import sys
import ROOT
import re
import os
from TopLJets2015.TopAnalysis.Plot import *

def showShapes(resultsDir,name,plotTitle,mass,boson,lumi,plotData=True,showSysts=True,showAllBkgs=True,plotpfix=''):

    """ show the shapes corresponding to a given analysis """

    colors=[ROOT.kGreen+1,ROOT.kAzure+3,ROOT.kOrange+2,ROOT.kGray,ROOT.kRed+2]
    
    shapeFiles=[]
    for f in os.listdir(resultsDir):
        if not '.root' in f : continue
        if len(f.split('_'))!=2: continue
        shapeFiles.append(os.path.join(resultsDir,f))
    print shapeFiles,boson

    for f in shapeFiles:

        fIn=ROOT.TFile.Open(f)

        #check consistency of the file with the boson to plot
        if 'shapes_22' in f and boson=='g':
            v='g'
            channel='#gamma'
        elif 'shapes_169' in f and boson=='z': 
            v='zmm'
            channel='Z(#mu#mu)'
        elif 'shapes_121' in f and boson=='z': 
            v='zee'
            channel='Z(ee)'
        else:
            continue
        extraText=plotTitle.format(channel)

        bkgH       = fIn.Get('bkg_%s'%v)
        fidsigH    = fIn.Get('fidsig_%s_m%d'%(v,mass))
        outfidsigH = fIn.Get('outfidsig_%s_m%d'%(v,mass))
        dataH      = fIn.Get('data_obs_%s'%(v))

        try:
            p=Plot('%s_%s_inc%s'%(name,v,plotpfix))
            p.plotformats = ['pdf','png','root']
            p.cmsLabel='#bf{CMS+Totem} #it{Preliminary}'
            p.xtit='Missing mass [GeV]'
            p.ytit='Events'
            if fidsigH and not plotData:
                p.add(fidsigH,            title='fiducial', color=ROOT.TColor.GetColor('#fdc086'), isData=False, spImpose=False, isSyst=False)
    
            if showAllBkgs:
                if outfidsigH and not plotData:
                    p.add(outfidsigH,     title='out fiducial',  color=ROOT.TColor.GetColor('#a6cee3'), isData=False, spImpose=False, isSyst=False)
                p.add(bkgH,               title='background',    color=ROOT.TColor.GetColor('#1f78b4'), isData=False, spImpose=False,    isSyst=False)
            else:
                allBkg=bkgH.Clone('allbkg')
                if outfidsigH and not plotData: 
                    allBkg.Add(outfidsigH)
                p.add(allBkg,               title='background',    color=ROOT.TColor.GetColor('#1f78b4'), isData=False, spImpose=False, isSyst=False)

            if plotData:
                dtitle='Data'
                p.add(dataH, title=dtitle,   color=1, isData=True, spImpose=False, isSyst=False)
                p.fitRatio=False

            if fidsigH:
                p.add(fidsigH.Clone(),    title='fiducial', color=ROOT.TColor.GetColor('#fdc086'), isData=False, spImpose=True,  isSyst=False)

            if showAllBkgs and outfidsigH:
                p.add(outfidsigH.Clone(), title='out fiducial',  color=ROOT.TColor.GetColor('#a6cee3'), isData=False, spImpose=True,  isSyst=False)
            
            p.ratiorange=[0.68,1.32]
            p.show('./',lumi*1000,extraText=extraText)


            if plotData:
                dataSubH=dataH.Clone(dataH.GetName()+'_sub')
                dataSubH.Add(bkgH,-1)
                p=Plot('%s_%s_sub_inc%s'%(name,v,plotpfix))
                p.plotformats = ['pdf','png','root']
                p.cmsLabel='#bf{CMS+Totem} #it{Preliminary}'
                p.frameMin=None
                p.frameMax=None
                p.doPoissonErrorBars=False
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                p.add(dataSubH, title='Data-bkg',   color=1, isData=True, spImpose=False, isSyst=False)
                if fidsigH:
                    p.add(fidsigH.Clone(),    title='fiducial', color=ROOT.TColor.GetColor('#fdc086'), isData=False, spImpose=True,  isSyst=False)
                if outfidsigH:
                    p.add(outfidsigH.Clone(), title='out fiducial',  color=ROOT.TColor.GetColor('#a6cee3'), isData=False, spImpose=True,  isSyst=False)
                p.show('./',lumi*1000,extraText=extraText,noRatio=True)

            if not showSysts: continue

            #background systs
            p=Plot('%s_%s_inc_bkgunc%s'%(name,v,plotpfix))
            p.plotformats = ['pdf','png','root']
            p.cmsLabel='#bf{CMS+Totem} #it{Preliminary}'
            p.noErrorsOnRatio=True
            p.doPoissonErrorBars=False
            p.xtit='Missing mass [GeV]'
            p.ytit='Events'
            p.add(bkgH, title='background', color=1, isData=True,spImpose=False, isSyst=False)
            ic=0
            for syst,title in [('bkgShapeEMUp',           'HP: up'),
                               ('bkgShapeEMDown',         'HP :down'),
                               ('bkgShapeSingleDiffUp',   'SD: pos arm mix'),
                               ('bkgShapeSingleDiffDown', 'SD: neg arm mix')]:
                h=fIn.Get('bkg_%s_%s'%(v,syst))
                p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                ic+=1
            p.ratiorange=[0.76,1.24]                
            p.show('./',lumi,noStack=True,extraText=plotTitle)

            #signal systs
            if fidsigH:
                p=Plot('%s_%s_sigunc%s'%(name,v,plotpfix))
                p.plotformats = ['pdf','png','root']
                p.cmsLabel='#bf{CMS+Totem} #it{Preliminary}'
                #fidsigH.Scale(1./5.)
                p.doPoissonErrorBars=False
                p.noErrorsOnRatio=True
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                p.add(fidsigH, title='signal', color=1, isData=True,spImpose=False, isSyst=False)
                ic=0
                for syst,title in [('ShapeEMUp',   'e#mu mix.'),
                                   ('CalibUp',   'Time dependence'),                     
                                   ('PPSEffUp',  'PPS efficiency'),
                                   ('PzModelUp', 'p_{z}(pp)')]:
                    h=fIn.Get('fidsig_%s_m%d_sig%s'%(v,mass,syst))
                    print ic,syst,title,h
                    p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                    ic+=1
                p.ratiorange=[0.5,1.46]
                p.show('./',lumi,noStack=True,extraText=plotTitle)

            #out fiducial signal systs
            if outfidsigH:
                p=Plot('%s_%s_outfidsigunc%s'%(name,v,plotpfix))
                p.plotformats = ['pdf','png','root']
                p.cmsLabel='#bf{CMS+Totem} #it{Preliminary}'
                #outfidsigH.Scale(1./5.)
                p.doPoissonErrorBars=False
                p.noErrorsOnRatio=True
                p.xtit='Missing mass [GeV]'
                p.ytit='Events'
                p.add(outfidsigH, title='out-fid. signal', color=1, isData=True,spImpose=False, isSyst=False)
                ic=0
                for syst,title in [('ShapeEMUp',   'e#mu mix.'),
                                   ('CalibUp',   'Time dependence'),
                                   ('PPSEffUp',  'PPS efficiency'),
                                   ('PzModelUp', 'p_{z}(pp)')]:
                    h=fIn.Get('outfidsig_%s_m%d_sig%s'%(v,mass,syst))
                    p.add(h, title=title, color=colors[ic], isData=False, spImpose=False, isSyst=False)
                    ic+=1

                p.ratiorange=[0.5,1.46]                    
                p.show('./',lumi,noStack=True,extraText=plotTitle)

        except Exception as e:
            print e
            pass

        fIn.Close()


def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    import optparse
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-b', '--boson',
                      dest='boson',   
                      default='zmm',
                      help='boson name [default: %default]')
    parser.add_option('-m', '--mass',
                      dest='mass',   
                      default=1000,
                      type=int,
                      help='mass [default: %default]')
    parser.add_option('-t', '--tag',
                      dest='tag', 
                      default='',
                      help='tag [default: %default]')
    parser.add_option('-u', '--unblind',
                      dest='unblind', 
                      default=False,
                      action='store_true',
                      help='unblind [default: %default]')
    parser.add_option('-s', '--showSysts',
                      dest='showSysts', 
                      default=False,
                      action='store_true',
                      help='showSysts [default: %default]')
    (options, args) = parser.parse_args()

    #from prepareOptimScanCards import OPTIMLIST
    #optimPt=int(re.search('optim_(\d+)',args[0]).group(1))
    #cuts=OPTIMLIST[optimPt][2].split(',')
    
    lumi=37.193
    bosonName=options.boson
    bosonName=bosonName.replace('mm','#mu#mu')
    bosonName=bosonName.replace('z','Z')
    if bosonName=='g' : 
        bosonName='#gamma'
        lumi=2.288
    plotpfix=options.tag
    mass=options.mass

    plotData=options.unblind
    if plotData:
        plotpfix += 'obs'

    title='pp#rightarrowpp%sX(%d)'%(bosonName,mass)
    if 'mm' in plotpfix:
        title+='\\multi-multi'
    if 'ms' in plotpfix:
        title+='\\multi-single'
    if 'sm' in plotpfix:
        title+='\\single-multi'
    if 'ss' in  plotpfix:
        title+='\\single-single'
    if 'exc' in plotpfix: title += ', N_{{jets}}=0'
    if '120' in plotpfix: title += ', 120 #murad'
    if '130' in plotpfix: title += ', 130 #murad'
    if '140' in plotpfix: title += ', 140 #murad'
    if '150' in plotpfix: title += ', 150 #murad'
    if 'lowpu' in plotpfix: title += ', N_{{vtx}}<20'
    if 'highpu' in plotpfix: title += ', N_{{vtx}}>20'

    showShapes(resultsDir=args[0],
               name='shapes',
               plotTitle=title,
               mass=mass,
               boson=options.boson,
               lumi=lumi,
               plotData=plotData,
               showSysts=options.showSysts,
               showAllBkgs=True, 
               plotpfix=plotpfix)

if __name__ == "__main__":
    main()
