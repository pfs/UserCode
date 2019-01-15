import ROOT
import optparse
import json
import sys
import os
from collections import OrderedDict
from TopLJets2015.TopAnalysis.Plot import fixExtremities
from TopLJets2015.TopAnalysis.gaussianFilterSmoother import GFSmoother

def showVariation(h,varList,warns,output):
    
    """a plot with the relative variations"""
    
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',500,300)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)

    frame=h.Clone('frame')
    frame.Reset('ICE')
    frame.Draw()
    frame.GetYaxis().SetTitle('Variation/Nominal')
    frame.GetYaxis().SetRangeUser(0.8,1.2)
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleSize(0.06)
    frame.GetXaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetLabelSize(0.05)

    leg=ROOT.TLegend(0.15,0.9,0.4,0.9-0.04*(len(varList)+1))
    leg.SetBorderSize(0)
    leg.SetTextSize(0.06)
    leg.SetFillStyle(0)  
    
    ratios=[]
    colors=[ROOT.kMagenta, ROOT.kAzure+7,ROOT.kMagenta+2, ROOT.kRed+1,ROOT.kMagenta-9,ROOT.kBlue-7]
    linestyles=[1,2,6]
    for i in xrange(0,len(varList)):
        ratios.append( varList[i].Clone('ratio%d'%i) )
        ratios[-1].Divide(h)
        ratios[-1].SetLineWidth(2)
        ratios[-1].SetLineStyle(linestyles[int(i/len(colors))])
        ratios[-1].SetLineColor(colors[i%len(colors)])
        ratios[-1].Draw('histsame')
        leg.AddEntry(ratios[-1],ratios[-1].GetTitle(),'l')
        
    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.06)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.93,'#bf{CMS} #it{simulation preliminary}')
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.93,'13 TeV')  

    if warns:
        tex.SetTextSize(0.04)
        tex.SetTextFont(52) 
        tex.SetTextColor(ROOT.kMagenta+3)
        for iw in range(len(warns)): 
            tex.DrawLatex(0.92,0.88-0.04*iw,warns[iw])
            
    c.Modified()
    c.Update()
    for ext in ['png','pdf']: c.SaveAs(output+'.'+ext)           


def getMirrored(var,nom,name) :

    """mirrors a shape"""

    mirrored_var=var.Clone(name)
    for xbin in range(nom.GetNbinsX()):
        mirrored_var.SetBinContent(xbin+1,2*nom.GetBinContent(xbin+1)-var.GetBinContent(xbin+1))
    return mirrored_var


def formatTemplate(h,hname,htit=None,ci=1,fs=0,norm=None):
    
    """ applies a bunch of standard formatting operations """

    h.SetDirectory(0)
    h.SetFillStyle(fs)
    h.SetLineColor(ci)
    h.SetMarkerColor(ci)
    h.SetName(hname)
    h.SetTitle(htit if htit else hname)
    if norm and norm>0: h.Scale(norm/h.Integral())
    
    #check bin contents
    for xbin in range(h.GetNbinsX()):
        if h.GetBinContent(xbin+1)>0: continue
        h.SetBinContent(xbin+1,1e-6)
    
    return h

def doFinalTemplateCheck(h,hup,hdn):

    """ checks if there are variations of yields or shape """

    keep=True
    report=''

    #check drastic changes of yields
    yields=h.Integral()
    yieldsUp=hup.Integral()
    yieldsDn=hdn.Integral()
    ivar=max(abs(yieldsUp/yields-1),abs(yieldsDn/yields-1))
    if ivar>0.3:
        report+='yields change by %d'%(100*ivar)
    
    #check negligible shape variations
    diffUp,diffDn=0,0
    for xbin in range(h.GetNbinsX()):
        diffUp += abs((hup.GetBinContent(xbin+1)/h.GetBinContent(xbin+1))*(yields/yieldsUp)-1)
        diffDn += abs((hdn.GetBinContent(xbin+1)/h.GetBinContent(xbin+1))*(yields/yieldsDn)-1)
    iShapeVar=max(diffUp,diffDn)
    if iShapeVar<0.01:
        report += 'shape variations change only by %3.1f%%\n'%(100*iShapeVar)

    #if there is <1% change can probably forget this one
    if iShapeVar<0.01 and ivar<0.01: keep=False
    
    return keep,report

def getBinByBinUncertainties(h):

    """loops over the bins of a template and builds the bin-by-bin uncertainties"""

    #init the up and down variations
    histos = [ [h.Clone('bin%dUp'%(xbin)),h.Clone('bin%dDown'%(xbin))] for xbin in xrange(1,h.GetNbinsX()+1) ]
    for hup,hdn in histos:
        formatTemplate(hup,hup.GetName())
        formatTemplate(hdn,hdn.GetName())

    for xbin in range(h.GetNbinsX()):
        val=h.GetBinContent(xbin+1)
        if val==0: continue
        valUnc=h.GetBinError(xbin+1)        
        histos[xbin][0].SetBinContent(xbin+1,val+valUnc)
        histos[xbin][1].SetBinContent(xbin+1,max(val-valUnc,1e-6))

    return histos


def getBinByBinUncertaintiesForSysts(h,hvars,byMax=False):
    
    """loops over the bins of a template and build the bin-by-bin stat unc associated to systs"""
    
    #init the up and down variations
    histos = [ [h.Clone('sysbin%dUp'%(xbin)),h.Clone('sysbin%dDown'%(xbin))] for xbin in xrange(1,h.GetNbinsX()+1) ]
    for hup,hdn in histos:
        formatTemplate(hup,hup.GetName())
        formatTemplate(hdn,hdn.GetName())

    for xbin in range(h.GetNbinsX()):

        #assign an uncertainty based on the max. found
        unc=0
        for ihvar in hvars:
            ival = ihvar.GetBinContent(xbin+1)
            if ival==0: continue
            relUnc = ihvar.GetBinError(xbin+1)/ival
            if byMax:
                unc = max(unc,relUnc)
            else:
                unc += relUnc**2

        #finalize computation of total unc. for this bin
        if not byMax: unc=ROOT.TMath.Sqrt(unc)

        #scale central yields up/down
        val=h.GetBinContent(xbin+1)
        histos[xbin][0].SetBinContent(xbin+1, (1+unc)*val)
        histos[xbin][1].SetBinContent(xbin+1, (1-unc)*val)

    return histos


def applySmoothing(h,ntoys=100,sigma=1):

    """ applies a Gaussian KDE smoothing to the histogram """

    print '[applySmoothing] with ',h.GetName(),'toys=',ntoys
    gfs=GFSmoother(h,ntoys=ntoys,sigma=sigma)
    h=gfs.smooth
    return h

def getUncertaintiesFromProjection(opt,fIn,d,proc_systs,hnom):

    """projects the _exp and _th 2D histograms to build the corresponding Up/Down uncertainty templates"""

    histos=[]
    errors=[]
    warns=[]

    #map all systematics available for projection
    allSysts={}
    for hname in ['{0}_exp'.format(d),'{0}_th'.format(d)]:
        h2d=fIn.Get('{0}/{0}_{1}'.format(hname,proc_systs['title']))
        try:
            for ybin in range(h2d.GetNbinsY()):
                allSysts[h2d.GetYaxis().GetBinLabel(ybin+1)]=(h2d,ybin+1)
        except:
            pass

    #project systematics
    for s in proc_systs['proj']:
        slist,norm,doEnvelope,smooth=proc_systs['proj'][s]
     
        try:            
            h2d,ybin=allSysts[ slist[0] ]
            varUp=h2d.ProjectionX('varup',ybin,ybin)
            fixExtremities(varUp) #add the overflow as for the main plot
            if smooth: applySmoothing(varUp)
        except:            
            errors.append('Failed to prepare %s'%slist[0])
            continue

        varDn=None

        #project down variation if available
        if not doEnvelope:
            try:
                h2d,ybin=allSysts[ slist[1] ]
                varDn=h2d.ProjectionX('vardn',ybin,ybin)
                fixExtremities(varDn) #add the overflow as for the main plot
                if smooth: applySmoothing(varDn)
            except:
                errors.append('Failed to prepare %s'%slist[1])
                continue

        #make an envelope if several are available
        elif len(slist)>2:

            try:
                varUp.Reset('ICE')
                varUp.Add(hnom) 
                varDn=hnom.Clone('vardn')
                for s_i in slist:

                    h2d,ybin=allSysts[ s_i ]
                    vartemp=h2d.ProjectionX('vartemp',ybin,ybin)
                    fixExtremities(vartemp) #add the overflow as for the main plot
                    if smooth: applySmoothing(vartemp)

                    #do max(var_i-nom,var_j-nom) per bin
                    for xbin in range(hnom.GetNbinsX()):
                        val         = vartemp.GetBinContent(xbin+1)
                        deltaVal    = val-hnom.GetBinContent(xbin+1)
                        if deltaVal>0:
                            curDeltaVal = varUp.GetBinContent(xbin+1)-hnom.GetBinContent(xbin+1)                 
                            if curDeltaVal<deltaVal:
                                varUp.SetBinContent(xbin+1,val)
                        else:
                            curDeltaVal = varDn.GetBinContent(xbin+1)-hnom.GetBinContent(xbin+1)                 
                            if curDeltaVal>deltaVal:
                                varDn.SetBinContent(xbin+1,val)
                            
                    #delete as no longer used
                    vartemp.Delete()

            except:
                errors.append('Failed to prepare %s'%s_i)
                continue


        #mirror the shape if down variation is not available yet
        if not varDn:
            varDn=getMirrored(varUp,hnom,'vardn')
            
        #check, format and add to the list if variation is to keep
        keep,checkReport=doFinalTemplateCheck(hnom,varUp,varDn)
        if len(checkReport):
            pfix='Keep but' if keep else 'Discard as'
            warns.append('%s %s for %s'%(pfix,checkReport.replace('\n',','),s))

        formatTemplate(varUp,s+'Up',norm=hnom.Integral() if norm else None)
        formatTemplate(varDn,s+'Down',norm=hnom.Integral() if norm else None)
        if keep:            
            histos.append(varUp)
            histos.append(varDn)

        #show plot with warnings found
        if opt.debug: 
            showVariation(hnom,
                          [varUp,varDn],
                          checkReport.split('\n'),
                          os.path.join(opt.output,'{0}_{1}_{2}'.format(hnom.GetTitle(),d,s)))

    #print errors/warnings found
    if len(errors) or len(warns):
        print '-'*50
        print 'Some errors/warnings found projecting syst templates for',proc_systs['title']
        print 'Please check the list below'
        for e in errors : print e
        for w in warns  : print w
        print '-'*50
            
    return histos


def getDirectUncertainties(opt,fIn,d,proc_systs,hnom):

    """ reads directly from the histograms and eventually normalizes, mirrors or takes an envelope """

    histos=[]
    errors=[]
    warns=[]

    #get systematics
    fIn.cd()
    for s in proc_systs['dir']:

        slist,norm,doEnvelope,smooth=proc_systs['dir'][s]

        varH=[]
        try:
            for s_i in slist:

                h=None

                #check if there is some placeholder to substitute
                if '{0}' in s_i:
                    hname=s_i.format(d)
                    hname='{0}/{0}_{1}'.format(hname,proc_systs["title"])
                    h=fIn.Get(hname)

                #otherwise look in the sub-directories
                else:                    
                    for k in fIn.GetListOfKeys():
                        if k.GetName()!=d: continue
                        for kk in k.ReadObj().GetListOfKeys():
                            if kk.GetName()!=d+'_'+s_i: continue
                            h=kk.ReadObj()
                            break

                if not h: continue
                h.GetName()
                if smooth: applySmoothing(varUp)
                varH.append(h)
        except:
            pass
        if len(varH)==0 : continue

        #mirror the shape if only one variation is available
        if len(varH)==1:
            varH.append( getMirrored(varH[0],hnom,'vardn') )

        #add to the histograms
        for i in xrange(0,2):
            pfix='Up' if i==0 else 'Down'
            formatTemplate(varH[i],s.format(d)+pfix,norm=hnom.Integral() if norm else None)
            histos.append(varH[i])

        #show plot with warnings found
        if opt.debug: 
            pfix=s.format(d) if '{0}' in s else s
            showVariation(hnom,
                          varH,
                          '',
                          os.path.join(opt.output,'{0}_{1}_{2}'.format(hnom.GetTitle(),d,pfix)))
        

    #if nothing found
    if len(histos)==0:
        raise Exception('No direct histograms were retrieved from {0} for {1} {2}'.format(fIn.GetName(),proc_systs["title"],d))

    return histos

            
def getTemplateHistos(opt,d,proc,proc_systs):

    """parses the input files for a specific distribution and builds the necessary templates for systematics"""
    
    histos=[]
    bbbUncHistos=[]
    systbbbUncHistos=[]
    
    #nominal histogram (use first found in inputs)
    for url in opt.input.split(','):        
        if not os.path.isfile(url) : continue
        fIn=ROOT.TFile.Open(url)
        h=fIn.Get('{0}/{0}_{1}'.format(d,proc_systs['title']))
        try:
            formatTemplate(h,'central',proc)
            if 'smooth' in proc_systs and proc_systs['smooth']: applySmoothing(h)
            histos.append(h)
            bbbUncHistos=getBinByBinUncertainties(h)
            break
        except:
            pass
        fIn.Close()

    if len(histos)==0:
        print 'Error: unable to find histogram',d,'for',proc
        return histos

    #associated experimental/weighted theory uncertainties (use first found in inputs)
    projFound=False
    dirFound=False
    for url in opt.input.split(','):        
        if not os.path.isfile(url) : continue
        fIn=ROOT.TFile.Open(url)
        try:

            if 'proj' in proc_systs and not projFound:
                histos+=getUncertaintiesFromProjection(opt,fIn,d,proc_systs,histos[0])
                projFound=True

            if 'dir' in proc_systs and not dirFound:
                dirHistos=getDirectUncertainties(opt,fIn,d,proc_systs,histos[0])
                histos+=dirHistos
                systbbbUncHistos=getBinByBinUncertaintiesForSysts(histos[0],dirHistos)
                dirFound=True

        except Exception as e:
            print e
            pass
        fIn.Close()


    #finalize bin-by bin uncertainties: if max. variation does not exceed threshold discard it
    for bbbColl,tag in [(bbbUncHistos,'cen'),(systbbbUncHistos,'sys')]:

        nHistos=len(histos)
        for xbin in range(0,len(bbbColl)):
            bbbUp,bbbDn=bbbColl[xbin]
            varUp=abs(bbbUp.GetBinContent(xbin+1)/histos[0].GetBinContent(xbin+1)-1)
            varDn=abs(bbbDn.GetBinContent(xbin+1)/histos[0].GetBinContent(xbin+1)-1)
            if max(varUp,varDn)<opt.bbbThr : continue
            bbbUp.SetName('{0}_{1}_{2}'.format(proc,d,bbbUp.GetName()))
            bbbDn.SetName('{0}_{1}_{2}'.format(proc,d,bbbDn.GetName()))
            histos.append(bbbUp)
            histos.append(bbbDn)


        #show overall sum for clarity, even if each bin has its own template
        if opt.debug and len(histos)>nHistos:
            totalbbb=[ histos[0].Clone('totalbbb%s'%(x)) for x in ['Up','Down'] ]
            for h in totalbbb: 
                h.SetTitle(h.GetName())
            nAdded=0
            for i in xrange(nHistos,len(histos),2):
                totalbbb[0].Add(histos[i])
                totalbbb[1].Add(histos[i+1])
                nAdded+=1        

            #subtract what has been added in excess
            for h in totalbbb:
                h.Add(histos[0],-nAdded)        
            showVariation(histos[0],
                          totalbbb,
                          None,
                          os.path.join(opt.output,'{0}_{1}_{2}bbb'.format(proc,d,tag)))
                      
    return histos

def prepareTemplateFile(opt,proc,proc_systs):

    """loops over the ROOT files and retrieves the relevant information """

    #read histos and variations
    histos={}
    for d in opt.distList.split(','): histos[d]=getTemplateHistos(opt,d,proc,proc_systs)

    #dump histograms to file
    url=os.path.join(opt.output,'templates_%s.root'%proc)
    fOut=ROOT.TFile.Open(url,'RECREATE')
    for d in histos:
        fOut.cd()
        outDir=fOut.mkdir(d)
        outDir.cd()
        for h in histos[d]: h.Write()
    fOut.Close()
    print 'Templates for',proc,'stored at',url

    return


def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='plotter files [%default]',  
                      default='test/analysis/top17010/0c522df/plots/plotter.root,test/analysis/top17010/0c522df/plots/syst_plotter.root,test/analysis/top17010/0c522df/plots/plotter_dydata.root',       
                      type='string')
    parser.add_option('-d', '--dist',          
                      dest='distList',       
                      help='distribution list [%default]',
                      default='emhighpt_mlb',
                      type='string')
    parser.add_option('-s', '--systs',          
                      dest='systs',       
                      help='description of the systematics [%default]',
                      default='test/analysis/top17010/systs_dict.json',
                      type='string')
    parser.add_option('--bbbThr',          
                      dest='bbbThr',
                      help='bin-by-bin relative uncertainty threshold to apply [%default]',
                      default=0.05,
                      type=float)
    parser.add_option('-o', '--out',          
                      dest='output',       
                      help='output directory',
                      default='test/analysis/top17010/0c522df/templates/',
                      type='string')
    parser.add_option('--debug',
                      dest='debug',
                      help='debug [%default]',
                      default=False,
                      action='store_true')
    (opt, args) = parser.parse_args()

    #prepare the output
    os.system('mkdir -p %s'%opt.output)

    #decode systematics map from the json file
    with open(opt.systs,'r') as cache:
        syst_dict=json.load( cache, encoding='utf-8', object_pairs_hook=OrderedDict ).items()

    #loop over processes
    for proc,proc_systs in syst_dict:
        prepareTemplateFile(opt,proc,proc_systs)



if __name__ == "__main__":
    sys.exit(main())
