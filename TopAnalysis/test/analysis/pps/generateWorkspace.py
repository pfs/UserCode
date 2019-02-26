import ROOT
import os
import sys
import argparse

def shushRooFit():

    """stop all the messaging from RooFit"""
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Integration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)


def addPDFwithParameterUncertainties(w,paramList,pdfName,newPDFName,addParamTag=None):

    """ adds a parametrization to affect parameters by their uncertainties and creates a new PDF using those """

    cmd='EDIT::%s(%s,'%(newPDFName,pdfName)
    for param in paramList:
        val=w.var(param).getVal()
        unc=w.var(param).getError()
        paramName=param+addParamTag if addParamTag else param
        w.factory("expr::{0}_param('{1}+{2}*@0',{{{0}_nuis[0,-5,5]}})".format(paramName,val,unc) )
        cmd+='{0}={1}_param,'.format(param,paramName)
    cmd=cmd[0:-1]+')'
    w.factory(cmd)


def parametrizeSignal(w,args,url,mu,sigmass,cuts,debug,outDir):

    """paraterizes the signal and adds to the workspace"""    

    #import signal events
    fIn=ROOT.TFile.Open(url)
    data=fIn.Get('data')
    wgtds=ROOT.RooDataSet("wgtsig_data","wgtsig_data",data,args,cuts,'wgt')        
    ds=ROOT.RooDataSet("sig_data","sig_data",args,ROOT.RooFit.Import(data),ROOT.RooFit.Cut(cuts))
    fIn.Close()

    #fit signal shape
    w.factory("RooCBShape::sig_core(mmiss, sig_mean_core[%3.0f,%3.0f],sig_sigma_core[20,200],  sig_alpha_core[0,10], sig_n_core[0,10])"%(sigmass-100,sigmass+100))
    w.factory("RooCBShape:sig_turnon(mmiss,sig_mean_turnon[600,900],sig_sigma_turnon[100,300], sig_alpha_turnon[1,4],  sig_n_turnon[4,6])")
    w.factory("SUM::sig_base(sig_core_frac[0.7,1.0]*sig_core,sig_turnon)")
    w.pdf('sig_base').fitTo(ds)

    sigExp={}
    baseSigExp=wgtds.sumEntries()*0.0143*37500*mu
    xangleFracs={120:1.85,130:1.83,140:1.0,150:1.67}
    for xangle in [120,130,140,150]:
        sigExp[xangle]=baseSigExp*xangleFracs[xangle]
        addPDFwithParameterUncertainties(w,
                                         ['sig_mean_core','sig_sigma_core','sig_alpha_core','sig_n_core',
                                          'sig_mean_turnon','sig_sigma_turnon','sig_alpha_turnon','sig_n_turnon',
                                          'sig_core_frac'],
                                         'sig_base',
                                         'sig%d'%xangle,
                                         '_%d'%xangle)

    if not debug: return sigExp

    #show fit result
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    frame=w.var('mmiss').frame()
    ds.plotOn(frame)
    w.pdf('sig_base').plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlue))
    w.pdf('sig_base').plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray),ROOT.RooFit.LineStyle(2),ROOT.RooFit.Components("*turnon*"))
    frame.Draw()
    frame.GetXaxis().SetTitle('Missing mass [GeV]')
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.17,0.9,'120,130,140,150 #murad')
    c.Modified()
    c.Update()   
    c.RedrawAxis()
    out_url=os.path.basename(url)
    for ext in ['png','pdf']:
        c.SaveAs(os.path.join(outDir,out_url.replace('.root','.'+ext)))

    return sigExp

def parametrizeBackground(w,args,url,cuts,debug,outDir):

    """paraterizes the signal and adds to the workspace"""    

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(url,x) for x in os.listdir(url) if 'Data13TeV' in x]:
        if 'Photon' in f : continue
        if 'MuonEG' in f : continue
        data.AddFile(f)

    ds=ROOT.RooDataSet("bkg_data","bkg_data",data,args,cuts)

    w.factory('RooCBShape::bkg_decay(mmiss,bkg_mean_decay[1000,2000],bkg_sigma_decay[100,2000],bkg_alpha_decay[0,100],bkg_n_decay[0,10])')
    w.factory('RooCBShape::bkg_turnon(mmiss,bkg_mean_turnon[600,1000],bkg_sigma_turnon[100,300],bkg_alpha_turnon[1,4],bkg_n_turnon[4,6])')
    w.factory("SUM::bkg_base(bkg_decay_frac[0.0,1.0]*bkg_decay,bkg_turnon)")

    totalBkg={}
    for xangle in [120,130,140,150]:
        xangleCut='xangle==%d'%xangle
        rds=ds.reduce(xangleCut)
        rds.SetName('data_obs%d'%xangle)
        rds.SetTitle('data_obs%d'%xangle)
        getattr(w,'import')(rds)
        totalBkg[xangle]=rds.sumEntries()

        #fit background shape
        w.pdf('bkg_base').fitTo(rds)

        addPDFwithParameterUncertainties(w,
                                         ['bkg_mean_decay','bkg_sigma_decay','bkg_alpha_decay','bkg_n_decay',
                                          'bkg_mean_turnon','bkg_sigma_turnon','bkg_alpha_turnon','bkg_n_turnon',
                                          'bkg_decay_frac'],
                                         'bkg_base',
                                         'bkg%d'%xangle,
                                         '_%d'%xangle)

        if not debug:  continue

        #show fit result
        c=ROOT.TCanvas('c','c',500,500)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.02)
        c.SetTopMargin(0.05)
        c.SetBottomMargin(0.1)
        frame=w.var('mmiss').frame()
        rds.plotOn(frame)
        w.pdf('bkg%d'%xangle).plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlue))
        w.pdf('bkg%d'%xangle).plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray),ROOT.RooFit.LineStyle(2),ROOT.RooFit.Components("*turnon*"))

        frame.Draw()
        frame.GetXaxis().SetTitle('Missing mass [GeV]')
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
        tex.DrawLatex(0.17,0.9,'%d #murad'%xangle)
        c.Modified()
        c.Update()   
        c.RedrawAxis()
        for ext in ['png','pdf']:
            c.SaveAs(os.path.join(outDir,'bkgpdf_{0}.{1}'.format(xangle,ext)))
            
    return totalBkg

def writeDataCard(w,outDir,sigExp,bkgExp):

    """writes the datacard and the workspace"""

    os.system('mkdir -p %s'%outDir)
    wsURL=os.path.join(outDir,'workspace.root')
    w.writeToFile(wsURL)

    for xangle in [120,130,140,150]:
        cat='a%d'%xangle
        with open(os.path.join(outDir,'shapes-parametric-%dmurad.datacard.dat'%xangle),'w') as dc:
            dc.write('imax *\n')
            dc.write('jmax *\n')
            dc.write('kmax *\n')
            dc.write('-'*50+'\n')
            dc.write('shapes sig      %10s workspace.root w:$PROCESS%d\n'%(cat,xangle))
            dc.write('shapes bkg      %10s workspace.root w:$PROCESS%d\n'%(cat,xangle))
            dc.write('shapes data_obs %10s workspace.root w:data_obs%d\n'%(cat,xangle))
            dc.write('-'*50+'\n')
            dc.write('bin %s\n'%cat)
            dc.write('observation -1\n')
            dc.write('-'*50+'\n')
            dc.write('%15s %15s %15s\n'%('bin',cat,cat))
            dc.write('%15s %15s %15s\n'%('process','sig','bkg'))
            dc.write('%15s %15s %15s\n'%('process','0', '1'))
            dc.write('%15s %15s %15s\n'%('rate','%3.2f'%sigExp[xangle], '%3.2f'%bkgExp[xangle]))
            dc.write('-'*50+'\n')
            
            #float the background normalization as well as the signal
            dc.write('mu_bkg%d rateParam %s bkg 1\n'%(xangle,cat))

            args=w.allVars()
            it = args.createIterator()
            while it:
                try:                    
                    obj=it.Next()                    
                    vname = obj.GetName()
                    if 'nuis' in vname:
                        addToDataCard=True
                        if ('bkg' in vname or 'sig' in vname) and not str(xangle) in vname:
                            addToDataCard=False            
                        if addToDataCard:
                            dc.write('%40s param 0 1\n'%vname)                    
                except:
                    break
            del it



def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='plots/analysis',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--sig',
                        dest='sig',
                        default=1000,
                        type=int,
                        help='signal point')
    parser.add_argument('-c', '--cuts',
                        dest='cuts', 
                        default='(cat==121 || cat==169) && l1pt>30 && l2pt>20 && bosonpt>50',
                        help='Output directory [default: %default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='analysis',
                        help='Output directory [default: %default]')
    parser.add_argument('-d', '--debug',
                        dest='debug', 
                        default=False,
                        action='store_true',
                        help='Save debug plots [default: %default]')
    parser.add_argument('--mu',
                        dest='mu',
                        default=1.0,
                        type=float,
                        help='signal strength [default: %default]')
    opt=parser.parse_args(args)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    shushRooFit()

    print '[generateWorkspace] with sig=%d and cuts: (%s)'%(opt.sig,opt.cuts)

    #start the workspace
    w=ROOT.RooWorkspace('w')
    w.factory('mmiss[500,2500]')
    w.factory('l1pt[20,6500]')
    w.factory('l2pt[20,6500]')
    w.factory('bosonpt[20,6500]')
    w.factory('xangle[100,200]')
    w.factory('cat[100,200]')
    w.factory('wgt[-1.0e9,1.0e9]')
    varSet=w.allVars()

    #parametrize the signal
    print '[generateWorkspace] signal parametrization'
    sigExp=parametrizeSignal(w,
                             varSet,
                             os.path.join(opt.input,'MC13TeV_2017_PPZX_140urad_%d.root'%opt.sig),
                             opt.mu,
                             float(opt.sig),
                             opt.cuts,
                             opt.debug,                             
                             opt.output)

    #parametrize the background
    print '[generateWorkspace] background parametrization'
    bkgExp=parametrizeBackground(w,
                                 varSet,
                                 opt.input,
                                 opt.cuts,
                                 opt.debug,
                                 opt.output)

    #write summary in datacards
    print '[generateWorkspace] writing datacard'
    writeDataCard(w,opt.output,sigExp,bkgExp)
                          
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



