import ROOT
import optparse
import sys
import os
from TopLJets2015.TopAnalysis.rounding import *


def doRinRout(nin,nout,ndyin,ndyout,ch1,ch2,refch='em'):

    """implementation of RinRout a la TOP-16-005"""    

    rinll    = nin[ch1][0]/nin[ch2][0]
    rinllUnc = rinll*ROOT.TMath.Sqrt(1./nin[ch1][0]+1./nin[ch2][0])

    kll    = ROOT.TMath.Sqrt(rinll)
    kllUnc = 0.5*kll*rinllUnc/rinll
    
    routin    = ndyout[ch1][0]/ndyin[ch1][0]
    routinUnc = routin*ROOT.TMath.Sqrt( (ndyout[ch1][1]/ndyout[ch1][0])**2 
                                        + (ndyin[ch1][1]/ndyin[ch1][0])**2 )

    nllout    = routin*(nin[ch1][0]-0.5*nin[refch][0]*kll)
    nlloutUnc = ROOT.TMath.Sqrt( (routinUnc*(nin[ch1][0]-0.5*nin[refch][0]*kll))**2
                                 + (routin**2)*nin[ch1][0]
                                 + nin[refch][0]*(0.5*routin*kll)**2
                                 + (0.5*routin*nin[refch][0]*kllUnc)**2 )
    
    sf    = nllout/ndyout[ch1][0]
    sfUnc = sf*ROOT.TMath.Sqrt( (nlloutUnc/nllout)**2+(ndyout[ch1][1]/ndyout[ch1][0])**2 )

    #add a small report with the numbers computed
    report  = 'R_{out/nin}(MC)=%s\n'%toLatexRounded(routin,routinUnc)
    report += 'k_{\\ell,\\ell}=%s\n'%toLatexRounded(kll,kllUnc)
    report += 'N_{\\rm in}(data)=%s\n'%nin[ch1][0]
    report += 'N_{\\rm out}=%s\n'%toLatexRounded(nllout,nlloutUnc)
    report += 'SF_{\\rm DY}=%s\n'%toLatexRounded(sf,sfUnc)

    return [sf,sfUnc,report]

def getPlots(name,fIn,chList=['ee','mm'],sig='DY',subtractBkg=False):

    """ gets the total and signal expectations in a region as well as the background-subtracted data """

    tag=''.join(chList)

    #sum up contributions from different channels to improve statistics
    dataH,sigH,bkgH=None,None,None
    for ch in chList:
        dirname='{0}{1}'.format(ch,name)               
        for key in fIn.Get(dirname).GetListOfKeys():
            objName=key.GetName()
            if 'Graph' in objName: continue
            if objName==dirname:
                if not dataH: dataH=key.ReadObj().Clone(name+sig+tag+'data')
                else        : dataH.Add(key.ReadObj())
            else:
                proc=objName.split('_')[-1]
                if proc==sig:
                    if not sigH: sigH=key.ReadObj().Clone(name+sig+tag)
                    else       : sigH.Add(key.ReadObj())
                else:
                    if not bkgH: bkgH=key.ReadObj().Clone(name+'bkg'+tag)
                    else       : bkgH.Add(key.ReadObj())

    #detach from current directory
    if dataH : dataH.SetDirectory(0)
    if bkgH  : bkgH.SetDirectory(0)
    if sigH  : sigH.SetDirectory(0)
    if dataH and bkgH and subtractBkg: 
        dataH.Add(bkgH,-1)
    
    #return results
    return {'sig':sigH,'data':dataH,'bkg':bkgH}
        
def transferDYestimate(h,tf):

    """apply the transfer factor and assign as shape envelope the effect of not applying it"""

    dy={}
    for key in ['nom','up','dn']:
        dy[key]=h.Clone(h.GetName()+'_'+key)
        dy[key].Reset('ICE')
        dy[key].SetDirectory(0)
        
    for xbin in xrange(1,h.GetNbinsX()+1):
        val_cr=h.GetBinContent(xbin)
        val_tf=tf.Eval(h.GetXaxis().GetBinCenter(xbin))

        dy['up'].SetBinContent(xbin,val_cr*(2*val_tf-1))
        dy['nom'].SetBinContent(xbin,val_cr*val_tf)
        dy['dn'].SetBinContent(xbin,val_cr)

    return dy


def estimateDY(srCat,crCat,dist,fIn,outDir):
    """sums up the ee and mm channels and compares with the results under the Z peak"""
                
    cr=getPlots(crCat+'_'+dist, fIn, chList=['zee','zmm'], sig='DY', subtractBkg=True)                
    sr=getPlots(srCat+'_'+dist, fIn, chList=['ee','mm'],   sig='DY', subtractBkg=True)

    
    #transfer factor (smoothed)
    tf_kfactor=cr['sig'].Integral()/sr['sig'].Integral()
    tf=sr['sig'].Clone('tf')
    tf.Divide(cr['sig'])
    tf.Scale(tf_kfactor)
    tf.SetLineColor(1)
    tf_func=ROOT.TF1('tf_func','[0]*x*x*x+[1]*x*x+[2]*x+[3]',tf.GetXaxis().GetXmin(),tf.GetXaxis().GetXmax())
    tf.Fit(tf_func,'MRQ0+')

    #data shape
    datady_sr=transferDYestimate(cr['data'],tf_func)
    for key in datady_sr: datady_sr[key].Scale(1./datady_sr[key].Integral())

    #closure test
    mcdy_sr=transferDYestimate(cr['sig'],tf_func)
    for key in mcdy_sr:
        mcdy_sr[key].Scale( sr['sig'].Integral()/mcdy_sr[key].Integral() )

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    drawOpt='hist'
    for key in mcdy_sr:
        mcdy_sr[key].GetYaxis().SetRangeUser(0,sr['sig'].GetMaximum()*1.3)
        mcdy_sr[key].SetLineWidth(2)
        mcdy_sr[key].SetFillStyle(0)
        mcdy_sr[key].SetLineColor(ROOT.kGray)
        mcdy_sr[key].Draw(drawOpt)
        drawOpt='histsame'
    sr['sig'].SetFillStyle(0)
    sr['sig'].SetLineColor(1)
    sr['sig'].SetMarkerStyle(20)
    sr['sig'].SetMarkerColor(1)
    sr['sig'].Draw('e1same')
    leg=ROOT.TLegend(0.15,0.9,0.4,0.8)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.SetFillStyle(0)
    leg.SetHeader(srCat)
    leg.AddEntry(mcdy_sr['nom'],'extrapolation from CR','l')
    leg.AddEntry(sr['sig'],'expectation in SR','l')
    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,'13 TeV')  
    for ext in ['png','pdf']:
        c.SaveAs(outDir+'/%s_dyclosure.%s'%(srCat+dist,ext))


    return datady_sr
            
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='plotter file [%default]',  
                      default='test/analysis/top17010/0c522df/plots/plotter.root',       
                      type='string')
    parser.add_option('-d', '--dist',          
                      dest='dist',       
                      help='distribution',
                      default='mlb',
                      type='string')
    parser.add_option('-o', '--out',          
                      dest='output',       
                      help='output directory',
                      default='test/analysis/top17010/0c522df/plots/',
                      type='string')
    (opt, args) = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    fIn=ROOT.TFile.Open(opt.input)

    #estimate the dy scale factors using Rin/Rout    
    nin,nout,ndyin,ndyout={},{},{},{}
    for ch in ['ee','mm','em']:

        chList=[ch]
        if ch!='em': chList += ['z'+ch]
        plots=getPlots('_mll', fIn, chList=chList, sig='DY')

        binMin = plots['data'].GetXaxis().FindBin(20)
        binMax = plots['data'].GetNbinsX()
        bin1   = plots['data'].GetXaxis().FindBin(76)
        bin2   = plots['data'].GetXaxis().FindBin(106)

        nin[ch]=(plots['data'].Integral(bin1,bin2),0.)
        nout[ch]=(plots['data'].Integral(binMin,binMax)-nin[ch][0],0.)

        dyInUnc=ROOT.Double(0)
        dyInCts=plots['sig'].IntegralAndError(bin1,bin2,dyInUnc)
        ndyin[ch]=(dyInCts,float(dyInUnc))

        dyOutUnc=ROOT.Double(0)
        dyOutCts=plots['sig'].IntegralAndError(binMin,binMax,dyOutUnc)
        ndyout[ch]=(dyOutCts-dyInCts,ROOT.TMath.Sqrt( float(dyOutUnc)**2 - float(dyInUnc)**2) )

    #a dict with the scale factors per channel
    dySF={'ee':doRinRout(nin,nout,ndyin,ndyout,'ee','mm'),
          'mm':doRinRout(nin,nout,ndyin,ndyout,'mm','ee')}
    sfProd=dySF['ee'][0]*dySF['mm'][0]
    sfProdUnc=ROOT.TMath.Sqrt((dySF['ee'][0]*dySF['mm'][1])**2+(dySF['mm'][0]*dySF['ee'][1])**2)        
    dySF['em']=[ROOT.TMath.Sqrt(sfProd),0.5*sfProdUnc/ROOT.TMath.Sqrt(sfProd)]
    dySF['em'].append('SF_{e\\mu}=%s\n'%toLatexRounded(dySF['em'][0],dySF['em'][1]))
    

    #estimate the DY shape in the different categories and scale to the final yields
    fOut=ROOT.TFile.Open(os.path.join(opt.output,'plotter_dydata.root'),'RECREATE')
    finalHistos={}
    for srCat in ['highpt', 'highpt1b', 'highpt2b', 'lowpt', 'lowpt1b', 'lowpt2b']:
        crCat = 'highpt1b' if 'highpt' in srCat else 'lowpt1b'
        datady_sr = estimateDY(srCat=srCat,crCat=crCat,dist=opt.dist,fIn=fIn,outDir=opt.output)
        
        for ch in ['ee','mm','em']:

            #scale yields for this channel
            dyExp=getPlots(srCat+'_'+opt.dist, fIn, chList=[ch], sig='DY')['sig'].Integral()
            dyExpUnc=dyExp*dySF[ch][1]
            dyExp*=dySF[ch][0]

            #scale the yields to match the expectations and write the shapes
            fOut.cd()
            for key in datady_sr:
                outName=ch+srCat+'_'+opt.dist
                if key!='nom' : outName += '_dyshape{0}'.format(key)
                outDir=fOut.mkdir(outName)
                outDir.cd()
                h=datady_sr[key].Clone(outName+'_DY (data)')
                h.SetTitle('DY (data)')
                h.Scale(dyExp/h.Integral())
                h.Write()
                finalHistos[outName]=h

    #add the inclusive shapes
    for ch in ['ee','mm','em']:

        for key in datady_sr.keys():
            incH=None
            for srCat in ['highpt','lowpt']:
                outName=ch+srCat+'_'+opt.dist
                if key!='nom' : outName += '_dyshape{0}'.format(key)
                ih=finalHistos[outName]
                if not incH:
                    incOutName=outName.replace(srCat,'')
                    incH=ih.Clone('{0}_DY (data)'.format(incOutName))
                    incH.SetDirectory(0)
                else:
                    incH.Add(ih)

            fOut.cd()
            outDir=fOut.mkdir(incOutName)
            outDir.cd()
            incH.Write()

    fOut.Close()
    
    #save the report
    with open('%s/dy_table.dat'%opt.output,'w') as cache:
        for key in dySF:
            cache.write('['+key+' channel]\n')
            cache.write(dySF[key][-1])
            cache.write('-'*50+'\n')
            
    print '{Shapes,report} have been saved in %s/{plotter_dydata.root,dy_report.dat}'%opt.output


if __name__ == "__main__":
    sys.exit(main())
