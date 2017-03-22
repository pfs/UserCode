#!/bin/python

import sys
import pprint
import ROOT
import optparse

ROOT.gROOT.SetBatch()

parser=optparse.OptionParser("Get initial datacards for TOP-16-019")
parser.add_option('--injWid', dest='injWid', help='injected width', type='string',
        default='')
parser.add_option('--injMas', dest='injMas', help='injected mass', type='string',
        default='')
parser.add_option('--lbcStr', dest='lbcStr', help='(csv) p_T categories', type='string',
        default='lowpt,highpt')
parser.add_option('--lfsStr', dest='lfsStr', help='(csv) lepton categories', type='string',
        default='EE,EM,MM')
parser.add_option('--btcStr', dest='btcStr', help='(csv) b multiplicity cats', type='string',
        default='1b,2b')

opt,args=parser.parse_args()

#####################################################################################
# settings
#####################################################################################
signals=['tW','tbart']
backgrounds=['Multiboson','tch','tbartV','W','DY']

wids=[x*0.2 for x in range(1,16)]
wids=[("%.1fw"%x).replace('.','p') for x in wids]
#wids=['0p2w']

masses=["169.5","172.5","175.5"]

lfs=opt.lfsStr.split(',')
ptC=opt.ptcStr.split(',')
nbC=opt.btcStr.split(',')
dist=['incmlb']

injMasMap={
            "": "",
            "169.5": "_MtopDown",
            "172.5": "",
            "175.5": "_MtopUp"
        }

channels=['%s%s%s_%s'%(a,b,c,d)
    for a in ptC
    for b in lfs
    for c in nbC
    for d in dist]

fIn = ROOT.TFile("shapes.root")
fOut= ROOT.TFile("MCasData_shapes.root","RECREATE")

#####################################################################################
# function to get template from shapes file
#####################################################################################
def makeMCTemplate(channel,injWid,injMas) :
    hist=None
    newHistName='data_obs_%s_%s_%s'%(channel,injWid,injMas)
    print newHistName

    for sig in signals :
        sigHist=fIn.Get('%s%s/%s%s'%(channel,injMasMap[injMas],sig,injWid)).Clone()
        if hist is None :
            hist=sigHist.Clone(newHistName);
        else :
            hist.Add(sigHist.Clone())
        del sigHist

    for bkg in backgrounds:
        bkgHist=fIn.Get('%s/%s'%(channel,bkg)).Clone()
        if hist is None :
            hist=bkgHist.Clone(newHistName);
        else :
            hist.Add(bkgHist.Clone())
        del bkgHist

    if hist is None : print "ERROR: data_obs is None for channel,",tChannel,',',injWid,',',injMas
    return hist

#####################################################################################
# produce templates and write to file
#####################################################################################
for wid,mas,ch in [(w,m,c)
        for w in wids
        for m in masses
        for c in channels] :
    fOut.cd()
    newHist=makeMCTemplate(ch,wid,mas);
    print wid,mas,ch,newHist.Integral();
    newHist.Write();

fIn.Close()
fOut.Close()
