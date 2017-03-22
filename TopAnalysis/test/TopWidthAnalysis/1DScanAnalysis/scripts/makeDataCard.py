#!/bin/python

import CombineHarvester.CombineTools.ch as ch
import ROOT
import os,sys
import optparse

ROOT.gROOT.SetBatch(True)

parser=optparse.OptionParser("Get initial datacards for TOP-16-019")
parser.add_option('--injWid', dest='injWid', help='injected width', type='string',
        default='')
parser.add_option('--injMas', dest='injMas', help='injected mass', type='string',
        default='')
parser.add_option('--ptcStr', dest='lbcStr', help='(csv) p_T categories', type='string',
        default='lowpt,highpt')
parser.add_option('--lfsStr', dest='lfsStr', help='(csv) lepton categories', type='string',
        default='EE,EM,MM')
parser.add_option('--btcStr', dest='btcStr', help='(csv) b multiplicity cats', type='string',
        default='1b,2b')

opt,args=parser.parse_args()

#####################################################################################
# set up the instance
#####################################################################################

injMas=opt.injMas
injWid=opt.injWid
if injMas != "" and injWid != "" :
    print "\n Injecting",injWid," width and ",injMas," mass."

signal = ['tbart','tW']
bkgs   = ['tbartV','tch','Multiboson','DY','W']

wids=[x*0.2 for x in range(1,16)]
wids=[("%.1fw"%x).replace('.','p') for x in wids]

lfs=opt.lfsStr.split(',')
ptC=opt.lbcStr.split(',')
nbC=opt.btcStr.split(',')
dist=['incmlb']

# Don't modify below this line
#####################################################################################
# produce categories / helpers
#####################################################################################
os.system('mkdir -p datacards')

cats=[]
iCat=0
for tpt,tlc,tnb,tds in [(a,b,c,d)
        for a in ptC
        for b in lfs
        for c in nbC
        for d in dist] :
    cats += [(iCat,'%s%s%s_%s'%(tpt,tlc,tnb,tds))]
    iCat+=1

channels=[b for (a,b) in cats]

#####################################################################################
# add processes
#####################################################################################
cb = ch.CombineHarvester()
cb.SetVerbosity(3)

cb.AddObservations(  wids, [''], [''], [''],                   cats)
cb.AddProcesses(     wids, [''], [''], [''], signal,           cats, True)
cb.AddProcesses(     wids, [''], [''], [''], ['Multiboson'],   cats, False)
cb.AddProcesses(     wids, [''], [''], [''], ['tbartV'],       cats, False)
cb.AddProcesses(     wids, [''], [''], [''], ['tch'],          cats, False)
cb.AddProcesses(     wids, [''], [''], [''], ['DY'],           cats, False)
cb.AddProcesses(     wids, [''], [''], [''], ['W'],            cats, False)

#####################################################################################
# implement systematics
#####################################################################################
#cb.ForEachProc(lambda x: x.set_shape())
#
#cb.cp().bin_id([0,1]).AddSyst(cb, 'eerate', 'lnN', ch.SystMap()(1.05))
#cb.cp().bin_id([7]).AddSyst(cb, 'lltrate', 'lnN', ch.SystMap()(1.05))
#cb.cp().bin_id(range(0,8)).AddSyst(cb, 'pileup', 'lnN', ch.SystMap()(1.05))
#cb.cp().bin_id([2,3,4,5]).AddSyst(cb, 'muEscale', 'lnN', ch.SystMap()(1.05))

# RATE SYSTEMATICS
cb.cp().bin_id([a for (a,b) in cats]).AddSyst(
        cb, 'lumi_13TeV', 'lnN', ch.SystMap()(1.062))
cb.cp().process(['DY']).AddSyst(
        cb, 'DYnorm_th', 'lnN', ch.SystMap()(1.30))
cb.cp().process(['W']).AddSyst(
        cb, 'Wnorm_th', 'lnN', ch.SystMap()(1.30))
cb.cp().process(['tW']).AddSyst(
        cb, 'tWnorm_th', 'lnN', ch.SystMap()(1.054))
cb.cp().process(['tch']).AddSyst(
        cb, 'tnorm_th', 'lnN', ch.SystMap()(1.044))
cb.cp().process(['Multiboson']).AddSyst(
        cb, 'VVnorm_th', 'lnN', ch.SystMap()(1.20))
cb.cp().process(['tbartV']).AddSyst(
        cb, 'tbartVnorm_th', 'lnN', ch.SystMap()(1.30))
cb.cp().AddSyst(
        cb, 'sel', 'lnN', ch.SystMap()(1.02))

cb.cp().process('tbart').AddSyst(
        cb, 'jer', 'shape', ch.SystMap()(1.00))
#cb.cp().AddSyst(
#        cb, 'jes', 'shape', ch.SystMap()(1.00))
#cb.cp().AddSyst(
#        cb, 'btag', 'shape', ch.SystMap()(1.00))
#cb.cp().AddSyst(
#        cb, 'PDF', 'shape', ch.SystMap()(1.00))
#cb.cp().AddSyst(
#        cb, 'pu', 'shape', ch.SystMap()(1.00))
#cb.cp().AddSyst(
#        cb, 'amcnloFxFx', 'shape', ch.SystMap()(1.00))
#cb.cp().AddSyst(
#        cb, 'Herwig', 'shape', ch.SystMap()(1.00))


#####################################################################################
# extract the shapes
#####################################################################################

# Note that order matters here, because of the way CombineHarvester
# Differentiates Obs shapes from other shapes.

# inject MC as data if desired
if injWid!="" and injMas!="" :
    cb.cp().process('data_obs').ExtractShapes(
            os.getcwd()+"/MCasData_shapes.root",
            "data_obs_$BIN_%s_%s"%(injWid,injMas),
            "data_obs_$BIN_%s_%s_%SYSTEMATIC");
else :
    cb.cp().data().ExtractShapes(
            os.getcwd()+"/MCasData_shapes.root",
            "$BIN/$PROCESS",
            "$BIN_$SYSTEMATIC/$PROCESS");

# now add in other shapes

#for sig in signal:
cb.cp().signals().ExtractShapes(
        os.getcwd()+"/shapes.root",
        "$BIN/$PROCESS$MASS",
        "$BIN_$SYSTEMATIC/$PROCESS$MASS");
#for bkg in bkgs:
cb.cp().backgrounds().ExtractShapes(
        os.getcwd()+"/shapes.root",
        "$BIN/$PROCESS",
        "$BIN_$SYSTEMATIC/$PROCESS");


# subdominant backgrounds cross section theoretical uncertainties
# cb.cp().process(['ttB4p']).AddSyst(cb, 'ttB4pxs', 'lnN', ch.SystMap()(1.0024))

#cb.WriteDatacard('datacards/test.txt')


#####################################################################################
# write the cards
#####################################################################################
cb.PrintAll()
writer = ch.CardWriter('$TAG/datacard.txt',
                       '$TAG/shapes.root')
#writer.SetWildcardMasses([])
writer.SetVerbosity(1)

for wid in wids:
    injDir="bias_%s_%s/"%(injWid,injMas) if (injWid!="" and injMas!="") else ""
    x = writer.WriteCards('datacards/%s%s'%(injDir,wid), cb.cp().mass([wid]))  # All cards combined
