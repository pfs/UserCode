#!/usr/bin/env python

import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombinePdfs.morphing as morphing
import ROOT
import glob
import optparse

ROOT.gROOT.SetBatch(True)


parser=optparse.OptionParser("Get morphing workspace & datacards for TOP-16-019")
parser.add_option('-i', '--input', dest='input', help='input datacards', type='string',
        default='datacards/bias_0p2w_172.5/*/*.txt')
parser.add_option('-o', '--outdir', dest='outdir', help='output directory', type='string',
        default='datacards/bias_0p2w_172.5/')
parser.add_option(      '--biasWid', dest='biasWid', help='injected width value', type='string',
        default='0p2w')

opt,args=parser.parse_args()


cmb = ch.CombineHarvester()

# Read all the cards.
# CH stores metadata about each object (Observation, Process, Systematic),
# this is extracted from the card names with some regex
for card in glob.glob(opt.input):
    #cmb.QuickParseDatacard(card, """dc_mt_mass(?<MASS>\d+)_$CHANNEL.txt""")
    print card
    cmb.QuickParseDatacard(card, "datacards/bias_%s_172.5/$MASS/datacard.txt"%opt.biasWid)

# The process name needs to be common between bins
cmb.ForEachObj(lambda obj: obj.set_mass(obj.mass().replace('p','.').replace('w','')))

# We'll have three copies of the observation, one for each mass point.
# Filter all but one copy.
cmb.FilterObs(lambda obj: obj.mass() != '0.2')

# Create workspace to hold the morphing pdfs and the mass
w = ROOT.RooWorkspace('morph', 'morph')
width = w.factory('gamma[0.2,3.0]')

# BuildRooMorphing will dump a load of debug plots here
debug = ROOT.TFile('%s/debug.root'%opt.outdir, 'RECREATE')

cmb.SetVerbosity(3)
# Run for each bin,process combination
for b in cmb.bin_set():
    for p in cmb.cp().bin([b]).process_set():
        morphing.BuildRooMorphing(w, cmb, b, p, width, verbose=True, file=debug)

# Just to be safe
width.setConstant(True)

# Now the workspace is copied into the CH instance and the pdfs attached to the processes
# (this relies on us knowing that BuildRooMorphing will name the pdfs in a particular way)
cmb.AddWorkspace(w, True)
cmb.cp().ExtractPdfs(cmb, 'morph', '$BIN_$PROCESS_morph', '')

# Print the contents of the model
cmb.PrintAll()


# Adjust the rateParams a bit - we currently have three for each bin (one for each mass),
# but we only want one. Easiest to drop the existing ones completely and create new ones
cmb.syst_type(['rateParam'], False)
cmb.cp().AddSyst(cmb, 'norm_$BIN', 'rateParam', ch.SystMap()(1.00))

# Have to set the range by hand
for sys in cmb.cp().syst_type(['rateParam']).syst_name_set():
    cmb.GetParameter(sys).set_range(0.5, 1.5)


# Write out the cards, one per bin
writer = ch.CardWriter('%s/morphDatacard.txt'%opt.outdir,
                       '%s/morphShapes.root'%opt.outdir)
writer.SetVerbosity(1)
writer.WriteCards('', cmb)
