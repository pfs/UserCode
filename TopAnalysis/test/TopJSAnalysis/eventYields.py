import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt

debug = True

"""
steer the script
"""
def main():
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input', dest='input', help='input directory', default='plots/plotter.root', type='string')
    parser.add_option('-o', '--obs', dest='obs', help='observable used for event yields', default='l0eta', type='string')
    (opt, args) = parser.parse_args()

    #read lists of samples
    processes = ['t#bar{t}', 'Single top', 'W', 'DY', 'Multiboson', 't#bar{t}+V', 'QCD']
    channels = ['E', 'M', 'L']
    stages = ['1_1l', '2_1l4j', '3_1l4j2b', '4_1l4j2b2w']
    
    infile = ROOT.TFile.Open('%s'%(opt.input))
    
    with open('yields.tex', 'w') as tex:
        for channel in channels:
            tex.write('\hline \n')
            tex.write('%s + jets channel \\\\\n'%(channel))
            tex.write('\hline \n')
            totalMC = len(stages)*[0.]
            for process in processes:
                istage = 0
                tex.write('%s'%(process.replace('t#bar{t}', '$t\\bar{t}$')))
                for stage in stages:
                    hist = infile.Get('%s%s_%s/%s%s_%s_%s'%(channel, stage, opt.obs, channel, stage, opt.obs, process))
                    value = hist.Integral()
                    totalMC[istage] += value
                    istage += 1
                    tex.write(' & %.1f'%(value))
                tex.write(' \\\\\n')
            # total MC unc
            tex.write('\hline \n')
            tex.write('Total MC')
            istage = 0
            for stage in stages:
                tex.write(' & %.1f'%(totalMC[istage]))
                istage += 1
            tex.write(' \\\\\n')
            tex.write('($t\\bar{t}$ uncertainty)')
            for stage in stages:
                hist = infile.Get('%s%s_%s/totalmcunc'%(channel, stage, opt.obs))
                tex.write(' & $\pm$ %.1f'%(getError(hist, 1)))
            tex.write(' \\\\\n')
            # data
            tex.write('\hline \n')
            tex.write('Data')
            for stage in stages:
                hist = infile.Get('%s%s_%s/%s%s_%s'%(channel, stage, opt.obs, channel, stage, opt.obs))
                tex.write(' & %.1f'%(hist.Integral()))
            tex.write(' \\\\\n')
            for stage in stages:
                hist = infile.Get('%s%s_%s/%s%s_%s'%(channel, stage, opt.obs, channel, stage, opt.obs))
                tex.write(' & $\pm$ %.1f'%(getError(hist, 2)))
            tex.write(' \\\\\n')

# power = 2 for stat uncertainties
# power = 1 for correlated syst uncertainties
def getError(hist, power = 2):
    error    = 0.
    for i in range(hist.GetNbinsX()+2):
        error    += hist.GetBinError(i)**power
    error = error**(1./power)
    return error

"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

