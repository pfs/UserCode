import ROOT
import sys
from array import array
import pandas as pd

df = pd.read_csv(sys.argv[1], sep=',', header = None, skiprows=1)
df = df.values.tolist()

bins=df[0][1:]+[df[1][-1]]
jerunc=ROOT.TH1F('jerunc',
                 ';Pseudo-rapidity;Relative SF unc.',
                 len(bins)-1,
                 array('f',bins))

fOut=ROOT.TFile.Open('jer_relunc.root','RECREATE')
for i in xrange(1,len(df)):
    name=df[i][0]
    jerunc.Reset('ICE')
    for xbin in range(jerunc.GetNbinsX()):
        jerunc.SetBinContent(xbin+1,df[i][xbin+1])
    jerunc.Clone().Write(name)
fOut.Close()
