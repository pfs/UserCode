import ROOT
import sys
import array
import numpy as np


dists=['HighVPtLowMJJA_vbfmva',
       'HighVPtHighMJJA_vbfmva',
       'LowVPtHighMJJA_vbfmva']
fOutName='inverse_cdfs.root'
url=sys.argv[1]

cdfInv=[]
origmva={}
finalmva={}
for d in dists:

    m='BDT_VBF0'+d.split('_')[0][0:-1]
    print 'Opening',url,'for dist=',d,'for m=',m
    
    origmva[d]={}

    fIn=ROOT.TFile.Open(url)
    for key in fIn.Get(d).GetListOfKeys():
        kname=key.GetName()
        if d==kname : continue
        if 'Graph' in kname : continue
        h=key.ReadObj()
        ptype='B'
        if 'EWK' in kname: ptype='S'
        if not ptype in origmva[d] : 
            origmva[d][ptype]=h.Clone()
            origmva[d][ptype].Sumw2()
        else:
            origmva[d][ptype].Add(h)

    probSum = array.array('d', np.arange(0.01,0.99,0.01))
    q = array.array('d', [0.0]*len(probSum))
    origmva[d]['B'].GetQuantiles(len(probSum), q, probSum)

    cdfInv.append(ROOT.TGraph())
    cdfInv[-1].SetName(m+'_cdfinv')
    cdfInv[-1].SetMarkerStyle(20)
    cdfInv[-1].SetPoint(0,-1,0)
    for i in range(len(probSum)):
        cdfInv[-1].SetPoint(i+1,q[i],probSum[i])
    cdfInv[-1].SetPoint(len(probSum)+1,1,1)

    fIn.Close()

#save results
print 'Results saved in',fOutName
fOut=ROOT.TFile(fOutName,'RECREATE')
for gr in cdfInv: gr.Write()
for m in finalmva:
    for c in finalmva[m]:
        origmva[m][c].Write()
fOut.Close()
