import ROOT

baseDir="/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis"
sfs={}
for i in 'BCDEF':
    fa=ROOT.TFile.Open("%s/Data13TeV_2017%s_SinglePhoton.root"%(baseDir,i))
    fmm=ROOT.TFile.Open("%s/Data13TeV_2017%s_DoubleMuon.root"%(baseDir,i))
    ha=fa.Get('xangle_a')
    ha.Scale(1./ha.Integral())
    hmm=fmm.Get('xangle_mm')
    hmm.Scale(1./hmm.Integral())    
    ha.Divide(hmm)
    sfs[i]=[ ha.GetBinContent(xbin+1) for xbin in range(ha.GetNbinsX())]
    fa.Close()
    fmm.Close()
    
print sfs
