import ROOT
import numpy as np

# methodList=[('TMVA_HighMJJ.root', 'vbf/TrainTree', 'BDT_VBF0HighMJJ'),
#             ('TMVA_LowMJJ.root',  'vbf/TrainTree', 'BDT_VBF0LowMJJ')
#             ]

methodList=[('TMVA_HighMJJ.root', 'vbf/TrainTree', 'BDT_VBF0HighMJJ'),
            ('TMVA_LowMJJ.root',  'vbf/TrainTree', 'BDT_VBF0LowMJJ'),
            ('TMVA_HighVPtLowMJJ.root',   'vbf/TrainTree', 'BDT_VBF0HighVPtLowMJJ'),
            ('TMVA_HighVPtHighMJJ.root',  'vbf/TrainTree', 'BDT_VBF0HighVPtHighMJJ'),
            ('TMVA_LowVPtHighMJJ.root',   'vbf/TrainTree', 'BDT_VBF0LowVPtHighMJJ')
            ]

fOutName='inverse_cdfs.root'

cdfInv=[]
origmva={}
finalmva={}
for f,t,m in methodList:

    print 'Opening',f,'for method=',m
    
    #build list of mva values for background
    origmva[m]={}
    finalmva[m]={}
    for c in 'SB':
        origmva[m][c]=ROOT.TH1F('origmva_%s_%s'%(m,c),';MVA;Events',50,-1,1)
        origmva[m][c].SetDirectory(0)
        origmva[m][c].Sumw2()
        finalmva[m][c]=ROOT.TH1F('finalmva_%s_%s'%(m,c),';MVA;Events',50,0,1)
        finalmva[m][c].SetDirectory(0)
        finalmva[m][c].Sumw2()

    mva=[]
    fIn=ROOT.TFile.Open(f)
    tree=fIn.Get(t)
    for i in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(i)
        mvaVal=getattr(tree,m)
        wgt=tree.weight
        if tree.classID==0: 
            origmva[m]['S'].Fill(mvaVal,wgt)
        else:
            mva.append(mvaVal)
            origmva[m]['B'].Fill(mvaVal,wgt)

    #compute quantiles
#    p=np.percentile(mva, range(1,99) )
	qArray = np.arange(5, 100, 5)
    p=np.percentile(mva, qArray )
    cdfInv.append(ROOT.TGraph())
    cdfInv[-1].SetName(m+'_cdfinv')
    cdfInv[-1].SetMarkerStyle(20)
    cdfInv[-1].SetPoint(0,-1,0)
    for i in xrange(0,len(p)):
#        cdfInv[-1].SetPoint(i+1,p[i],(i+1.)/100.)
        cdfInv[-1].SetPoint(i+1,p[i],(i+1.)/len(p))
    cdfInv[-1].SetPoint(len(p)+1,1,1)


    for i in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(i)
        mvaVal=cdfInv[-1].Eval(getattr(tree,m))
        wgt=tree.weight
        if tree.classID==0: 
            finalmva[m]['S'].Fill(mvaVal,wgt)
        else:
            finalmva[m]['B'].Fill(mvaVal,wgt)

    fIn.Close()

#save results
print 'Results saved in',fOutName
fOut=ROOT.TFile(fOutName,'RECREATE')
for gr in cdfInv: gr.Write()
for m in finalmva:
    for c in finalmva[m]:
        origmva[m][c].Write()
        finalmva[m][c].Write()
fOut.Close()
