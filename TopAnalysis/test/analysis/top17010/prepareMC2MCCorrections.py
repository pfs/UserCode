import ROOT
import os

def getRatio(num,den,name):
    r=ROOT.TGraphErrors()
    r.SetName(name)
    r.SetTitle(num.GetTitle())
    r.SetMarkerStyle(num.GetMarkerStyle())
    r.SetMarkerColor(num.GetMarkerColor())
    r.SetLineColor(num.GetLineColor())
    r.SetLineWidth(num.GetLineWidth())
    r.SetLineStyle(num.GetLineStyle())

    x,y1,y2=ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,num.GetN()):
        num.GetPoint(i,x,y1)
        ey1=num.GetErrorY(i)
        den.GetPoint(i,x,y2)
        ey2=den.GetErrorY(i)
        
        if y1==0: continue
        newy=float(y2)/float(y1)
        if newy==0: continue
        newyUnc=ROOT.TMath.Sqrt((float(y2)*ey1)**2+(float(y1)*ey2)**2)/(float(y1)**2)
                
        npt=r.GetN()
        r.SetPoint(npt,float(x),newy)
        r.SetPointError(npt,0,newyUnc)
    r.Sort()
    return r

PROCLIST=['MC13TeV_2016_TTJets',
          'MC13TeV_2016_TTJets_fsrdn',   'MC13TeV_2016_TTJets_fsrup',
          'MC13TeV_2016_TTJets_hdampup', 'MC13TeV_2016_TTJets_hdampdn',
          'MC13TeV_2016_TTJets_uedn',
          'MC13TeV_2016_TTJets_erdon'
          ]
PROCNOM=PROCLIST[0]
relResponses={}
for proc in PROCLIST:
    inDir   = '/store/cmst3/group/top/RunIIReReco/2016/0c522df/{0}'.format(proc)
    outFile = 'test/analysis/top17010/expectedBJEC_{0}.root'.format(proc)
    os.system('python scripts/saveExpectedBtagEff.py -i {0} -o {1} --taggers DeepCSV:medium:0.4941'.format(inDir,outFile))
    
    fIn=ROOT.TFile.Open(outFile)
    relResponses[proc]={}
    for item,name in [('jresp_b',          'Rb'),
                      ('jresp_c',          'Rc'),
                      ('jresp_udsg',       'Rudsg'),
                      ('DeepCSV_medium/b', 'eb'),
                      ('DeepCSV_medium/c', 'ec'),
                      ('DeepCSV_medium/udsg', 'eudsg'),
                      ]:
        if proc==PROCNOM: 
            relResponses[proc][name]=fIn.Get('%s'%(item))
            relResponses[proc][name].SetName(name)
        else:
            relResponses[proc][name]=getRatio(fIn.Get('%s'%(item)),relResponses[PROCNOM][name],name)

#save MC two MC corrections
outFile = 'test/analysis/top17010/mc2mc_corrections.root'
fOut=ROOT.TFile.Open(outFile,'RECREATE')
for proc in PROCLIST:
    if proc==PROCNOM:
        continue
    fOut.cd()
    outDir=fOut.mkdir(proc)
    outDir.cd()
    for key in relResponses[proc]:
        relResponses[proc][key].Write(key)
    fOut.cd()
fOut.Close()
