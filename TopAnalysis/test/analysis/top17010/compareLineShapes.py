import ROOT
from TopLJets2015.TopAnalysis.Plot import *

def getScenario(mt,gt):
    mask=int('0xffff',16)
    scenario=(int((gt-0.7)/0.01) & mask)
    scenario |= ((int((mt-169)/0.25) & mask) << 16)
    return scenario


simScenarios={    
    (169.5,1.23) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets_m1695.root',
    (171.5,1.28) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets_m1715.root',
    (172.5,1.31) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets.root',
    (173.5,1.34) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets_m1735.root',
    (175.5,1.39) : '/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/MC13TeV_2016_TTJets_m1755.root',
    }
rwScenarios=[ (169.5,1.23), (171.5,1.28), (172.5,1.0), (172.5,2.0), (173.5,1.34), (175.5,1.39) ]



def getDist(url,dname='mlb'):
    fIn=ROOT.TFile.Open(url)

    h=None
    if dname=='mlb':
        h=fIn.Get('em_mlb').Clone()
        h.Add(fIn.Get('ee_mlb'))
        h.Add(fIn.Get('mm_mlb'))
        h.SetDirectory(0)
        h.Scale(1./h.Integral())
    elif 'genmass' in dname:
        try:
            h=fIn.Get(dname).Clone()
            h.SetDirectory(0)
        except:
            pass
    fIn.Close()
    return h


simScenariosH={}
for key in simScenarios:
    simScenariosH[key]=getDist(simScenarios[key])
    simScenariosH[key].SetName('sim%d'%len(simScenariosH))
    simScenariosH[key].SetTitle('(%3.1f,%3.2f)'%key)


rwScenariosH={}
rwScenariosMassH={}
for key in rwScenarios:
    url='/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/scenario%d/MC13TeV_2016_TTJets.root'%getScenario(key[0],key[1])
    print url
    rwScenariosH[key]=getDist(url)
    rwScenariosH[key].SetName('rw%d'%len(rwScenariosH))
    rwScenariosH[key].SetTitle('(%3.1f,%3.2f)'%key)


    rwScenariosMassH[key]=(getDist(url,'genmass'),getDist(url,'rwgt_genmass'))
    rwScenariosMassH[key][0].SetName('mass%d'%len(rwScenariosH))
    rwScenariosMassH[key][1].SetName('rwmass%d'%len(rwScenariosH))
    rwScenariosMassH[key][0].SetTitle('(172.5,1.32)')
    rwScenariosMassH[key][1].SetTitle('(%3.1f,%3.2f) wgt'%key)
    norig=rwScenariosMassH[key][0].Integral()
    rwScenariosMassH[key][0].Scale(1./norig)
    rwScenariosMassH[key][1].Scale(1./norig)


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

p=Plot('mlbwidth_closure','13 TeV')
p.ratiorange=[0.958,1.042]
p.savelog=True
p.doPoissonErrorBars=False
p.range=(1e-4,0.014)
p.ytit='PDF'
p.add(h=simScenariosH[(172.5,1.31)].Clone(),title=simScenariosH[(172.5,1.31)].GetTitle(),color=1,isData=True,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.add(h=rwScenariosH[(172.5,1.0)].Clone(), title=rwScenariosH[(172.5,1.0)].GetTitle(),color=2,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.add(h=rwScenariosH[(172.5,2.0)].Clone(), title=rwScenariosH[(172.5,2.0)].GetTitle(),color=3,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.show('./',lumi=35900,noStack=True)


p=Plot('mlbwidth_m_closure','13 TeV')
p.ratiorange=[0.858,1.142]
p.savelog=True
p.doPoissonErrorBars=False
p.range=(1e-4,0.014)
p.ytit='PDF'
p.add(h=simScenariosH[(172.5,1.31)].Clone(),title=simScenariosH[(172.5,1.31)].GetTitle(),color=1,isData=True,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.add(h=rwScenariosH[(171.5,1.28)].Clone(), title=rwScenariosH[(171.5,1.28)].GetTitle(),color=2,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.add(h=rwScenariosH[(173.5,1.34)].Clone(), title=rwScenariosH[(173.5,1.34)].GetTitle(),color=3,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
p.show('./',lumi=35900,noStack=True)


for name,key in [('m1695',(169.5,1.23)),
                 ('m1715',(171.5,1.28)),
                 ('m1735',(173.5,1.34)),
                 ('m1755',(175.5,1.39)),
                 ]:
    p=Plot('%s_closure'%name,'13 TeV')
    p.ratiorange=[0.95,1.05]
    p.savelog=True
    p.range=(1e-4,0.02)
    p.ytit='PDF'
    p.doPoissonErrorBars=False
    pval=simScenariosH[key].Chi2Test(rwScenariosH[key],'WW')
    pval_ks=simScenariosH[key].KolmogorovTest(rwScenariosH[key],'')
    p.add(h=simScenariosH[key],title='Simulation',color=1,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=True)
    p.add(h=rwScenariosH[key], title='Reweighted',color=2,isData=True,spImpose=False,isSyst=False,doDivideByBinWidth=True)
    p.show('./',lumi=35900,noStack=True,extraText='(m_{t},#Gamma_{t})=(%3.1f,%3.2f)\\p-val(#chi^{2})=%3.3f\\p-val(KS)=%3.3f'%(key[0],key[1],pval,pval_ks))


p=Plot('genmass_closure','13 TeV')
p.ratiorange=[0.,2]
p.savelog=True
p.range=(1e-4,0.3)
p.ytit='PDF'
p.doPoissonErrorBars=False
key=rwScenariosMassH.keys()[0]
p.add(h=rwScenariosMassH[key][0], title=rwScenariosMassH[key][0].GetTitle(),color=1,isData=True,spImpose=False,isSyst=False,doDivideByBinWidth=False)
colors=[ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]
bwigners=[]
for name,key in [('m1695',(169.5,1.23)),
                 ('m1715',(171.5,1.28)),
                 ('m1735',(173.5,1.34)),
                 ('m1755',(175.5,1.39)),
                 ]:
    iplot=len(bwigners)
    ci=colors[iplot]
    bwigners.append( ROOT.TF1("bwigner"+name,
                     "[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))",
                     169.,176.) )
    bwigners[-1].SetParLimits(1,169,176)
    bwigners[-1].SetParLimits(2,1.0,3.0)
    bwigners[-1].SetLineColor(ci)
    rwScenariosMassH[key][1].Fit(bwigners[-1],'MR+','',169.1,175.9)
    p.add(h=rwScenariosMassH[key][1], title=rwScenariosMassH[key][1].GetTitle(),color=ci,isData=False,spImpose=False,isSyst=False,doDivideByBinWidth=False)

p.show('./',lumi=35900,noStack=True,overlay=bwigners)

