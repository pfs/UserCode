import os
import sys
import json
import ROOT

"""
"""
def main():

    impactsUrl=sys.argv[1].split(',')

    systsToMerge={
        'QCD multijets normalization':['MultiJetsNorm'],
        'QCD multijets shape':['MultiJetsShape'],
        '\\tbar QCD scale':['ttRenScale','ttCombScale','ttFactScale'],
        'Parton shower QCD scale':['ttPartonShower'],
        'W QCD scale':['wCombScale','wRenScale','wFactScale'],
        'tW QCD scale':['tWscale'],
        '\\cPqb-tagging':['BtagEff','CtagEff','LtagEff'],
        'Luminosity':['lumi_13TeV'],
        '$m_{\\cPqt}$':['Mtop'],
        '\\cPqt~\\pt':['topPt'],
        'Trigger efficiency':['Trigger'],
        'Lepton efficiency':['EleEfficiency','MuEfficiency'],
        'Lepton energy scale':['EleScale','MuScale'],
        'Other backgrounds':['DYnorm_th','VVnorm_th','tnorm_th','tWnorm_th','Wnorm_th'],
        'JES - pileup':['PileUpDataMC','PileUpPtHF','PileUpPtRef','PileUpPtEC1','PileUpPtEC2','PileUpPtBB'],
        'JES - relative':['RelativeJEREC2','RelativeJEREC1','RelativeStatHF','RelativeStatEC','RelativePtBB','RelativePtHF','RelativeJERHF','RelativePtEC1','RelativePtEC2','RelativeFSR'],
        'JES - \\pt':['Fragmentation','SinglePionHCAL','SinglePionECAL'],
        'JES - scale':['AbsoluteStat','AbsoluteMPFBias','AbsoluteFlavMap','AbsoluteScale'],
        'JES - flavour':['FlavorPureCharm','FlavorPureBottom','FlavorPureGluon','FlavorPureQuark'],
        'JES - time':['TimePt','TimeEta']
        }

    for url in impactsUrl:
        jsonF=open(sys.argv[1])
        results=json.load(jsonF,encoding='utf-8')
        jsonF.close()

        totalUnc=0
        systList={}
        for p in results['params']:
            systName,impact=p['name'],p['impact_r']

            for s in systsToMerge:
                for tag in systsToMerge[s]:
                    if tag in systName:
                        systName=s
            if not systName in systList: systList[systName]=0
            systList[systName] += impact**2
            totalUnc           += impact**2
        

        for s in systList:
            systList[s]=ROOT.TMath.Sqrt(systList[s])
            print '%30s & %3.3f \\\\'%(s,systList[s])
        totalUnc=ROOT.TMath.Sqrt(totalUnc)
        print '%30s & %3.3f \\\\'%('Total',totalUnc)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
