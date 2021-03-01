import ROOT
from array import array
import os
import optparse
import sys

def saveExpectedBtagEff(opt):
    """loops over events and computes expected tagging efficiency per flavour"""

    #prepare data chain
    inputDir=opt.input
    data=ROOT.TChain('analysis/data') 
    input_list=[os.path.join(inputDir,f) for f in os.listdir( inputDir ) ]
    for i in xrange(0,min(5,len(input_list))):
        data.Add(input_list[i])
    print 'Projecting tagging efficiency from ',data.GetEntries(),' events'

    #prepare histograms
    ptBins = [0,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,800,1000]
    preTagH=ROOT.TH1F('preTagH',';p_{T} [GeV];',len(ptBins)-1,array('d',ptBins))
    preTagH.Sumw2()
    tagH=preTagH.Clone('tagH')


    #loop over taggers and flavours
    taggerList=[ x for x in opt.taggers.split(',')]
    flavConds=[('b',"abs(j_hadflav)==5"),
               ('c',"abs(j_hadflav)==4"),
               ('udsg','abs(j_hadflav)!=5 && abs(j_hadflav)!=4')]
    effgrs={}
    for x in taggerList:

        print '\t Starting',x
        tagger,op,opval=x.split(':')        
        cut='j_%s>%s'%('deepcsv',opval)
        
        #count number of tagged jets
        effgrs[(tagger,op)]={}
        for flav,cond in flavConds:
            print '\t computing for',flav,cond
            preTagH.Reset('ICE')
            tagH.Reset('ICE')
            data.Draw("j_pt >> preTagH",cond,'goff')
            data.Draw('j_pt >> tagH',cond + ' && ' + cut,'goff')
            effgrs[(tagger,op)][flav]=ROOT.TGraphAsymmErrors()
            effgrs[(tagger,op)][flav].Divide(tagH,preTagH)

    #jet response
    print 'Projecting out jet energy response per flavour'
    ptVar='j_pt' 
    genPtVar='g_pt[j_g]' 
    respgrs=[]
    respH=ROOT.TH2F('respH',';p_{T} [GeV];p_{T,rec}/p_{T,gen};',len(ptBins)-1,array('d',ptBins),1000,0.,5)
    respH.Sumw2() 
    for flav,cond in flavConds:
        respH.Reset('ICE') 
        data.Draw('{0}/{1}:{0} >> respH'.format(ptVar,genPtVar),cond + ' && j_g>=0','goff') 
        respgrs.append( ROOT.TGraphAsymmErrors(respH.ProfileX()) )  
        respgrs[-1].SetName('jresp_%s'%flav)

    #save in file
    fOut=ROOT.TFile.Open(opt.output,'RECREATE')
    for key in effgrs:
        fOut.cd()
        outDir=fOut.mkdir('%s_%s'%key)
        outDir.cd()
        for flav in effgrs[key]: effgrs[key][flav].Write(flav)
        fOut.cd()
    for gr in respgrs: 
        gr.Write()
    fOut.Close()

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',        dest='input',    help='input directory with files [%default]', default='/store/cmst3/group/top/RunIIReReco/f93b8d8/MC13TeV_2017_TTJets')
    parser.add_option('-o', '--out',       dest='output',   help='output file [%default]',           default='data/era2017/expectedBtagEff.root')
    parser.add_option('--taggers',   dest='taggers',  help='tagger:cut,tagger:cut [%default]',      
                      default='DeepCSV:loose:0.1522,DeepCSV:medium:0.4941,DeepCSV:tight:0.8001')
    (opt, args) = parser.parse_args()

    saveExpectedBtagEff(opt)

if __name__ == "__main__":
    sys.exit(main())
