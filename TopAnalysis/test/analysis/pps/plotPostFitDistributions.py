import sys, os, optparse, copy, math
import ROOT
from array import array
import re
import pprint
from TopLJets2015.TopAnalysis.Plot import convertToPoissonErrorGr

## usage:
## python plotPostFitDistributions.py /eos/cms/store/cmst3/group/hintt/fits/datacards_2019-09-19_elPt25muPt20_ttbarBRFix/sphericity/allFlavors.card/fitDiagnosticsobs.root --outdir <directoryForPlots> (--inclusive)

mcStack=[('bkg'       , 'background',   'postfit',  ROOT.kGray) ,
         ('outfidsig' , 'non-fiducial', 'postfit',  ROOT.kAzure+6 ) ,         
         ('fidsig'    , 'fiducial' ,    'prefit',   ROOT.kRed+2), ]

def convertHisto(raw_h, target_h):
    
    """fills a copy of the target_h with the raw_h entries"""
    
    h=target_h.Clone( raw_h.GetName()+'_target')
    h.SetDirectory(0)
    h.Reset('ICE')

    for ibin in range(1,target_h.GetXaxis().GetNbins()+1):
        h.SetBinContent(ibin, raw_h.GetBinContent(ibin))
        h.SetBinError  (ibin, raw_h.GetBinError(ibin)  )
                
    return h
    
def plotShapes(categs_dict,categ_group,lumi,doPostfitOverPrefit=True):

    """does the plot of a collection of shapes"""

    title=categ_group['title']
    name=categ_group['name']

    c = ROOT.TCanvas("c","c",800*4,1200)
    c.SetTopMargin(0)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)

    pads=[]
    sub_pads=[]
    frames=[]
    data=[]
    postpre_ratios=[]        
    data2fitGr=[]
    relUncGr=[]
    leg=None
    leg_ratio=None

    #find the maximum in y-axis
    #maxY=0.
    #for (cat,xangle) in categ_group['group']:
    #    maxY=max(maxY,categs_dict[cat]['data'].GetMaximum()*1.5)

    for icat,(cat,xangle) in enumerate(categ_group['group']):

        c.cd()

        maxY=categs_dict[cat]['data'].GetMaximum()*1.5

        ix=icat*0.25
        pads.append( ROOT.TPad("p{}".format(icat), "p{}".format(icat), ix, 0., ix+0.25, 1.) )
        pads[icat].SetTopMargin(0)
        pads[icat].SetLeftMargin(0)
        pads[icat].SetRightMargin(0)
        pads[icat].SetBottomMargin(0)
        pads[icat].Draw()
        pads[icat].cd()

        #main pad: data and post-fit expectations
        lm=0.15 if icat==0 else 0.1
        rm=0.05
        sub_pads.append([])
        if doPostfitOverPrefit:
            sub_pads[icat].append( ROOT.TPad("p{}_1".format(icat), "p{}_1".format(icat),0., 0.5, 1, 1.) )
        else:
            sub_pads[icat].append( ROOT.TPad("p{}_1".format(icat), "p{}_1".format(icat),0., 0.4, 1, 1.) )
        sub_pads[icat][0].SetTopMargin(0.1)
        sub_pads[icat][0].SetBottomMargin(0.03)
        sub_pads[icat][0].SetLeftMargin(lm)
        sub_pads[icat][0].SetRightMargin(rm)
        sub_pads[icat][0].SetGridx()        
        sub_pads[icat][0].SetGridy()
        sub_pads[icat][0].SetFillStyle(0)
        sub_pads[icat][0].Draw()
        sub_pads[icat][0].cd()

        postfit=categs_dict[cat]['postfit']
        prefit=categs_dict[cat]['prefit']

        for m in [800,1600]:
            ibin=categs_dict[cat]['data'].GetXaxis().FindBin(m-50)
            jbin=categs_dict[cat]['data'].GetXaxis().FindBin(m+50)
            nobs=categs_dict[cat]['data'].Integral(ibin,jbin)
            print(m,ibin,jbin,nobs,icat,name,cat,categs_dict[cat]['data'].Integral())


        data.append(  convertToPoissonErrorGr( categs_dict[cat]['data'] ) )
        for i in range(data[icat].GetN()):
            data[icat].SetPointEXhigh(i,0.)
            data[icat].SetPointEXlow(i,0.)
        data[icat].SetMarkerStyle(20)
        data[icat].SetMarkerSize(1.)
        data[icat].SetMarkerColor(1)
        data[icat].SetLineColor(1)
        data[icat].SetLineWidth(3)

        frames.append([])
        frames[icat].append( postfit['total'].Clone('frame') )
        frames[icat][0].Reset('ICE')
        frames[icat][0].GetYaxis().SetTitle('Events')
        frame_yfont=0.06 #if icat==0 else 0.
        if not doPostfitOverPrefit: frame_yfont=0.052
        frame_ytitlefont=frame_yfont if icat==0 else 0.

        frames[icat][0].GetYaxis().SetTitleOffset(1.2 if doPostfitOverPrefit else 1.4)
        frames[icat][0].GetYaxis().SetNdivisions(505)
        frames[icat][0].GetYaxis().SetTitleSize(frame_ytitlefont)
        frames[icat][0].GetYaxis().SetLabelSize(frame_yfont)
        frames[icat][0].GetXaxis().SetTitleSize(0)
        frames[icat][0].GetXaxis().SetLabelSize(0)
        frames[icat][0].GetYaxis().SetRangeUser(0,maxY)
        frames[icat][0].GetXaxis().SetNdivisions(505)
        frames[icat][0].Draw()

        if icat==3:
            leg = ROOT.TLegend(0.6,0.5,0.94,0.8)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.06 if doPostfitOverPrefit else 0.045)
            leg.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
            leg.AddEntry(data[icat],'Data','ep')

        for proc,ptitle,fit_type,ci in mcStack:
            ispostfit=(fit_type=='postfit')
            h=postfit[proc] if ispostfit else prefit[proc]
            h.SetTitle(ptitle)
            if ispostfit:
                h.SetLineColor(1)
                h.SetFillStyle(1001)
                h.SetFillColor(ci)
                if icat==3:
                    leg.AddEntry(h,ptitle,'f')
            else:
                h.SetLineWidth(2)
                h.SetFillStyle(0)
                h.SetLineColor(ci)
                if icat==3:
                    leg.AddEntry(h,ptitle, 'l')

            h.Draw('histsame')


        if not leg is None:
            leg.Draw()

        #label this plot
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.07)
        tex.SetNDC()
        try:
            sub_categ_txt='Angle=%d#murad'%xangle
        except:
            sub_categ_txt=xangle
        if icat==0:
            tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
            tex.DrawLatex(0.14,0.96,'#bf{CMS-Totem} #it{Preliminary}')            
            tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
            tex.SetTextSize(0.06)
            for it,txt in enumerate(title):
                tex.DrawLatex(0.16,0.86-it*0.06,'%s'%txt)
        tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        tex.SetTextSize(0.06)
        tex.DrawLatex(0.93,0.86,sub_categ_txt)
        if icat==3:
            tex.SetTextSize(0.07)
            tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
            tex.DrawLatex(0.96,0.96,'%3.1f fb^{-1} (#sqrt{s} = 13 TeV)'%(lumi*1e-3))

        data[icat].Draw('PZ0')
        sub_pads[icat][0].RedrawAxis()

        #post/pre-fit expectations for the background
        if doPostfitOverPrefit:
            pads[icat].cd()
            sub_pads[icat].append( ROOT.TPad("p{}_1".format(icat), "p{}_1".format(icat),0., 0.3, 1, 0.49) )
            sub_pads[icat][1].SetTopMargin(0.04)
            sub_pads[icat][1].SetBottomMargin(0.04)
            sub_pads[icat][1].SetLeftMargin(lm)
            sub_pads[icat][1].SetRightMargin(rm)
            sub_pads[icat][1].SetGridx()
            sub_pads[icat][1].SetGridy()
            sub_pads[icat][1].Draw()
            sub_pads[icat][1].cd()

            frames[icat].append( frames[icat][-1].Clone('frame_{}_2'.format(icat) ) )
            frames[icat][1].GetYaxis().SetTitle('Postfit/Prefit')
            frames[icat][1].GetYaxis().SetTitleOffset(0.4)
            frame_yfont=0.18 if icat==0 else 0.
            frames[icat][1].GetYaxis().SetTitleSize(frame_yfont)
            frames[icat][1].GetYaxis().SetLabelSize(frame_yfont)
            frames[icat][1].GetYaxis().SetNdivisions(5)
            frames[icat][1].GetXaxis().SetTitleSize(0.0)
            frames[icat][1].GetXaxis().SetLabelSize(0.0)
            frames[icat][1].GetYaxis().SetRangeUser(0.66,1.304)
            frames[icat][1].GetXaxis().SetNdivisions(505)
            frames[icat][1].Draw()
            
            postpre_ratios.append( postfit['bkg'].Clone('bkg_ratio_{}'.format(icat) ) )
            postpre_ratios[-1].Divide( prefit['bkg'] )
            postpre_ratios[-1].SetLineWidth(2)
            postpre_ratios[-1].SetFillStyle(0)
            postpre_ratios[-1].Draw('histsame')
            for ibin in range(1,postfit['bkg'].GetNbinsX()+1):
                pos=postfit['bkg'].GetBinContent(ibin)
                pre=prefit['bkg'].GetBinContent(ibin)
                if pos<0.1 and pre<0.1:
                    postpre_ratios[icat].SetBinContent(ibin,1)

            sub_pads[icat][1].RedrawAxis()


        #data/MC
        pads[icat].cd()
        if doPostfitOverPrefit:
            sub_pads[icat].append( ROOT.TPad("p{}_2".format(icat), "p{}_2".format(icat),0., 0., 1, 0.29) )
        else:
            sub_pads[icat].append( ROOT.TPad("p{}_2".format(icat), "p{}_2".format(icat),0., 0., 1, 0.39) )
        sub_pads[icat][-1].SetTopMargin(0.04)
        sub_pads[icat][-1].SetBottomMargin(0.4)
        sub_pads[icat][-1].SetLeftMargin(lm)
        sub_pads[icat][-1].SetRightMargin(rm)
        sub_pads[icat][-1].SetGridx()
        sub_pads[icat][-1].SetGridy()
        sub_pads[icat][-1].Draw()
        sub_pads[icat][-1].cd()

        frame_xfont=0.12
        frame_yfont=0.12 if icat==0 else 0
        if not doPostfitOverPrefit:
            frame_xfont=0.085
            frame_yfont=0.085 if icat==0 else 0
        frames[icat].append( frames[icat][-1].Clone('frame_{}_3'.format(icat) ) )
        frames[icat][-1].GetYaxis().SetTitle('Data/Pred.')
        frames[icat][-1].GetYaxis().SetTitleOffset(0.6 if doPostfitOverPrefit else 0.9)
        frames[icat][-1].GetYaxis().SetTitleSize(frame_yfont)
        frames[icat][-1].GetYaxis().SetLabelSize(frame_yfont)
        frames[icat][-1].GetYaxis().SetNdivisions(5)
        frames[icat][-1].GetXaxis().SetTitleSize(frame_xfont)
        frames[icat][-1].GetXaxis().SetLabelOffset(0.01)
        frames[icat][-1].GetYaxis().SetLabelOffset(0.03)
        frames[icat][-1].GetXaxis().SetTitleOffset(1.2)
        frames[icat][-1].GetXaxis().SetLabelSize(frame_xfont)
        frames[icat][-1].GetXaxis().SetTitle('Missing mass [GeV]')
        frames[icat][-1].GetXaxis().SetNdivisions(505)
        frames[icat][-1].GetYaxis().SetRangeUser(0.66,1.304)
        frames[icat][-1].Draw()

        relUnc=postfit['total'].Clone()
        data2fitGr.append( data[icat].Clone('data2fit') )
        data2fitGr[icat].Set(0)
        for i in range(data[icat].GetN()):

            den=postfit['total'].GetBinContent(i+1)
            denUnc=postfit['total'].GetBinError(i+1)
            if float(den)>0 : 

                relUnc.SetBinContent(i+1,1)
                relUnc.SetBinError(i+1,denUnc/den)

                x,y=ROOT.Double(0),ROOT.Double(0)

                data[icat].GetPoint(i,x,y)
                eyhi=data[icat].GetErrorYhigh(i)
                eylo=data[icat].GetErrorYlow(i)

                np=data2fitGr[icat].GetN()
                data2fitGr[icat].SetPoint(np,x,y/den)
                data2fitGr[icat].SetPointEYhigh(np,eyhi/den)
                data2fitGr[icat].SetPointEYlow(np,eylo/den)

            else:
                relUnc.SetBinContent(i+1,0)
                relUnc.SetBinError(i+1,0)

        relUncGr.append( ROOT.TGraphErrors(relUnc) )
        relUnc.Delete()
        relUncGr[icat].SetFillStyle(3244)
        relUncGr[icat].SetFillColor(1)
        relUncGr[icat].SetMarkerSize(2.0)
        relUncGr[icat].Draw('E20')

        data2fitGr[icat].Draw('PZ0')

        if icat==3:
            if doPostfitOverPrefit:
                leg_ratio = ROOT.TLegend(0.08,0.14,0.5,0.23)
            else:
                leg_ratio = ROOT.TLegend(0.12,0.43,0.6,0.5)
            leg_ratio.SetTextAlign(ROOT.kHAlignCenter+ROOT.kVAlignCenter)
            leg_ratio.SetNColumns(2)
            leg_ratio.SetFillStyle(0)
            leg_ratio.SetBorderSize(0)
            leg_ratio.SetTextSize(0.1 if doPostfitOverPrefit else 0.08)
            leg_ratio.AddEntry(data2fitGr[icat],'Data','ep')
            leg_ratio.AddEntry(relUncGr[icat],'Postfit unc.','f')
            leg_ratio.Draw()

        sub_pads[icat][-1].RedrawAxis()

    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        name=name if doPostfitOverPrefit else name+'_simple'
        c.SaveAs('{}.{}'.format(name,ext))



def parseCategoriesFromDatacard(url,sources_url):

    """parses the necessary information from the datacard"""

    categ_dict={}
    print('[parseCategoriesFromDatacard] with {}'.format(url))

    #start by identifying the categories
    cats=[]
    with open(url,'r') as fin:
        for l in fin.readlines():            
            if l.find('Combination of')!=0 : continue
            cats=dict([ tuple(x.split('=')) for x in l.split() if x.find('cat')==0])            
            break

    #identify the cuts and the base histogram in each category
    base_dir=os.path.dirname(url)
    for c,f in cats.items():

        categ_dict[c]={}

        source_card=os.path.join(sources_url,os.path.basename(f))
        with open( source_card,'r' ) as fin:
            for l in fin.readlines():

                #decode preselection
                if l.find('# Namespace(')==0 :

                    #convert the options string to a dict
                    opts=dict([tuple(x.split('=',1)) for x in l.split(',') if '=' in x])

                    #decode preselection
                    presel=opts[' presel']
                    try:
                        nvtx=int(re.findall("nvtx>=([0-9]+)",presel)[0])
                    except:
                        try:
                            nvtx=-int(re.findall("nvtx<([0-9]+)",presel)[0])
                        except:
                            nvtx=None
                    categ_dict[c]={'cat':int(re.findall("cat==([0-9]+)",presel)[0]),
                                   'protonCat':int(re.findall("protonCat==([0-9]+)",presel)[0]),
                                   'xangle':int(re.findall("xangle==([0-9]+)",presel)[0]),
                                   'nvtx':nvtx}

                #get baseline histogram
                elif l.find('shapes data_obs')==0:
                    rshapes=l.split()[3]
                    fIn=ROOT.TFile.Open(rshapes)
                    for k in fIn.GetListOfKeys():
                        if not 'data_obs' in k.GetName() : continue
                        categ_dict[c]['data']=k.ReadObj()
                        categ_dict[c]['data'].SetDirectory(0)
                        categ_dict[c]['shape']=categ_dict[c]['data'].Clone('{}_template'.format(c))
                        categ_dict[c]['shape'].Reset('ICE')
                        categ_dict[c]['shape'].SetDirectory(0)
                        break
                    fIn.Close()

    print('Parsed info on {} categories'.format(len(categ_dict)))
    return categ_dict

def fillFitShapes(url,outDir,categ_dict,xtitle='Missing mass [GeV]'):

    """parses the fitDiagnostics file and reads the histograms"""

    print('[fillFitShapes with {} categories expected'.format(len(categ_dict)))

    inF = ROOT.TFile(url)
    print(url)
    for shapes_dir in ['shapes_prefit','shapes_fit_s']:

        fit_type='prefit' if 'prefit' in shapes_dir else 'postfit'
        for c in categ_dict: categ_dict[c][fit_type]={}

        #iterate over categories, read and convert shapes
        for ch in inF.Get(shapes_dir).GetListOfKeys():
            chName=ch.GetName()
            chDir=ch.ReadObj()
            print(chDir.GetName())
            for proc in chDir.GetListOfKeys():
                pname=proc.GetName()
                if pname=='total_covar' : continue
                raw_h=proc.ReadObj()
                if not raw_h.InheritsFrom('TH1'): continue
                categ_dict[chName][fit_type][pname]=convertHisto(raw_h=proc.ReadObj(),
                                                                 target_h=categ_dict[chName]['shape'])

    return categ_dict

def groupCategories(categ_dict):

    """returns the categories to plot in the same canvas"""

    categ_groups={}
    for c in categ_dict:
        boson=categ_dict[c]['cat']
        nvtx=categ_dict[c]['nvtx']
        protonCat=categ_dict[c]['protonCat']
        xangle=categ_dict[c]['xangle']

        key=(boson,protonCat,nvtx)
        if not key in categ_groups:
            title=[]
            if boson==121:    title.append('pp#rightarrow ppZeeX')
            if boson==169:    title.append('pp#rightarrow ppZ#mu#muX')
            if boson==22:     title.append('pp#rightarrow pp#gammaX')
            if protonCat==1 : title.append('multi-multi')
            if protonCat==2 : title.append('multi-single')
            if protonCat==3 : title.append('single-multi')
            if protonCat==4 : title.append('single-single')
            if not nvtx is None:
                title.append( 'N_{{vtx}}{}{}'.format('<' if nvtx<0 else '#geq',abs(nvtx)))
                pname='mmiss_{}_{}_{}{}'.format(boson,protonCat,'lt' if nvtx<0 else 'gt',abs(nvtx))
            else:
                pname='mmiss_{}_{}'.format(boson,protonCat)
            categ_groups[key]={'title':title, 'name':pname, 'group':[]}
        categ_groups[key]['group'].append( (c,xangle) )
        
    for key in categ_groups:
        categ_groups[key]['group']=sorted(categ_groups[key]['group'], key=lambda x: x[1])

    print('[groupCategories] {} groups found'.format(len(categ_groups)))

    return categ_groups

def mergeCategoriesFromList( categ_dict, categ_list,categ_name ):

    merged_categ={'cat':categ_name}
    for c in categ_list:
        if not 'data' in merged_categ:
            merged_categ['nvtx']=[]
            merged_categ['xangle']=[]
            data=categ_dict[c]['data']
            merged_categ['data']=data.Clone('{}_{}'.format(data.GetName(),categ_name))
            merged_categ['data'].SetDirectory(0)

            for fit in ['postfit','prefit']:
                merged_categ[fit]={}
                for proc,h in categ_dict[c][fit].items():
                    merged_categ[fit][proc]=h.Clone('{}_{}'.format(h.GetName(),categ_name))

            merged_categ['protonCat']=categ_dict[c]['protonCat']
            merged_categ['shape']=categ_dict[c]['shape'].Clone('{}_template'.format(categ_name))
        else:

            merged_categ['data'].Add( categ_dict[c]['data'] )
            for fit in ['postfit','prefit']:
                for proc,h in categ_dict[c][fit].items():
                    merged_categ[fit][proc].Add(h)

        merged_categ['nvtx'].append(categ_dict[c]['nvtx'])
        merged_categ['xangle'].append(categ_dict[c]['xangle'])

    return merged_categ


def mergeCategories(categ_dict):

    """merges by boson"""

    categ_groups={}
    categs_to_merge={}
    for c in categ_dict:
        boson=categ_dict[c]['cat']
        nvtx=categ_dict[c]['nvtx']
        protonCat=categ_dict[c]['protonCat']
        xangle=categ_dict[c]['xangle']

        key=boson
        if not key in categ_groups:
            title=[]
            if boson==121:    title.append('pp#rightarrow ppZeeX')
            if boson==169:    title.append('pp#rightarrow ppZ#mu#muX')
            if boson==22:     title.append('pp#rightarrow pp#gammaX')
            pname='mmiss_{}'.format(boson)
            categ_groups[key]={'title':title, 'name':pname, 'group':[]}
            categs_to_merge[key]={}

        protonCatStr='multi-multi'
        if protonCat==2 : protonCatStr='multi-single'
        if protonCat==3 : protonCatStr='single-multi'
        if protonCat==4 : protonCatStr='single-single'

        #this is done to keep the same format as used in groupCategories method
        mergedcat='mergedcat{}_{}'.format(boson,protonCat)
        categ_groups[key]['group'].append( (mergedcat,protonCatStr) )
        if not mergedcat in categs_to_merge[key]: categs_to_merge[key][mergedcat]=[]
        categs_to_merge[key][mergedcat].append(c)

    #now sum up the histos        
    merged_categ_dict={}
    for key in categ_groups:
        categ_groups[key]['group']=sorted(list(set(categ_groups[key]['group'])), key=lambda x: x[1] )        
        for g,_ in categ_groups[key]['group']:
            merged_categ_dict[g] = mergeCategoriesFromList( categ_dict, categs_to_merge[key][g], g )
    pprint.pprint(merged_categ_dict)
    return merged_categ_dict,categ_groups


if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog datacard [opts] ', version='%prog 1.0')
    parser.add_option('--outdir',    type='string'       , default=''    , help='output directory where the postfit plots are saved.')
    parser.add_option('--prefit',    action='store_true' , default=False , help='make prefit distributions')
    (options, args) = parser.parse_args()

    url=sys.argv[1]
    base_dir=os.path.dirname(url)
    boson='z'
    if 'PPzmmX' in url: boson='zmm'
    if 'PPzeeX' in url: boson='zee'
    if 'PPgX' in url : boson='g'
    from generateBinnedWorkspace import LUMI
    lumi=LUMI['169'] if 'z' in boson else LUMI['22']
    datacard_url='{}/{}_datacard.dat'.format(base_dir,boson)
    datacard_sources_url='{}/{}_cards/'.format(base_dir,boson)
    
    categs_dict=parseCategoriesFromDatacard(datacard_url,datacard_sources_url)
    categs_dict=fillFitShapes(url,options.outdir,categs_dict)
    #pprint.pprint(categs_dict)

    categ_groups=groupCategories(categs_dict)
    #pprint.pprint(categ_groups)

    merged_categs_dict, merged_categ_groups=mergeCategories(categs_dict)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    for cg in categ_groups:
        plotShapes(categs_dict,categ_groups[cg],lumi)

    for cg in merged_categ_groups:
        plotShapes( merged_categs_dict, merged_categ_groups[cg], lumi )
        plotShapes( merged_categs_dict, merged_categ_groups[cg], lumi, False)
