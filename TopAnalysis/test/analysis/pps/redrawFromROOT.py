import ROOT
import numpy as np

def redraw(url,zoom_xmin,zoom_xmax,zoom_sym,txtmoves,extrapolate_first0,newName=None,legtextsize=None,
           axistextsize=None,txtsizes={},startxleg=None,startyleg=None):
    inf=ROOT.TFile.Open(url)
    c=inf.Get('c')
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    c.Draw()
    p1=c.GetPrimitive('p1')
    frame=p1.GetPrimitive('frame')
    frame.GetYaxis().SetTitle('Events / bin')
    if not axistextsize is None:
        frame.GetYaxis().SetTitleSize(axistextsize[0])
        frame.GetXaxis().SetTitleSize(axistextsize[0])
        frame.GetYaxis().SetLabelSize(axistextsize[1])
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.GetXaxis().SetLabelSize(axistextsize[1])
    extrapolated_pts=[]
    xmin,xmax=None,None
    for x in p1.GetListOfPrimitives():

        if x.GetName()=='totalmcuncband':
            
            grx=ROOT.TGraph()
            for xbin in range(x.GetNbinsX()):
                ycts=x.GetBinContent(x.GetNbinsX()-xbin-1)
                if ycts==0: continue
                for ix in range(xbin,xbin+3):
                    grx.SetPoint(grx.GetN(),
                                 x.GetXaxis().GetBinCenter(x.GetNbinsX()-ix-1),
                                 x.GetBinContent(x.GetNbinsX()-ix-1))
                break

            fname='expo'
            grx.Fit(fname,'','N')
            gaus=grx.GetFunction(fname)
            first0 = x.GetXaxis().GetBinCenter(x.GetNbinsX()-xbin)
            extrapolated_pts.append( (x.GetNbinsX()-xbin,first0, gaus.Eval(first0), 2*x.GetBinError(x.GetNbinsX()-xbin-1) ) )
            if extrapolate_first0 :
                x.SetBinContent(x.GetNbinsX()-xbin,extrapolated_pts[-1][2])
                x.SetBinError(x.GetNbinsX()-xbin,extrapolated_pts[-1][3])

        if x.InheritsFrom('TLatex'):
            pos=(x.GetX(),x.GetY())
            print(pos,x.GetTitle())
            if pos in txtsizes:
                x.SetTextSize(txtsizes[pos])
            if pos in txtmoves:
                pos=txtmoves[pos]
            if x.GetTitle()=='#scale[0.9]{m_{X}=1000 GeV}':
                x.SetTitle('pp#rightarrowppZX(1000) - 1 pb')
            x.SetX(pos[0])
            x.SetY(pos[1])


        #update labels
        if x.InheritsFrom('TLegend'):
            x.SetTextSize(0.05 if legtextsize is None else legtextsize)
            if not startxleg is None:
                x.SetX1NDC(startxleg)
            if not startyleg is None:
                dy=abs(x.GetY1NDC()-x.GetY2NDC())
                x.SetY1NDC(startyleg)
                x.SetY2NDC(startyleg+dy)
            for ientry in x.GetListOfPrimitives():
                newTitle = ientry.GetLabel()
                print(newTitle)
                for tkn,newtkn in [('ZX(','pp#rightarrowppZX('),                                   
                                   ('#gammaX(','pp#rightarrowpp#gammaX('),
                                   ('fiducial (1 pb)','fiducial'),
                        #('#scale[0.8]{ - 100pb}','- 100pb')
                               ]:
                    newTitle=newTitle.replace(tkn,newtkn)
                ientry.SetLabel(newTitle)
        
        #determine the appropriate zoom from the data
        if x.GetTitle()=='Data':
            n=x.GetN()            
            xarr=x.GetX()
            xarr=np.array( [xarr[i] for i in range(n)] )
            yarr=x.GetY()
            yarr=np.array( [yarr[i] for i in range(n)] )
            xmin=xarr[yarr>0][0]
            xmax=xarr[yarr>0][-1]    

    #set the zoom ranges
    if not zoom_xmin is None: xmin=zoom_xmin
    if not zoom_xmax is None: xmax=zoom_xmax
    if zoom_sym:
        xabs=max([abs(xmin),abs(xmax)])
        xmin=-xabs
        xmax=xabs

    #zoom
    if not xmin is None and not xmin is None:
        frame.GetXaxis().SetRangeUser(xmin,xmax)
    p1.cd()
    p1.Modified()
    p1.Update()
    p1.RedrawAxis()
    ROOT.gPad.RedrawAxis();
    ROOT.gPad.RedrawAxis("G");


    #zoom also ratio
    try:
        p2=c.GetPrimitive('p2')
        ratioframe=p2.GetPrimitive('ratioframe')
        ratioframe.GetXaxis().SetRangeUser(xmin,xmax)
        grout=[]
        if extrapolate_first0 :

            lastx=xarr[yarr>0][-1]
            lasty=yarr[yarr>0][-1]
            rbin=extrapolated_pts[-1][0]
            rval=lasty/extrapolated_pts[-1][2]
            
            grout.append( ROOT.TGraph() )
            grout[-1].SetMarkerStyle(22)
            grout[-1].SetMarkerColor(9)
            grout[-1].SetMarkerSize(1.6)
            grout[-1].SetLineColor(9)
            grout[-1].SetFillStyle(0)
            grout.append( grout[-1].Clone() )
            grout[-1].SetMarkerStyle(23)

            for xp2 in p2.GetListOfPrimitives():
                if xp2.GetName()=='ratioframe':
                    xp2.SetBinContent(rbin,1)
                    xp2.SetBinError(rbin,1)
                
                if xp2.InheritsFrom('TGraphAsymmErrors'):
                    lastx=xarr[yarr>0][-1]
                    lasty=yarr[yarr>0][-1]
                    np_xp2=xp2.GetN()
                    rval=lasty/extrapolated_pts[-1][0]
        
                    if rval>1.6:
                        grout[-2].SetPoint(grout[-2].GetN(),lastx,1.6-0.04)
                    elif rval<0.4:
                        grout[-1].SetPoint(grout[-1].GetN(),lastx,0.4+0.04)
                    else:
                        xp2.SetPoint(np_xp2,lastx,rval)
                        # xp2.SetPointEYhigh(np_xp2,rval)
                        # xp2.SetPointEYlow(np_xp2,rval)
            p2.cd()
            for g in grout:
                g.Draw('p')
            p2.Modified()
            p2.Update()
            p2.RedrawAxis()
            ROOT.gPad.RedrawAxis();
            ROOT.gPad.RedrawAxis("G");

    except:
        print('No second pad')
        pass


    c.cd()
    c.Modified()
    c.Update()
    c.RedrawAxis()
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gPad.RedrawAxis();
    ROOT.gPad.RedrawAxis("G");
    if newName is None: newName=ROOT.gSystem.BaseName(url)
    newurl='/eos/user/p/psilva/www/EXO-19-009/paper_fr/'+newName
    c.SaveAs(newurl.replace('.root','.pdf'))
    c.SaveAs(newurl)
    raw_input()
    inf.Close()


ROOT.gROOT.SetBatch(True)

#figure 3
for rp in range(1,5):
    redraw('mmass_1000_sigacc_{}exc.root'.format(rp),
           0,1800,False,
           {(0.18, 0.86):(0.18, 0.84),(0.18, 0.41000000000000003):(0.18, 0.71),(0.18,0.46):(0.18,0.76)},
           True,
           legtextsize=0.04,
           startyleg=0.73,
           axistextsize=(0.05,0.04),
           txtsizes={(0.18, 0.86):0.05,(0.18, 0.41000000000000003):0.04,(0.18, 0.46):0.04,(0.95, 0.965):0.04}
    )

#figure 4
#for fs in ['mm','ee','a'] :                             
#    redraw('/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/plots/csi_{}rpinhpur0neg.root'.format(fs),0.035,0.22,False,{(0.2, 0.89):(0.18, 0.87),(0.2, 0.54):(0.18, 0.78),(0.2, 0.44000000000000006):(0.18, 0.69),},True,startxleg=0.5,legtextsize=0.045)
#    redraw('/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/plots/csi_{}rpinhpur0pos.root'.format(fs),0.035,0.22,False,{(0.2, 0.89):(0.18, 0.87),(0.2, 0.54):(0.18, 0.78),(0.2, 0.44000000000000006):(0.18, 0.69),},True,startxleg=0.5,legtextsize=0.045)
#    redraw('/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/plots/ypp_{}rpinhpur.root'.format(fs),-0.8,1.0,False,{(0.2, 0.89):(0.18, 0.87),(0.2, 0.54):(0.18,0.8)},False if fs!='mm' else True,startxleg=0.5,legtextsize=0.045)

#figure 5
#for fig in ['csi1_0','csi2_0']:
#    redraw('/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/bkg_ptll40/{}.root'.format(fig),0.035,0.18,False,{(0.2, 0.89):(0.18,0.87),(0.2, 0.49000000000000005):(0.18,0.8)},False,fig+'_emu_ptll40.root')
#for fig in ['mpp_0','mmiss_0']:
#    redraw('/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/bkg_ptll40/{}.root'.format(fig),None,None,False,{(0.2, 0.89):(0.18,0.87)},False,fig+'_emu_ptll40.root')
