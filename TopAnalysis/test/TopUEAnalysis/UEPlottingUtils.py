#!/usr/bin/env python                                                                                                                                               
import ROOT
from collections import OrderedDict

#cf. http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
PLOTFORMATS=[
    ('#1f78b4',20,1001),
    ('#fb9a99',20,1001),
    ('#e31a1c',20,1001),
    ('#b15928',20,1001),
    ('#33a02c',20,1001),
    ('#b2df8a',20,1001),
    ('#fdbf6f',20,1001),
    ('#ff7f00',20,1001),
    ('#ffff99',20,1001),
    ]

#ffff99 PW+HW++ EE5C
#e31a1c ERD on CR
#b2df8a ISR dn
#fdbf6f FSR up
#fb9a99 QCD based CR
#1f78b4 PW+PY8 CUETP8M2T4
#33a02c ISR up
#b15928 Gluon move CR
#ff7f00 FSR dn

COMPARISONSETS=[
    ('PW+PY8 CUETP8M2T4', [ ('nominal',         ['t#bar{t}']), 
                            ('#deltaCUET8P2MT4',['t#bar{t} UEup',     't#bar{t} UEdn']),
                            ('FSR',             ['t#bar{t} fsr up',   't#bar{t} fsr dn']),
                            ('ISR',             ['t#bar{t} isr up',   't#bar{t} isr dn']),
                            ('hdamp',           ['t#bar{t} hdamp up', 't#bar{t} hdamp dn']),
                            ('CR',              ['t#bar{t} QCDbased', 't#bar{t} ERDon', 't#bar{t} gluon move']) ] 
     ),
    ('FSR up', [ ('nominal',         ['t#bar{t} fsr up'])]),
    ('FSR dn', [ ('nominal',         ['t#bar{t} fsr dn'])]),
    ('ISR up', [ ('nominal',         ['t#bar{t} isr up'])]),
    ('ISR dn', [ ('nominal',         ['t#bar{t} isr dn'])]),
    ('QCD based CR', [ ('nominal',         ['t#bar{t} QCDbased'])]),
    ('ERD on CR',    [ ('nominal',         ['t#bar{t} ERDon'])]),
    ('Gluon move CR',[ ('nominal',         ['t#bar{t} gluon move'])]),
    ('aMC@NLO+PY8 CUETP8M2T4', [ ('nominal', ['t#bar{t} aMC@NLO']) ]),
    ('PW+HW++ EE5C'      , [ ('nominal', ['t#bar{t} Herwig++']) ]),
    ]

PLOTTINGSETS=(
    ['PW+PY8 CUETP8M2T4'],
    ['QCD based CR','ERD on CR','Gluon move CR'],
    ['ISR up','ISR dn'],
    ['FSR up','FSR dn']
#    ['PW+HW++ EE5C']
)




"""
"""
def showProfile(plotColl,sliceAxis,xtitle,ytitle,outName):

    #check ranges and prepare frames
    imax,meanMin,meanMax=0,1.0e9,-1.0e9
    for i in xrange(0,len(plotColl)):
        p=plotColl[i]["Data"]
        meanMin=min(meanMin,p.mean[0])
        if meanMax<p.mean[0]:
            meanMax=p.mean[0]
            imax=i
        
    #check how models should be plotted (order and format)
    modelIndices=OrderedDict()
    for idx in xrange(0,len(PLOTTINGSETS)):
        for key in PLOTTINGSETS[idx]:
            modelIndices[key]=(idx,len(modelIndices))


    #make graphs
    incProfiles={}
    profiles={}
    for i in xrange(0,len(plotColl)):
        for key in plotColl[i]:
            
            #filter only for needed keys
            if key!='Data' and not key in modelIndices : continue

            p=plotColl[i][key]
            mean=p.mean[0]
            statUnc,expUnc,thUnc=p.mean[1]
            totalUnc=ROOT.TMath.Sqrt(statUnc**2+expUnc**2+thUnc**2)

            #start graphs
            if not key in incProfiles:
                incProfiles[key]=ROOT.TGraphErrors()
                incProfiles[key].SetName("prof_%s"%key)

                ci,fill,marker=1,0,20
                if key!='Data': 
                    ci,marker,fill=PLOTFORMATS[ modelIndices[key][1] ]
                    ci=ROOT.TColor.GetColor(ci)

                #print key,ci,marker,fill,modelIndices[key][1]
                incProfiles[key].SetMarkerStyle(marker)
                incProfiles[key].SetFillStyle(fill)
                incProfiles[key].SetMarkerColor(ci)
                incProfiles[key].SetLineColor(ci)
                incProfiles[key].SetFillColor(ci)                
                profiles[key]=incProfiles[key].Clone('excprof_%s'%key)
            
            #fill plots
            if i==0:          
                xcen,xwid=1,0
                if key!='Data':
                    xwid=1.0/float(len(PLOTTINGSETS))
                    xcen=0.5+(modelIndices[key][0]+0.5)*xwid
                    xwid/=2
                incProfiles[key].SetPoint(0,xcen,mean)
                incProfiles[key].SetPointError(0,xwid,totalUnc)
            else:
                xcen=sliceAxis.GetBinCenter(i)
                xwid=sliceAxis.GetBinWidth(i)*0.5
                if key=='Data' : 
                    xwid=0
                else:
                    xcen-=xwid
                    xwid=2*xwid/float(len(PLOTTINGSETS))
                    xcen+=(modelIndices[key][0]+0.5)*xwid
                    xwid/=2
                np=profiles[key].GetN()
                profiles[key].SetPoint(np,xcen,mean)
                profiles[key].SetPointError(np,xwid,totalUnc)

    showFormattedUEPlots(incProfiles,profiles,modelIndices,meanMax*1.5,meanMin*0.5,imax,sliceAxis,xtitle,ytitle,outName)
    incProfilesRatio=getRatiosWithRespectTo(incProfiles,'Data')
    profilesRatio=getRatiosWithRespectTo(profiles,'Data')
    showFormattedUEPlots(incProfilesRatio,profilesRatio,modelIndices,1.5,0.7,1,sliceAxis,xtitle,'K^{%s}_{MC/data}'%ytitle,outName+'_ratio')

    for key in incProfiles:
        incProfiles[key].Delete()
        profiles[key].Delete()
        incProfilesRatio[key].Delete()
        profilesRatio[key].Delete()


def showDifferential(plotColl,obsAxis,sliceAxis,sliceTitle,xtitle,ytitle,outName):
    """show differential distributions"""

    def sliceScale(i): return ROOT.TMath.Power(10,-i)

    #check ranges and prepare frames
    imax,minVal,maxVal=0,1.0e9,-1.0e9
    x,y=ROOT.Double(0),ROOT.Double(0)
    allProfiles=[]
    for i in xrange(0,len(plotColl)):
        p=plotColl[i]["Data"]
        key=('Data',None)
        allProfiles.append( {'Data':p.plot[0].Clone('diffdata_%d'%i)} )
        allProfiles[-1]['Data'].SetMarkerStyle(20+i)
        allProfiles[-1]['Data'].SetFillStyle(0)
        allProfiles[-1]['Data'].SetLineColor(1)
        allProfiles[-1]['Data'].SetMarkerColor(1)
        dataTitle='Data'
        if sliceAxis and i>0:
            dataTitle='%d<%s<%d #scale[0.8]{(x10^{-%d})}'%(sliceAxis.GetBinLowEdge(i),sliceTitle,sliceAxis.GetBinUpEdge(i),i)
        allProfiles[-1]['Data'].SetTitle(dataTitle)
        for np in xrange(0,allProfiles[-1]['Data'].GetN()):
            allProfiles[-1]['Data'].GetPoint(np,x,y)
            ey=allProfiles[-1]['Data'].GetErrorY(np)
            if y==0 : continue
            newy=y*sliceScale(i)
            newey=ey*sliceScale(i)            
            allProfiles[-1]['Data'].SetPoint(np,x,newy)
            allProfiles[-1]['Data'].SetPointError(np,0,newey)
            minVal=min(newy*0.1,minVal)
            maxVal=max(newy*100,maxVal)

    #check how models should be plotted (order and format)
    modelIndices=OrderedDict()
    for idx in xrange(0,len(PLOTTINGSETS)):
        for key in PLOTTINGSETS[idx]:
            modelIndices[key]=(idx,len(modelIndices))

    #make graphs
    for i in xrange(0,len(plotColl)):
        for key in plotColl[i]:
            
            #filter only for needed keys
            if key=='Data' : continue
            if not key in modelIndices : continue
            
            p=plotColl[i][key]
            allProfiles[i][key]=p.plot[0].Clone('diff%s_%d'%(key,i))
            ci,marker,fill=PLOTFORMATS[ modelIndices[key][1] ]
            ci=ROOT.TColor.GetColor(ci)
            allProfiles[i][key].SetMarkerStyle(marker)
            allProfiles[i][key].SetFillStyle(fill)
            allProfiles[i][key].SetMarkerColor(ci)
            allProfiles[i][key].SetLineColor(ci)
            allProfiles[i][key].SetFillColor(ci)                
            
            for np in xrange(0,allProfiles[i][key].GetN()):

                allProfiles[i][key].GetPoint(np,x,y)
                ex=allProfiles[i][key].GetErrorX(np)
                ey=allProfiles[i][key].GetErrorY(np)
            
                newy=y*sliceScale(i)
                newey=ey*sliceScale(i)
                allProfiles[i][key].SetPoint(np,x,newy)
                allProfiles[i][key].SetPointError(np,ex,newey)

    showFormattedDifferentialUEPlots(allProfiles,modelIndices,maxVal,minVal,imax,obsAxis,xtitle,ytitle,outName)

    #compute profile ratios and shift in x
    allProfileRatios=[]
    for profList in allProfiles:
        allProfileRatios.append( getRatiosWithRespectTo(profList,"Data") )
        for key in allProfileRatios[-1]:
            if key=='Data' : continue
            if not key in modelIndices: continue
            for np in xrange(0,allProfileRatios[-1][key].GetN()):
                allProfileRatios[-1][key].GetPoint(np,x,y)
                ex=allProfileRatios[-1][key].GetErrorX(np)
                ey=allProfileRatios[-1][key].GetErrorY(np)
                newex=2*ex/float(len(PLOTTINGSETS)) 
                newx=(x-ex)+(modelIndices[key][0]+0.5)*newex  
                newex/=2.                         
                allProfileRatios[-1][key].SetPoint(np,newx,y)
                allProfileRatios[-1][key].SetPointError(np,newex,ey)
    showFormattedDifferentialRatiosUEPlots(allProfileRatios,modelIndices,obsAxis,sliceAxis,sliceTitle,xtitle,ytitle,outName)
    
    #free mem
    for i in xrange(0,len(allProfiles)):
        for key in allProfiles[i]:
            allProfiles[i][key].Delete()
            allProfileRatios[i][key].Delete()

def getRatiosWithRespectTo(profiles,refKey):
    """
    Compute the ratio using a given key as reference
    """
    profilesRatio={}
    x,y=ROOT.Double(0),ROOT.Double(0)
    xref,yref=ROOT.Double(0),ROOT.Double(0)    
    for key in profiles:
        profilesRatio[key]=profiles[key].Clone('ratio2%s_%s'%(refKey,profiles[key].GetName()))
        for np in xrange(0,profiles[key].GetN()):

            profiles[key].GetPoint(np,x,y)
            profiles[refKey].GetPoint(np,xref,yref)
            ratio=-1 if yref==0 else y/yref

            ex=profiles[key].GetErrorX(np)
            ey=profiles[key].GetErrorY(np)
            eyref=profiles[refKey].GetErrorY(np)
                
            ratioUnc=(eyref*y)**2
            if key!=refKey:ratioUnc+=(ey*yref)**2
            ratioUnc=0 if yref==0 else ROOT.TMath.Sqrt(ratioUnc)/(yref**2)

            profilesRatio[key].SetPoint(np,x,ratio)
            profilesRatio[key].SetPointError(np,ex,ratioUnc)

    return profilesRatio

def showFormattedUEPlots(incProfiles,profiles,modelIndices,meanMax,meanMin,imax,sliceAxis,xtitle,ytitle,outName):
    """This method dumps the formatted plots to the canvas"""

    legAtLeft=True if sliceAxis.GetBinCenter(imax)>(sliceAxis.GetXmax()+sliceAxis.GetXmin())*0.5 else False
    frame=ROOT.TH1F('frame','frame',sliceAxis.GetNbins(),sliceAxis.GetXbins().GetArray())
    frame.GetYaxis().SetRangeUser(meanMin,meanMax)
    incframe=ROOT.TH1F('incframe','incframe',1,0.5,1.5)
    incframe.GetYaxis().SetRangeUser(meanMin,meanMax)

    #start the canvas
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetBottomMargin(0.0)

    c.cd()
    p1=ROOT.TPad('p1','p1',0,0,0.8,0.95)
    p1.SetTopMargin(0.01)
    p1.SetRightMargin(0.01)
    p1.SetLeftMargin(0.12)
    p1.SetBottomMargin(0.1)
    p1.SetGridy()
    p1.Draw()

    c.cd()
    p2=ROOT.TPad('p2','p2',0.8,0,1,0.95)
    p2.SetTopMargin(0.01)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.01)
    p2.SetBottomMargin(0.1)
    p2.SetGridy()
    p2.Draw()

    #inclusive frame
    p2.cd()
    incframe.Draw()
    incframe.GetXaxis().SetNdivisions(0)
    incframe.GetXaxis().SetTitleSize(0.15)
    incframe.GetXaxis().SetTitle('Inclusive')
    incframe.GetXaxis().SetTitleOffset(0.25)
    for key in incProfiles:
        if key =='Data': continue
        incProfiles[key].Draw('2')
    incProfiles['Data'].Draw('ep')
    p2.RedrawAxis()
                
    #differential frame
    p1.cd()

    frame.Draw()
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.04)
    for key in profiles:
        if key =='Data': continue
        profiles[key].Draw('2')        
    profiles['Data'].Draw('ep')
        
    #draw the legend
    startLegX=0.15 if legAtLeft else 0.50
    startLegY=0.9  if legAtLeft else 0.95
    leg=ROOT.TLegend(startLegX,startLegY,startLegX+0.4,startLegY-len(PLOTTINGSETS)*0.06) 
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry( incProfiles['Data'], 'Data','ep' )
    for key in modelIndices: leg.AddEntry(profiles[key],key,'f')
    leg.Draw()
        
    #the header
    tex2=ROOT.TLatex()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.045)
    tex2.SetNDC()
    tex2.DrawLatex(0.16,0.93,'#bf{CMS} #it{preliminary}')

    p1.RedrawAxis()

    #the lumi/sqrts
    c.cd()
    tex3=ROOT.TLatex()
    tex3.SetTextFont(42)
    tex3.SetTextSize(0.04)
    tex3.SetNDC()
    tex3.DrawLatex(0.7,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

    # all done
    c.Modified()
    c.Update()
    for ext in ['pdf','png']: c.SaveAs('%s.%s'%(outName,ext))

    p1.Delete()
    p2.Delete()
    c.Delete()

def showFormattedDifferentialUEPlots(allProfiles,modelIndices,maxVal,minVal,imax,sliceAxis,xtitle,ytitle,outName):
    """Differential distributions"""
    
    legAtLeft=True if sliceAxis.GetBinCenter(imax)>(sliceAxis.GetXmax()+sliceAxis.GetXmin())*0.5 else False
    frame=ROOT.TH1F('frame','frame',sliceAxis.GetNbins(),sliceAxis.GetXbins().GetArray())
    frame.GetYaxis().SetRangeUser(minVal,maxVal)
    
    #start the canvas
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetLogy()

    frame.Draw()
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.04)

    #draw the legend
    startLegX=0.15 if legAtLeft else 0.50
    startLegY=0.9  if legAtLeft else 0.93
    leg=ROOT.TLegend(startLegX,startLegY,startLegX+0.4,startLegY-(len(allProfiles)+1)*0.04) 
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)

    for profList in allProfiles:
        profList[COMPARISONSETS[0][0]].Draw('20')
        profList['Data'].Draw('ep')
        leg.AddEntry(profList['Data'],profList['Data'].GetTitle(),'ep')
    leg.AddEntry(allProfiles[0][COMPARISONSETS[0][0]],COMPARISONSETS[0][0],'f')
    leg.Draw()
        
    #the header
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.045)
    tex.SetNDC()
    tex.DrawLatex(0.16,0.9,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.65,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

    # all done
    c.RedrawAxis()
    c.Modified()
    c.Update()
    for ext in ['pdf','png']: c.SaveAs('%s.%s'%(outName,ext))

    c.Delete()

def showFormattedDifferentialRatiosUEPlots(allProfileRatios,modelIndices,obsAxis,sliceAxis,sliceTitle,xtitle,ytitle,outName):
    """Ratios of differential distributions"""

    npads=len(allProfileRatios)
    height=500 if npads==1 else 1000
    dyBottom,dyTop=50,120
    padHeight=(1-float(dyBottom+dyTop)/float(height))/float(npads)
    
    fontSize=0.04 if npads==1 else 0.04/(2*padHeight)

    frame=ROOT.TH1F('frame','frame',obsAxis.GetNbins(),obsAxis.GetXbins().GetArray())
    frame.GetYaxis().SetNdivisions(10)
    frame.GetYaxis().SetRangeUser(0.44,1.56)
    frame.GetYaxis().SetTitle('K^{%s}_{MC/data}'%ytitle)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(fontSize)
    frame.GetYaxis().SetLabelSize(fontSize)
    frame.GetXaxis().SetTitleSize(0.0)
    frame.GetXaxis().SetLabelSize(0.0)
    finalFrame=frame.Clone('finalframe')
    finalFontSize=fontSize*padHeight*height/(padHeight*height+dyBottom)
    finalFrame.GetXaxis().SetTitleSize(finalFontSize)
    finalFrame.GetXaxis().SetLabelSize(finalFontSize)
    finalFrame.GetXaxis().SetTitle(xtitle)
    finalFrame.SetTitleOffset(1.1)
    finalFrame.GetYaxis().SetTitleSize(finalFontSize)
    finalFrame.GetYaxis().SetLabelSize(finalFontSize)


    #start the canvas
    c=ROOT.TCanvas('c','c',500,height)
    c.SetTopMargin(float(dyTop)/float(height))

    c.SetLeftMargin(0.)
    c.SetRightMargin(0.0)
    c.SetBottomMargin(0.)
    allPads=[]
    allCaptions=[]
    for i in xrange(0,npads):
        c.cd()
        padAnchor=(1-float(dyTop)/float(height))-i*padHeight
        yup,ydn=padAnchor,padAnchor-padHeight
        if i==npads-1 : ydn=0
        allPads.append( ROOT.TPad('p%d'%i,'p%d'%i,0,yup,1,ydn) )
        allPads[-1].SetTopMargin(0.01)
        allPads[-1].SetLeftMargin(0.2)
        allPads[-1].SetRightMargin(0.03)
        allPads[-1].SetBottomMargin(0.01 if i!=npads-1 else float(dyBottom)/(height*padHeight) )
        allPads[-1].Draw()
        allPads[-1].SetGridy()
        allPads[-1].cd()
        if i==npads-1: 
            finalFrame.Draw()
        else : frame.Draw()
        for key in allProfileRatios[i]:
            if key =='Data': continue
            allProfileRatios[i][key].Draw('2')
        allProfileRatios[i]['Data'].Draw('ep')
        
        if npads==1: continue
        allCaptions.append( ROOT.TLatex() )
        allCaptions[-1].SetTextFont(42)
        allCaptions[-1].SetTextSize( 0.8*(fontSize if i!=npads-1 else finalFontSize ) )
        allCaptions[-1].SetNDC()
        caption='inclusive'
        if i>0:
            caption='%d<%s<%d'%(sliceAxis.GetBinLowEdge(i),sliceTitle,sliceAxis.GetBinUpEdge(i))
        allCaptions[-1].DrawLatex(0.8,0.92,caption)

        
    #the header
    c.cd()
    allPads.append( ROOT.TPad('ph','ph',0,1,1,1-float(dyTop)/float(height)) )
    allPads[-1].SetTopMargin(0)
    allPads[-1].SetBottomMargin(0)
    allPads[-1].SetLeftMargin(0)
    allPads[-1].SetRightMargin(0)
    allPads[-1].Draw()
    allPads[-1].cd()

    leg=ROOT.TLegend(0.1,0.05,1,0.7)
    leg.SetTextFont(42)
    leg.SetTextSize(fontSize*1.2 if npads!=1 else fontSize*2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(allProfileRatios[0]['Data'],'Data','ep')
    for key in modelIndices: leg.AddEntry(allProfileRatios[0][key],key,'f')
    leg.SetNColumns((len(modelIndices)+1)/3)
    leg.Draw()

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(fontSize*1.8 if npads!=1 else fontSize*2.5)
    tex.SetNDC()
    tex.DrawLatex(0.1,0.8,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.67,0.8,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

    # all done    
    c.Modified()
    c.Update()
    for ext in ['pdf','png']: c.SaveAs('%s_ratio.%s'%(outName,ext))
    for p in allPads: p.Delete()
    c.Delete()
