


"""
"""
def showMatrices(opt,varIdx=0,wgtIdx=0):

    #prepare canvas
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(55)
    c=ROOT.TCanvas('c','c',1000,1000)
    c.SetTopMargin(0.01)
    c.SetRightMargin(0.1)
    c.SetBottomMargin(0.15)

    
    #open ROOT file with the definitions
    url='%s/UEanalysis.root'%opt.out
    sliceVarsQ,obsVarsQ=readSlicesAndObservables(url=url)
    fIn=ROOT.TFile.Open(url)
    mmDir='migmatrices_%d_%d'%(varIdx,wgtIdx)

    #build the fully combined migration matrices
    bigMatrix={}
    for s in sliceVarsQ:
        nRecBins=sliceVarsQ[s].GetXaxis().GetNbins()
        nGenBins=nRecBins if 'nj' in s else nRecBins/2
        for o in obsVarsQ:
            if 'nch' in o and 'nch' in s: continue
            obsNRecBins=obsVarsQ[o].GetXaxis().GetNbins()
            obsNGenBins=obsVarsQ[o].GetXaxis().GetNbins()/2

            key=(s,o)
            nbinsX=1+3*nGenBins*obsNGenBins
            nbinsY=1+3*nRecBins*obsNRecBins
            bigMatrix[key]=ROOT.TH2F('%s_%s'%(s,o),'',nbinsX,0,nbinsX,nbinsY,0,nbinsY)
            bigMatrix[key].SetDirectory(0)

            #fill the matrix
            genBinCtr=0
            for xbin in xrange(1,nbinsX+1):
                if xbin!=1 :
                    genBinCtr+=1
                    if genBinCtr>obsNGenBins: genBinCtr=1

                recBinCtr=0
                for ybin in xrange(1,nbinsY+1):
                    if ybin!=1:
                        recBinCtr+=1
                        if recBinCtr>obsNRecBins: recBinCtr=1                

                    sliceGenBin=(xbin-2)/(3*obsNGenBins)+1
                    obsGenRegionBin=((xbin-2)/obsNGenBins)%3
                    obsGenRegion=getRegionName(obsGenRegionBin)

                    sliceRecBin=(ybin-2)/(3*obsNRecBins)+1
                    obsRecRegionBin=((ybin-2)/obsNRecBins)%3
                    obsRecRegion=getRegionName(obsRecRegionBin)   

                    if xbin==1 and ybin==1:
                        mmH=fIn.Get('%s/%s00_%s'%(mmDir,s,o))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(1,1))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(1,1))
                    elif xbin==1:                                                
                        mmH=fIn.Get('%s/%s0%d_%s%s'%(mmDir,s,sliceRecBin,o,obsRecRegion))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(1,recBinCtr))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(1,recBinCtr))
                    elif ybin==1:
                        mmH=fIn.Get('%s/%s%d0_%s%s'%(mmDir,s,sliceGenBin,o,obsGenRegion))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(genBinCtr,1))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(genBinCtr,1))
                    else:
                        mmH=fIn.Get('%s/%s%d%d_%s%s%s'%(mmDir,s,sliceGenBin,sliceRecBin,o,obsGenRegion,obsRecRegion))
                        bigMatrix[key].SetBinContent(xbin,ybin,mmH.GetBinContent(genBinCtr,recBinCtr))
                        bigMatrix[key].SetBinError(xbin,ybin,mmH.GetBinError(genBinCtr,recBinCtr))
                        
                        
            #label x-axis
            sliceBin,obsBin=0,0
            for xbin in xrange(1,nbinsX+1):
                label=''
                if xbin==1: label='fail'
                else:
                    if (xbin+1-obsNGenBins/2)%(3*obsNGenBins)==0:
                        label='#color[38]{[%3.1f,%3.1f[}'%(sliceVarsQ[s].GetXaxis().GetBinLowEdge(2*sliceBin+1),sliceVarsQ[s].GetXaxis().GetBinUpEdge(2*sliceBin+2))
                        sliceBin+=1
                    elif (xbin+4)%(obsNGenBins)==0:
                        obsBin+=1
                        label='#it{%s}'%getRegionName((obsBin-1)%3)
                if len(label)==0: continue    
                bigMatrix[key].GetXaxis().SetBinLabel(xbin,label)

            #label y-axis
            sliceBin,obsBin=0,0
            for ybin in xrange(1,nbinsY+1):
                label=''
                if ybin==1: label='fail'
                else:
                    if (ybin+1-obsNRecBins/2)%(3*obsNRecBins)==0:
                        label='#color[38]{[%3.1f,%3.1f[}'%(sliceVarsQ[s].GetXaxis().GetBinLowEdge(sliceBin+1),sliceVarsQ[s].GetXaxis().GetBinUpEdge(sliceBin+1))
                        sliceBin+=1
                    elif (ybin+4)%(obsNRecBins)==0:
                        obsBin+=1
                        label='#it{%s}'%getRegionName((obsBin-1)%3)
                if len(label)==0: continue    
                bigMatrix[key].GetYaxis().SetBinLabel(ybin,label)
            
            #display the matrix
            c.Clear()
            #c.SetLogz()
            bigMatrix[key].Draw('colz')
            bigMatrix[key].GetZaxis().SetLabelSize(0.03)
            bigMatrix[key].GetYaxis().SetLabelSize(0.03)
            bigMatrix[key].GetXaxis().SetLabelSize(0.03)
            bigMatrix[key].GetXaxis().SetTickLength(0)
            bigMatrix[key].GetYaxis().SetTickLength(0)
            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.025)
            tex.SetNDC()
            tex.DrawLatex(0.01,0.96,'#bf{#splitline{Reco.}{level}}')
            tex.DrawLatex(0.95,0.1,'#bf{#splitline{Gen.}{level}}')
            stit='p_{T}(t#bar{t})'
            if 'nj' in s : stit='N(jets)'
            if 'nch' in s : stit='N(ch)'                
            if 'mll' in s : stit='M(l,l)'                
            if 'dphill' in s : stit='#Delta#phi(l,l)'
            otit='N(ch)'
            if 'ptsum' in o : otit='#Sigmap_{T}(ch)'
            if 'avgpt' in o : otit='#bar{p}_{T}(ch)'
            tex.DrawLatex(0.10,0.03,'#scale[1.2]{#bf{CMS}} #it{simulation preliminary} %s vs %s'%(stit,otit))

            line=ROOT.TLine()           
            line.SetLineWidth(1)
            for xbin in xrange(1,nbinsX+1):
                if xbin==1 or (xbin-1)%(obsNGenBins)==0:
                    ls=1 if (xbin-1)%(3*obsNGenBins)==0 else 9
                    lc=1 if (xbin-1)%(3*obsNGenBins)==0 else ROOT.kGray
                    line.SetLineStyle(ls)
                    line.SetLineColor(lc)  
                    line.DrawLine(bigMatrix[key].GetXaxis().GetBinUpEdge(xbin),bigMatrix[key].GetYaxis().GetXmin(),bigMatrix[key].GetXaxis().GetBinUpEdge(xbin),bigMatrix[key].GetYaxis().GetXmax())
            for ybin in xrange(1,nbinsY+1):
                if ybin==1 or (ybin-1)%(obsNRecBins)==0:
                    ls=1 if (ybin-1)%(3*obsNRecBins)==0 else 9
                    lc=1 if (ybin-1)%(3*obsNRecBins)==0 else ROOT.kGray
                    line.SetLineStyle(ls)
                    line.SetLineColor(lc)                      
                    line.DrawLine(bigMatrix[key].GetXaxis().GetXmin(),bigMatrix[key].GetYaxis().GetBinUpEdge(ybin),bigMatrix[key].GetXaxis().GetXmax(),bigMatrix[key].GetYaxis().GetBinUpEdge(ybin))

                    
            c.Modified()
            c.Update()
            raw_input()
            #for ext in ['png','pdf']:
            #    c.SaveAs('%s/%s_%s_migration.%s'%(opt.out,s,o,ext))

    #save to ROOT file
    fOut=ROOT.TFile.Open('%s/UEanalysis.root'%opt.out,'UPDATE')
    outDir=fOut.Get('final_'+mmDir)
    try:
        outDir.cd()
    except:
        outDir=fOut.mkdir('final_'+mmDir)
        outDir.cd()
    for k in bigMatrix:
        bigMatrix[k].SetDirectory(outDir)
        bigMatrix[k].Write(bigMatrix[k].GetName(),ROOT.TObject.kOverwrite)
    fOut.Close()
