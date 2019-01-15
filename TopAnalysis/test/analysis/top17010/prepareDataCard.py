#!/usr/bin/env python

import os
import sys
import optparse
import json
from time import time
from datetime import datetime
from collections import OrderedDict
from TopLJets2015.TopAnalysis.Plot import fixExtremities
import ROOT

def customizeSignalShapes(binName,sigName,sigF,baseSigF,outDir):

    """ 
    takes specific signal and scales the base signal simulation + systematics by the ratio sig/orig_signal
    the result is stored in {0}_{1}.shapes.root {0}=signal name {1}=binName
    """
    
    #get signal
    fIn=ROOT.TFile.Open(sigF)
    sigH=fIn.Get('{0}_mlb'.format(binName))
    sigH.SetDirectory(0)
    fixExtremities(sigH) # this needs to be done as the signal is the raw output of TOP-17-010.cc
    fIn.Close()

    #get original signal
    fIn=ROOT.TFile.Open(baseSigF)
    origSigH=fIn.Get('{0}_mlb/central'.format(binName))
    origSigH.SetDirectory(0)

    #force the integral to be the same as this is a shape analysis
    sf=origSigH.Integral()/sigH.Integral()
    sigH.Scale(sf)
    
    #transfer factor
    tfH=sigH.Clone('tf')
    tfH.Divide(origSigH)
    tfH.SetDirectory(0)
    
    #prepare the output applying the transfer factor to each histogram found
    fOut=ROOT.TFile.Open(os.path.join(outDir,'{0}.shapes.root'.format(sigName)),'RECREATE')
    fdir=fOut.mkdir('{0}_mlb'.format(binName))
    for k in fIn.Get('{0}_mlb'.format(binName)).GetListOfKeys():
        h=k.ReadObj()
        h.SetDirectory(fdir)
        for xbin in range(h.GetNbinsX()):
            h.SetBinContent(xbin+1,h.GetBinContent(xbin+1)*tfH.GetBinContent(xbin+1))
        fdir.cd()
        h.Write()
    fOut.Close()

    #free mem
    tfH.Delete()
    fIn.Close()

def customizeData(binName,dataDef,bkgList,templDir,outDir):

    """ customizes data to use as observation """
    dataDef=dataDef.replace("\"","")
    dataType,dataF,dataHisto=dataDef.split(',')
    
    #get the data
    inF=ROOT.TFile.Open(dataF)
    dataH=inF.Get(dataHisto)
    dataH.SetDirectory(0)
    inF.Close()
    
    #specific for pseudo-data
    #sum up signal to background expectations when dataType is 'sig' type
    #round each bin to an integer
    if dataType=='sig':
        for proc in bkgList:
            inF=ROOT.TFile.Open('{0}/templates_{1}.root'.format(templDir,proc))
            h=inF.Get('{0}_mlb/central'.format(binName))
            dataH.Add(h)
            inF.Close()

        for xbin in range(dataH.GetNbinsX()):
            dataH.SetBinContent(xbin+1,int(dataH.GetBinContent(xbin+1)))
    
        dataType=dataH.GetTitle()
        for tkn in [' ','#','_','^','{','}',',',':','+','-']:
            dataType=dataType.replace(tkn,'')
        dataType=dataType.replace('.','p')
        dataType='pseudodata.%s'%dataType

    #save to file
    dataShapesURL=os.path.join(outDir,'%s.shapes.root'%dataType)
    fOut=ROOT.TFile.Open(dataShapesURL,'RECREATE')
    fdir=fOut.mkdir('{0}_mlb'.format(binName))
    fdir.cd()
    dataH.SetTitle('')
    dataH.SetDirectory(fdir)
    dataH.Write('data_obs')
    fOut.Close()

    return dataShapesURL

def shapeIsPresent(url,sname,binName):

    """checks if syst shape is present"""

    isPresent=False
    fIn=ROOT.TFile.Open(url)
    varUp=fIn.Get('{0}_mlb/{1}Up'.format(binName,sname))
    varDn=fIn.Get('{0}_mlb/{1}Down'.format(binName,sname))
    try:
        isPresent=True if varUp.Integral()>0 and varDn.Integral()>0 else False
    except:
        pass
    fIn.Close()

    return isPresent


def systIsConsistent(sname,binName,proc):

    """a series of rules to veto inconsistent systematics"""

    if sname=='toppt'       and proc!='ttbar'         : return False
    if sname=='mmtrig'      and binName.find('mm')!=0 : return False
    if sname=='eetrig'      and binName.find('ee')!=0 : return False
    if sname=='emtrig'      and binName.find('em')!=0 : return False
    if sname=='esel'        and binName.find('mm')==0 : return False
    if sname=='msel'        and binName.find('ee')==0 : return False
    if sname.find('ees')==0 and binName.find('mm')==0 : return False
    if sname.find('mes')==0 and binName.find('ee')==0 : return False

    return True


def printHeader(dc,binName,procList,templDir,outDir):

    """ create datacard header """

    shapeFiles={}

    nProcs=len(procList)
    ts=time()

    dc.write('# TOP-17-010 datacard for category %s\n'%binName)
    dc.write('# generated @ %s by %s\n'%(datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'),os.environ.get("USER")))
    dc.write('-'*50 + '\n')
    dc.write('imax *\n')
    dc.write('jmax *\n')
    dc.write('kmax *\n')
    dc.write('-'*50 + '\n')
    for i in range(nProcs):
        templFile='{0}/templates_{1}.root'.format(templDir,procList[i]) if i>0 else os.path.join(outDir,'{0}.shapes.root'.format(procList[0]))
        dc.write('shapes {0} {1} {2} $CHANNEL_mlb/central $CHANNEL_mlb/$SYSTEMATIC\n'.format(procList[i],binName,templFile))
        shapeFiles[procList[i]]=templFile
    dc.write('shapes data_obs {0} {1}/_DATAOBSSHAPES_ $CHANNEL_mlb/$PROCESS $CHANNEL_mlb/$PROCESS_$SYSTEMATIC \n'.format(binName,outDir))
    dc.write('-'*50 + '\n')
    dc.write('bin %s\n'%binName)
    dc.write('observation -1\n')
    dc.write('-'*50 + '\n')
    dc.write('%30s'%'bin' + ('%15s'%binName)*nProcs + '\n')
    dc.write('%30s'%'process' + ''.join(['%15s'%p for p in procList]) + '\n')
    dc.write('%30s'%'process' + ''.join(['%15s'%i for i in range(nProcs)]) + '\n')
    dc.write('%30s'%'rate'    + ('%15s'%'-1')*nProcs + '\n')
    dc.write('-'*50 + '\n')

    return shapeFiles


def printShapeSysts(dc,syst_dict,binName,shapeFiles):

    """maps which systs apply to each process and then dumps the datacard accordingly"""

    procList=[x for x,_ in syst_dict]

    allSysts={}
    for proc,proc_desc in syst_dict:
        for item in proc_desc:
            try:
                for syst in proc_desc[item].keys():
                    syst=syst.format(binName)
                    if not syst in allSysts: allSysts[syst]=[]
                    allSysts[syst].append(proc)
            except Exception as e:
                pass

    for s in allSysts:
        validShape=False
        sline='%30s  shape   '%s
        for p in procList:
            if p in allSysts[s] and shapeIsPresent(shapeFiles[p],s,binName) and systIsConsistent(s,binName,p):
                sline+='%15s'%'1'
                validShape=True
            else : 
                sline+='%15s'%'-'
        if not validShape: continue
        dc.write(sline+'\n')


def printBinByBinUncs(dc,binName,procList,shapeFiles):

    """
    loops over the template files and looks for bin-by-bin uncertainty histograms
    if found they are added to the datacard
    """

    #check if there are bin-by-bin uncertainty histograms
    dirName='%s_mlb'%binName
    nprocs=len(procList)
    for i in range(nprocs):

        p=procList[i]

        #loop over all template histos and match bin-by-bin uncertainty name
        matchTags=['{0}_{1}_mlb_bin'.format(p,binName),
                   '{0}_{1}_mlb_sysbin'.format(p,binName)]
        fIn=ROOT.TFile.Open(shapeFiles[p])
        for obj in fIn.Get(dirName).GetListOfKeys():
            objname=obj.GetName()
            if not 'Up' in objname: continue

            #dump uncertainty to datacard if it matches a tag
            for m in matchTags:
                if not m in objname: continue                        
                uncName=objname[0:-2]
                sline='%30s  shape   '%uncName
                for j in range(nprocs):
                    sline+='%15s'%('1' if j==i else '-')
                dc.write(sline+'\n')
                break

def printRateSysts(dc,binName,procList):
    """ this dumps to the datacard a series of hardcoded rate systematics """
 
    rateSysts=[
        ('lumi_13TeV',                  1.025,  None,     ['dy']),
        ('dynorm_{0}'.format(binName),  1.30,   ['dy'],   None),
        ('twnorm_th',                   1.15,   ['tw'],   None),
        ('vvnorm_th',                   1.20,   ['vv'],   None),
        ]

    for uncName,uncVal,whiteList,blackList in rateSysts:
        sline='%30s  lnN     '%uncName
        for p in procList:
            hasSyst=True
            if whiteList and not p in whiteList: hasSyst=False
            if blackList and p in blackList: hasSyst=False
            sline+='%15s'%('%3.3f'%uncVal if hasSyst else '-')
        dc.write(sline+'\n')
            


def main():

    usage = 'usage: %prog dataDef1="type1,plotter1,dir1" ... [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--dist',          
                      dest='dist',       
                      help='distribution',
                      default=None,
                      type='string')
    parser.add_option('-t', '--templ',          
                      dest='templ',
                      help='template directory [%default]',  
                      default='/eos/cms/store/cmst3/group/top/TOP17010/0c522df/templates',
                      type='string')
    parser.add_option('-s', '--sig',
                      dest='signal',
                      help='signal to use in the datacard to fit the data [%default]',  
                      default=None,
                      type='string')
    parser.add_option('-o', '--outdir',          
                      dest='outDir',
                      help='output directory [%default]',  
                      default='datacard',
                      type='string')
    parser.add_option('--systs',          
                      dest='systs',       
                      help='description of the systematics [%default]',
                      default='test/analysis/top17010/systs_dict.json',
                      type='string')
    (opt, args) = parser.parse_args()

    #category
    cat=opt.dist.split('_')[0]

    #signal to use in the fit
    sigName,sigFile=opt.signal.split(',')

    #prepare output directory (append category and signal name
    opt.outDir='%s/%s/%s'%(opt.outDir,cat,sigName)
    os.system('mkdir -p %s'%opt.outDir)

    #decode process and systs needed from the json file
    with open(opt.systs,'r') as cache:
        syst_dict=json.load( cache, encoding='utf-8', object_pairs_hook=OrderedDict ).items()
    procList=[x for x,_ in syst_dict]


    #customize signal shapes
    sigName=syst_dict[0][0]
    customizeSignalShapes(binName=cat,
                          sigName=sigName,
                          sigF=sigFile,
                          baseSigF=os.path.join(opt.templ,'templates_%s.root'%sigName),
                          outDir=opt.outDir)

    #customize data
    dataFiles=[customizeData(binName=cat,
                             dataDef=dataDef.split("=")[1],
                             bkgList=procList[1:],
                             templDir=opt.templ,
                             outDir=opt.outDir)
               for dataDef in args]
        
    #dump the template datacard
    dcTemplURL=os.path.join(opt.outDir,'datacard.dat.templ')
    with open(dcTemplURL,'w') as dc:
        shapeFiles=printHeader(dc=dc,
                               binName=cat,
                               procList=procList,
                               templDir=opt.templ,
                               outDir=os.path.abspath(opt.outDir))
        printShapeSysts(dc=dc,
                        syst_dict=syst_dict,
                        binName=cat,
                        shapeFiles=shapeFiles)
        printBinByBinUncs(dc=dc,
                          procList=procList,
                          binName=cat,
                          shapeFiles=shapeFiles)    
        printRateSysts(dc=dc,
                       binName=cat,
                       procList=procList)

    #customize for the different data
    for dataURL in dataFiles: 
        fname=os.path.basename(dataURL)
        tag=fname.split('.')[1]
        if tag=='shapes' : tag='data'
        dcURL=os.path.join(opt.outDir,'%s.datacard.dat'%tag)
        os.system("sed s/_DATAOBSSHAPES_/%s/g < %s > %s"%(fname,dcTemplURL,dcURL))

if __name__ == "__main__":
    sys.exit(main())
