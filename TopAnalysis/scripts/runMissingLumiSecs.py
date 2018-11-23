#!/usr/bin/env python

import das_client
import sys
import os
import json

WORKAREA='grid'

def getListOfFiles(inputPD,runSel):
      das_query_string = 'das_client.py --query="file dataset=%s  run in [%s] | grep file.name" --limit=0' % (inputPD,runSel)
      das_output = os.popen(das_query_string).read() #query DAS via bash
      das_output = das_output.replace("\"\"","'")
      das_output = das_output.replace("root","root',")
      das_output = das_output.replace("\n","") #don't want newline in string
      das_output = das_output[:-1] #remove last comma in list of secondary files
      return das_output

def main():
      job=sys.argv[1]
      #check if any lumi section is missing
      jsonF='%s/results/notFinishedLumis.json'%job
      runSel=[]
      try:
            with open(jsonF) as missingLumis :
                  data = json.load(missingLumis)
                  runSel=[int(x) for x in data.keys()]
      except:
            sys.exit(0) 
      if len(runSel)==0 : sys.exit(0)

      print 'Starting',job
      print '\t %d runs with missing lumi sections'%len(runSel)
      cfg=job.replace('crab_','')+'_cfg.py'      
      newCfg='%s/%s'%(WORKAREA,os.path.basename(cfg.replace('_cfg','_ext_cfg')))
      os.system('mkdir -p grid_new')
      newCfgFile=open(newCfg,'w')
      lumiMaskSet=False
      for line in open(cfg,'r'):
            newLine=line      
            if 'workArea' in newLine      : newLine='config.General.workArea = \"%s\"\n'%WORKAREA
            if 'requestName' in newLine   : newLine=newLine[:-2]+"_ext\"\n"
            if 'lumiMask' in newLine      : 
                  newLine='config.Data.lumiMask = \"%s\"\n'%os.path.abspath(jsonF)
                  lumiMaskSet=True
            if 'unitsPerJob' in newLine   : newLine='config.Data.unitsPerJob = 25\n'
            if 'outLFNDirBase' in newLine : newLine=newLine[:-3]+"_ext\"\n"
            newCfgFile.write( newLine )

      if not lumiMaskSet:
            print 'Setting lumi mask now'
            newCfgFile.write('config.Data.lumiMask = \"%s\"\n'%os.path.abspath(jsonF))

      newCfgFile.close()
      print '\t new cfg to process missing lumis @',newCfg
         
      os.system('crab submit -c %s' % newCfg)
      print '\t jobs have been submitted'

if __name__ == "__main__":
    main()
