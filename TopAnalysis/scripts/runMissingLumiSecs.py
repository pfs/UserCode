#!/usr/bin/env python

import sys
import os
import json

def main():

      job=sys.argv[1]

      #parse the work area from the location of the directory
      WORKAREA='/'.join(job.split('/')[0:-1])
      if WORKAREA=='': WORKAREA='./'
      print 'Base work area is:',WORKAREA

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
      newWorkArea='%s_new'%WORKAREA
      newCfg='%s/%s'%(newWorkArea,os.path.basename(cfg.replace('_cfg','_ext_cfg')))
      os.system('mkdir -p %s'%newWorkArea)
      newCfgFile=open(newCfg,'w')
      lumiMaskSet=False
      for line in open(cfg,'r'):
            newLine=line      
            if 'workArea' in newLine      : newLine='config.General.workArea = \"%s\"\n'%newWorkArea
            if 'requestName' in newLine   : newLine=newLine[:-2]+"_ext\"\n"
            if 'lumiMask' in newLine      : 
                  newLine='config.Data.lumiMask = \"%s\"\n'%os.path.abspath(jsonF)
                  lumiMaskSet=True
            if 'unitsPerJob' in newLine   : newLine='config.Data.unitsPerJob = 5\n'
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
