import pickle
import os
import sys
from random import shuffle
from collections import defaultdict

rpData = defaultdict(list)

baseDir=sys.argv[1]
fList=[f for f in os.listdir(baseDir+'/mixing/Chunks') if 'Data' in f]

toCheck=[]
print 'Checking and merging',len(fList),'files'
for f in fList:

    if 'ZeroBias' in f : continue
    if 'DoubleEG' in f : continue
    if 'Photon' in f : continue

    fullf=os.path.join(baseDir+'/mixing/Chunks',f)

    try:
        with open(fullf,'r') as cache:
            a=pickle.load(cache)
        for key in a: 
            newKey=(key[0],int(key[1]),int(key[2]))
            rpData[newKey] += a[key]
    except Exception as e:
        toCheck.append(f)

print 'Total number of events for mixing'
for key in rpData:
    print key,len(rpData[key])

print 'Writing mixing bank @ %s/mixing/mixbank.pck'%baseDir
with open('mixbank.pck','w') as cache:
    pickle.dump(rpData,cache)
os.system('mv -v mixbank.pck %s/mixing/mixbank.pck'%baseDir)


if len(toCheck)>0:
    print '-'*50
    print 'Found these files with errors'
    print toCheck
    print '-'*50
