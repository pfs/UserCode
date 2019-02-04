import pickle
import os
import sys
from random import shuffle
from collections import defaultdict

rpData = defaultdict(list)
scanDir=sys.argv[1]
for f in os.listdir(scanDir):
    if not '.pck' in f: continue
    if not 'Data' in f: continue
    with open(os.path.join(scanDir,f),'r') as cache:
        a=pickle.load(cache)
        for key in a: rpData[key] += a[key]

print 'Total number of events for mixing'
for key in rpData:
    print key,len(rpData[key])

with open(os.path.join(scanDir,'mixbank.pck'),'w') as cache:
    pickle.dump(rpData,cache)
