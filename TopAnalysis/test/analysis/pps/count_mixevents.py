from MixedEventSummary import MixedEventSummary
import os
import pickle

d='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis/mixing/'
flist=[os.path.join(d,f) for f in os.listdir(d) if '.pck' in f]

nperangle={120:0,130:0,140:0,150:0}
for f in flist:
    print '\t',f
    with open(f,'r') as cachefile:
        rpData=pickle.load(cachefile)
        for key in rpData:
            fs=key[2]
            if fs!=169: continue
            xangle=key[1]
            nperangle[xangle] += len(rpData[key])

print(nperangle)
