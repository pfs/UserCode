import optparse
import json
import sys
from collections import OrderedDict
	
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default=None,              type='string')
    (opt, args) = parser.parse_args()

    #read lists of samples
    with open(opt.json,'r') as cachefile:
        samples=json.load(cachefile, encoding='utf-8', object_pairs_hook=OrderedDict).items()

    table={'data':[],'mc':[]}
    for tag,sample in samples:
        xsec,isData,dset=sample[0:3]
        table['data' if isData==1 else 'mc'].append((tag,xsec,dset.split('/')[1]))
        
    for tag in table:
        print '\n'
        print '\\multicolumn{%d}{c}{%s}\\\\'%(2 if tag=='data' else 3,tag)
        for name,xsec,dset in table[tag]:
            if tag=='data':
                print '%30s & %80s \\\\'%(name,dset.replace('_','\_'))
            else:
                print '%30s & %15s & %80s \\\\'%(name,str(xsec),dset.replace('_','\_'))


        
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

