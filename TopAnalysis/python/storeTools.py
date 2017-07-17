import ROOT
import pickle

"""
Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
"""
def getEOSlslist(directory, mask='', prepend='root://eoscms//eos/cms/'):
    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'

    eos_cmd = 'eos'
    data = Popen([eos_cmd, 'ls', directory],stdout=PIPE)
    out,err = data.communicate()

    full_list = []

    ## if input file was single root file:
    if directory.endswith('.root'):
        if len(out.split('\n')[0]) > 0:
            return [prepend + directory]

    ## instead of only the file name append the string to open the file in ROOT
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        full_list.append(prepend + directory + '/' + line)

    ## strip the list of files if required
    if mask != '':
        stripped_list = [x for x in full_list if mask in x]
        return stripped_list

    ## return 
    return full_list
