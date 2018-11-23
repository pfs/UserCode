import ROOT
import pickle
import os

def getEOSlslist(directory, mask='', prepend='root://eoscms//eos/cms/'):

    """
    Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
    """


    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'

    #eos should be mounted
    data = Popen(['ls', '/eos/cms/%s'%directory],stdout=PIPE)
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

def getChunksInSizeOf(chunkSize,directoryList,mask='',prepend='root://eoscms//eos/cms/'):
    
    """groups files in directory in chunks of a given size"""

    chunkList=[[]]
    cursize=0.
    for directory in directoryList:

        fList=getEOSlslist(directory=directory, mask=mask, prepend='/eos/cms')

        for f in fList:
            fsize=float(os.path.getsize(f))/(1024.e+6)
            if cursize+fsize>chunkSize: 
                chunkList.append([])
                cursize=fsize
            else:
                cursize += fsize
            chunkList[-1].append(f.replace('/eos/cms/',prepend))

    return chunkList
