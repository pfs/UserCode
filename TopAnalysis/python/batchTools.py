import commands
from time import sleep

"""
check periodically if a list of jobs is still running in the batch
"""
def babySitBatchJobs(jobNumbers,waitTime=10):
    while True:
        command_out = commands.getstatusoutput("bjobs | awk '{print $1}'")[1]
        runningJobs=[int(s) for s in command_out.split() if s.isdigit()]
        nJobsPending=0
        for j in jobNumbers:
            if j in runningJobs:
                nJobsPending+=1
        if nJobsPending==0:
            print 'All done for this step'
	    break
        else:
            print 'Querying batch in %ds (%d jobs missing)'%(waitTime,nJobsPending)
            sleep(waitTime)
