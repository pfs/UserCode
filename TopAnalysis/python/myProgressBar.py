import sys

"""
simple progress bar
"""
def drawProgressBar(percent, barLen = 100):
    sys.stdout.write("\r")
    progress = ""
    for i in range(barLen):
        if i < int(barLen * percent):
            progress += "="
        else:
            progress += " "
    sys.stdout.write("[ %s ] %.1f%%" % (progress, percent * 100))
    sys.stdout.flush()
