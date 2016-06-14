import os
import sys
import optparse
import ROOT
from createDataCard import getDistsFrom

"""
extrapolate linearly the yields at a given value
"""
def linearYieldExtrapolation(node0,node1,val):
    #extrapolate linearly the yields
    a=(node1[1].Integral()-node0[1].Integral())/(node1[0]-node0[0])
    b=node0[1].Integral()-a*node0[0]
    return (a*val+b)


"""
gets a morphed shape at a given value
"""
def getMorphedShape(name,node0,node1,val):


    # use a linear interpolation NIM A 425 (1999) 357-360
    # for 2D distributions do this for each projection in Y
    if node0[1].InheritsFrom('TH2'):
        h=node0[1].Clone(name)
        maxbins=ROOT.TMath.Min(node0[1].GetYaxis().GetNbins(),node1[1].GetYaxis().GetNbins())+1
        for ybin in xrange(1,maxbins):
            px0=node0[1].ProjectionX('px0',ybin,ybin)
            px1=node1[1].ProjectionX('px1',ybin,ybin)
            newNorm=linearYieldExtrapolation( (node0[0],px0), (node1[0],px1), val)
            if newNorm>0 :
                pxint=ROOT.th1fmorph('pxint','pxint',px0,px1,node0[0],node1[0],val,newNorm)
                for xbin in xrange(1,pxint.GetXaxis().GetNbins()+1):
                    h.SetBinContent(xbin,ybin,pxint.GetBinContent(xbin))
                    h.SetBinError(xbin,ybin,pxint.GetBinError(xbin))
                pxint.Delete()
            px0.Delete()
            px1.Delete()
    else:
        newNorm=linearYieldExtrapolation(node0,node1,val)
        h=ROOT.th1fmorph(name,name,node0[1],node1[1],node0[0],node1[0],val,newNorm)
    h.SetDirectory(0)

    #all done
    return h


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputs',         dest='inputs',      help='input plotters and masses (csv)', default=None,                        type='string')
    parser.add_option(      '--ttxsec',         dest='ttxsec',      help='ttbar xsec used by default',      default=831.76,                      type=float)
    parser.add_option(      '--scan',           dest='scan',        help='scan values (csv)',               default=None,                        type='string')
    parser.add_option('-d', '--dist',           dest='dist',        help='distribution',                    default='nbtags',                    type='string')
    parser.add_option(      '--cats',           dest='cats',        help='categories',                      default='1j,2j,3j,4j',               type='string')
    parser.add_option('-o', '--out',            dest='out',         help='output file',                     default='stealth_stop_plotter.root', type='string')
    (opt, args) = parser.parse_args()

    #parse the inputs
    nodes=[]
    nodesList=opt.inputs.split(',')
    for nodeInfo in nodesList:
        mass,fname=nodeInfo.split(':')
        nodes.append( (float(mass),ROOT.TFile.Open(fname)) )
    nodes=sorted(nodes, key=lambda node: node[0])


    #mass scan
    minMass, maxMass, massStep = [ float(x) for x in opt.scan.split(',') ]
    isteps=int((maxMass-minMass)/massStep)+1
    
    #categories
    cats = [ x for x in opt.cats.split(',') ]

    #load combine library and interpolate shapes
    ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit.so')
    shapesToSave={}
    for i in xrange(0,isteps):
        mass=minMass+i*massStep
        dists=[]
        newName='stealthstopm%3.0f'%(10*mass)
        for j in xrange(0,len(nodes)):

            #mass coincides with node
            if mass==nodes[j][0]:
                print mass,'will be taken from node',j
                keyFilter='t#bar{t}'
                if mass==169.5 : keyFilter += ' m=169.5'
                if mass==175.5 : keyFilter += ' m=175.5'
                for cat in cats:
                    for histo in ['%s_%s'%(opt.dist,cat),
                                  '%sshapes_%s_exp'%(opt.dist,cat),
                                  '%sshapes_%s_gen'%(opt.dist,cat)]:

                        if not histo in shapesToSave: shapesToSave[histo]=[]
                        _,exp=getDistsFrom(directory=nodes[j][1].Get(histo),keyFilter=keyFilter )
                        shapesToSave[histo].append( exp[ exp.keys()[0] ].Clone(newName) )
                        #shapesToSave[histo][-1].Scale(1./opt.ttxsec)
                        shapesToSave[histo][-1].SetDirectory(0)

            #mass has to be interpolated
            elif mass>=nodes[j][0] and mass<nodes[j+1][0]: 
                print mass,'will be interpolated from',nodes[j][0],'-',nodes[j+1][0]
                keyFilter0,keyFilter1='t#bar{t}','t#bar{t}'
                if nodes[j][0]==169.5 : keyFilter0 += ' m=169.5'
                if nodes[j][0]==175.5 : keyFilter0 += ' m=175.5'
                if nodes[j+1][0]==169.5 : keyFilter1 += ' m=169.5'
                if nodes[j+1][0]==175.5 : keyFilter1 += ' m=175.5'
                for cat in cats:
                    for histo in ['%s_%s'%(opt.dist,cat),
                                  '%sshapes_%s_exp'%(opt.dist,cat),
                                  '%sshapes_%s_gen'%(opt.dist,cat)
                                  ]:
                        
                        if not histo in shapesToSave: shapesToSave[histo]=[]
                        _,exp0=getDistsFrom(directory=nodes[j][1].Get(histo),   keyFilter=keyFilter0 )
                        _,exp1=getDistsFrom(directory=nodes[j+1][1].Get(histo), keyFilter=keyFilter1 )
                        shapesToSave[histo].append( getMorphedShape(newName,
                                                                    (nodes[j][0],   exp0[ exp0.keys()[0] ] ),
                                                                    (nodes[j+1][0], exp1[ exp1.keys()[0] ] ),
                                                                    mass) )
                        #shapesToSave[histo][-1].Scale(1./opt.ttxsec)
                        

    #save to file
    fOut=ROOT.TFile.Open(opt.out,'RECREATE')
    for histo in shapesToSave:
        fOut.cd()
        outDir=fOut.mkdir(histo)
        outDir.cd()
        for h in shapesToSave[histo]: h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    fOut.Close()
    


"""                                                                                                                                                                                                               
for execution from another script                                                                                                                                                                           
"""
if __name__ == "__main__":
    sys.exit(main())

