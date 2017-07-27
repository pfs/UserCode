import ROOT
import sys

from roofitTools import *

ALLPDFS={
('S_cor','mjj'):
    [
        #'RooBifurGauss::S_cor1mjj_{0}(mjj,mu_scormjj_{0}[20,100],sigmaL_scormjj_{0}[1,14],sigmaR_scormjj_{0}[1,14])',
        
        'RooCBShape::S_cor1mjj_{0}(mjj,mu_scormjj_{0}[75,90],sigma_scormjj_{0}[5,13],alpha_scormjj_{0}[0.01,20],n_scormjj_{0}[0])',
        'RooGamma::S_cor2mjj_{0}(mjj,gamma_scormjj_{0}[50.,0.1,100],beta_scormjj_{0}[5,0.01,100],lambda_scormjj_{0}[0])',
        'SUM::S_cormjj_{0}( f_scormjj_{0}[0.1,0.98]*S_cor1mjj_{0},S_cor2mjj_{0} )'
    ],
('S_wro','mjj'):
    [
        #'RooCBShape::S_wro1mjj_{0}(mjj,mu_swromjj_{0}[25,250],sigma_swromjj_{0}[5,13],alpha_swromjj_{0}[0.01,20],n_swromjj_{0}[0])',
        #'RooGamma::S_wro2mjj_{0}(mjj,gamma_swromjj_{0}[50.,0.1,100],beta_swromjj_{0}[5,0.01,100],lambda_swromjj_{0}[1])',
        'RooBifurGauss::S_wro1mjj_{0}(mjj,mu_swromjj_{0}[25,75],sigmaL_swromjj_{0}[15,100],sigmaR_swromjj_{0}[15,100])',
        'RooLandau::S_wro2mjj_{0}(mjj,mpv_swromjj_{0}[25,200],width_swromjj_{0}[50,5,100])',
        'SUM::S_wromjj_{0}( f_swromjj_{0}[0,1]*S_wro1mjj_{0}, S_wro2mjj_{0} )'
    ],
('B','mjj'):
    [
    'RooLandau::B_mjj_{0}(mjj,mpv_bmjj_{0}[20,300],width_bmjj_{0}[20,0.01,100])'
    ],
('S_cor','mthad'):
    [
     'RooBifurGauss::S_cor1mthad_{0}(mthad,mu_scormthad_{0}[100,200],sigmaL_scormthad_{0}[10,20],sigmaR_scormthad_{0}[10,20])',
     'RooGamma::S_cor2mthad_{0}(mthad,gamma_scormthad_{0}[50,0.1,100],beta_scormthad_{0}[0.5,0.01,100],lambda_scormthad_{0}[0])',
     'SUM::S_cormthad_{0}( f_scormthad_{0}[0,1]*S_cor1mthad_{0}, S_cor2mthad_{0} )'
    ],
('S_wro','mthad'):
    [
     'RooBifurGauss::S_wro1mthad_{0}(mthad,mu_swromthad_{0}[100,200],sigmaL_swromthad_{0}[20,100],sigmaR_swromthad_{0}[20,100])',
     'RooLandau::S_wro2mthad_{0}(mthad,mpv_swromthad_{0}[100,200],width_swromthad_{0}[5,100])',
     'SUM::S_wromthad_{0}( f_swromthad_{0}[0,1]*S_wro1mthad_{0}, S_wro2mthad_{0} )'
    ],
('B','mthad'):
    [
      'RooLandau::B_mthad_{0}(mthad,mpv_bmthad_{0}[20,300],width_bmthad_{0}[20,0.01,100])'
    ],
('S_cor','mtlep'):
    [
     'RooBifurGauss::S_cor1mtlep_{0}(mtlep,mu_scormtlep_{0}[160,180],sigmaL_scormtlep_{0}[10,30],sigmaR_scormtlep_{0}[10,40])',
     'RooGamma::S_cor2mtlep_{0}(mtlep,gamma_scormtlep_{0}[50,0.1,100],beta_scormtlep_{0}[0.5,0.01,100],lambda_scormtlep_{0}[0])',
     'SUM::S_cormtlep_{0}( f_scormtlep_{0}[0,1]*S_cor1mtlep_{0}, S_cor2mtlep_{0} )'
    ],
('S_wro','mtlep'):
        [
         'RooBifurGauss::S_wro1mtlep_{0}(mtlep,mu_swromtlep_{0}[20,180],sigmaL_swromtlep_{0}[5,100],sigmaR_swromtlep_{0}[5,100])',
         'RooLandau::S_wro2mtlep_{0}(mtlep,mpv_swromtlep_{0}[20,300],width_swromtlep_{0}[20,10,100])',
         'SUM::S_wromtlep_{0}( f_swromtlep_{0}[0,1]*S_wro1mtlep_{0}, S_wro2mtlep_{0} )'
        ],
('B','mtlep'):
        [
          'RooLandau::B_mtlep_{0}(mtlep,mpv_bmtlep_{0}[20,300],width_bmtlep_{0}[20,0.01,100])'
        ],
}


"""
"""
def parametrize(data,varName,w,pdfDef):

    #parameterize baseShape
    for p in pdfDef: pdf=w.factory(p)

    #fit PDF
    pdf.fitTo(data)

    #showFitResult(fitVar=varName,
    #              data=data,
    #              pdf=pdf,
    #              categs=[''],
    #              w=w,
    #             tagTitle=data.GetName(),
    #              rangeX=(0,500),
    #              showComponents=[('*1*','1'),('*1*,*2*','2')],
    #              outDir='/dev/null',
    #              paramList=[])
    #raw_input()

    #freeze parameters
    iter = pdf.getParameters(data).createIterator()
    iparam = iter.Next()
    while iparam :
        pname=iparam.GetName()
        title='f'
        if 'alpha' in pname : title='#alpha'
        if 'beta' in pname  : title='#beta'
        if 'gamma' in pname : title='#gamma'
        if 'lambda' in pname : title='#lambda'
        if 'mu' in pname    : title='#mu'
        if 'sigma' in pname :
            if 'sigmaL' in pname   : title='#sigma^{L}'
            elif 'sigmaR' in pname : title='#sigma^{R}'
            else                   : title='#sigma'
        if 'mpv' in pname          : title='MPV'
        if 'width' in pname        : title='w'
        if 'cor' in pname   : title += '_{cor}'
        if 'wro' in pname   : title += '_{wro}'
        iparam.SetTitle(title)
        iparam.setConstant(True)
        iparam = iter.Next()

    return pdf

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #get the workspace from the file
    url=sys.argv[1]
    fIn=ROOT.TFile.Open(url)
    w=fIn.Get('w')
    fIn.Close()

    doSignal=True
    if not 'TTbar' in url : doSignal=False

    waux=None
    try:
        fIn=ROOT.TFile.Open(sys.argv[2])
        waux=fIn.Get('w')
        fIn.Close()
    except:
        pass


    data=w.data('data')
    corSummary={}
    for varName,varTitle in [('mjj','M(jj)'),
                             ('mthad','M(t_{had})'),
                             ('mtlep','M(t_{lep})')
                             ]:
        w.var(varName).SetTitle(varTitle)
        for ch in ['e','mu']:
            for jcount in ['1l4j2b','1l4j2q','1l4j1b1q']:
                cat=ch+jcount
                catCut="sample==sample::%s"%cat

                if not doSignal:

                    #data
                    redData=data.reduce(ROOT.RooFit.Cut(catCut))
                    w.cat('sample').setLabel(cat)

                    #define the PDF
                    pdfDef=[x.format(cat) for x in  ALLPDFS[('B',varName)] ]
                    pdf=parametrize(redData,varName,w,pdfDef)

                    pdfAux=None
                    try:
                        pdfAux=waux.pdf('S_wro%s_%s'%(varName,cat))
                    except:
                        pass

                    showFitResult(fitVar=varName,
                            data=redData,
                            pdf=pdf,
                            categs=[''],
                            w=w,
                            tagTitle=cat,
                            rangeX=(0,500),
                            showComponents=[],
                            outDir='./',
                            paramList=[],
                            pdfAux=pdfAux)
                else:

                    for c in ['cor','wro']:

                        #project the data
                        cut=catCut
                        if c=='cor':
                            if varName=='mjj'   : cut += ' && corwjj==1'
                            if varName=='mthad' : cut += ' && corthad==1'
                            if varName=='mtlep' : cut += ' && cortlep==1'
                        if c=='wro':
                            if varName=='mjj'   : cut += ' && corwjj==0'
                            if varName=='mthad' : cut += ' && corthad==0'
                            if varName=='mtlep' : cut += ' && cortlep==0'
                        redData=data.reduce(ROOT.RooFit.Cut(cut))
                        w.cat('sample').setLabel(cat)

                        #define the PDF
                        pdfDef=[x.format(cat) for x in  ALLPDFS[('S_%s'%c,varName)] ]
                        pdf=parametrize(redData,varName,w,pdfDef)

                        if c=='cor':
                            mean,median,mode,effWid=getEffSigma(w.var(varName),pdf,0,200)
                            corSummary[(cat,varName)]=(mean,median,mode,effWid)

                    #combine correct+wrong combinations
                    redData=data.reduce(catCut)
                    pdf=w.factory('SUM::S_{0}_{1}( f_s{0}_{1}[0,1]*S_cor{0}_{1}, S_wro{0}_{1} )'.format(varName,cat))
                    pdf.fitTo(redData)
                    showFitResult(fitVar=varName,
                                data=redData,
                                pdf=pdf,
                                categs=[''],
                                w=w,
                                tagTitle=cat,
                                rangeX=(0,500),
                                showComponents=[('S_wro*','wrong permutations'),
                                                ('S_cor*,S_wro*','correct permutations')],                                
                                outDir='./',
                                paramList=[])

    if doSignal:
        with open(url.replace('.root','.dat'),'w') as fOut:
            fOut.write('Category | Variable | Eff. width | Mean | Median\n')
            for cat,varName in corSummary:
                fOut.write('%s %s '%(cat,varName))
                for x in corSummary[(cat,varName)]:
                    fOut.write('%3.1f '%x)
                fOut.write('\n')

    w.writeToFile('pdf_'+url,True)

"""
"""
if __name__ == "__main__":
    main()
