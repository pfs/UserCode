#!/usr/bin/env python

"""
specifications for systematic uncertainties
"""
def xsecSystSpecs(analysis='TOP-16-006'):

    rateSysts,sampleSysts=[],[]

    #TOP-16-006 specific
    if analysis=='TOP-16-006':
        rateSysts=[
            ('lumi_13TeV',       1.027,    'lnN',    []                   ,['Multijetsdata']),
            #('DYnorm_th',        1.038,    'lnN',    ['DYl','DYc','DYb']  ,[]),
            #('Wnorm_th',         1.037,    'lnN',    ['Wl' ,'Wc','Wb']    ,[]),
            ('DYnorm_th',        1.038,    'lnN',    ['DY']  ,[]),
            ('Wnorm_th',         1.037,    'lnN',    ['W']   ,[]),
            ('tWnorm_th',        1.054,    'lnN',    ['tW']               ,[]),
            ('tnorm_th',         1.044,    'lnN',    ['tch']              ,[]),
            ('VVnorm_th',        1.20,     'lnN',    ['Multiboson']       ,[]),
            ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV']           ,[]),
            ]

        sampleSysts=[
            #ttbar modelling
            ('Mtop',            {'tbart'         : ['tbartm=169.5','tbartm=175.5'],  'tW':['tWm=169.5','tWm=175.5'] },                True ,  True, False),
            ('ttPartonShower',  {'tbart'         : ['tbartscaledown','tbartscaleup']},                                                False , True, False),            
            ('NLOgenerator',    {'tbart'         : ['tbartaMCNLO']},                                                                  False,  True, False),
            ('Hadronizer',      {'tbart'         : ['tbartHerwig']},                                                                  False , True, True),

            #tWinterference
            ('tWttinterf',       {'tW'            : ['tWDS']},                                                                        False , True, True),            

            #QCD SCALES
            ('tWscale',         {'tW'            : ['tWscaledown','tWscaleup']},                                                      False , True, False),            

            #Madgraph W+jets
            #('wFactScale',           { 'Wl': ['id3mur1muf0.5','id2mur1muf2'], 
            #                           'Wc': ['id3mur1muf0.5','id2mur1muf2'], 
            #                           'Wb': ['id3mur1muf0.5','id2mur1muf2'] },  False, False),
            #('wRenScale',            { 'Wl': ['id7mur0.5muf1','id4mur2muf1'],  
            #                           'Wc': ['id7mur0.5muf1','id4mur2muf1'],  
            #                           'Wb': ['id7mur0.5muf1','id4mur2muf1'] },  False, False),
            #('wCombScale',           { 'Wl': ['id9mur0.5muf0.5','id5mur2muf2'], 
            #                           'Wc': ['id9mur0.5muf0.5','id5mur2muf2'],
            #                           'Wb': ['id9mur0.5muf0.5','id5mur2muf2'] },  False, False),

            #amc@NLO W+jets
            ('wFactScale',           { 'W': ['id1003muR0.10000E+01muF0.50000E+00','id1002muR0.10000E+01muF0.20000E+01'] },  False, False, False),
            ('wRenScale',            { 'W': ['id1007muR0.50000E+00muF0.10000E+01','id1004muR0.20000E+01muF0.10000E+01'] },  False, False, False),
            ('wCombScale',           { 'W': ['id1009muR0.50000E+00muF0.50000E+00','id1005muR0.20000E+01muF0.20000E+01'] },  False, False, False),
            
            #ttbar Powheg
            ('ttFactScale',          { 'tbart': ['muR1muF0.5hdampmt172.5',  'muR1muF2hdampmt172.5'] },     True , False, False),
            ('ttRenScale',           { 'tbart': ['muR0.5muF1hdampmt172.5',  'muR2muF1hdampmt172.5'] },     True , False, False),
            ('ttCombScale',          { 'tbart': ['muR0.5muF0.5hdampmt172.5','muR2muF2hdampmt172.5'] },     True , False, False),
            ]

    if analysis=='TOP-16-015':

        rateSysts=[
            ('lumi_5TeV',        1.023,    'lnN',    []                   ,['Multijetsdata']),
            ('DYnorm',           1.30,     'lnN',    ['DY']  ,[]),
            ('Wnorm_th',         1.037,    'lnN',    ['W']   ,[]),
            ('tWnorm',           1.30,     'lnN',    ['tW']               ,[]),
            ('VVnorm',           1.30,     'lnN',    ['Multiboson']       ,[])
            ]

        sampleSysts=[

            #ttbar modelling
            ('ttPartonShower',  {'tbart': ['tbartPSup','tbartPSdown']}, False , True, False),            
            ('Hadronizer',      {'tbart': ['tbartHW']},                 False , True, True),

            #QCD scales for ttbar
            ('ttFactScale',          { 'tbart': ['genUnc1','genUnc2'] }, True, False, False),
            ('ttRenScale',           { 'tbart': ['genUnc3','genUnc6'] }, True, False, False),
            ('ttCombScale',          { 'tbart': ['genUnc4','genUnc8'] }, True, False, False),

            #QCD scales for W+jets
            ('wFactScale',           { 'W': ['genUnc1','genUnc2'] },   False, False, False),
            ('wRenScale',            { 'W': ['genUnc3','genUnc6'] },   False, False, False),
            ('wCombScale',           { 'W': ['genUnc4','genUnc8'] },   False, False, False),
            ]

    if analysis=='TopRadius':
        rateSysts=[]
        sampleSysts=[]


    #return specifications
    return rateSysts,sampleSysts

