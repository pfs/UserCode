#!/usr/bin/env python

import numpy as np
import optparse
import os
import sys
import pickle
import json
from collections import OrderedDict
from keras.models import load_model
import pandas as pd
import root_pandas as rp

def fitModels(data,features,opt,alwaysOptim=False,doGBC=False,doDNN=False):

    """ fits different models to the data """
    
    from sklearn import preprocessing
    #from sklearn.decomposition import PCA
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.ensemble import RandomForestClassifier
    from xgboost import XGBClassifier 
    from sklearn.model_selection import GridSearchCV
    from keras.models import Sequential
    from keras.layers.core import Dense, Dropout
    from keras.layers import BatchNormalization
    from keras.layers.advanced_activations import LeakyReLU
    from keras.optimizers import Adam
    from keras.callbacks import EarlyStopping,ModelCheckpoint, ReduceLROnPlateau

    #best_models={'lda':{}, 'rfc':{}, 'rfc10':{}, 'dnn':{}, 'gbc':{}, 'gbc10':{}}
    best_models={'lda':{}, 'rfc':{}, 'rfc10':{}}
    
    #preprocess by scaling to a zero mean/unit variance distribution hypothesis
    scaler = preprocessing.StandardScaler()
    data['X']=scaler.fit_transform(data['X'])
    #pca=PCA()
    #data['X']=pca.fit_transform(data['X'])

    xangle_list=data['s'][:,0]
    for xangle in [120,130,140,150]:

        filt=(xangle_list==xangle)
        X=data['X'][filt]
        y=data['y'][filt]
        print '[fitModels] @ crossing angle={0} murad has {1} events ({2} signal events)'.format(xangle,len(y),len(y[y==1]))

        #prepare the train (and test) dataset    
        df=pd.DataFrame(X)
        df.columns=features
        df['class']=y.tolist()
        train=df.sample(frac=opt.trainFrac,random_state=200)
        test=df.drop(train.index)
        
        #optimize a random forest classifier
        rfc = RandomForestClassifier()
        if alwaysOptim or xangle==120:
            rfc_params = {
                'bootstrap': [True],
                'n_estimators': [100,200,300,400],
                'max_depth': [5,10,15,20],
                #'min_samples_split':[0.25,0.5],
                #'min_samples_leaf':[0.25,0.5],
                'max_features':['sqrt'],
                }
            rfc_optim = GridSearchCV(estimator = rfc, param_grid = rfc_params, cv = 3, n_jobs = -1, verbose = 2)
            rfc_optim.fit(train[features], train['class'])
            rankedFeatures=sorted(zip(features,rfc_optim.best_estimator_.feature_importances_), key=lambda x:x[1],reverse=True)
            best_models['rfc'][xangle]=(rfc_optim.best_estimator_,rfc_optim.best_params_,rankedFeatures,features)
        else:
            print 'Re-using best parameters found for 120murad'
            rfc.set_params( **(best_models['rfc'][120][1]) )
            rfc.fit(train[features], train['class'])
            rankedFeatures=sorted(zip(features,rfc.feature_importances_), key=lambda x:x[1],reverse=True)
            best_models['rfc'][xangle]=(rfc,best_models['rfc'][120][1],rankedFeatures,features)

        #re-train using only best 10 features
        best10=[rankedFeatures[ix][0] for ix in range(10)]
        rfc10=RandomForestClassifier()
        rfc10.set_params(**rfc_optim.best_params_)
        rfc10.fit(train[best10],train['class'])
        ranked10Features=sorted(zip(best10,rfc10.feature_importances_), key=lambda x:x[1],reverse=True)
        best_models['rfc10'][xangle]=(rfc10,rfc_optim.best_params_,ranked10Features,best10)
    
        #train a linear discriminant using the two best ranked variables
        best2=[ranked10Features[ix][0] for ix in range(2)]
        lda=LinearDiscriminantAnalysis()
        lda_params={
            'solver':['lsqr'],
            'shrinkage':np.arange(0,1,0.1),
            }
        lda_optim = GridSearchCV(estimator = lda, param_grid = lda_params, cv = 3, n_jobs = -1, verbose = 2)
        lda_optim.fit(train[best2], train['class'])        
        best_models['lda'][xangle]=(lda_optim.best_estimator_,lda_optim.best_params_,None,best2)


        #optimize a gradient boost classifier
        if doGBC:
            gbc = XGBClassifier()
            gbc_params = {
                'learning_rate':[0.05,0.15],
                'n_estimators': [10,50,100,200],
                'max_depth': [2,5,10,15], #20,50],
                #'min_samples_split':[0.25,0.5],
                #'min_samples_leaf':[0.25,0.5],
                }
            gbc_optim = GridSearchCV(estimator = gbc, param_grid = gbc_params, cv = 3, n_jobs = -1, verbose = 2)
            gbc_optim.fit(train[features], train['class'])
            rankedFeatures=sorted(zip(features,gbc_optim.best_estimator_.feature_importances_), key=lambda x:x[1],reverse=True)
            best_models['gbc'][xangle]=(gbc_optim.best_estimator_,gbc_optim.best_params_,rankedFeatures,features)
        
            #re-train using only best 10 features
            best10=[rankedFeatures[ix][0] for ix in range(10)]
            gbc10 = XGBClassifier()
            gbc10.set_params(**gbc_optim.best_params_)
            gbc10.fit(train[best10],train['class'])
            ranked10Features=sorted(zip(best10,gbc10.feature_importances_), key=lambda x:x[1],reverse=True)
            best_models['gbc10'][xangle]=(gbc10,gbc_optim.best_params_,ranked10Features,best10)
            
        
        #train a DNN
        if doDNN:
            dnn = Sequential()

            def addCommonStructureBetweenDense(dnn):
                dnn.add(BatchNormalization())
                dnn.add(Dropout(0.1))
                dnn.add(LeakyReLU(0.2))

            dnn.add(Dense(512, kernel_initializer='glorot_normal',  bias_initializer='glorot_uniform', input_dim=len(features)))
            addCommonStructureBetweenDense(dnn)
            for i in [512,256,128]:
                dnn.add(Dense(i, kernel_initializer='glorot_normal',  bias_initializer='glorot_uniform'))
                addCommonStructureBetweenDense(dnn)
            dnn.add(Dense(1, kernel_initializer='glorot_normal',  bias_initializer='glorot_uniform', activation='sigmoid'))
        
            dnn.compile(loss='binary_crossentropy',optimizer=Adam(lr=0.01),metrics=['accuracy'])
            reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=3, min_lr=0.0001, cooldown=3, verbose=1)
            es = EarlyStopping(monitor='val_loss', mode='auto', verbose=1, patience=6)
            best_url=os.path.join(opt.output,'bestdnn_%d.hd5'%xangle)
            mc = ModelCheckpoint(best_url, monitor='val_loss', verbose=1, save_best_only=True)
            history = dnn.fit(train[features], train['class'], 
                              validation_data=(test[features], test['class']),
                              callbacks=[reduce_lr,es,mc], 
                              epochs=100,
                              batch_size=128
                              )
            best_models['dnn'][xangle]=(best_url,history.history,None,features)

    #save results  to pickle file
    pfix='ZeroBias' if opt.zeroBiasTrain else  'Zmm'
    out_url=os.path.join(opt.output,'pu_models_%s.pck'%pfix)
    with open(out_url,'w') as cache:    
        pickle.dump(best_models, cache, pickle.HIGHEST_PROTOCOL)
        pickle.dump(scaler,      cache, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(pca,         cache, pickle.HIGHEST_PROTOCOL)
    print 'Best fit models and preprocessing scaler have been stored in',out_url

def predict(data,features,baseName,opt):

    """ runs the prediction for the trained models and dumps a tree """
    
    print '[predict] with',baseName,'with',len(data),'events'

    #load models and standard scaler
    with open(opt.model,'r') as cache:
        best_models=pickle.load(cache)
        scaler=pickle.load(cache)

    #scale data and switch to pandas DataFrame
    df=pd.DataFrame(scaler.transform(data))
    df.columns=features

    #run all predictions
    pred=pd.DataFrame()
    for key in best_models:
        if key!='rfc': continue
        for xangle in best_models[key]:
            tag='%s_%d'%(key,xangle)
            print tag,'for',baseName
            clf=best_models[key][xangle][0]
            features=best_models[key][xangle][-1]
            y_prob=clf.predict_proba(df[features])[:,0]
            pred[tag]=y_prob

    #write to output
    rp.to_root(pred, baseName, key='pudiscr',store_index=False)
    if opt.output:        
        os.system('xrdcp -f {0} root://eoscms//{1}/{0}'.format(baseName,opt.output.replace('/eos/cms/','')))
        os.system('rm {0}'.format(baseName))


def runTrainJob(url,features,spectators,categs,onlyThis,opt):


    """ 
    converts the ROOT trees to numpy arrays for the 
    - features to use in the training
    - the spectator variables 
    - the classes to predict
    a selection string is used to filter the original events in the trees
    """

    #build the chain
    from ROOT import TChain
    
    t=TChain('tree')
    for f in os.listdir(url):
        if onlyThis:
            if f!=onlyThis: continue
        else:
            if opt.zeroBiasTrain:
                if not 'ZeroBias' in f: 
                    continue
            else:
                if not 'DoubleMuon' in f:
                    continue           
        t.AddFile(os.path.join(url,f))
    print 'Data chain has {0} events'.format(t.GetEntries())

    #convert to numpy arrays
    from root_numpy import tree2array

    data={}
    cut=opt.selection if opt.selection else ''
    if len(cut) : print cut
    step=100 if opt.zeroBiasTrain else None
    data['X']=tree2array(t, branches=features,   selection=cut, step=step)
    data['s']=tree2array(t, branches=spectators, selection=cut, step=step)
    data['y']=tree2array(t, branches=categs,     selection=cut, step=step)

    #convert features and spectators to array of floats
    for tag,branches in [('X',features),('s',spectators)]:
        data[tag]=data[tag].astype([(b,'<f4') for b in branches])
        data[tag]=data[tag].view(np.float32).reshape(data[tag].shape + (-1,))

    #filter out runs in which the RP were out
    nEntries=len(data['s'])
    print nEntries,'events available at start'
    if opt.RPout:
        with open(opt.RPout,'r') as cache:

            #read RP out json file
            runLumi=json.load(cache,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
            runLumiList={int(x[0]):x[1] for x in runLumi}

            #build a filter
            filt=[]
            nFilt=0
            print 'Filtering out runs in which the RP were out'
            for i in range(nEntries):

                if i%10000==0 : 
                    fracDone=int(float(100.*i)/float(nEntries))
                    fracFilt=int(float(100.*nFilt)/float(nEntries))
                    sys.stdout.write('\r [ %d/100 ] done [ %d/100 ] to be filtered' %(fracDone,fracFilt))

                run,lumi=int(data['s'][i][1]),int(data['s'][i][2])
                rpInFlag=True
                if run in runLumiList:                 
                    for lran in runLumiList[run]:
                        if lumi<lran[0]: continue
                        if lumi>lran[1]: continue
                        rpInFlag=False
                        nFilt+=1
                        break
            
                filt.append(rpInFlag)

            #apply filter
            filt=np.array(filt)
            for key in data: data[key]=data[key][filt]

            print '%d%% events removed as RP were out of the run'%fracFilt

    print 'Converted to numpy array' 

    #using imputer to fill NaN with the median
    print 'Checking/substituting for NaN in features with the Imputer'
    from sklearn.preprocessing import Imputer
    imputer = Imputer(strategy="median",verbose=1)
    data['X']=imputer.fit_transform(data['X'])

    if opt.model is None:

        pfix='ZeroBias' if opt.zeroBiasTrain else  'Zmm'
        out_url=os.path.join(opt.output,'train_data_%s.pck'%pfix)
        with open(out_url,'w') as cache:    
            pickle.dump(data,          cache, pickle.HIGHEST_PROTOCOL)
            pickle.dump(features,      cache, pickle.HIGHEST_PROTOCOL)
            pickle.dump(spectators,    cache, pickle.HIGHEST_PROTOCOL)
            pickle.dump(categs,        cache, pickle.HIGHEST_PROTOCOL)
            pickle.dump(opt.selection, cache, pickle.HIGHEST_PROTOCOL)
        print 'A dump of the data in numpy format is saved @',out_url

        fitModels(data,features,opt)
        
    else:

        predict(data['X'],features,onlyThis,opt)


def runTrainJobPacked(args):
    
    """wrapper for parallel threads"""

    try:
        runTrainJob(*args)
    except Exception as e:
        print '>'*50
        print e
        print '>'*50
        return False

def main():


    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                      dest='input',   
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/Chunks/',
                      help='input directory with the files [default: %default]')
    parser.add_option('--RPout',
                      dest='RPout', 
                      default=None,
                      type='string',
                      help='json with the runs/lumi sections in which RP are out')
    parser.add_option('--trainFrac',
                      dest='trainFrac',
                      default=0.5,
                      type=float,
                      help='fraction to use for training [default: %default]')
    parser.add_option('-o', '--output',
                      dest='output',
                      default=None,
                      help='output directory [default: %default]')
    parser.add_option('--onlyMissing',
                      dest='onlyMissing',
                      default=False,
                      action='store_true',
                      help='run on only missing [default: %default]')
    parser.add_option('--zeroBiasTrain',
                      dest='zeroBiasTrain',
                      default=False,
                      action='store_true',
                      help='use zero bias events for the training [default: %default]')
    parser.add_option('-s', '--selection',
                      dest='selection',   
                      default=None,
                      help='selection [default: %default]')
    parser.add_option('-m', '--model',
                      dest='model',   
                      default=None,
                      help='pickle file with trained models [default: %default]')
    (opt, args) = parser.parse_args()

    if opt.output: os.system('mkdir -p %s'%opt.output)

    #build the feature,spectator, category
    features=['nvtx','rho']
    for v in ['PFMult','PFPz','PFHt','PFChMult','PFChPz','PFChHt']:
        for r in ['EB','EE','HE','HF']:
            if 'Ch' in v and 'H' in r : 
                continue
            for o in ['Sum','Diff']:
                features.append(v+o+r)
    #keep order for beamXangle,run,lumi in the first position as it will be used to filter events...
    spectators=['beamXangle','run','lumi','nrawmu-2']  
    categs='trainCat'

    #fit directly the data if pickle file is available
    if '.pck' in opt.input:
        with open(opt.input,'r') as cache:
            data=pickle.load(cache)
        fitModels(data,features,opt)
        return 

    #build the tasks
    task_list=[]
    buildTrainData=True
    if opt.model is None:
        task_list.append( (opt.input, features, spectators, categs, None, opt) )    
    else:
        if os.path.isfile(opt.input): 
            task_list.append( (os.path.dirname(opt.input), features, spectators, categs, os.path.basename(opt.input), opt) )
        else:
            for f in os.listdir(opt.input):
                if not os.path.isfile(os.path.join(opt.input,f)) : continue
                outLoc=os.path.join(opt.output,f)
                if not '/eos/cms' in outLoc : outLoc='/eos/cms/'+outLoc
                if opt.onlyMissing and os.path.isfile(outLoc): continue
                task_list.append( (os.path.dirname(f), features, spectators, categs, os.path.basename(f), opt) )

    #run it
    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map(runTrainJobPacked, task_list)

if __name__ == "__main__":
    sys.exit(main())
