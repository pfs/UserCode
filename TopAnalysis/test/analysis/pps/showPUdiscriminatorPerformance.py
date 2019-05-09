#!/usr/bin/env python

import sys
import pickle
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
from keras.models import load_model
from matplotlib import pyplot as plt
from importanceViaColumnPermutation import *

def plotFeatureImportance(models,outName):
    
    """ plots the feature importance """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    features=[]
    for _,m_feats in models:
        for f,_ in m_feats:
            features.append(f)
    features=list(set(features))

    feature_importance={}
    for m,m_feats in models:
        feature_importance[m]=[0.]*len(features)
        for f,fval in m_feats:
            idx=features.index(f)
            feature_importance[m][idx]=fval

    y_pos=np.arange(len(features))
    for m in feature_importance:
        ax.barh(y_pos,feature_importance[m],align='center',alpha=0.5,label=m)
    plt.yticks(y_pos,features)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tick_params(axis='both', which='minor', labelsize=8)

    plt.xlabel('Feature importance')
    fig.text(0.20, 0.96, r'CMS', fontweight='bold')
    fig.text(0.26, 0.96, r'Preliminary')
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.1)
    plt.legend(loc='upper right')
    plt.grid()
    for ext in ['.png','.pdf']:
        plt.savefig(outName+ext)    
    plt.clf()

def showFeatureCorrelation(data,features,truth,outfile):

    #compute correlation
    corr_bkg=data.loc[df[truth]==0][features].corr(method='pearson')
    corr_sig=data.loc[df[truth]==1][features].corr(method='pearson')
    corr_ratio=corr_bkg/corr_sig

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(corr_ratio,cmap=plt.cm.jet,vmin=-2,vmax=2)
    ax.xaxis.tick_bottom()
    clb=fig.colorbar(cax)
    clb.set_label(r'$\rho$(B) / $\rho$(S)',rotation=270.,labelpad=15)
    plt.xlabel('Variable index')
    plt.ylabel('Variable index')
    fig.text(0.15, 0.9, r'CMS', fontweight='bold')
    fig.text(0.22, 0.9, r'Preliminary')
    plt.grid()
    for ext in ['.png','.pdf']:
        plt.savefig(outfile+ext)    
    plt.clf()


def showConfusionMatrix(clfList,data,truth,outfile):
    
    """displays the confusion matrix"""

    from sklearn.metrics import confusion_matrix

    for m,clf,features in clfList:
        if isinstance(clf,basestring) : clf=load_model(clf)     
        y_pred = clf.predict(data[features])
        conf = confusion_matrix(data[truth], y_pred)


def showROCs(clfList,data,truth,outfile):

    """ show final ROC curves"""

    from sklearn.metrics import auc, roc_curve

    bestThr=-1

    fig=plt.figure(111)
    plt.plot([0,1], [0,1], ':', color='gray')

    roc_curve_list=[]

    for m,clf,features in clfList:
        try:
            invert=False

            #special case for DNN which comes stored in hd5 file
            if isinstance(clf,basestring) : 
                clf=load_model(clf)     
                invert=True
            y_prob = clf.predict_proba(data[features])[:,0]
            if invert:
                y_prob=1-y_prob

        except Exception as e:
            print '<'*50
            print m
            print e
            print '<'*50
            continue
        fpr, tpr, thresholds = roc_curve(data[truth],y_prob,pos_label=0)
        aucVal = 1-auc(tpr, fpr)
        roc_curve_list.append( (m,aucVal,fpr, tpr, thresholds) )
        if m=='rfc':
            data['rfc']=y_prob

    for v in ['nvtx','rho','nrawmu']:
        fpr,tpr,thresholds=roc_curve(data['class'],data[v],pos_label=0)
        fpr=1-fpr # these variables are cut by upper value v<max_val => signal-like
        tpr=1-tpr
        aucVal=1-auc(tpr,fpr)
        roc_curve_list.append( (v,aucVal,fpr, tpr, thresholds) )
  
    roc_curve_list.sort(key=lambda x: x[1], reverse=True)        
    for m, aucVal,fpr, tpr, thresholds in roc_curve_list:
        idx=min(range(len(fpr)), key=lambda i: abs(fpr[i]-0.10))
        ls='-'
        if m in ['nvtx','rho','nrawmu'] : ls ='--'
        plt.plot(tpr, fpr, linestyle=ls, label='{0} auc = {1:.3f}'.format(m,aucVal))
        if m=='rfc':
            #plt.plot([0.25,0.25,0], [0,fpr[idx],fpr[idx]], linewidth=2, color='gray')
            plt.plot([tpr[idx],tpr[idx],0], [0,fpr[idx],fpr[idx]], linewidth=2, color='gray')
            bestThr=thresholds[idx]

        print '%s %3.3f %3.3f %3.3f %3.3f'%(m,tpr[idx],fpr[idx],thresholds[idx],aucVal)
 

    plt.xlabel('Signal efficiency')
    plt.ylabel('Background efficiency')
    plt.legend(loc='upper left')
    plt.grid()
    fig.text(0.12, 0.9, r'CMS', fontweight='bold')
    fig.text(0.18, 0.9, r'Preliminary')
    for ext in ['.png','.pdf']:
        plt.savefig(outfile+ext)    
    plt.clf()

    return bestThr,data

with open(sys.argv[1],'r') as cache:
    best_models=pickle.load(cache)
    scaler=pickle.load(cache)
    #pca=pickle.load(cache)

with open(sys.argv[2],'r') as cache:
    data=pickle.load(cache)
    features=pickle.load(cache)
    spectators=pickle.load(cache)

data['X']=scaler.transform(data['X'])
#data['X']=pca.transform(data['X'])
xangle_list=data['s'][:,0]
for xangle in [120,130,140,150]:

    clfList=[]
    rfcList=[]
    for key in best_models:
        if not xangle in best_models[key] : continue
        clfList.append( (key,best_models[key][xangle][0],best_models[key][xangle][-1]) )
        if not 'rfc' in key: continue
        rfcList.append( (key,best_models[key][xangle][2]) )

    #feature importance (for random forest classifier only)
    plotFeatureImportance(rfcList,'featimportance_%d'%xangle)

    #oob = oob_classifier_accuracy(rfcList[0][1], data['X'], data['y'])
    #imp = permutation_importances(rfcList[0][1], data['X'], data['y'], oob_classifier_accuracy)

    
    #filter the data for this crossing angle
    filt=(xangle_list==xangle)
    X=data['X'][filt]
    s=data['s'][filt][:,-1]
    y=data['y'][filt]
    df=pd.DataFrame(X)
    df.columns=features
    df['nrawmu']=s.tolist()
    df['class']=y.tolist()
    
    #showFeatureCorrelation(data=df,features=features,truth='class',outfile='featurecorr_%d'%xangle)
    #showConfusionMatrix(clfList=clfList,data=df,truth='class',outfile='confmatrix_%d'%xangle)
    bestThr,df=showROCs(clfList=clfList,data=df,truth='class',outfile='rocs_%d'%xangle)
    print bestThr
