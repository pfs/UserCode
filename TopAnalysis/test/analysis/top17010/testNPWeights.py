import pickle
import ROOT

with open('test/analysis/top17010/TT2l_npsysts.pck','rb') as fin:
    np_weights=pickle.load(fin)


#print(np_weights['CRerd'].keys())

cat=('lpt',)
cat2=('hpt',)
#cat=('1b',)
#cat2=('2b',)

var_list=['CRerd','CRgmove','CRqcd','UEup','UEdown','hdampup','hdampdown']
print('Weights for cat=',cat)
print('m(l,b)\t {}'.format('\t\t'.join(var_list)))
print('='*110)

for i,mlb in enumerate([22.5496,52.9061,66.6138,82.2203,101.453,124.2,168.71]):
    line='{:3.1f}\t& '.format(mlb)
    for var in var_list:
        cat_wgt=np_weights[var][cat]['Up'].Eval(mlb)
        cat2_wgt=np_weights[var][cat2]['Up'].Eval(mlb)
        line += '%3.2f (%3.2f) \t &'%(cat_wgt,cat2_wgt)
    print line[:-1] + '\\\\'


ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
c=ROOT.TCanvas('c','c',600,600)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
mg1=ROOT.TMultiGraph()
mg2=ROOT.TMultiGraph()
colors=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]
markers=[20,20,20,22,23,22,23]
for i,var in enumerate(var_list):
    np_weights[var][cat]['Up'].SetTitle(var)
    np_weights[var][cat]['Up'].SetMarkerStyle( markers[i] )
    np_weights[var][cat]['Up'].SetMarkerColor( colors[i] )
    np_weights[var][cat]['Up'].SetLineColor( colors[i] )
    mg1.Add(np_weights[var][cat]['Up'],'p')

    np_weights[var][cat2]['Up'].SetTitle(var)
    np_weights[var][cat2]['Up'].SetMarkerStyle( markers[i] )
    np_weights[var][cat2]['Up'].SetMarkerColor( colors[i] )
    np_weights[var][cat2]['Up'].SetLineColor( colors[i] )
    mg2.Add(np_weights[var][cat2]['Up'],'p')

mg1.Draw('a')
mg1.GetXaxis().SetTitle('m(l,b) [GeV]')
mg1.GetYaxis().SetTitle('Weight')
mg1.GetYaxis().SetRangeUser(0.94,1.06)
leg=c.BuildLegend()
leg.SetBorderSize(0)
leg.SetHeader(' '.join(cat))
c.SetGridx()
c.SetGridy()
c.SaveAs('{}_rivet_Weights.png'.format('_'.join(cat)))

mg2.Draw('a')
mg2.GetXaxis().SetTitle('m(l,b) [GeV]')
mg2.GetYaxis().SetTitle('Weight')
mg2.GetYaxis().SetRangeUser(0.94,1.06)
leg=c.BuildLegend()
leg.SetBorderSize(0)
leg.SetHeader(' '.join(cat2))
c.SetGridx()
c.SetGridy()
c.SaveAs('{}_rivet_Weights.png'.format('_'.join(cat2)))
