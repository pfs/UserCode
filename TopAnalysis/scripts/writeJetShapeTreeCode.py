#! /usr/bin/env python
import sys

def main():
  jets = ["j", "gj"]
  observables = ["mult", "width", "ptd", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
  reco        = ["charged", "puppi", "all"]
  
  ##cc
  #create branches
  for j in jets:
    for o in observables:
      for r in reco:
        print('t->Branch("'+j+'_'+o+'_'+r+'",   tjsev.'+j+'_'+o+'_'+r+',   "'+j+'_'+o+'_'+r+'[n'+j+']/F");')
    print('\n')
  
  #reset values
  for j in jets:
    for o in observables:
      for r in reco:
        print('tjsev.'+j+'_'+o+'_'+r+'[i]=-99;')
    print('\n')
  
  #header
  for j in jets:
    for o in observables:
      for r in reco:
        print('Float_t '+j+'_'+o+'_'+r+'[50];')
    print('\n')
  
  #C observable filling
  ns    = ["1", "2", "3"]
  betas = ['02', '05', '10', '20']
  for n in ns:
    for b in betas:
      print('tjsev.gj_c'+n+'_'+b+'_charged[i] = getC('+n+', '+str(float(b)/10.)+', genJets[i]);')
      print('tjsev.gj_c'+n+'_'+b+'_all[i]     = getC('+n+', '+str(float(b)/10.)+', genJets[i], true);')
      print('tjsev.gj_c'+n+'_'+b+'_puppi[i]   = getC('+n+', '+str(float(b)/10.)+', genJets[i], true, true);')
  print('\n')
  for n in ns:
    for b in betas:
      print('tjsev.j_c'+n+'_'+b+'_charged[ij] = getC('+n+', '+str(float(b)/10.)+', jets[ij]);')
      print('tjsev.j_c'+n+'_'+b+'_all[ij]     = getC('+n+', '+str(float(b)/10.)+', jets[ij], true);')
      print('tjsev.j_c'+n+'_'+b+'_puppi[ij]   = getC('+n+', '+str(float(b)/10.)+', jets[ij], true, true);')
  print('\n')
  
  for n in ns:
    for b in betas:
      print('if (obs == "c'+n+'_'+b+'"):')
      print('    label = "C_{'+n+'}^{('+str(float(b)/10.)+')}"')
      print('    highbin = 0.7')
      print('    nbins = 35')

if __name__ == "__main__":
	sys.exit(main())
