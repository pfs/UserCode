from ROOT import TLorentzVector, TVector2, TVector3
from math import *
from copy import copy

# parameter used to set lower bounds while searching
mc=0.0

"""
Calculate the M_{T2} variable given two Four-Vector inputs and the MET four-vector
"""
def calcMt2(vis1in, vis2in, childin):
	if vis1in.M()>vis2in.M():
		vis1=copy(vis2in)
		vis2=copy(vis1in)
	else:
		vis1=copy(vis1in)
		vis2=copy(vis2in)

	child=copy(childin)

	mag=mt2Sqrt(vis1.Px()*vis1.Px()+vis1.Py()*vis1.Py())
	cospart=vis1.Px()/mag
	sinpart=vis1.Py()/mag

	ma,pxa,pya=vis1.M(),vis1.Px(),vis1.Py()
	mb,pxb,pyb=vis2.M(),vis2.Px(),vis2.Py()
	mc,pxc,pyc=child.M(),child.Px(),child.Py()

	vis1.SetXYZM(mag,0.0,0.0,ma)
	vis2.SetXYZM(pxb*cospart+pyb*sinpart, pyb*cospart-pxb*sinpart, 0.0, mb)
	child.SetXYZM(pxc*cospart+pyc*sinpart, pyc*cospart-pxc*sinpart, 0.0, mc)

	outputmt2=0.0
	solved=False

	vis1M2=vis1.M()*vis1.M()
	vis2M2=vis2.M()*vis2.M()
	mc2=mc*mc
	vis1Px2=vis1.Px()*vis1.Px()
	vis2Px2=vis2.Px()*vis2.Px()
	childPx2=child.Px()*child.Px()
	vis2Py2=vis2.Py()*vis2.Py()
	childPy2=child.Py()*child.Py()
	vis1Pt2=vis1Px2
	vis2Pt2=vis2Px2+vis2Py2
	childPt2=childPx2+childPy2
	vis1Et2=vis1M2+vis1Pt2
	vis2Et2=vis2M2+vis2Pt2
	childEt2=mc2+childPt2
	vis1Et=mt2Sqrt(vis1Et2)
	vis2Et=mt2Sqrt(vis2Et2)

	Mtmin=0.0
	Mtmax=0.0

	if (not(vis1.M()<=0.0 or vis2.M()<=0.0)):
		xlmin=vis1.Px()*mc/vis1.M()
		xrmin=vis2.Px()*mc/vis2.M()
		yrmin=vis2.Py()*mc/vis2.M()

		altxlmin=child.Px()-xlmin
		altxrmin=child.Px()-xrmin
		altyrmin=child.Py()-yrmin

		Mtlmin=vis1.M()+mc
		Mtrmin=vis2.M()+mc

		Mtratlmin=mt2Sqrt(vis2M2+mc2+2.0*(vis2Et*mt2Sqrt(mc2+altxlmin*altxlmin+childPy2)-vis2.Px()*altxlmin-vis2.Py()*child.Py()))
		Mtlatrmin=mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*mt2Sqrt(mc2+altxrmin*altxrmin+altyrmin*altyrmin)-vis1.Px()*altxrmin))

		if (Mtlmin>=Mtratlmin):
			solved=True
			outputmt2=Mtlmin
		elif (Mtrmin>=Mtlatrmin):
			solved=True
			outputmt2=Mtrmin
		else:
			if (Mtlmin>Mtrmin):
				Mtmin=Mtlmin
			else:
				Mtmin=Mtrmin

			if (Mtlatrmin<Mtratlmin):
				Mtmax=Mtlatrmin
			else:
				Mtmax=Mtratlmin

		backupmid=mt2Sqrt(mc2+(child.Px()*child.Px()+child.Py()*child.Py())*0.25)
		backup1=mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*backupmid-0.5*vis1.Px()*child.Px()))
		backup2=mt2Sqrt(vis2M2+mc2+2.0*(vis2Et*backupmid-0.5*(vis2.Px()*child.Px()+vis2.Py()*child.Py())))
		if (backup1>backup2):
			if (backup1<Mtmax):
				Mtmax=backup1
		else:
			if (backup2<Mtmax):
				Mtmax=backup2

	elif (not(vis2.M()<=0.0)):
		xrmin=vis2.Px()*mc/vis2.M()
		yrmin=vis2.Py()*mc/vis2.M()

		altxrmin=child.Px()-xrmin
		altyrmin=child.Py()-yrmin

		Mtlmin=mc
		Mtrmin=vis2.M()+mc

		Mtlatrmin=mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*mt2Sqrt(mc2+altxrmin*altxrmin+altyrmin*altyrmin)-vis1.Px()*altxrmin))

		if (Mtrmin>=Mtlatrmin):
			solved=True
			outputmt2=Mtrmin
		else:
			if (Mtlmin>Mtrmin):
				Mtmin=Mtlmin
			else:
				Mtmin=Mtrmin

			Mtmax=Mtlatrmin

		backupmid=mt2Sqrt(mc2+(child.Px()*child.Px()+child.Py()*child.Py())*0.25)
		backup1=mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*backupmid-0.5*vis1.Px()*child.Px()))
		backup2=mt2Sqrt(vis2M2+mc2+2.0*(vis2Et*backupmid-0.5*(vis2.Px()*child.Px()+vis2.Py()*child.Py())))
		if (backup1>backup2):
			if (backup1<Mtmax):
				Mtmax=backup1
		else:
			if (backup2<Mtmax):
				Mtmax=backup2

	else:
		Mtmin=mc
		trialmid=mt2Sqrt(mc2+0.25*childPt2)
		trial1=mt2Sqrt(mc2+2.0*(vis1Et*trialmid-vis1.Px()*child.Px()*0.5))
		trial2=mt2Sqrt(mc2+2.0*(vis2Et*trialmid-0.5*(vis2.Px()*child.Px()+vis2.Py()*child.Py())))
		if (trial1>trial2):
			Mtmax=trial1
		else:
			Mtmax=trial2


	if (not solved):
		solved=True
		outputmt2=Mtmin+(Mtmax-Mtmin)*0.5

		C1=1.0/vis1Et2

		A2=vis2M2+vis2Py2
		B2=-2.0*vis2.Px()*vis2.Py()
		C2=1.0/(vis2M2+vis2Px2)

		preF1=vis1Et2*mc2

		preD2=-2.0*child.Px()*A2-B2*child.Py()
		preE2=-2.0*child.Py()/C2-B2*child.Px()
		preF2=vis2Et2*childEt2-childPx2*vis2Px2-childPy2*vis2Py2+B2*child.Px()*child.Py()

		G=B2*0.5*C2
		J1=-vis1M2*C1
		J2=(B2*B2*0.25*C2-A2)*C2

		alpha=G*G-J1-J2
		p0_4=alpha*alpha-4.0*J1*J2

		while (outputmt2>Mtmin and outputmt2<Mtmax):
			q1=mc2+vis1M2-outputmt2*outputmt2
			D1=q1*vis1.Px()
			F1=preF1-q1*q1*0.25

			q2=outputmt2*outputmt2-mc2-vis2M2
			D2=preD2+q2*vis2.Px()
			E2=preE2+q2*vis2.Py()
			F2=preF2-q2*(q2*0.25+vis2.Px()*child.Px()+vis2.Py()*child.Py())

			H=E2*0.5*C2

			K1=-D1*C1
			L1=-F1*C1

			K2=(B2*E2*0.5*C2-D2)*C2
			L2=(E2*E2*0.25*C2-F2)*C2

			beta=2.0*G*H-K1-K2
			gamma=H*H-L1-L2

			p0_4nonzero=p0_4
			if fabs(p0_4)<1e-9:
				p0_4nonzero=1e-9

			p0_3=(2.0*alpha*beta-4.0*(J1*K2+J2*K1))/p0_4nonzero
			p0_2=(2.0*alpha*gamma+beta*beta-4.0*(J1*L2+J2*L1+K1*K2))/p0_4nonzero
			p0_1=(2.0*beta*gamma-4.0*(K1*L2+K2*L1))/p0_4nonzero
			p0_0=(gamma*gamma-4.0*L1*L2)/p0_4nonzero

			p2_2=0.1875*p0_3*p0_3-p0_2*0.5
			p2_1=p0_3*p0_2*0.125-0.75*p0_1
			p2_0=p0_3*p0_1*0.0625-p0_0

			p2_2nonzero=p2_2
			if fabs(p2_2)<1e-9:
				p2_2nonzero=1e-9

			p3_1=(4.0*p2_0+3.0*p0_3*p2_1)/p2_2nonzero-4.0*p2_1*p2_1/(p2_2nonzero*p2_2nonzero)-2.0*p0_2
			p3_0=3.0*p0_3*p2_0/p2_2nonzero-4.0*p2_1*p2_0/(p2_2nonzero*p2_2nonzero)-p0_1

			p4_0=0
			if fabs(p3_1)<1e-9:
				if p2_2>0:
					p4_0=-1
				elif p2_2==0:
					p4_0=0
				else:
					p4_0=1
			else:
				p4_0=p2_1*p3_0/p3_1-p2_2*p3_0*p3_0/(p3_1*p3_1)-p2_0

			negroots=1
			posroots=0

			if ((p0_4<0.0 and p2_2<0.0) or (p0_4>0.0 and p2_2>0.0)):
				negroots+=1
			elif((p0_4<0.0 and p2_2>0.0) or (p0_4>0.0 and p2_2<0.0)):
				posroots+=1

			if ((p2_2<0.0 and p3_1<0.0) or (p2_2>0.0 and p3_1>0.0)):
				negroots+=1
			elif((p2_2<0.0 and p3_1>0.0) or (p2_2>0.0 and p3_1<0.0)):
				posroots+=1

			if ((p3_1<0.0 and p4_0<0.0) or (p3_1>0.0 and p4_0>0.0)):
				negroots+=1
			elif((p3_1<0.0 and p4_0>0.0) or (p3_1>0.0 and p4_0<0.0)):
				posroots+=1

			if (posroots==negroots):
				Mtmin=outputmt2
			else:
				Mtmax=outputmt2

			outputmt2=Mtmin+(Mtmax-Mtmin)*0.5
	return outputmt2

"""
Special version of sqrt() to handle extreme inputs
"""
def mt2Sqrt(x):
	if isnan(x):
		return 0.0
	elif (x<=0.0):
		return 0.0
	elif (isinf(x)):
		return 1e99999999999999999999999999999999
	else:
		prev2root=-1.0
		prevroot=-1.0
		root=1.0
		while((root!=prevroot) and (root!=prev2root)):
			prev2root=prevroot
			prevroot=root
			root=(root+x/root)*0.5

		return root
