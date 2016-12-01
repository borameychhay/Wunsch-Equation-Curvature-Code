import pandas as pd
import numpy as np
import sys
import math


def computecurvature(U1, U2, V1, V2):
	nu1 = len(U1)
	nu2 = len(U2)
	nv1 = len(V1)
	nv2 = len(V2)
	n = max(nu1, nu2, nv1, nv2)

	for i in range(n-nu1):
		U1.append(0)
	for i in range(n-nu2):
		U2.append(0)
	for i in range(n-nv1):
		V1.append(0)
	for i in range(n-nv2):
		V2.append(0)

	u1 = np.zeros(n)
	for i in range(n):
		u1[i] = U1[i]
	u2 = np.zeros(n)
	for i in range(n):
		u2[i] = U2[i]
	v1 = np.zeros(n)
	for i in range(n):
		v1[i] = V1[i]
	v2 = np.zeros(n)
	for i in range(n):
		v2[i] = V2[i]




	##derivatives
	#cosine
	u1x = np.zeros(n)
	for i in range(n):
		u1x[i] = (i)*u1[i]
	#sine
	u2x = np.zeros(n)
	for i in range(n):
		u2x[i] = (-1)*(i)*u2[i]

	#cosine
	v1x = np.zeros(n)
	for i in range(n):
		v1x[i] = (i)*v1[i]
	#sine
	v2x = np.zeros(n)
	for i in range(n):
		v2x[i] = (-1)*(i)*v2[i]

	##2nd derivatives
	#sine
	u1xx = np.zeros(n)
	for i in range(n):
		u1xx[i] = (-1)*(i)**2*u1[i]
	#cosine
	u2xx = np.zeros(n)
	for i in range(n):
		u2xx[i] = (-1)*(i)**2*u2[i]

	#sine
	v1xx = np.zeros(n)
	for i in range(n):
		v1xx[i] = (-1)*(i)**2*v1[i]
	#cosine
	v2xx = np.zeros(n)
	for i in range(n):
		v2xx[i] = (-1)*(i)**2*v2[i]


	##Hilbert transforms
	#cosine
	Ju1 = np.zeros(n)
	for i in range(n):
		Ju1[i] = (-1)*u1[i]
	#sine
	Ju2 = np.zeros(n)
	for i in range(n):
		Ju2[i] = u2[i]

	#cosine
	Jv1 = np.zeros(n)
	for i in range(n):
		Jv1[i] = (-1)*v1[i]
	#sine
	Jv2 = np.zeros(n)
	for i in range(n):
		Jv2[i] = v2[i]

	#sine
	Ju1x = np.zeros(n)
	for i in range(n):
		Ju1x[i] = u1x[i]
	#cosine
	Ju2x = np.zeros(n)
	for i in range(n):
		Ju2x[i] = (-1)*u2x[i]

	#sine
	Jv1x = np.zeros(n)
	for i in range(n):
		Jv1x[i] = v1x[i]
	#cosine
	Jv2x = np.zeros(n)
	for i in range(n):
		Jv2x[i] = (-1)*v2x[i]

	#cosine
	Ju1xx = np.zeros(n)
	for i in range(n):
		Ju1xx[i] = (-1)*u1xx[i]
	#sine
	Ju2xx = np.zeros(n)
	for i in range(n):
		Ju2xx[i] = u2xx[i]

	#cosine
	Jv1xx = np.zeros(n)
	for i in range(n):
		Jv1xx[i] = (-1)*v1xx[i]
	#sine
	Jv2xx = np.zeros(n)
	for i in range(n):
		Jv2xx[i] = v2xx[i]

	##Computing alpha term
	#cosinecosine
	alphacc = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			alphacc[i][j] = -u1x[i]*v2[j]/2 + u2[i]*v1x[j]/2
	#cosinesine
	alphacs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			alphacs[i][j] = -u1x[i]*v1[j]/2 -v2[i]*u2x[j]/2 + v1x[i]*u1[j]/2 + u2[i]*v2x[j]/2
	#sinesine
	alphass = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			alphass[i][j] = -u2x[i]*v1[j]/2 + u1[i]*v2x[j]/2
	##Putting alpha in an array
	##Cosine term and Sine term
	alphacosine = np.zeros(2*n)
	alphasine = np.zeros(2*n)
	for i in range(n):
		for j in range(n):
			alphacosine[abs(i-j)] += alphacc[i][j]/2 + alphass[i][j]/2
			alphacosine[i+j] += alphacc[i][j]/2 - alphass[i][j]/2
			if i>j:
				alphasine[i+j] += alphacs[i][j]/2
				alphasine[i-j] += -1*alphacs[i][j]/2
			if i<j:
				alphasine[i+j] += alphacs[i][j]/2
				alphasine[j-i] += alphacs[i][j]/2
			if i==j:
				alphasine[i+j] += alphacs[i][j]/2
	##Computing Lambda alpha
	Lalphacosine = np.zeros(2*n)
	Lalphasine = np.zeros(2*n)
	for i in range(2*n):
		Lalphacosine[i] = i*alphacosine[i]
		Lalphasine[i] = i*alphasine[i]


	##Computing Lambda beta term
	#cosinecosine
	Lbetacc = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			Lbetacc[i][j] = u1x[i]*Jv2x[j] + u2[i]*Jv1xx[j]/2 -v1x[i]*Ju2x[j] -v2[i]*Ju1xx[j]/2
	#sinesine
	Lbetass = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			Lbetass[i][j] = u2x[i]*Jv1x[j] + u1[i]*Jv2xx[j]/2 -v2x[i]*Ju1x[j] -v1[i]*Ju2xx[j]/2
	#cosinesine
	Lbetacs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			Lbetacs[i][j] = u1x[i]*Jv1x[j] + Jv2x[i]*u2x[j] + Jv1xx[i]*u1[j]/2 + u2[i]*Jv2xx[j]/2 -v1x[i]*Ju1x[j] -Ju2x[i]*v2x[j] -Ju1xx[i]*v1[j]/2 -v2[i]*Ju2xx[j]/2
	##Putting Lambda beta in an array
	##Cosine term and Sine term
	Lbetacosine = np.zeros(2*n)
	Lbetasine = np.zeros(2*n)
	for i in range(n):
		for j in range(n):
			Lbetacosine[abs(i-j)] += Lbetacc[i][j]/2 + Lbetass[i][j]/2
			Lbetacosine[i+j] += Lbetacc[i][j]/2 - Lbetass[i][j]/2
			if i>j:
				Lbetasine[i+j] += Lbetacs[i][j]/2
				Lbetasine[i-j] += -1*Lbetacs[i][j]/2
			if i<j:
				Lbetasine[i+j] += Lbetacs[i][j]/2
				Lbetasine[j-i] += Lbetacs[i][j]/2
			if i==j:
				Lbetasine[i+j] += Lbetacs[i][j]/2


	##Computing Lambda delta term
	#cosinecosine
	Ldeltacc = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			Ldeltacc[i][j] = u1x[i]*Jv2x[j] + u2[i]*Jv1xx[j]/2 + v1x[i]*Ju2x[j] + v2[i]*Ju1xx[j]/2
	#sinesine
	Ldeltass = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			Ldeltass[i][j] = u2x[i]*Jv1x[j] + u1[i]*Jv2xx[j]/2 + v2x[i]*Ju1x[j] + v1[i]*Ju2xx[j]/2
	#cosinesine
	Ldeltacs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			Ldeltacs[i][j] = u1x[i]*Jv1x[j] + Jv2x[i]*u2x[j] + Jv1xx[i]*u1[j]/2 + u2[i]*Jv2xx[j]/2 + v1x[i]*Ju1x[j] + Ju2x[i]*v2x[j] + Ju1xx[i]*v1[j]/2 + v2[i]*Ju2xx[j]/2
	##Putting Lambda delta in an array
	##Cosine term and Sine term
	Ldeltacosine = np.zeros(2*n)
	Ldeltasine = np.zeros(2*n)
	for i in range(n):
		for j in range(n):
			Ldeltacosine[abs(i-j)] += Ldeltacc[i][j]/2 + Ldeltass[i][j]/2
			Ldeltacosine[i+j] += Ldeltacc[i][j]/2 - Ldeltass[i][j]/2
			if i>j:
				Ldeltasine[i+j] += Ldeltacs[i][j]/2
				Ldeltasine[i-j] += -1*Ldeltacs[i][j]/2
			if i<j:
				Ldeltasine[i+j] += Ldeltacs[i][j]/2
				Ldeltasine[j-i] += Ldeltacs[i][j]/2
			if i==j:
				Ldeltasine[i+j] += Ldeltacs[i][j]/2
	##Computing the delta term
	deltacosine = np.zeros(2*n)
	deltasine = np.zeros(2*n)
	for i in range(2*n):
		if i != 0:
			deltacosine[i] = Ldeltacosine[i]/i
			deltasine[i] = Ldeltasine[i]/i


	#Computing the LBu term
	#cosinecosine
	LBucc = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			LBucc[i][j] = u1x[i]*Ju2x[j] + u2[i]*Ju1xx[j]/2
	#sinesine
	LBuss = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			LBuss[i][j] = u2x[i]*Ju1x[j] + u1[i]*Ju2xx[j]/2
	#cosinesine
	LBucs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			LBucs[i][j] = u1x[i]*Ju1x[j] + Ju2x[i]*u2x[j] + Ju1xx[i]*u1[j]/2 + u2[i]*Ju2xx[j]/2
	##Putting Lambda Bu in an array
	##Cosine term and Sine term
	LBucosine = np.zeros(2*n)
	LBusine = np.zeros(2*n)
	for i in range(n):
		for j in range(n):
			LBucosine[abs(i-j)] += LBucc[i][j]/2 + LBuss[i][j]/2
			LBucosine[i+j] += LBucc[i][j]/2 - LBuss[i][j]/2
			if i>j:
				LBusine[i+j] += LBucs[i][j]/2
				LBusine[i-j] += (-1)*LBucs[i][j]/2
			if i<j:
				LBusine[i+j] += LBucs[i][j]/2
				LBusine[j-i] += LBucs[i][j]/2
			if i==j:
				LBusine[i+j] += LBucs[i][j]/2


	#Computing the LBv term
	#cosinecosine
	LBvcc = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			LBvcc[i][j] = v1x[i]*Jv2x[j] + v2[i]*Jv1xx[j]/2
	#sinesine
	LBvss = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			LBvss[i][j] = v2x[i]*Jv1x[j] + v1[i]*Jv2xx[j]/2
	#cosinesine
	LBvcs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			LBvcs[i][j] = v1x[i]*Jv1x[j] + Jv2x[i]*v2x[j] + Jv1xx[i]*v1[j]/2 + v2[i]*Jv2xx[j]/2
	##Putting Lambda Bv in an array
	##Cosine term and Sine term
	LBvcosine = np.zeros(2*n)
	LBvsine = np.zeros(2*n)
	for i in range(n):
		for j in range(n):
			LBvcosine[abs(i-j)] += LBvcc[i][j]/2 + LBvss[i][j]/2
			LBvcosine[i+j] += LBvcc[i][j]/2 - LBvss[i][j]/2
			if i>j:
				LBvsine[i+j] += LBvcs[i][j]/2
				LBvsine[i-j] += (-1)*LBvcs[i][j]/2
			if i<j:
				LBvsine[i+j] += LBvcs[i][j]/2
				LBvsine[j-i] += LBvcs[i][j]/2
			if i==j:
				LBvsine[i+j] += LBvcs[i][j]/2
	##Computing the Bv term
	Bvcosine = np.zeros(2*n)
	Bvsine = np.zeros(2*n)
	for i in range(2*n):
		if i != 0:
			Bvcosine[i] = LBvcosine[i]/i
			Bvsine[i] = LBvsine[i]/i


	##Now adding up the terms
	One = 0
	for i in range(2*n):
		One += Ldeltacosine[i]*deltacosine[i]/2 + Ldeltasine[i]*deltasine[i]/2

	Two = 0
	for i in range(2*n):
		if i == 0:
			Two += (2)*Lbetacosine[i]*alphacosine[i] + (2)*Lbetasine[i]*alphasine[i]
		else:
			Two += (1)*Lbetacosine[i]*alphacosine[i] + (1)*Lbetasine[i]*alphasine[i]

	Three = 0
	for i in range(2*n):
		Three += (-3)*Lalphacosine[i]*alphacosine[i]/2 + (-3)*Lalphasine[i]*alphasine[i]/2

	Four = 0
	for i in range(2*n):
		Four += (-2)*LBucosine[i]*Bvcosine[i] + (-2)*LBusine[i]*Bvsine[i]

	#Curvature without dividing by the denominator
	Curvature = 0
	Curvature = One + Two + Three + Four




	print Curvature
	if Curvature < 0:
		print One
		print Two
		print Three
		print Four


		print alphacosine
		print alphasine
		print Lalphacosine
		print Lalphasine




		sys.exit()





run = 1
if run == 0:
	usine = [0,16,2]
	ucos = [0,0,10]
	vsine = [0,3,7]
	vcos = [0,12,1]

	computecurvature(usine, ucos, vsine, vcos)



#Arrays with random entries
while run == 1:
	#A are the sine terms and B are the cosine terms
	#The first term must be zero since we are dealing with mean zero vfs
	#This is the length of the arrays we are scanning
	n = 3
	Usine = np.random.uniform(low =-1, high =1, size =n)
	Usine[0] = 0
	Ucos = np.random.uniform(low =-1, high =1, size =n)
	Ucos[0] = 0
	#The first term must be zero since we are dealing with mean zero vfs
	#C are the sine terms and D are the cosine terms								
	Vsine = np.random.uniform(low =-1, high =1, size =n)
	Vsine[0] = 0
	Vcos = np.random.uniform(low =-1, high =1, size =n)
	Vcos[0] = 0
	computecurvature(Usine, Ucos, Vsine, Vcos)
	


#Arrays with random entries with Gaussian decay
while run == 2:
	#A are the sine terms and B are the cosine terms
	#The first term must be zero since we are dealing with mean zero vfs
	#This is the length of the arrays we are scanning
	n = 30
	#We can adjust the variance of the Gaussian distribution here
	var = 1
	Usine = np.random.uniform(low =-1, high =1, size =n)
	Usine[0] = 0
	Ucos = np.random.uniform(low =-1, high =1, size =n)
	Ucos[0] = 0
	#The first term must be zero since we are dealing with mean zero vfs
	#C are the sine terms and D are the cosine terms								
	Vsine = np.random.uniform(low =-1, high =1, size =n)
	Vsine[0] = 0
	Vcos = np.random.uniform(low =-1, high =1, size =n)
	Vcos[0] = 0
	for i in range(n):
		Usine[i] = Usine[i]*math.exp(-i*i/(2*var))
		Ucos[i] = Ucos[i]*math.exp(-i*i/(2*var))
		Vsine[i] = Vsine[i]*math.exp(-i*i/(2*var))
		Vcos[i] = Vcos[i]*math.exp(-i*i/(2*var))
	computecurvature(Usine, Ucos, Vsine, Vcos)






