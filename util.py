

def computecurvature(U1, U2, V1, V2):
	import numpy as np
	n = max(len(U1), len(U2), len(V1), len(V2))

	u1 = np.zeros(n)
	u2 = np.zeros(n)
	v1 = np.zeros(n)
	v2 = np.zeros(n)
	for i in range(n):
		u1[i] = U1[i]
		u2[i] = U2[i]
		v1[i] = V1[i]
		v2[i] = V2[i]

	u1x = np.zeros(n)
	u2x = np.zeros(n)
	v1x = np.zeros(n)
	v2x = np.zeros(n)
	u1xx = np.zeros(n)
	u2xx = np.zeros(n)
	v1xx = np.zeros(n)
	v2xx = np.zeros(n)
	##1st derivatives
	for i in range(n):
		#cosine
		u1x[i] = (i)*u1[i]
		#sine
		u2x[i] = (-1)*(i)*u2[i]
		#cosine
		v1x[i] = (i)*v1[i]
		#sine
		v2x[i] = (-1)*(i)*v2[i]

	##2nd derivatives
		#sine
		u1xx[i] = (-1)*(i)**2*u1[i]
		#cosine
		u2xx[i] = (-1)*(i)**2*u2[i]
		#sine
		v1xx[i] = (-1)*(i)**2*v1[i]
		#cosine
		v2xx[i] = (-1)*(i)**2*v2[i]

	Ju1 = np.zeros(n)
	Ju2 = np.zeros(n)
	Jv1 = np.zeros(n)
	Jv2 = np.zeros(n)
	Ju1x = np.zeros(n)
	Ju2x = np.zeros(n)
	Jv1x = np.zeros(n)
	Jv2x = np.zeros(n)
	Ju1xx = np.zeros(n)
	Ju2xx = np.zeros(n)
	Jv1xx = np.zeros(n)
	Jv2xx = np.zeros(n)
	##Hilbert transforms of basic functions
	for i in range(n):
		#cosine
		Ju1[i] = (-1)*u1[i]
		#sine
		Ju2[i] = u2[i]
		#cosine
		Jv1[i] = (-1)*v1[i]
		#sine
		Jv2[i] = v2[i]

	##Hilbert transforms of 1st derivatives
		#sine
		Ju1x[i] = u1x[i]
		#cosine
		Ju2x[i] = (-1)*u2x[i]
		#sine
		Jv1x[i] = v1x[i]
		#cosine
		Jv2x[i] = (-1)*v2x[i]

	##Hilbert transforms of 2nd derivatives
		#cosine
		Ju1xx[i] = (-1)*u1xx[i]
		#sine
		Ju2xx[i] = u2xx[i]
		#cosine
		Jv1xx[i] = (-1)*v1xx[i]
		#sine
		Jv2xx[i] = v2xx[i]

	##Computing alpha term
	alphacc = np.empty([n,n])
	alphacs = np.empty([n,n])
	alphass = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			#cosinecosine
			alphacc[i][j] = -u1x[i]*v2[j]/2 + u2[i]*v1x[j]/2
			#cosinesine
			alphacs[i][j] = -u1x[i]*v1[j]/2 -v2[i]*u2x[j]/2 + v1x[i]*u1[j]/2 + u2[i]*v2x[j]/2
			#sinesine
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
	Lbetacc = np.empty([n,n])
	Lbetass = np.empty([n,n])
	Lbetacs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			#cosinecosine
			Lbetacc[i][j] = u1x[i]*Jv2x[j] + u2[i]*Jv1xx[j]/2 -v1x[i]*Ju2x[j] -v2[i]*Ju1xx[j]/2
			#sinesine
			Lbetass[i][j] = u2x[i]*Jv1x[j] + u1[i]*Jv2xx[j]/2 -v2x[i]*Ju1x[j] -v1[i]*Ju2xx[j]/2
			#cosinesine
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
	Ldeltacc = np.empty([n,n])
	Ldeltass = np.empty([n,n])
	Ldeltacs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			#cosinecosine
			Ldeltacc[i][j] = u1x[i]*Jv2x[j] + u2[i]*Jv1xx[j]/2 + v1x[i]*Ju2x[j] + v2[i]*Ju1xx[j]/2
			#sinesine
			Ldeltass[i][j] = u2x[i]*Jv1x[j] + u1[i]*Jv2xx[j]/2 + v2x[i]*Ju1x[j] + v1[i]*Ju2xx[j]/2
			#cosinesine
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
	LBucc = np.empty([n,n])
	LBuss = np.empty([n,n])
	LBucs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			#cosinecosine
			LBucc[i][j] = u1x[i]*Ju2x[j] + u2[i]*Ju1xx[j]/2
			#sinesine
			LBuss[i][j] = u2x[i]*Ju1x[j] + u1[i]*Ju2xx[j]/2
			#cosinesine
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
	LBvss = np.empty([n,n])
	LBvcs = np.empty([n,n])
	for i in range(n):
		for j in range(n):
			#cosinecosine
			LBvcc[i][j] = v1x[i]*Jv2x[j] + v2[i]*Jv1xx[j]/2
			#sinesine
			LBvss[i][j] = v2x[i]*Jv1x[j] + v1[i]*Jv2xx[j]/2
			#cosinesine
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
	Two = 0
	Three = 0
	Four = 0
	for i in range(2*n):
		One += Ldeltacosine[i]*deltacosine[i]/2 + Ldeltasine[i]*deltasine[i]/2
		if i == 0:
			Two += (2)*Lbetacosine[i]*alphacosine[i] + (2)*Lbetasine[i]*alphasine[i]
		else:
			Two += (1)*Lbetacosine[i]*alphacosine[i] + (1)*Lbetasine[i]*alphasine[i]
		Three += (-3)*Lalphacosine[i]*alphacosine[i]/2 + (-3)*Lalphasine[i]*alphasine[i]/2
		Four += (-2)*LBucosine[i]*Bvcosine[i] + (-2)*LBusine[i]*Bvsine[i]

	#Curvature without dividing by the denominator
	Curvature = 0
	Curvature = One + Two + Three + Four




	return Curvature
	if Curvature < 0:
		print U1
		print U2
		print V1
		print V2
