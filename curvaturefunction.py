from util import computecurvature

def runonetime(usine, ucos, vsine, vcos):
	return computecurvature(usine, ucos, vsine, vcos)

def runrand():
	import numpy as np
	alwayspos = True
	#Arrays with random entries
	while alwayspos:
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
		print computecurvature(Usine, Ucos, Vsine, Vcos)
		if computecurvature(Usine, Ucos, Vsine, Vcos) < 0:
			alwayspos = False

def runrandgauss():
	import numpy as np
	import math
	alwayspos = True
	#Arrays with random entries with Gaussian decay
	while alwayspos:
		#A are the sine terms and B are the cosine terms
		#The first term must be zero since we are dealing with mean zero vfs
		#This is the length of the arrays we are scanning
		n = 3
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
		print computecurvature(Usine, Ucos, Vsine, Vcos)
		if computecurvature(Usine, Ucos, Vsine, Vcos) < 0:
			alwayspos = False


if __name__ == '__main__':
	runrandgauss()




























