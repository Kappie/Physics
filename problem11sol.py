from __future__ import division
import numpy as np
import scipy.linalg as spl
import numpy.linalg as npl
import matplotlib.pyplot as plt

# Parameters
J = -1
d = 2
betaList = np.array([0.1, 0.2, 0.3, 0.4, 0.43, 0.44, 0.45, \
	 0.46, 0.47, 0.48, 0.49, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
betaZoom = np.linspace(0.42, 0.45, 20)
ChiList = [2,3,5,10,15]
MaxIter = 10000
accuracy = 1e-8
np.random.seed(1)


def ExactMag(x):
	betaCrit = 0.5*np.log(1 + np.sqrt(2))
	if x > betaCrit:
		return ( 1 - np.sinh(2*x)**(-4) )**(1/8)
	else:
		return 0


def ConstructTensors(J, beta):
	Hbond = J*np.array([[1,-1],[-1,1]])
	d = Hbond.shape[0]
	Qs = spl.sqrtm( np.exp(-beta*Hbond) )

	# Lattice tensors
	a = np.zeros([d,d,d,d])
	b = np.zeros([d,d,d,d])
	for i in xrange(d):
		for j in xrange(d):
			for k in xrange(d):
				for l in xrange(d):
					for s in xrange(d):
						a[i,j,k,l] = a[i,j,k,l] +  Qs[i,s]*Qs[j,s]*Qs[k,s]*Qs[l,s]
						b[i,j,k,l] = b[i,j,k,l] +  (-2*s+1)*Qs[i,s]*Qs[j,s]*Qs[k,s]*Qs[l,s]

	return a, b


def SymmetrizeTC(T, C):
	C = C + np.transpose(C)
	for i in xrange(d):
		T[i] = T[i] + np.transpose(T[i])
	return T,C


def GrowC(a, T, C, d, Chi):
	Tens = np.tensordot(T, a, axes=[0,1])
	Tens = np.tensordot(Tens, C, axes=[1,0])
	Tens = np.tensordot(Tens, T, axes=([4,3],[1,0]))
	Tens = np.transpose(Tens, (0,1,3,2))
	return np.reshape(Tens, [Chi*d, Chi*d])


def GrowT(a, T, d, Chi):
	Tens = np.tensordot(T, a, axes=[0,0])
	Tens = np.transpose(Tens, (3,1,2,0,4))
	return np.reshape(Tens, [d, Chi*d, Chi*d])


def LatticeGrowthStep(a, T, C, d, Chi):
	newC = GrowC(a,T,C,d,Chi)
	U, s, V = npl.svd(newC, full_matrices=0)
	V = V[:Chi,:]
	newC = np.dot( np.dot(V, newC), np.transpose(V))
	newT = GrowT(a,T,d,Chi)
	newT = np.tensordot(newT, V, axes=[1,1])
	newT = np.tensordot(newT, np.transpose(V), axes=[1,0])
	svals = s[:Chi]
	return newT/np.amax(newT), newC/np.amax(newC), svals/np.sum(svals)


def doCTM(T, C, Chi, J, beta, MaxIter, accuracy):
	a, b = ConstructTensors(J, beta)
	svals = np.ones(Chi) / Chi
	for it in xrange(MaxIter):
		svals_old = svals
		T, C = SymmetrizeTC(T,C)
		T, C, svals = LatticeGrowthStep(a,T,C,d,Chi)
		diff = np.sum( np.abs(svals - svals_old) )
		if diff < accuracy:
			break
	return a, b, T, C, it


def ComputeMagnetisation(a, b, T ,C):
	env = np.tensordot(C, T, axes =[1,1])
	env = np.tensordot(env, env, axes=[2,0])
	env = np.tensordot(env, env, axes=([3,0],[0,3]))
	return np.tensordot(env, b, axes=([0,1,2,3],[0,1,2,3])) \
		/ np.tensordot(env, a, axes=([0,1,2,3],[0,1,2,3]))





# ===========================================
# Do CTM for different values of beta and Chi
def PlotMag(betaVals, ChiVals, plotName):
	mag = np.zeros([len(ChiList), len(betaVals)])

	for i, Chi in enumerate(ChiList):
		C = np.random.rand(Chi,Chi)
		T = np.random.rand(d,Chi,Chi)
		print '\nSet Chi = %2d \n' %Chi
		for j, beta in enumerate(reversed(betaVals)): # reversed because otherwise stuck in mag=0 phase
			a, b, T, C, it = doCTM(T, C, Chi, J, beta, MaxIter, accuracy)
			mag[i,j] = np.abs( ComputeMagnetisation(a, b, T, C) )
			diff_with_mag_exact = mag[i,j] - ExactMag(beta)
			print 'Chi: %2d, beta: %.3f, #iterations: %4d, magnetisation: %.8f, diff: %0.8f' \
				 %(Chi, beta, it, mag[i,j], diff_with_mag_exact)

	mag_exact = [ExactMag(x) for x in np.linspace(betaVals[0],betaVals[-1],1000)]

	plt.figure()
	for i, Chi in enumerate(ChiList):
		plt.plot(betaVals[::-1], np.abs(mag[i]), 'o:', ms=4, label='Chi =%d' %Chi)
	plt.plot(np.linspace(betaVals[0],betaVals[-1],1000), mag_exact, '-', label='Exact result')
	plt.xlabel('beta = 1/T')
	plt.ylabel('Magnetisation')
	plt.title('Magnetisation(beta,Chi)')
	plt.legend(loc='lower right')
	plt.savefig(plotName)


PlotMag(betaList, ChiList, 'Magnetisation_Chi')
PlotMag(betaZoom, ChiList, 'Magnetisation_Chi_PT')
