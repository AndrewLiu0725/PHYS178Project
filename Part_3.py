'''
	Is it necessary that neurons subserving integration during simulation also
	show persistent activity during working memory, Article asks this question explicitly 
	on page 2

	It is necessary that the neurons continue to fire in the attractor sites
	during working memory

'''
import matplotlib.pyplot as plt
import numpy as np
from sympy.solvers import solve
from scipy.optimize import fsolve 

#=================================
#			Constants
#=================================
a = 270
b = 108
d = 0.154
gamma = 0.641
taus = 100*10**-3
tauampa = 2*10**-3
J11, J22 = 0.2609, 0.2609
J12, J21 = 0.0497, 0.0497
Jext = 5.2*(10**-4)
I0 = 0.3255
stdnoise = 0.02
mu0 = 30
dt = 0.1*10**-3

#=================================
# 			Phase Plane
#=================================
'''
S1, S2 = experiment()
plt.plot(S2,S1)
plt.xlabel('S2')
plt.ylabel('S1')
plt.title('Phase Plane')
plt.show()

'''
# Plotting Nuclines for phase plane analysis
# dS1/dt == 0
# dS2/dt == 0
# Since S1 and S2 cant be solved directly need numerical solutions

# Not an efficient method will try solving simultaneous equations
# and using sympy to solve
def nullclines(cprime):

	def gate(S):
		# Input Currents 
		I = [Jext*mu0*(1+(cprime/100)), 
			 Jext*mu0*(1-(cprime/100)),0]
		# x variable
		x = [J11*S[0] - J12*S[1] + I0 + I[0],
			 J22*S[1] - J21*S[0] + I0 + I[1]]

		# H values	 
		H = [(a*x[i] - b)/(1 - np.exp(-d*(a*x[i] - b))) for i in range(2)] 

		# dS1/dt and dS2/dt
		expr = [-(S[i]/taus) + (1-S[i])*gamma*H[i] for i in range(2)]
		return expr


	# Brute Force approach

	iterr = np.arange(0,1,0.001)
	ds1 = []
	ds2 = []

	for i in iterr:
	
		for j in iterr:
			x,y = gate([i,j])
		
			if 0.01 > x > -0.01:
				ds1 += [[i,j]]
			if 0.01 > y > -0.01:
				ds2 += [[i,j]]

		
	# Stack the results for easy plotting		
	ds1=np.vstack(ds1)
	ds2=np.vstack(ds2)

	plt.scatter(ds1[:,0],ds1[:,1],s=1)
	plt.scatter(ds2[:,0],ds2[:,1],s=1)
	plt.legend(['dS1/dt=0','dS2/dt=0'])
	plt.title('Nullclines Cprime: {}'.format(cprime))
	plt.show()


nullclines(100)



	 





	