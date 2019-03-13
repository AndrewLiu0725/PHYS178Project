'''
	Is it necessary that neurons subserving integration during simulation also
	show persistent activity during working memory, Article asks this question explicitly 
	on page 2

	It is necessary that the neurons continue to fire in the attractor sites
	during working memory

'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve 
from twovar import *
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
cprime = 0
#=================================
# 			Phase Plane
#=================================
S1, S2 = experiment()
plt.plot(S2,S1)
plt.xlabel('S2')
plt.ylabel('S1')
plt.title('Phase Plane')
plt.show()


# Plotting Nuclines for phase plane analysis
# dS1/dt == 0
# dS2/dt == 0
# Since S1 and S2 cant be solved directly need numerical solutions

# Not an efficient method will try solving simultaneous equations
# and using sympy to solve

def S_one(S1,S2):
	# I_noise_1 How to implement in this code? Suggested Noiseless version
	I1 = Jext*mu0*(1+cprime/100)
	x1 = J11*S1 - J12*S2 + I0 + I1 
	H1 = (a*x1 - b)/(1-np.exp(-d*(a*x1 - b)))
	return -(S1/taus) + (1-S1)*gamma*H1

def S_two(S2,S1):
	I2 = Jext*mu0*(1-cprime/100)
	x2 = J22 * S2 - J21*S1 + I0 + I2
	H2 = (a*x2 - b)/(1-np.exp(-d*(a*x2 - b)))
	return -(S2/taus) + (1-S2)*gamma*H2 	

S1_S2 = []
itera = np.arange(0,1,0.01)
for i in itera:
	S1_S2+=[np.roots(S_one,1/i,args=(i))[0]]	

S2_S1 = []
for i in itera:
	S2_S1 += [fsolve(S_two,1/i,args=(i))[0]]

plt.plot(itera,S1_S2,itera,S2_S1)		 





	