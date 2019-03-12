import matplotlib.pyplot as plt
import numpy as np

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
cprime = 52

def H(xi):
	return (a*xi-b)/(1-np.exp(-d*(a*xi-b)))

interval = 3
steps = int(interval/dt)
time = np.linspace(0, interval, steps)

def experiment():
	H1, H2, S1, S2 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
	H1[0] = 2.5
	H2[0] = 2.5
	Inoise1, Inoise2 = np.zeros(steps+1), np.zeros(steps+1)
	for index, t in enumerate(time):
		Inoise1[index+1] = Inoise1[index] + dt*(-Inoise1[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		Inoise2[index+1] = Inoise2[index] + dt*(-Inoise2[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		x1 = J11*S1[index] - J12*S2[index] + I0 + Jext*mu0*(1+cprime/100) + Inoise1[index]
		x2 = J22*S2[index] - J21*S1[index] + I0 + Jext*mu0*(1-cprime/100) + Inoise2[index]
		H1[index+1] = H(x1)
		H2[index+1] = H(x2)
		S1[index+1] = S1[index] + dt*(-S1[index]/taus + (1 - S1[index])*gamma*H1[index])
		S2[index+1] = S2[index] + dt*(-S2[index]/taus + (1 - S2[index])*gamma*H2[index])
	return H1[1:], H2[1:]

plt.figure
for i in range(10):
	result = experiment()
	plt.plot(time, result[0], 'r', time, result[1], 'c')
plt.plot(time, 15*np.ones(steps))
plt.xlabel('Time')
plt.ylabel('Firing rate')
plt.show()