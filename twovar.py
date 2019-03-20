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
stdnoise = 0.5
mu0 = 30
dt = 0.1*10**-3
cprime = 0

def H(xi):
	return (a*xi-b)/(1-np.exp(-d*(a*xi-b)))
starttime = -0.5
endtime = 1.5
steps = int(abs(starttime - endtime)/dt)
time = np.linspace(starttime, endtime, steps)

def experiment():
	H1, H2, S1, S2 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
	H1[0] = 0
	H2[0] = 0
	Inoise1, Inoise2 = np.zeros(steps+1), np.zeros(steps+1)
	for index, t in enumerate(time):
		Inoise1[index+1] = Inoise1[index] + dt*(-Inoise1[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		Inoise2[index+1] = Inoise2[index] + dt*(-Inoise2[index] + np.random.normal(0, 1, 1)[0]*np.sqrt(tauampa*stdnoise**2))/tauampa
		if t > 0:
			x1 = J11*S1[index] - J12*S2[index] + I0 + Jext*mu0*(1+cprime/100) + Inoise1[index]
			x2 = J22*S2[index] - J21*S1[index] + I0 + Jext*mu0*(1-cprime/100) + Inoise2[index]
		else:
			x1 = J11*S1[index] - J12*S2[index] + I0 + Inoise1[index]
			x2 = J22*S2[index] - J21*S1[index] + I0 + Inoise2[index]
		H1[index+1] = H(x1)
		H2[index+1] = H(x2)
		S1[index+1] = S1[index] + dt*(-S1[index]/taus + (1 - S1[index])*gamma*H1[index])
		S2[index+1] = S2[index] + dt*(-S2[index]/taus + (1 - S2[index])*gamma*H2[index])
	return H1[1:], H2[1:]

def slided(data):
	timestep = 5*10**-3
	slided_data = []
	for index, value in enumerate(data):
		if index % int(timestep/dt) == 0:
			slided_data.append(value)
	return slided_data


def smoothing(data):
	length = len(data)
	smoothed_data = np.zeros(length)
	width = int(10*10**-3/dt)
	for i in range(length):
		if length - (i+1) < width:
			smoothed_data[i] = np.average(data[i:])
		else: smoothed_data[i] = np.average(data[i: i+width])
	return smoothed_data

'''plt.figure
cprime = 20
starttime = -0.1
endtime = 0.1
time = np.linspace(starttime, endtime, steps)
result = experiment()
plt.plot(time*1000, smoothing(result[0]), color = 'red')
plt.plot(time*1000, smoothing(result[1]), '--', color = 'blue')
plt.plot(time*1000, 15*np.ones(steps))
plt.xlabel('Time(ms)')
plt.ylabel('Firing rate(Hz)')
plt.ylim(top = 20)
plt.text(-75, 16, 'Threshold')
plt.show()'''

plt.figure

for cprime in [0, 10]:
	for i in range(10):
		result = experiment()
		if cprime == 0: hue = 'red'
		else: hue = 'blue'
		plt.plot(time*1000, smoothing(result[0]), color = hue)
		plt.plot(time*1000, smoothing(result[1]), color = hue)
plt.plot(time*1000, 15*np.ones(steps))
plt.xlabel('Time(ms)')
plt.ylabel('Firing rate(Hz)')
plt.ylim(top = 20)
plt.text(-125, 16, 'Threshold')
plt.show()