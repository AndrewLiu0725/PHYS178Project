import matplotlib.pyplot as plt
import numpy as np

gL_pyramidal = 25
taum_pyramidal = 20
gL_interneuron = 20
taum_interneuron = 10

VL = -70
Vth = -50
Vrest = -55
VE = 0
VI = -70
absrefractoryperiod_pyramidal = 2
absrefractoryperiod_interneuron = 1

gAMPA_etoe = 0.000237
gAMPA_etoi = 0.000289
gGABA_itoe = 0.0013
gGABA_itoi = 0.001

wplus = 2.4
w = 1
f = 0.15
wminus = 1 - f*(wplus - 1)/(1 - f) 

tauAMPA = 2
tauGABA = 2
taur = 2

dt = 0.1
interval = 50
steps = int(interval/dt)
time = np.linspace(0, interval, steps)
mu0 = 30
Jext = 0.2243*10**-3
cprime = 0

r1, r2, r3, r4 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
#r1[0], r2[0], r3[0], r4 [0] = 1, 1, 1, 1
S1, S2, S3, S4 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
V1, V2, V3, V4 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)
V1[0], V2[0], V3[0], V4[0] = Vrest, Vrest, Vrest, Vrest
t1, t2, t3, t4 = 0, 0, 0, 0
I1, I2, I3, I4 = np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1), np.zeros(steps+1)

def phie(Isyn):
	return (310*Isyn - 125)/(1 - np.exp(-0.16*(310*Isyn - 125)))
def phii(Isyn):
	return (615*Isyn - 177)/(1 - np.exp(-0.087*(615*Isyn - 177)))

def LIFe(Vi, ti, Ii):
	if Vi > Vth:
		print('ya')
		Vi = Vrest
		ti += dt
	else:
		if 0 < t < 2:
			Vi = Vrest
			ti += dt
		else:
			Vi += (dt/(taum_pyramidal*gL_pyramidal)*(-gL_pyramidal*(Vi - VL) + Ii*10**3))
	return Vi

def LIFi(Vi, ti, Ii):
	if Vi > Vth:
		print('ya')
		Vi = Vrest
		ti += dt
	else:
		if 0 < t < 1:
			Vi = Vrest
			ti += dt
		else:
			Vi += (dt/(taum_interneuron*gL_interneuron)*(-gL_interneuron*(Vi - VL) + Ii*10**3))
	return Vi

for index, t in enumerate(time):
	# 11-variable dynamics system (without NMDA)
	I1[index + 1] = Jext*(1 + cprime/100) + gAMPA_etoe*(V1[index] - VE)*(wplus*S1[index] + wminus*S2[index] + wminus*S3[index]) + gGABA_itoe*(V1[index] - VI)*S4[index]
	I2[index + 1] = Jext*(1 - cprime/100) + gAMPA_etoe*(V2[index] - VE)*(wminus*S1[index] + wplus*S2[index] + wminus*S3[index]) + gGABA_itoe*(V2[index] - VI)*S4[index]
	I3[index + 1] = gAMPA_etoe*(V3[index] - VE)*(wminus*S1[index] + wminus*S2[index] + w*S3[index]) + gGABA_itoe*(V3[index] - VI)*S4[index]
	I4[index + 1] = gAMPA_etoi*(V4[index] - VE)*(w*S1[index] + w*S2[index] + w*S3[index]) + gGABA_itoi*(V4[index] - VI)*S4[index]
	r1[index + 1] = r1[index] + (dt/taur)*(-r1[index] + phie(I1[index]))
	r2[index + 1] = r2[index] + (dt/taur)*(-r2[index] + phie(I2[index]))
	r3[index + 1] = r3[index] + (dt/taur)*(-r3[index] + phie(I3[index]))
	r4[index + 1] = r4[index] + (dt/taur)*(-r4[index] + phii(I4[index]))
	S1[index + 1] = S1[index] + dt*(-S1[index]/tauAMPA + r1[index]*10**-3)
	S2[index + 1] = S2[index] + dt*(-S2[index]/tauAMPA + r2[index]*10**-3)
	S3[index + 1] = S3[index] + dt*(-S3[index]/tauAMPA + r3[index]*10**-3)
	S4[index + 1] = S4[index] + dt*(-S4[index]/tauGABA + r4[index]*10**-3)
	# LIF model 
	V1[index + 1] = LIFe(V1[index], t1, I1[index])
	V2[index + 1] = LIFe(V2[index], t2, I2[index])
	V3[index + 1] = LIFe(V3[index], t3, I3[index])
	V4[index + 1] = LIFi(V4[index], t4, I4[index])
plt.figure(1)
plt.plot(time, r1[1:], 'r', time, r2[1:], 'b-')
#plt.plot(time, 15*np.ones(steps))
plt.xlabel('Time')
plt.ylabel('Firing rate')

plt.figure(2)
plt.plot(time, S1[1:], 'r', time, S2[1:], 'b-')
plt.show()