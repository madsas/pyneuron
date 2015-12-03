import pytestlib as ptl
import pyneurlib as pdl
import numpy as np
import matplotlib.pyplot as plt

cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon.csv',ex_flag = True)

#stimulation parameters
mydt = 0.01
mysimtime = 7
mydelay = 1
myrho = 7900

#seeking parameters
init_amp = 1e3
init_step = .2e3
x_soma = 15
x_axon = 270


#looking at waveforms directly dpp
dur1 = .4
dur2 = .1
amp = 500
dppamp = 300
y = 10

stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]

[t,i] = pdl.make_dppbal(mydelay, dur1, dur2, dppamp, amp, mysimtime, mydt)
sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
sim.set_exstim([t,i],x_dist = x_soma, y_dist = y, rho = myrho)
sim.go()
sim.show(showAx=0)
plt.figure()
plt.plot(t,i)

