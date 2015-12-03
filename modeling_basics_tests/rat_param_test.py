import pyneurlib as pdl
import numpy as np
import matplotlib.pyplot as plt

#cell = pdl.RGC_Neuron('fohlmeister_geo_params_rat.csv',ex_flag = True)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params.csv',ex_flag = True)
cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultralight.csv',ex_flag = True)
mydt = 0.005
mysimtime = 4
mydelay = 1
myrho = 7900 #ec medium linear resistance in ohm*cm
#myrho = 294 #ec medium linear resistance in ohm*cm
mydur = 1
myamp = 1500 #nA
x = 1300 #above the soma
#x = 40
y = 5 #z/sqrt(2) <--Rattay
[t,i] = pdl.make_square(mydelay,mydur,myamp,mysimtime,mydt)
sim = pdl.Simulation(cell,mydt,sim_time=mysimtime)
sim.set_exstim([t,i],x,y,myrho)
sim.go()
t,soma,axon = sim.get_recording()

plt.close('all')
plt.plot(t,axon)
plt.plot(t,soma)
