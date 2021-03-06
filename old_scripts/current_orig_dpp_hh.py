import matplotlib.pyplot as plt
import sqlite3 as sqlite
import pyneurlib as pdl
import numpy as np

cell = pdl.HH_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultralight.csv',ex_flag = True, ex_rand = 70)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_light.csv',ex_flag = True)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params.csv',ex_flag = True)

mydt = .01
#mysimtime = 500
mysimtime = 8
myrho = 7900
mydelay = 1.5
#dur1 = .5
dur1 = 4
dur2 = .1
dpAmp = .1e3
stimAmp = .6e3

x_dist_val = 55


for x_dist_val in [0,15,55,250]:
	[t,i]=pdl.make_dpp(mydelay,dur1,dur2,mydt,mysimtime,dpAmp,stimAmp)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist=x_dist_val,y_dist=3,rho = myrho) 
	sim.go()
	plt.figure()
	sim.show(showAx=False)

	plt.plot(t,i/max(i)*20)

	[t,i]=pdl.make_square(mydelay+dur1,stimAmp,dur2,mysimtime,mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist=x_dist_val,y_dist=3,rho = myrho) 
	sim.go()
	sim.show(showAx=False)
	plt.plot(t,i/max(i)*20)


