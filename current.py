import matplotlib.pyplot as plt
import sqlite3 as sqlite
import pyneurlib as pdl
import numpy as np

cell = pdl.RGC_Neuron('fohlmeister_geo_params_light.csv',True,True)

mydt = .01
mydur = .3
mysimtime = 500
mydelay = 100
bigv = []
reps = range(10)
for i_amp in [.2e3, .4e3, .5e3, .7e3, 1e3, 1.2e3, 1.5e3, 2e3, 2.5e3]: #microamps
	smallv = []
	for rep in reps: 
		print i_amp

		triamp = i_amp
		[t,i]=pdl.make_triphasic(mydelay,mydur,mydt,mysimtime,triamp)
		i = i + .1*(np.random.rand(len(t))-0.5)
		sim = pdl.Simulation(cell,mydt)
#sim.set_IClamp()
		myrho = 4e6
		sim.set_exstim([t,i],x_dist=10,y_dist=7,rho = myrho) #this rho isn't the best!
		sim.go()
		smallv.append(sim.get_recording()[1])
	bigv.append(smallv)
#	sim.show()

