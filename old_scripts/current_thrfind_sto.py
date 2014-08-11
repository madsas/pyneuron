import matplotlib.pyplot as plt
import sqlite3 as sqlite
import pyneurlib as pdl
import numpy as np

def get_max_sq(cell, amp, x_dist_val,y): #square wave
	[t,i]=pdl.make_square(mydelay,amp,mydur,mysimtime,mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist=x_dist_val,y_dist = y,rho = myrho)
	sim.go()
#	plt.figure()
#	sim.show()
	_,som_rec,ax_rec = sim.get_recording()
	return som_rec,ax_rec

#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultralight.csv',ex_flag = True, ex_rand = 70)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_light.csv',ex_flag = True)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params.csv',ex_flag = True)

"""
global mydt, mydur, mysimtime, myrho, mydelay
mydt = .01
mydur = .3
mysimtime = 200
myrho = 7900
mydelay = 100
"""

global mydt, mydur, mysimtime, myrho, mydelay
#mydt = .01
mydt = .005
#mydt = .04
mydur = .1
mysimtime = 5
#myrho = 7900
#myrho = 79000
myrho = 7900
mydelay = 2

#x_dist_vals = [0, 15, 55, 270]
#x_dist_vals = [15,270]
#x_dist_vals = [15]
x_dist_vals = [270]
#triamps = [np.linspace(0,.15e3,10), np.linspace(0,.1e3,10), np.linspace(0,.6e3,10), np.linspace(0,.6e3,10)]
#triamps = [np.linspace(.02e3,.1e3,40), np.linspace(.2e3,.6e3,40)]
#triamps = [np.linspace(.043e3,.09e3,70), np.linspace(.32e3,.4e3,70)]
triamps = [np.linspace(.43e3,.53e3,150)]
y = 3 #fix this

bigbigout = []
#for randamp in [10,50,100,200]:
for xind,x_dist_val in enumerate(x_dist_vals): 
#	bigout = [[],[],[]]
	bigout = []
	for triamp in triamps[xind]:
#		for rind,randamp in enumerate([10,20,30]):
		for randamp in [22]:
#		for rind,randamp in enumerate([20,25,30]):
			#print 'RandAmp:',randamp
			cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True, ex_rand = randamp)
			pr = []
			for i in range(20): 
				#print
				#print 'Round:',i
				msom_rec = max(get_max_sq(cell, triamp, x_dist_val,y)[0])
				max_rec = max(get_max_sq(cell, triamp, x_dist_val,y)[1])
				if (msom_rec > 0 or max_rec > 0): out=1
				else: out=0
				#print out
				pr.append(out)
			bigout.append(np.average(pr))

	plt.figure()
	plt.plot(triamps, bigout)
#	plt.plot(bigout[0])
#	plt.plot(bigout[1])
#	plt.plot(bigout[2])
	plt.ylim(-.5,1.5)
	bigbigout.append(bigout)
