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
cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True, ex_rand = 10)
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
mydt = .001
mydur = 1
mysimtime = 10
#myrho = 7900
#myrho = 79000
myrho = 7900
mydelay = 2

#threshold finding parameters
init_triamp = .5e3
init_ampstep = .2e3
amp_thr = .025e3
thr = 0 #mV

x_dist_vals = [0, 15, 55, 270]

for y in [3,30,100]:
	out = [] #to store output
	for x_dist_val in x_dist_vals: 
		print
		print
		print x_dist_val
		triamp = init_triamp
		ampstep = init_ampstep
		downflg = 0
		cnt = 0
		while ampstep >= amp_thr: 
			print triamp
			msom_rec = max(get_max_sq(cell, triamp, x_dist_val,y)[0])
			max_rec = max(get_max_sq(cell, triamp, x_dist_val,y)[1])

			if (msom_rec > thr or max_rec > thr): 
				triamp -= ampstep
				downflg = 1
			elif downflg:
				ampstep *= 0.5
				triamp += ampstep
			else: triamp +=ampstep #don't reduce ampstep if the thing has never gone down

			cnt += 1
			if cnt > 15: break #make sure things don't go overboard
		out.append(triamp)

	plt.plot(x_dist_vals,out)



