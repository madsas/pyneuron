import matplotlib.pyplot as plt
import sqlite3 as sqlite
import pyneurlib as pdl
import numpy as np

def get_max_tp(cell, amp, x_dist_val,y): #square wave
	[t,i]=pdl.make_triphasic(mydelay,mydur,amp,mysimtime,mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist=x_dist_val,y_dist = y,rho = myrho)
	sim.go()
	_,som_rec,ax_rec = sim.get_recording()
	return som_rec,ax_rec

def get_thresh(cell, init_amp, init_step, x, y, spike_thr = 0, amp_thr = .025e3, cnt_lim = 15):
	amp = init_amp
	step = init_step
	downflg = 0
	cnt = 0
	while step >= amp_thr: 
		msom_rec = max(get_max_tp(cell, amp, x,y)[0])
		max_rec = max(get_max_tp(cell, amp, x,y)[1])
		if (msom_rec > spike_thr or max_rec > spike_thr): 
			amp -= step
			downflg = 1
		elif downflg:
			step *= 0.5
			amp += step
		else: amp +=step #don't reduce ampstep if the thing has never gone down
		cnt += 1
		if cnt > cnt_lim: break #make sure things don't go overboard
	return amp

#--------------------------------------------------
cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True, ex_rand = 10)

global mydt, mydur, mysimtime, myrho, mydelay
mydt = .01
mydur = .3
mysimtime = 6
mydelay = 1
myrho = 7900
#--------------------------------------------------

#threshold finding parameters
init_amp = .5e3
init_ampstep = .2e3

#distance finding parameters
y0 = 10
x_soma = 15
x_axon = 270

#Calculate threshold for soma
baseline_thr = get_thresh(cell, init_amp, init_ampstep, x_soma, y0)
print 'baseline:',baseline_thr

#adjust distance to axon until soma threshold is matched
for dely in range(0,50,2):
	y_new = y0 - dely
	thr = get_thresh(cell, init_amp, init_ampstep, x_axon, y_new,cnt_lim=20)
	print thr
	if thr == baseline_thr:
		print "Now we found it:",y_new
		


