#import pyneurlib as pdl
import pyneurlib_temp as pdl
import numpy as np
import matplotlib.pyplot as plt
import pytestlib as ptl
import pickle

#Define Functions
def save_obj(obj, name): 
    with open(name, 'wb') as f: pickle.dump(obj,f)

plt.close('all')

#Universal Parameters

TEST_TYPE = 1 #<--------- CHANGE THIS

x_soma = 15
x_axon = 270
mydelay = 1
mysimtime = 17
#myrho = 7900
myrho = 125 #from abramaniam .8 S/m
mydt = 0.01
init_amp = -1e3
init_step = .2e3
y = 200 #effective distance required to get threshold over soma as in figure 5 of abramaniam

cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_rat.csv',ex_flag = True)

#looking at waveforms directly square
if TEST_TYPE == 1:
    dur = .1
    amp = -50
    y = 200 

    stim_params = [mydelay, dur, mysimtime, mydt, myrho]
    [t,i] = pdl.make_square(mydelay, dur, amp, mysimtime, mydt)
    sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
#    sim.set_exstim([t,i],x_dist = 300, y_dist = y, rho = myrho)
    sim.set_exstim([t,i],x_dist = 20, y_dist = y, rho = myrho)
    sim.go()
    sim.show()
    #plt.figure()
    #plt.plot(t,i)
    #plt.show()


#make 3d square wave plot duration vs. distance
if TEST_TYPE == 2:
	cntlim = 30
	durlist = [np.round(i,1) for i in np.linspace(.1,3,20)]
	distlist = np.linspace(1,400,30)
	outsoma = np.zeros((len(durlist),len(distlist)))
	for durind, dur in enumerate(durlist):
		stim_params = [mydelay, dur, mysimtime, mydt, myrho]
		for distind, dist in enumerate(distlist):
			print 'dur:',dur,'dist:',dist
			thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, dist, y, stim_params, cnt_lim = cntlim)
			print thr
			outsoma[durind, distind] = thr

	fig, ax = plt.subplots()
	im = ax.imshow(outsoma, cmap=plt.cm.RdBu, vmin=0, vmax=8, extent=[min(distlist), max(distlist), min(durlist), max(durlist)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('Electrode Distance along cell [microns]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. Square Pulse Location and Duration (Soma)')
	plt.show()
	
	save_obj(outsoma, 'outsoma.p')

#make 3d DPP Balanced plot prepulse duration vs. distance
if TEST_TYPE == 3:
	cntlim = 30
	durlist = [np.round(i,1) for i in np.linspace(.1,3,10)]
	distlist = np.linspace(1,400,10)
	outsoma = np.zeros((len(durlist),len(distlist)))
	for durind, dur in enumerate(durlist):
		dur1 = dur - .1
		dur2 = .1
		stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
		for distind, dist in enumerate(distlist):
			print 'dur:',dur,'dist:',dist
			thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, dist, y, stim_params, cnt_lim = cntlim, dpp_amp = -50)
			print thr
			outsoma[durind, distind] = thr

	fig, ax = plt.subplots()
	im = ax.imshow(outsoma, cmap=plt.cm.RdBu, vmin=0, vmax=8, extent=[min(distlist), max(distlist), min(durlist), max(durlist)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('Electrode Distance along cell [microns]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. DPP Location and Duration (Soma)')
	plt.show()
	
	save_obj(outsoma, 'outsoma_dppbal.p')
