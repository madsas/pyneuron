import pytestlib as ptl
import pyneurlib as pdl
import numpy as np
import matplotlib.pyplot as plt
import pickle

#Helper Functions--------------------------------------------------
import pickle
def save_obj(obj, name): 
	with open(name, 'wb') as f: pickle.dump(obj, f)

def load_obj(name):
	with open(name,'rb') as f: return pickle.load(f)

def plot_3d(data, list1, list2, title, list1name, list2name, cbarlabel = 'Stimulation Threshold [uA]', myvmin = 0, myvmax = 1):
	fig, ax = plt.subplots()
	im = ax.imshow(np.flipud(data), cmap=plt.cm.RdBu, vmin=myvmin, vmax=myvmax, extent=[min(list1),max(list1),min(list2),max(list2)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	cb.set_label(cbarlabel, rotation = 270)
	plt.xlabel(list1name)
	plt.ylabel(list2name)
	plt.title(title)
	plt.show()
#End Helper Functions--------------------------------------------------

#define cell
#ADD RANDOMNESS LATER
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True,ex_rand = 24)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params.csv',ex_flag = True,ex_rand = 24)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon.csv',ex_flag = True,ex_rand = 24)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_switched.csv',ex_flag = True,ex_rand = 24)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon.csv',ex_flag = True,ex_rand = 24, dendrite_flag = False)
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon.csv',ex_flag = True)
cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_rat.csv',ex_flag = True)

#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_somaparams.csv',ex_flag = True,ex_rand = 24, dendrite_flag = True)

#stimulation parameters
mydt = 0.01
mysimtime = 4
mydelay = 0
myrho = 7900

#seeking parameters
init_amp = 1e3
init_step = .2e3
x_soma = 15
x_axon = 270

y_soma_sq = 95
y_axon_sq = 63
y_soma_tp = 28
y_axon_tp = 1
y_soma_tp_rand = 44
y_axon_tp_rand = 30
y_soma_tp_rand_norm = 30
y_axon_tp_rand_norm = 26

#Just add types as you go
TEST_TYPE = 25

#Catalog


#Square wave distance finder
if TEST_TYPE == 1:
	y0 = 90
	mydur = .3
	stim_params = [mydelay, mydur, mysimtime, mydt, myrho]
	baseline_thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y0, stim_params)
	print 'baseline:',baseline_thr
	for dely in range(0,50,2):
		y_new = y0 - dely
		print 'y:',y_new
		thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, y_new, stim_params)
		print thr

#Triphasic wave distance finder
if TEST_TYPE == 2:
	y0 = 30
	mydur = .3
	stim_params = [mydelay, mydur, mysimtime, mydt, myrho]
	baseline_thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y0, stim_params)
	print 'baseline:',baseline_thr
	for dely in range(0,y0,2):
		y_new = y0 - dely
		print 'y:',y_new
		thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, y_new, stim_params)
		print thr

#DPP simple at adjusted distances
elif TEST_TYPE == 3:
	dur1 = .5
	dur2 = .1
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
	mydpp_amp = .9e3 #90% of threshold of 1 uA

	thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y_soma, stim_params,dpp_amp = mydpp_amp)
	print thr
	thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y_axon, stim_params,dpp_amp = mydpp_amp)
	print thr

#Loop of DPP over range of values
elif TEST_TYPE == 4:
	dur1 = .5
	dur2 = .1
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]

	for mydpp_amp in np.linspace(0.1e3,.9e3,10):
		print 'DPP_AMP:',mydpp_amp
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y_soma_sq, stim_params,dpp_amp = mydpp_amp)
		print thr
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y_axon_sq, stim_params,dpp_amp = mydpp_amp)
		print thr

#Loop of DPP over range of values at TP equi-potential values
elif TEST_TYPE == 5:
	dur1 = .5
	dur2 = .1
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]

	for mydpp_amp in np.linspace(0.1e3,.9e3,10):
		print 'DPP_AMP:',mydpp_amp
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params,dpp_amp = mydpp_amp)
		print thr
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params,dpp_amp = mydpp_amp)
		print thr

elif TEST_TYPE == 6:
	dur = .6
	stim_params = [mydelay, dur, mysimtime, mydt, myrho]

	thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y_soma, stim_params)
	print thr
	thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, y_axon, stim_params)
	print thr

#Jepson curve triphasic
elif TEST_TYPE == 7:
	dur = .3
	stim_params = [mydelay, dur, mysimtime, mydt, myrho]
	stimrange = np.linspace(.8e3,1.7e3,25)
	adj_dist = [y_soma_tp_rand_norm,y_axon_tp_rand_norm]
	bigbigout = []
	for xind,x in enumerate([x_soma, x_axon]):
		bigout = []
		for triamp in stimrange:
			print triamp
			pr = []
			for i in range(50):
				[msom_rec,max_rec] = ptl.get_max(cell, 'tp',triamp, x, adj_dist[xind], stim_params)
				if (max(msom_rec) > 0 or max(max_rec) > 0): pr.append(1)
				else: pr.append(0)
			bigout.append(np.average(pr))
		bigbigout.append(bigout)
		plt.plot(stimrange,bigout)
		plt.xlabel('Stimulation Amplitude [nA]')
		plt.ylabel('Spike Probability')
		plt.title('Spiking vs. Triphasic Stimulus Amplitude')
	plt.show()

#Jepson curve DPP Balanced
elif TEST_TYPE == 8:
	dur1 = .4
	dur2 = .1
	#dpp = .25e3
	dpp = .1e3
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
	stimrange = np.linspace(.2e3,1.6e3,50)
	#stimrange = np.linspace(.01e3,.6e3,50)
	adj_dist = [y_soma_tp_rand_norm,y_axon_tp_rand_norm]
	bigbigout = []
	for xind,x in enumerate([x_soma, x_axon]):
		bigout = []
		for triamp in stimrange:
			pr = []
			#for i in range(150):
			for i in range(50):
				[msom_rec,max_rec] = ptl.get_max(cell, 'dppbal',triamp, x, adj_dist[xind], stim_params, dppamp = dpp)
				if (max(msom_rec) > 0 or max(max_rec) > 0): pr.append(1)
				else: pr.append(0)
			bigout.append(np.average(pr))
		bigbigout.append(bigout)
		plt.plot(stimrange,bigout)
		plt.xlabel('Stimulation Amplitude [nA]')
		plt.ylabel('Spike Probability')
		plt.ylim(-.5,1.3)
		plt.title('Spiking vs. DPP (Balanced) Stimulus Amplitude')
	plt.show()
	save_obj(bigbigout, 'bigbigout-tt8') 

				
#plot distance vs. threshold for triphasic
if TEST_TYPE == 9:
	y0 = 60
	mydur = .3
	stim_params = [mydelay, mydur, mysimtime, mydt, myrho]
	yrange = np.linspace(0,y0,y0)

	#soma
	somap = []
	print 'soma:'
	for y in yrange:
		print y
		thr =ptl.get_thresh(cell, 'tp', init_amp, init_step, x_soma, y, stim_params,cnt_lim = 30) 
		print thr
		somap.append(thr)

	axonp = []
	print 'axon:'
	for y in yrange:
		print y
		thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_axon, y, stim_params,cnt_lim = 30)
		print thr
		axonp.append(thr)
	
	plt.figure()
	plt.plot(yrange, somap, label='soma')
	plt.plot(yrange, axonp, label='axon')
	plt.ylabel('Stimulation Threshold (nA)')
	plt.xlabel('Distance of electrode (microns)')
	plt.title('Threshold vs. Distance (Triphasic)')
	plt.legend()
	plt.show()

#plot distance vs. threshold for DPP
if TEST_TYPE == 10:
	y0 = 60
	dur1 = .4
	dur2 = .1
	dpp = .4e3
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
	yrange = np.linspace(0,y0,y0)

	#soma
	somap = []
	for y in yrange:
		print y
		thr =ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y, stim_params,cnt_lim = 30, dpp_amp = dpp) 
		somap.append(thr)

	axonp = []
	for y in yrange:
		print y
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y, stim_params,cnt_lim = 30, dpp_amp = dpp)
		axonp.append(thr)
	
	plt.figure()
	plt.plot(yrange, somap, label='soma')
	plt.plot(yrange, axonp, label='axon')
	plt.ylabel('Stimulation Threshold (nA)')
	plt.xlabel('Distance of electrode (microns)')
	plt.title('Threshold vs. Distance (DPP)')
	plt.legend()
	plt.show()

#loop over duration
if TEST_TYPE == 10:
	mysimtime = 10
	durlist = [.05, .1, .3, .6, 1, 2, 3]
	tpsoma = []
	tpaxon = []
	for dur in durlist:
		stim_params = [mydelay, dur, mysimtime, mydt, myrho]
		thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = 30)
		tpsoma.append(thr)
		thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = 30)
		tpaxon.append(thr)
	plt.figure()
	plt.plot(durlist, tpsoma, label='soma')
	plt.plot(durlist, tpaxon, label = 'axon')
	plt.show()
	
if TEST_TYPE == 12:
	mydur = .3
	stim_params = [mydelay, mydur, mysimtime, mydt, myrho]
	thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = 30)
	print thr
	thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = 30)
	print thr

	thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = 30)
	print thr
	thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = 30)
	print thr

	print 'now dpp:'
	dur1 = .2
	dur2 = .1
	dpp = .9e3
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
	thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = 30, dpp_amp = dpp)
	print thr
	thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = 30,dpp_amp = dpp)
	print thr

#Duration vs. Threshold
if TEST_TYPE == 13:
	durlist = [.2, .3, .6, 1, 2, 3]
	cntlim = 50
	tpsoma = []
	tpaxon = []
	tpsomatp = []
	tpaxontp = []
	tpsomadpp = []
	tpaxondpp = []
	for dur in durlist:
		print 'dur',dur
		stim_params = [mydelay, dur, mysimtime, mydt, myrho]
		thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = cntlim)
		tpsoma.append(thr)
		thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = cntlim)
		tpaxon.append(thr)

		thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = cntlim)
		tpsomatp.append(thr)
		thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = cntlim)
		tpaxontp.append(thr)

		dur1 = dur - .1
		dur2 = .1
		dpp = 2e3
		stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = cntlim, dpp_amp = dpp)
		tpsomadpp.append(thr)
		thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = cntlim,dpp_amp = dpp)
		tpaxondpp.append(thr)
	fig = plt.figure(figsize=(8,9))
	ax = fig.add_subplot(111)
	ax.set_position([.1,.1,.5,.8])
	ax.plot(durlist, tpsoma,label = 'Square Pulse Soma')
	ax.plot(durlist, tpaxon, label = 'Square Pulse Axon')
	ax.plot(durlist, tpsomatp, label = 'Triphasic Pulse Soma')
	ax.plot(durlist, tpaxontp, label = 'Triphasic Pulse Axon')
	ax.plot(durlist, tpsomadpp, label = 'DPP Soma')
	ax.plot(durlist, tpaxondpp, label = 'DPP Axon')
	leg = ax.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5))
	plt.show()
	plt.xlabel('Duration [ms]')
	plt.ylabel('Stimulation Threshold [nA]')
	plt.title('Duration vs. Threshold')

#make 3d dpp plot (duration vs. dpp amp)
if TEST_TYPE == 14:
	cntlim = 50
	durlist = [.2, .3, .6, 1, 2, 3]
	dpplist = np.linspace(.01e3,3e3,20)
	outsoma = np.zeros((len(durlist),len(dpplist)))
	outaxon = np.zeros((len(durlist),len(dpplist)))
	for durind, dur in enumerate(durlist):
		dur1 = dur - .1
		dur2 = .1
		stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
		for dppind, dpp in enumerate(dpplist):
			thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = cntlim, dpp_amp = dpp)
			outsoma[durind, dppind] = thr
			thr = ptl.get_thresh(cell, 'dpp', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = cntlim,dpp_amp = dpp)
			outaxon[durind, dppind] = thr

	fig, ax = plt.subplots()
	im = ax.imshow(outsoma, cmap=cm.RdBu, vmin=abs(outsoma).min(), vmax=abs(outsoma).max(), extent=[0, 1, 0, 1])
	im.set_interpolation('bilinear')
	cb = fig.colorbar(im, ax=ax)

	fig, ax = plt.subplots()
	im = ax.imshow(outaxon, cmap=cm.RdBu, vmin=abs(outaxon).min(), vmax=abs(outaxon).max(), extent=[0, 1, 0, 1])
	im.set_interpolation('bilinear')
	cb = fig.colorbar(im, ax=ax)
	
#make 3d dppbal plot (duration vs. dpp amp)
if TEST_TYPE == 15:
	cntlim = 30
	#durlist = [.2, .3, .6, 1, 2, 3]
	#durlist = [np.round(i,1) for i in np.linspace(.2,4,25)]
	durlist = [np.round(i,1) for i in np.linspace(.3,4,25)]
	dpplist = np.linspace(.01e3,3e3,20)
	outsoma = np.zeros((len(durlist),len(dpplist)))
	outaxon = np.zeros((len(durlist),len(dpplist)))
	for durind, dur in enumerate(durlist):
		dur1 = dur - .1
		dur2 = .1
		stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
		for dppind, dpp in enumerate(dpplist):
			print dur, dpp
			thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = cntlim, dpp_amp = dpp)
			print 'soma',thr
			outsoma[durind, dppind] = thr
			thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = cntlim,dpp_amp = dpp)
			print 'axon',thr
			outaxon[durind, dppind] = thr

	plot_3d(outsoma/1000, dpplist, distlist, 'Threshold vs. DPP (Balanced) Amplitude and Duration (Soma)', 'DPP Amplitude [nA]', 'Pulse Duration [ms]')
	plot_3d(outaxon/1000, dpplist, distlist, 'Threshold vs. DPP (Balanced) Amplitude and Duration (Axon)', 'DPP Amplitude [nA]', 'Pulse Duration [ms]')
	plot_3d(abs(outaxon/outsoma), dpplist, distlist, 'Ratio of Threshold vs. DPP (Axon vs. Soma)', 'DPP Amplitude [nA]', 'Pulse Duration [ms]', 'Stimulation Threshold Ratio', myvmax = 6)

	save_obj(outsoma,'outsoma-minaxon-tt15.p')
	save_obj(outaxon,'outaxon-minaxon-tt15.p')


#make 3d square wave plot duration vs. distance
if TEST_TYPE == 16:
	cntlim = 30
	durlist = [np.round(i,1) for i in np.linspace(0,4,50)]
	distlist = np.linspace(1,100,20)
	outsoma = np.zeros((len(durlist),len(distlist)))
	outaxon = np.zeros((len(durlist),len(distlist)))
	for durind, dur in enumerate(durlist):
		stim_params = [mydelay, dur, mysimtime, mydt, myrho]
		for distind, dist in enumerate(distlist):
			print 'dur:',dur,'dist:',dist
			thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, dist, stim_params, cnt_lim = cntlim)
			outsoma[durind, distind] = thr
			thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, dist, stim_params, cnt_lim = cntlim)
			outaxon[durind, distind] = thr

	fig, ax = plt.subplots()
	im = ax.imshow(outsoma, cmap=cm.RdBu, vmin=0, vmax=8, extent=[min(distlist), max(distlist), min(durlist), max(durlist)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('Electrode Distance [microns]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. Square Pulse Amplitude and Duration (Soma)')

	fig, ax = plt.subplots()
	im = ax.imshow(outaxon, cmap=cm.RdBu, vmin=0, vmax=8, extent=[min(distlist), max(distlist), min(durlist), max(durlist)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('Electrode Distance [microns]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. Square Pulse Amplitude and Duration (Soma)')

#make 3d dppbal plot (distance vs. dpp amp)
if TEST_TYPE == 17:
	cntlim = 30
#	durlist = [.2, .3, .6, 1, 2, 3]
	#durations will be 200 50 100
	dur1 = .2 #fix this
	dur2 = .05
	dur3 = .1
	stim_params = [mydelay, dur1, dur2, dur3, mysimtime, mydt, myrho]

	dpplist = np.linspace(.01e3,3e3,20)
	distlist = np.linspace(1,100,20)
	outsoma = np.zeros((len(distlist),len(dpplist)))
	outaxon = np.zeros((len(distlist),len(dpplist)))
	for distind, dist in enumerate(distlist):
		for dppind, dpp in enumerate(dpplist):
			print dist,dpp
			thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_soma, dist, stim_params, cnt_lim = cntlim, dpp_amp = dpp)
			outsoma[distind, dppind] = thr
			thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_axon, dist, stim_params, cnt_lim = cntlim,dpp_amp = dpp)
			outaxon[distind, dppind] = thr

	plot_3d(outsoma/1000, dpplist, distlist, 'Threshold vs. DPP (Balanced) Distance and Duration (Soma)', 'DPP Amplitude [nA]', 'Electrode Distance [microns]')
	plot_3d(outaxon/1000, dpplist, distlist, 'Threshold vs. DPP (Balanced) Distance and Duration (Axon)', 'DPP Amplitude [nA]', 'Electrode Distance [microns]')
	plot_3d(abs(outaxon/outsoma), dpplist, distlist, 'Ratio of Threshold vs. DPP (Axon vs. Soma)', 'DPP Amplitude [nA]', 'Electrode Distance [microns]', 'Stimulation Threshold Ratio', myvmax = 100)

	save_obj(outsoma,'outsoma-tt17.p')
	save_obj(outaxon,'outaxon-tt17.p')

#looking at waveforms directly dpp
if TEST_TYPE == 18:
	dur1 = .2
	dur2 = .05
	dur3 = .2
	amp = 100
	dppamp = 80
	#amp = 0
	#dppamp = 0
	y = 10

	stim_params = [mydelay, dur1, dur2, dur3, mysimtime, mydt, myrho]

	[t,i] = pdl.make_dppbal(mydelay, dur1, dur2, dur3, dppamp, amp, mysimtime, mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist = x_soma, y_dist = y, rho = myrho)
	sim.go()
	sim.show(showAx=0)
	plt.figure()
	plt.plot(t,i)
	plt.show()

#	thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_soma, y, stim_params, cnt_lim = cntlim, dpp_amp = dppamp)
#	print thr

#looking at waveforms directly square
if TEST_TYPE == 19:
	dur = .1
	amp = 300
	y = 1
	cntlim = 30

	stim_params = [mydelay, dur, mysimtime, mydt, myrho]

	[t,i] = pdl.make_square(mydelay, dur, amp, mysimtime, mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist = x_axon, y_dist = y, rho = myrho)
	sim.go()
	sim.show()
	#plt.figure()
	#plt.plot(t,i)

	thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y, stim_params, cnt_lim = cntlim)
	print thr

#looking at waveforms directly triphasic
if TEST_TYPE == 20:
	dur = .3
	#amp = 200
	amp = 2000
	y = y_soma_tp_rand_norm
	y = 30
	x_soma = 0
	#mydt = .1
	cntlim = 30

	mydelay = 0
	mysimtime = 4
	#mysimtime = 15
	stim_params = [mydelay, dur, mysimtime, mydt, myrho]

	cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_somaparams.csv',ex_flag = True, dendrite_flag = False)
	#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_somaparams.csv',ex_flag = True, dendrite_flag = True)

	[t,i] = pdl.make_triphasic(mydelay, dur, amp, mysimtime, mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist = x_soma, y_dist = y, rho = myrho)
	sim.go()
	sim.show(showAx = 0)
	#sim.show(showAx = 2)
	#plt.figure()
	#plt.plot(t,i)

	#thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_soma, y, stim_params, cnt_lim = cntlim)
	#print thr

#Jepson Analysis Extended - soma only
##because we're doing comparison between two conditions, all parameters are specified here
elif TEST_TYPE == 21:

	cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_somaparams.csv',ex_flag = True,ex_rand = 24, dendrite_flag = False)

	x_soma = 0 #very sensitive parameter for some reason
	mydt = 0.01
	mysimtime = 4
	mydelay = 0
	myrho = 7900
	dur = .3
	stim_params = [mydelay, dur, mysimtime, mydt, myrho]

	stimrange = np.linspace(1e3,2.8e3,25)
	adj_dist = y_soma_tp_rand_norm

	bigbigout = []
	for densval in np.linspace(0, 0.003, 10):
	#for diamval in np.linspace(0, .003, 10):
		cell.set_channel_density('gcabar',densval)
		#cell.change_geometry('diam',diamval)
		bigout = []
		for triamp in stimrange:
			print triamp
			pr = []
			for i in range(50):
				msom_rec = ptl.get_max(cell, 'tp',triamp, x_soma, adj_dist, stim_params)
				if (max(msom_rec) > 0): pr.append(1)
				else: pr.append(0)
			bigout.append(np.average(pr))
		bigbigout.append(bigout)
		plt.plot(stimrange,bigout,label=str(densval))
		#plt.plot(stimrange,bigout,label=str(diamval))

	"""
	cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_axonparams.csv',ex_flag = True,ex_rand = 24, dendrite_flag = False)
	bigout = []
	for triamp in stimrange:
		print triamp
		pr = []
		for i in range(50):
			msom_rec = ptl.get_max(cell, 'tp',triamp, x_soma, adj_dist, stim_params)
			if (max(msom_rec) > 0): pr.append(1)
			else: pr.append(0)
		bigout.append(np.average(pr))
	plt.plot(stimrange,bigout)
	"""

	plt.xlabel('Stimulation Amplitude [nA]')
	plt.ylabel('Spike Probability')
	plt.title('Spiking vs. Triphasic Stimulus Amplitude')
	plt.legend(loc=0)

	plt.show()

	with open('vary_ca.p', 'wb') as f: pickle.dump(bigbigout, f)

#looking at waveforms directly triphasic just soma
if TEST_TYPE == 22:
	dur = .3
	amp = 1400
	y = y_soma_tp_rand_norm
	#cntlim = 30

	mydt = 0.01
	mysimtime = 4
	mydelay = 0
	myrho = 7900
	stim_params = [mydelay, dur, mysimtime, mydt, myrho]

	#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_somaparams.csv',ex_flag = True,ex_rand = 24, dendrite_flag = False)
	#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight_minaxon_somaparams.csv',ex_flag = True,ex_rand = 24, dendrite_flag = True)

	[t,i] = pdl.make_triphasic(mydelay, dur, amp, mysimtime, mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist = x_soma, y_dist = y, rho = myrho)
	sim.go()
	sim.show(showAx = 0)
	#plt.figure()
	#plt.plot(t,i)

	#thr = ptl.get_thresh(cell, 'tp', init_amp, init_step, x_soma, y, stim_params, cnt_lim = cntlim)
	#print thr

#make 3d dppbal plot (distance vs. dpp perc)
if TEST_TYPE == 23:
	cntlim = 40
#	durlist = [.2, .3, .6, 1, 2, 3]
	#durations will be 200 50 100
	dur1 = .2 #fix this
	dur2 = .05
	dur3 = .1
	init_amp = -1.5e3
	stim_params = [mydelay, dur1, dur2, dur3, mysimtime, mydt, myrho]
	spk_thr = 10

	dppperclist = np.linspace(.5,1,20)
	distlist = np.linspace(3,20,20)
	outsoma = np.zeros((len(distlist),len(dppperclist)))
	outaxon = np.zeros((len(distlist),len(dppperclist)))
	for distind, dist in enumerate(distlist):
		for dppind, dpp in enumerate(dppperclist):
			print dist,dpp
			thr = ptl.get_thresh_dppperc(cell, 'dppbal', init_amp, init_step, x_soma, dist, stim_params, cnt_lim = cntlim, dpp_perc = dpp, spike_thr = spk_thr)
			outsoma[distind, dppind] = thr
			print thr
			thr = ptl.get_thresh_dppperc(cell, 'dppbal', init_amp, init_step, x_axon, dist - 2, stim_params, cnt_lim = cntlim,dpp_perc = dpp, spike_thr = spk_thr)
			#try subtraction to get ratio 1 for control values
			outaxon[distind, dppind] = thr

#	plot_3d(outsoma/1000, dppperclist, distlist, 'Threshold vs. DPP (Balanced) Electrode Distance and Amplitude (Soma)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]')
#	plot_3d(outaxon/1000, dppperclist, distlist, 'Threshold vs. DPP (Balanced) Electrode Distance and Amplitude (Axon)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]')
#	plot_3d(abs(outaxon/outsoma), dppperclist, distlist, 'Ratio of Threshold vs. DPP (Axon vs. Soma)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]', 'Stimulation Threshold Ratio', myvmax = (outaxon/outsoma).max())

	save_obj(outsoma,'outsoma-tt23.p')
	save_obj(outaxon,'outaxon-tt23.p')

#looking at waveforms directly dpp (unbalanced)
if TEST_TYPE == 24:
	dur1 = .2
	dur2 = .05
#	dur3 = .2
	amp = -100
	dppamp = -80
	#amp = 0
	#dppamp = 0
	y = 10

	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]

	[t,i] = pdl.make_dpp(mydelay, dur1, dur2, dppamp, amp, mysimtime, mydt)
	sim = pdl.Simulation(cell,mydt,sim_time = mysimtime)
	sim.set_exstim([t,i],x_dist = x_soma, y_dist = y, rho = myrho)
	sim.go()
	sim.show(showAx=0)
	plt.figure()
	plt.plot(t,i)
	plt.show()

#	thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_soma, y, stim_params, cnt_lim = cntlim, dpp_amp = dppamp)
#	print thr

#make 3d dpp plot (distance vs. dpp perc) UNBALANCED
if TEST_TYPE == 25:
	cntlim = 30
#	durlist = [.2, .3, .6, 1, 2, 3]
	#durations will be 200 50 100
	dur1 = .2 #fix this
	dur2 = .05
#	dur3 = .1
	init_amp = -1.5e3
	stim_params = [mydelay, dur1, dur2, mysimtime, mydt, myrho]
	spk_thr = 1

	dppperclist = np.linspace(.5,1,20)
	distlist = np.linspace(3,20,20)
	outsoma = np.zeros((len(distlist),len(dppperclist)))
	outaxon = np.zeros((len(distlist),len(dppperclist)))
	for distind, dist in enumerate(distlist):
		for dppind, dpp in enumerate(dppperclist):
			print dist,dpp
			thr = ptl.get_thresh_dppperc(cell, 'dpp', init_amp, init_step, x_soma, dist, stim_params, cnt_lim = cntlim, dpp_perc = dpp, spike_thr = spk_thr)
			outsoma[distind, dppind] = thr
			print thr
			thr = ptl.get_thresh_dppperc(cell, 'dpp', init_amp, init_step, x_axon, dist - 2, stim_params, cnt_lim = cntlim,dpp_perc = dpp, spike_thr = spk_thr)
			#try subtraction to get ratio 1 for control values
			outaxon[distind, dppind] = thr

#	plot_3d(outsoma/1000, dppperclist, distlist, 'Threshold vs. DPP (Balanced) Electrode Distance and Amplitude (Soma)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]')
#	plot_3d(outaxon/1000, dppperclist, distlist, 'Threshold vs. DPP (Balanced) Electrode Distance and Amplitude (Axon)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]')
#	plot_3d(abs(outaxon/outsoma), dppperclist, distlist, 'Ratio of Threshold vs. DPP (Axon vs. Soma)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]', 'Stimulation Threshold Ratio', myvmax = (outaxon/outsoma).max())

	save_obj(outsoma,'outsoma-tt25.p')
	save_obj(outaxon,'outaxon-tt25.p')
