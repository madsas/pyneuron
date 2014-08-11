import pytestlib as ptl
import pyneurlib as pdl
import numpy as np
import matplotlib.pyplot as plt

#define cell
#ADD RANDOMNESS LATER
#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True)
cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = True,ex_rand = 24)

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

y_soma_sq = 95
y_axon_sq = 63
y_soma_tp = 28
y_axon_tp = 1
y_soma_tp_rand = 44
y_axon_tp_rand = 30
y_soma_tp_rand_norm = 30
y_axon_tp_rand_norm = 26

#Just add types as you go
TEST_TYPE = 7

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
	stimrange = np.linspace(.6e3,1.6e3,50)
	adj_dist = [y_soma_tp_rand_norm,y_axon_tp_rand_norm]
	bigbigout = []
	for xind,x in enumerate([x_soma, x_axon]):
		bigout = []
		for triamp in stimrange:
			pr = []
			for i in range(150):
				[msom_rec,max_rec] = ptl.get_max(cell, 'tp',triamp, x, adj_dist[xind], stim_params)
				if (max(msom_rec) > 0 or max(max_rec) > 0): pr.append(1)
				else: pr.append(0)
			bigout.append(np.average(pr))
		bigbigout.append(bigout)
		plt.plot(stimrange,bigout)
		plt.xlabel('Stimulation Amplitude [nA]')
		plt.ylabel('Spike Probability')
		plt.title('Spiking vs. Triphasic Stimulus Amplitude')

				
#plot distance vs. threshold
if TEST_TYPE == 8:
	y0 = 60
	mydur = .3
	stim_params = [mydelay, mydur, mysimtime, mydt, myrho]
	yrange = np.linspace(0,y0,y0)

	#soma
	somap = []
	for y in yrange:
		print y
		thr =ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, y, stim_params,cnt_lim = 30) 
		somap.append(thr)

	axonp = []
	for y in yrange:
		print y
		thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, y, stim_params,cnt_lim = 30)
		axonp.append(thr)
	
	plt.figure()
	plt.plot(yrange, somap)
	plt.plot(yrange, axonp)
	plt.show()

#loop over duration
if TEST_TYPE == 9:
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
	
if TEST_TYPE == 10:
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
if TEST_TYPE == 11:
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

#make 3d dpp plot
if TEST_TYPE == 12:
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
	
#make 3d dppbal plot
if TEST_TYPE == 13:
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
			thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_soma, y_soma_tp_rand_norm, stim_params, cnt_lim = cntlim, dpp_amp = dpp)
			outsoma[durind, dppind] = thr
			thr = ptl.get_thresh(cell, 'dppbal', init_amp, init_step, x_axon, y_axon_tp_rand_norm, stim_params, cnt_lim = cntlim,dpp_amp = dpp)
			outaxon[durind, dppind] = thr

	fig, ax = plt.subplots()
	im = ax.imshow(outsoma, cmap=cm.RdBu, vmin=-9000, vmax=2000, extent=[.01e3, 3e3, .2, 3], aspect = 'auto')
#	im.set_interpolation('bilinear')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('DPP Amplitude [nA]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. DPP (Balanced) Amplitude and Duration (Soma)')

	fig, ax = plt.subplots()
	im = ax.imshow(outaxon, cmap=cm.RdBu, vmin=(outaxon).min(), vmax=(outaxon).max(), extent=[.01e3, 3e3, .2, 3], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)

	plt.xlabel('DPP Amplitude [nA]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. DPP (Balanced) Amplitude and Duration (Axon)')

#make 3d square wave plot duration vs. distance
if TEST_TYPE == 14:
	cntlim = 30
	durlist = [np.round(i,1) for i in np.linspace(0,4,20)]
	distlist = np.linspace(1,100,20)
	outsoma = np.zeros((len(durlist),len(distlist)))
	outaxon = np.zeros((len(durlist),len(distlist)))
	for durind, dur in enumerate(durlist):
		dur1 = dur - .1
		dur2 = .1
		stim_params = [mydelay, dur, mysimtime, mydt, myrho]
		for distind, dist in enumerate(distlist):
			print 'dur:',dur,'dist:',dist
			thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_soma, dist, stim_params, cnt_lim = cntlim)
			outsoma[durind, distind] = thr
			thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x_axon, dist, stim_params, cnt_lim = cntlim)
			outaxon[durind, distind] = thr

	fig, ax = plt.subplots()
	im = ax.imshow(outsoma, cmap=cm.RdBu, vmin=-0, vmax=1000, extent=[min(distlist), max(distlist), min(durlist), max(durlist)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('Electrode Distance [microns]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. Square Pulse Amplitude and Duration (Soma)')

	fig, ax = plt.subplots()
	im = ax.imshow(outaxon, cmap=cm.RdBu, vmin=0, vmax=1000, extent=[min(distlist), max(distlist), min(durlist), max(durlist)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	plt.xlabel('Electrode Distance [microns]')
	plt.ylabel('Pulse duration [ms]')
	plt.title('Threshold vs. Square Pulse Amplitude and Duration (Soma)')
	
