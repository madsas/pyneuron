import pyneurlib as pdl
import pytestlib as ptl
import numpy as np
import matplotlib.pyplot as plt

cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultralight.csv',ex_flag = True)

#stimulation parameters
mydt = 0.01
mysimtime = 15
mydelay = 1
myrho = 7900

#seeking parameters
#init_amp = -1e3
init_amp = -100
init_step = .2e3
#these are based on ultralight params
x_soma = 15
x_ah = 50
x_ais = 120
x_axon = 500
x_list = [x_soma, x_ah, x_ais, x_axon]
x_list_label = ['Soma','AIS','Narrow Segment','Distal Axon']
y_list = [5,10,20,40,60,100] #microns

durange = np.linspace(.01,10,20)
for y in y_list:
	data = np.zeros((len(x_list), len(durange)))
	for ix,x in enumerate(x_list):
		for imydur,mydur in enumerate(durange):
			print 'y=',y,x, mydur
			stim_params = [mydelay, mydur, mysimtime, mydt, myrho]
			thr = ptl.get_thresh(cell, 'sq', init_amp, init_step, x, y, stim_params,cnt_lim = 30)
			print thr
			data[ix, imydur] = thr

#plot
	plt.figure()
	hlist = []
	#for ii,i in enumerate(np.transpose(data)):
	for ii,i in enumerate(data):
		grot = plt.plot(durange,i,label=x_list_label[ii])
		hlist.append(grot)
	legend(hlist)
	plt.title('Strength Duration Curves at '+str(y)+' microns')
	plt.xlabel('Time [ms]')
	plt.ylabel('Threshold [uA]')
