import pyneurlib as pdl
import numpy as np
import matplotlib.pyplot as plt

#cell = pdl.RGC_Neuron('fohlmeister_geo_params_ultralight.csv',ex_flag = True)
cell = pdl.RGC_Neuron('fohlmeister_geo_params_light.csv',ex_flag = True)
mydt = 0.05
mysimtime = 10
mydelay = 1
myrho = 7900 #ec medium linear resistance
mydur = 1.5
myamp = -1500
#x = 600
x = 700
y = 30 
[t,i] = pdl.make_triphasic(mydelay,mydur,myamp,mysimtime,mydt)
sim = pdl.Simulation(cell,mydt,sim_time=mysimtime)
sim.set_exstim([t,i],x,y,myrho)
sim.go(spaceflg=True)
t,soma,axon = sim.get_recording()

#plot
seg_len = 15
savstr = 'how_long_to_normalize_axonal_stim'
plt.close('all')
plt.plot(arange(len(axon))*15,axon)
plt.title('Superposition of voltage traces over axonal space')
plt.ylabel('Amplitude [mV]')
plt.xlabel('Position on axon [um]')
savefig(savstr+'_curves.png',bbox_inches='tight')
	#heatmap
plt.figure()
plt.pcolor(arange(len(axon))*seg_len,arange(axon.shape[1])*mydt,np.transpose(axon))
#plt.axis([0,len(axon)*seg_len,0,axon.shape[1]*mydt])
plt.axis('tight')
plt.title('Axonal voltage over time')
plt.ylabel('Time [ms]')
plt.xlabel('Position on axon [um]')
cb = plt.colorbar()
cb.set_label('Amplitude [mv]')
savefig(savstr+'_colormap.png',bbox_inches='tight')
