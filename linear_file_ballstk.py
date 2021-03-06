import neuron
from neuron import h
import numpy as np
import csv

#Parameters
#matrix of format L, diam, nseg, gkbar, gabar, gcabar, gkcbar, gnabar
#pas.e, pas.g, cad.depth, and cad.taur are all assumed to be fixed
#the row names are soma, axon hillock, axon narrow part, and distal axon

#these params are based off of EC2.5
#paramFile='fohlmeister_geo_params.csv'
paramFile='fohlmeister_geo_params.csv'
grot=[]
with open(paramFile,'rb') as f:
	wri=csv.reader(f,delimiter=',')
	for row in wri:
		grot.append(row)
#shave and flip and make int
grot2=np.transpose(np.array([i[1:] for i in grot[1:]]))
params=[[] for i in range(grot2.shape[0])] #this has to be a list to contain mixed type
for i in range(grot2.shape[0]):
	for j in range(grot2.shape[1]):
		try: 
			params[i].append(int(grot2[i,j]))
		except ValueError:
			params[i].append(float(grot2[i,j]))


#make dict of parts
nams=['soma','ah','narrow','distal']
grot=[neuron.h.Section() for i in range(4)]
parts=dict(zip(nams,grot))
params=dict(zip(nams,params))

#define geometry
for i in nams: 
	[parts[i].L,parts[i].diam,parts[i].nseg]=list((params[i][:3]))
	
#insert passive
	parts[i].insert('pas')
#	parts[i].Ra=110
	parts[i].Ra=110e3
	for seg in parts[i]:
		seg.pas.e=-65
		seg.pas.g=.00005

	parts[i].insert('spike')
	h('ena=35')
	h('ek=-75')
	for seg in parts[i]:
		[seg.spike.gkbar,seg.spike.gabar,seg.spike.gcabar,seg.spike.gkcbar,seg.spike.gnabar]=list(params[i][3:])

	parts[i].insert('cad')
	for seg in parts[i]:
		seg.cad.depth=3
		seg.cad.taur=10 #are these true still? 
	
	parts[i].insert('extracellular')
	h('xraxial=1e9')
	h('xg=1e9')
	h('xc=0')

#connect them all
parts['soma'].connect(parts['ah'],0,0)
parts['ah'].connect(parts['narrow'],0,0)
parts['narrow'].connect(parts['distal'],0,0)
h('celsius=22')

"""
stim=neuron.h.IClamp(cell(.5))
stim.delay=100
stim.amp=0.02
stim.dur=400
"""

neuron.h.dt=.001

"""
t_ext = neuron.h.Vector(np.arange(20 / neuron.h.dt + 1) * neuron.h.dt)
for seg in cell: 
	v_ext=neuron.h.Vector(np.random.rand(len(t_ext)-.5))
	v_ext.play(seg._ref_v,t_ext)
smoke you
"""
"""
stim=neuron.h.IClamp(parts['soma'](.5))
stim.delay=100
stim.dur=300
stim.amp=.02
"""
#t_ext = neuron.h.Vector(np.arange(400 / neuron.h.dt) * neuron.h.dt)
#v_ext=neuron.h.Vector(len(t_ext))
#for j,i in enumerate(t_ext): 
#	v_ext.x[j]=np.sin(i)
#v_ext=neuron.h.Vector(np.sin(t_ext))
#v_ext.play(stim._ref_amp,t_ext)
#neuron.h.Vector(np.sin(t_ext)).play(stim._ref_amp,t_ext)

#EXTRACELLULAR
#convention is that soma is 0 x and 0 y
#pos=[parts['soma'].L/2, 50] #in microns
pos=[0, 50] #in microns
t_ext=neuron.h.Vector(np.arange(500/neuron.h.dt)*neuron.h.dt)
i_ext=neuron.h.Vector(np.concatenate((np.zeros(100/neuron.h.dt),(np.ones(300/neuron.h.dt)*250*5),np.zeros(100/neuron.h.dt)))) #100 pA
#rext=.2 #ohms times meter TOO LOW
rext=7900
cumLen=0 #counter
for i in nams: 
	for seg in parts[i]:
		dist=np.sqrt(pow(pos[0]-cumLen,2)+pow(pos[1],2))
		v_ext=neuron.h.Vector(1000*np.array(i_ext)*1e-9*rext/(np.pi*4*dist*1e-6))
		v_ext.play(seg._ref_e_extracellular,t_ext)
		cumLen+=15


rec_v=neuron.h.Vector()
rec_v.record(parts['soma'](.5)._ref_v)
rec_va=neuron.h.Vector()
rec_va.record(parts['distal'](.5)._ref_v)
rec_d=neuron.h.Vector()
rec_d.record(parts['distal'](1)._ref_v)
rec_db=neuron.h.Vector()
rec_db.record(parts['distal'](0)._ref_v)
rec_t=neuron.h.Vector()
rec_t.record(neuron.h._ref_t)

neuron.h.finitialize(-65) #self.cell.E
neuron.init()
neuron.run(500)

import numpy as np
time=np.array(rec_t)
voltage_soma=np.array(rec_v)
voltage_axon=np.array(rec_va)
voltage_distal=np.array(rec_d)
voltage_distal_backup=np.array(rec_db)

import matplotlib.pyplot as plt
plt.plot(time,voltage_soma,label="soma")
plt.plot(time,voltage_axon,label="middle axon")
plt.plot(time,voltage_distal,label="distal axon")
plt.plot(time,voltage_distal_backup,label="distal axon backup")
plt.legend(loc="upper right")
plt.show()

