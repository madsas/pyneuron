import neuron
from neuron import h
import numpy as np

#Parameters
#matrix of format L, diam, nseg, gkbar, gabar, gcabar, gkcbar, gnabar
#pas.e, pas.g, cad.depth, and cad.taur are all assumed to be fixed
#the row names are soma, axon hillock, axon narrow part, and distal axon

#these params are based off of EC2.5
paramFile='fohlmeister_geo_params.csv'
grot=[]
with open(paramFile,'rb') as f:
	wri=csv.reader(f,delimiter=',')
	for row in wri:
		grot.append(row)
#shave and flip
params=np.transpose(np.array([i[1:] for i in grot[1:]]))

#make dict of parts
nams=['soma','ah','narrow','distal']

soma=neuron.h.Section()
soma.L=25
soma.diam=25
soma.nseg=1

cell.insert('pas')
cell.Ra=110
for seg in cell:
	seg.pas.e=-65
	seg.pas.g=.00005

cell.insert('spike')
for seg in cell:
	seg.spike.gkbar=.012
	seg.spike.gabar=.036
	seg.spike.gcabar=.0022
	seg.spike.gkcbar=.00005
	seg.spike.gnabar=.05 #added later

h('celsius=22')
h('ena=35')
h('ek=-75')

cell.insert('cad')
for seg in cell:
	seg.cad.depth=3 #working values
	seg.cad.taur=10
#	seg.cad.taur=1.5 #default paper values
#	seg.cad.depth=2.5

"""
stim=neuron.h.IClamp(cell(.5))
stim.delay=100
stim.amp=0.02
stim.dur=400
"""

#neuron.h.dt=.001
neuron.h.dt=.01

"""
t_ext = neuron.h.Vector(np.arange(20 / neuron.h.dt + 1) * neuron.h.dt)
for seg in cell: 
	v_ext=neuron.h.Vector(np.random.rand(len(t_ext)-.5))
	v_ext.play(seg._ref_v,t_ext)
smoke you
"""
stim=neuron.h.IClamp(cell(.5))
stim.delay=0
stim.dur=400
stim.amp=1e9
t_ext = neuron.h.Vector(np.arange(400 / neuron.h.dt) * neuron.h.dt)
#v_ext=neuron.h.Vector(len(t_ext))
#for j,i in enumerate(t_ext): 
#	v_ext.x[j]=np.sin(i)
v_ext=neuron.h.Vector(np.sin(t_ext))
v_ext.play(stim._ref_amp,t_ext)
#neuron.h.Vector(np.sin(t_ext)).play(stim._ref_amp,t_ext)

rec_v=neuron.h.Vector()
rec_v.record(cell(.5)._ref_v)
rec_t=neuron.h.Vector()
rec_t.record(neuron.h._ref_t)

#neuron.h.dt=.001
neuron.h.finitialize(-65) #self.cell.E
neuron.init()
neuron.run(500)

import numpy as np
time=np.array(rec_t)
voltage=np.array(rec_v)

import matplotlib.pyplot as plt
plt.plot(time,voltage)
plt.show()

