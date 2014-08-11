import neuron
from neuron import h
import numpy as np


cell=neuron.h.Section()
cell.L=25
cell.diam=25
cell.nseg=2

cell.insert('pas')
cell.Ra=110
for seg in cell:
	seg.pas.e=-65
	seg.pas.g=.00005
#	seg.pas.g=.005

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

cell.insert('extracellular')

"""
stim=neuron.h.IClamp(cell(.5))
stim.delay=100
stim.amp=0.02
stim.dur=400
"""

#neuron.h.dt=.001
dt = .01
neuron.h.dt=dt


"""
t_ext = neuron.h.Vector(np.arange(300.3 / neuron.h.dt) * neuron.h.dt)
for seg in cell: 
	unit=np.ones(10/neuron.h.dt)
	amp=.1
	resist=7900000
	rho=resist/(4*np.pi*2e-6)
	volt=unit*amp*1e-9*rho*1000
	v_ext=neuron.h.Vector(np.concatenate((np.zeros(300/neuron.h.dt)*-65,volt*2/3.0,volt*-1.0,volt*1/3.0)))
	v_ext.play(seg._ref_e_extracellular,t_ext)
"""
amp = 140
t_ext = neuron.h.Vector(np.arange(400.3 / neuron.h.dt) * neuron.h.dt)
v_ext = neuron.h.Vector(np.concatenate((np.zeros(200/dt),np.ones(.1/dt)*1/3*amp,np.ones(.1/dt)*-1*amp,np.ones(.1/dt)*2/3*amp,np.zeros(200/dt))))
#for seg in cell:
v_ext.play(cell(0)._ref_e_extracellular,t_ext)

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
"""

rec_v=neuron.h.Vector()
rec_v.record(cell(.5)._ref_v)
rec_t=neuron.h.Vector()
rec_t.record(neuron.h._ref_t)
"""
imem = neuron.h.List()
memirec = neuron.h.Vector()
memirec.record(seg._ref_i_membrane,neuron.h.dt)
imem.append(memirec)
"""

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

