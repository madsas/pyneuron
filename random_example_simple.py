#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import neuron
from neuron import h

#initialize
dt=0.001
soma=neuron.h.Section(name='soma')
soma.L=25
soma.diam=25
soma.nseg=1
soma.insert('pas')
soma.insert('hh')
soma.insert('extracellular')
soma.xc[0]=0
soma.xg[0]=1e9
soma.xraxial[0]=1e9

#extracellular vector
amp=-140
t_ext = neuron.h.Vector(np.arange(50/dt)*dt)
v_ext=neuron.h.Vector(np.concatenate((np.zeros(1/dt),np.ones(1/dt)*amp,np.zeros(48/dt))))
for seg in soma:
	v_ext.play(seg._ref_e_extracellular,t_ext)

#stimluate and record
rec_v=neuron.h.Vector()
rec_v.record(soma(.5)._ref_v)
rec_t=neuron.h.Vector()
rec_t.record(neuron.h._ref_t)

neuron.h.finitialize(-65) 
neuron.init()
neuron.run(50)

time=np.array(rec_t)
voltage=np.array(rec_v)
plt.plot(time,voltage)
plt.show()

