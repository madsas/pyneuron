import matplotlib.pyplot as plt
import sqlite3 as sqlite
import pyneurlib as pdl
import numpy as np

cell = pdl.RGC_Neuron()

[mydt,mydur,triAmp]=[.01,400,.005]

t=np.arange(mydur/mydt)*mydt
#v=pdl.make_triphasic(t,triAmp)
v=np.sin(t/10)
sim = pdl.Simulation(cell,dt=mydt,amp=.04,dur=mydur,delay=0)
sim.set_IClamp([t,v])
#sim.set_IClamp()
sim.go()
sim.show(titl="Sine",showStim=True)

