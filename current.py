import matplotlib.pyplot as plt
import sqlite3 as sqlite
import pyneurlib as pdl
import numpy as np

cell = pdl.RGC_Neuron()

[mydt,mydur,triAmp]=[.01,400,.005]

#t=np.arange(mydur/mydt)*mydt
#v=pdl.make_triphasic(t,triAmp)
sim = pdl.Simulation(cell,dt=mydt,amp=.01,dur=mydur)
#sim.set_IClamp([t,v])
sim.set_IClamp()
sim.go()
sim.show(showStim=True)
