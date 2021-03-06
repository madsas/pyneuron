import pyneurlib as pdl
import matplotlib.pyplot as plt

plt.close('all')
#create cell
#cell = pdl.RGC_Neuron_nav('fohlmeister_geo_params_ultraultralight-nona.csv',ex_flag = False)
cell = pdl.RGC_Neuron_nav('fohlmeister_geo_params_ultraultralight-nona.csv',ex_rand = 24)
#cell2 = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_flag = False)
cell2 = pdl.RGC_Neuron('fohlmeister_geo_params_ultraultralight.csv',ex_rand=24)

#Parameters
mydt = 0.01
mysimtime = 500
mydelay = 100
myrho = 7900

#stimluate
sim = pdl.Simulation(cell, dur = 300, amp = 10)
sim.set_IClamp()
sim.go()
sim.show()

sim2 = pdl.Simulation(cell2, dur = 300)
sim2.set_IClamp()
sim2.go()
sim2.show()
