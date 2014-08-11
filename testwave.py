import matplotlib.pyplot as plt
import pyneurlib as pdl

x,y = pdl.make_dppbal(1,.5,.1,.9e3,1.2e3,10,.01)
plt.plot(x,y)
plt.show()
