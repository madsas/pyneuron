import neuron
from neuron import h


soma = neuron.h.Section()

#passive
soma.Ra = 110
soma.insert('na12')
soma.insert('pas')
h('celsius = 22')

for seg in soma:
	seg.pas.g = 0.00005
	seg.pas.e = -65
	h('ena=35')
	h('ek=-75')


	
