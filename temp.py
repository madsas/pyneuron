x=neuron.h.Section()
x.L=20
x.diam=20
x.nseg=1

x.insert('extracellular')
for seg in x: 
	seg.extracellular.xc[0]=0

