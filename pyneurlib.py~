import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import neuron
from neuron import h

"""
There are weird rules here about variables. It seems that somehow variables are kept
in permament spots in some workspace as long as they are not destroyed. 
Must be some hoc thing
"""


def make_compartment(length=150, diameter=3, nseg=1):
	"""
	Returns a compartment.

	comp = make_compartment(120, 4) # comp.L: 120; comp.diam: 4; comp.nsg: 1
	comp = make_compartment()       # comp.L: 150; comp.diam: 3; comp.nsg: 1
	"""
	compartment = neuron.h.Section()
	compartment.L = length
	compartment.diam = diameter
	compartment.nseg = nseg
	return compartment


class RGC_Neuron(object):
	"""
	This class will produce RGC_Neuron objects with a standard soma (L=25 um,
	diam=25 um) with dendrites connected on specified sites of the
	soma. For the dendrites the following parameters can be changed:

	* d_sites: can be 'ipsi' or 'contra' or others (tbd)
	* d_length:   length of each dendrite
	* d_diameter: diameter of each dendrite

	To check the morphology with NEURON gui:
	>>> from neuron import gui
	"""
	def __init__(self, d_sites=[],d_length=150, d_diam=3):
		"""
		This method will be executed when you run
		>>> mycell = RGC_Neuron()
		"""
# passive properties
		self.gp = 0.00005 # [S/cm^2]
		self.E = -65 # []
		self.Ra = 110
# active properties
		self.gkbar_spike = 0.012
		self.gabar_spike = 0.036
		self.gcabar_spike = 0.0022
		self.gkcbar_spike = 0.00005
		self.gnabar_spike=.05 #added later and changed from original 
		self.depth_cad=3
		self.taur_cad=10
# set environment parameters
		h('celsius = 22')
#		h('ena=35')
#		h('ek=-75')
# creating compartments
		self.soma = make_compartment(25, 25) #default soma
		self.dendrites = {}
		if d_sites:
			for site in sites:
				dendrite = make_compartment(d_length, d_diam)
				if site == 'contra':
# connecting the contralateral dendrite to the soma
					dendrite.connect(self.soma,1,0)
				elif site == 'ipsi':
# connecting the ipsilateral dendrite to the soma
					dendrite.connect(self.soma, 0,0)
				else:
					print("Here is something wrong...")
					print("site: %s" %(site))
				self.dendrites.update({site: dendrite})

		# initialize parameters
		self.set_passive_parameters(self.gp, self.E, self.Ra)
		self.set_active_parameters(self.gkbar_spike,self.gabar_spike,self.gcabar_spike,self.gkcbar_spike,self.gnabar_spike,self.taur_cad,self.depth_cad)

	def set_passive_parameters(self, gp=0.004, E=-60, rho=200):
		for sec in neuron.h.allsec():
			sec.Ra = rho
			sec.insert("pas")
			for seg in sec:
				seg.pas.g = gp
				seg.pas.e = E

	def set_active_parameters(self,gkbar_spike,gabar_spike,gcabar_spike,gkcbar_spike,gnabar_spike,taur_cad,depth_cad):
		for sec in neuron.h.allsec():
			sec.insert("spike")
			h('ena=35')
			h('ek=-75')
			for seg in sec: 
				seg.spike.gkbar=gkbar_spike
				seg.spike.gabar=gabar_spike
				seg.spike.gcabar=gcabar_spike
				seg.spike.gkcbar=gkcbar_spike
				seg.spike.gnabar=gnabar_spike
		for sec in neuron.h.allsec():
			sec.insert("cad")
			h('ena=35')
			h('ek=-75')
			for seg in sec:
				seg.cad.taur=taur_cad
				seg.cad.depth=depth_cad

	def change_diameter(self, diam):
		for dendrite in self.dendrites.itervalues():
			dendrite.diam = diam

	def change_length(self, length):
		for dendrite in self.dendrites.itervalues():
			dendrite.L = length
      
    
class Simulation(object):
	"""
	Objects of this class control a current clamp simulation. Example of use:
	>>> cell = Cell()
	>>> sim = Simulation(cell)
	>>> sim.go()
	>>> sim.show()

	Default initialization values are  delay=100, amp=.01, dur=300, sim_time=500, dt=0.01
	Make sure to specificy dt yourself when playing a vector in
	"""
	def __init__(self, cell, dt=0.01, delay=100, amp=.01, dur=300, sim_time=500):
		self.cell = cell
		self.sim_time = sim_time
		self.dt = dt
		self.go_already = False
		self.hasCV=False
		self.delay=delay
		self.amp=amp
		self.dur=dur

	def set_IClamp(self,customVector=[]):
		"""
		Initializes values for current clamp.

		If custom vector is played into amplitude, default or given is overriden
		Custom vector is input as 2-element list with first element being delay, and scond being stimulus
		"""
		stim = neuron.h.IClamp(self.cell.soma(0.5))

		if not customVector: 
			stim.delay = self.delay
			stim.amp = self.amp
			stim.dur = self.dur
		else: 
			stim.delay=0
			stim.amp=1e9
			stim.dur=len(customVector[0])*self.dt
			t_ext=neuron.h.Vector(customVector[0])
			v_ext=neuron.h.Vector(customVector[1])
			v_ext.play(stim._ref_amp,t_ext)
			#be careful to keep the Vectors you are playing in existence
			self.v_ext=v_ext
			self.t_ext=t_ext
			self.hasCV=True

		self.stim = stim

	def show(self,titl="Voltage vs. Time",showStim=False):
		if self.go_already:
			x = np.array(self.rec_t)
			y = np.array(self.rec_v)
			plt.plot(x, y)
			plt.title(titl)
			plt.xlabel("Time [ms]")
			plt.ylabel("Voltage [mV]")
			plt.axhline(0,color="black",ls="--")
			if showStim and self.hasCV:
				v = np.array(self.t_ext)
				w = np.array(self.v_ext)
				pltScale=np.abs(max(y))/np.abs(max(w))
				plt.plot(v,w*pltScale,color='r')
		else:
			print("""First you have to `go()` the simulation.""")
		plt.show()
    
	def set_recording(self):
# Record Time
		self.rec_t = neuron.h.Vector()
		self.rec_t.record(neuron.h._ref_t)
# Record Voltage
		self.rec_v = neuron.h.Vector()
		self.rec_v.record(self.cell.soma(0.5)._ref_v)

	def get_recording(self):
		time = np.array(self.rec_t)
		voltage = np.array(self.rec_v)
		return time, voltage

	def go(self, sim_time=None):
		self.set_recording()
		neuron.h.dt = self.dt
		neuron.h.finitialize(self.cell.E)
		neuron.init()
		if sim_time:
			neuron.run(sim_time)
		else:
			neuron.run(self.sim_time)
		self.go_already = True

	def get_tau_eff(self, ip_flag=False, ip_resol=0.01):
		time, voltage = self.get_recording()
		vsa = np.abs(voltage-voltage[0]) #vsa: voltage shifted and absolut
		v_max = np.max(vsa)
		exp_val = (1-1/np.exp(1)) * v_max # 0.6321 * v_max
		ix_tau = np.where(vsa > ( exp_val ))[0][0] 
		tau = time[ix_tau] - self.stim.delay 
		return tau
  
	def get_Rin(self):
		"""
		This function returnes the input resistance.
		"""
		_, voltage = self.get_recording()
		volt_diff = max(voltage) - min(voltage)
		Rin = np.abs(float(volt_diff / self.stim.amp))
		return Rin

#Random functions
def make_triphasic(t,maxAmp):
	part=np.floor(len(t)/4.0)
	unit=np.ones(part)
	return np.concatenate((unit*maxAmp*(2.0/3),unit*-maxAmp,unit*-maxAmp,unit*maxAmp*(1.0/3)))

