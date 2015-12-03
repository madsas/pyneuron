import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import neuron
from neuron import h
import csv
import sys

"""
There are weird rules here about variables. It seems that somehow variables are kept
in permament spots in some workspace as long as they are not destroyed. 
Must be some hoc thing
"""


"""
def make_compartment(length=150, diameter=3, nseg=1):
	#Returns a compartment.

	#comp = make_compartment(120, 4) # comp.L: 120; comp.diam: 4; comp.nsg: 1
	#comp = make_compartment()       # comp.L: 150; comp.diam: 3; comp.nsg: 1
	compartment = neuron.h.Section()
	compartment.L = length
	compartment.diam = diameter
	compartment.nseg = nseg
	return compartment
	"""


class RGC_Neuron(object):
	"""
	This class will produce RGC_Neuron objects with a standard soma (L=25 um,
	diam=25 um) and with an axon consisting of an axon hillock, narrow region,
	and distal region. 

	To check the morphology with NEURON gui:
	>>> from neuron import gui
	"""
	def __init__(self, param_file, dendrite_flag=True, ex_flag=False, ex_rand=0):
		"""
		This method will be executed when you run
		>>> mycell = RGC_Neuron()
		"""
		self.dflag = dendrite_flag
		self.ex_rand = ex_rand
		if ex_rand !=0: ex_flag = True
#Load parameters
		params = read_param_file(param_file)	
		self.seg_len = params[-1][-1] #can't handle different length segments yet
#make dict of parts
		self.nams = ['soma','ah','narrow','distal'] if self.dflag else ['soma']
		self.parts=dict(zip(self.nams,[neuron.h.Section() for i in range(len(self.nams))]))
		self.params=dict(zip(self.nams,params))

#define geometry
		self.set_geometry()

#store normalized sodium channel densities for random noise scaling
		self.nchan_den = [i[-2] for i in params]
		self.nchan_den = dict(zip(self.nams,np.array(self.nchan_den)/max(self.nchan_den)))
# passive properties
		self.gp = 0.00005 # [S/cm^2]
		self.E = -65 # []
		self.Ra = 110
# some active properties
		self.depth_cad=3
		self.taur_cad=10
# set environment parameters
		h('celsius = 22')

# initialize parameters
		self.set_passive_parameters()
		self.set_active_parameters()
		if ex_flag: self.set_ex()

	def set_geometry(self):
		for i in self.nams:
			[self.parts[i].L,self.parts[i].diam,self.parts[i].nseg]=list((self.params[i][:3]))
		if self.dflag:
			self.parts['soma'].connect(self.parts['ah'],0,0)
			self.parts['ah'].connect(self.parts['narrow'],0,0)
			self.parts['narrow'].connect(self.parts['distal'],0,0)

	def set_passive_parameters(self):
		for sec in neuron.h.allsec():
			sec.Ra = self.Ra
			sec.insert("pas")
			for seg in sec:
				seg.pas.g = self.gp
				seg.pas.e = self.E

	def set_active_parameters(self):
#		for sec in neuron.h.allsec():
		for i in self.nams: 
			self.parts[i].insert("spike")
			h('ena=35')
			h('ek=-75')
			for seg in self.parts[i]: 
				[seg.spike.gkbar,seg.spike.gabar,seg.spike.gcabar,seg.spike.gkcbar,seg.spike.gnabar]=list(self.params[i][3:-1])
		for sec in neuron.h.allsec():
			sec.insert("cad")
			h('ena=35')
			h('ek=-75')
			for seg in sec:
				seg.cad.taur=self.taur_cad
				seg.cad.depth=self.depth_cad

	def set_channel_density(self, channame, val, partnams='all'):
		#'partnams' can be a list of names
		if partnams == 'all': secnams = self.nams
		else: secnams = partnams
		for i in secnams:
			for seg in self.parts[i]:
				setattr(seg.spike, channame, val)
	
	def change_geometry(self, dimname, val, partnams='all'):
		if partnams == 'all': secnams = self.nams
		else: secnams = partnams
		for i in secnams:
			setattr(self.parts[i], dimname, val)

		
	def set_ex(self):
		for sec in neuron.h.allsec():
			sec.insert('extracellular')
			for i in range(2):
				sec.xg[i] = 1e9
				sec.xraxial[i] = 1e9
				sec.xc[i] = 0

#This is the same as the object above except with different SODIUM channel types
class RGC_Neuron_nav(object):
	"""
	This class will produce RGC_Neuron objects with a standard soma (L=25 um,
	diam=25 um) and with an axon consisting of an axon hillock, narrow region,
	and distal region. 

	To check the morphology with NEURON gui:
	>>> from neuron import gui
	"""
	def __init__(self, param_file, dendrite_flag=True, ex_flag=False, ex_rand=0):
		"""
		This method will be executed when you run
		>>> mycell = RGC_Neuron()
		"""
		self.dflag = dendrite_flag
		self.ex_rand = ex_rand
		if ex_rand !=0: ex_flag = True
#Load parameters
		params = read_param_file(param_file)	
		self.seg_len = params[-1][-1] #can't handle different length segments yet
#make dict of parts
		self.nams = ['soma','ah','narrow','distal'] if self.dflag else ['soma']
		self.parts=dict(zip(self.nams,[neuron.h.Section() for i in range(len(self.nams))]))
		self.params=dict(zip(self.nams,params))

#define geometry
		self.set_geometry()

#store normalized sodium channel densities for random noise scaling
		self.nchan_den = [i[-2] for i in params]
		self.nchan_den = dict(zip(self.nams,np.array(self.nchan_den)/max(self.nchan_den)))
# passive properties
		self.gp = 0.00005 # [S/cm^2]
		self.E = -65 # []
		self.Ra = 110
# some active properties
		self.depth_cad=3
		self.taur_cad=10
# set environment parameters
		h('celsius = 22')

# initialize parameters
		self.set_passive_parameters()
		self.set_active_parameters(['soma','ah'])
		if ex_flag: self.set_ex()

	def set_geometry(self):
		for i in self.nams:
			[self.parts[i].L,self.parts[i].diam,self.parts[i].nseg]=list((self.params[i][:3]))
		if self.dflag:
			self.parts['soma'].connect(self.parts['ah'],0,0)
			self.parts['ah'].connect(self.parts['narrow'],0,0)
			self.parts['narrow'].connect(self.parts['distal'],0,0)

	def set_passive_parameters(self):
		for sec in neuron.h.allsec():
			sec.Ra = self.Ra
			sec.insert("pas")
			for seg in sec:
				seg.pas.g = self.gp
				seg.pas.e = self.E

	def set_active_parameters(self,nav16list):
#		for sec in neuron.h.allsec():
		for i in self.nams: 
			self.parts[i].insert("spike_nona")
			if i in nav16list:
				self.parts[i].insert("na16")
				for seg in self.parts[i]: 
					[seg.spike_nona.gkbar,seg.spike_nona.gabar,seg.spike_nona.gcabar,seg.spike_nona.gkcbar,seg.na16.gnabar]=list(self.params[i][3:-1])
			else:
				self.parts[i].insert("na12")
				for seg in self.parts[i]: 
					[seg.spike_nona.gkbar,seg.spike_nona.gabar,seg.spike_nona.gcabar,seg.spike_nona.gkcbar,seg.na12.gnabar]=list(self.params[i][3:-1])
			h('ena=35')
			h('ek=-75')

		for sec in neuron.h.allsec():
			sec.insert("cad")
			h('ena=35')
			h('ek=-75')
			for seg in sec:
				seg.cad.taur=self.taur_cad
				seg.cad.depth=self.depth_cad

	def set_ex(self):
		for sec in neuron.h.allsec():
			sec.insert('extracellular')
			for i in range(2):
				sec.xg[i] = 1e9
				sec.xraxial[i] = 1e9
				sec.xc[i] = 0

class HH_Neuron(object):
	"""
	This class will produce regular HH_Neuron objects with a standard soma (L=25 um,
	diam=25 um) and with an axon consisting of an axon hillock, narrow region,
	and distal region. Channel properties and densities are uniform througout

	To check the morphology with NEURON gui:
	>>> from neuron import gui
	"""
	def __init__(self, param_file, dendrite_flag=True, ex_flag=False, ex_rand=0):
		"""
		This method will be executed when you run
		>>> mycell = RGC_Neuron()
		"""
		self.dflag = dendrite_flag
		self.ex_rand = ex_rand
		if ex_rand !=0: ex_flag = True
#Load parameters
		self.params = read_param_file(param_file)	
		self.seg_len = self.params[-1][-1] #can't handle different length segments yet
#make dict of parts
		self.nams = ['soma','ah','narrow','distal'] if self.dflag else ['soma']
		self.parts=dict(zip(self.nams,[neuron.h.Section() for i in range(len(self.nams))]))
		self.params=dict(zip(self.nams,self.params))

#define geometry
		self.set_geometry()

# passive properties
		self.gp = 0.00005 # [S/cm^2]
		self.E = -65 # []
		self.Ra = 110
# some active properties
#		self.depth_cad=3
#		self.taur_cad=10
# set environment parameters
#		h('celsius = 22')

# initialize parameters
		self.set_passive_parameters()
		self.set_active_parameters()
		if ex_flag: self.set_ex()

	def set_geometry(self):
		for i in self.nams:
			[self.parts[i].L,self.parts[i].diam,self.parts[i].nseg]=list((self.params[i][:3]))
		if self.dflag:
			self.parts['soma'].connect(self.parts['ah'],0,0)
			self.parts['ah'].connect(self.parts['narrow'],0,0)
			self.parts['narrow'].connect(self.parts['distal'],0,0)

	def set_passive_parameters(self):
		for sec in neuron.h.allsec():
			sec.Ra = self.Ra
			sec.insert("pas")
			for seg in sec:
				seg.pas.g = self.gp
				seg.pas.e = self.E

	def set_active_parameters(self):
#		for sec in neuron.h.allsec():
		for i in self.nams: 
			self.parts[i].insert("hh")

	def set_ex(self):
		for sec in neuron.h.allsec():
			sec.insert('extracellular')
			for i in range(2):
				sec.xg[i] = 1e9
				sec.xraxial[i] = 1e9
				sec.xc[i] = 0

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


	def set_IClamp(self,customVector=[],sect = 'soma', region = 0.5):
		"""
		Initializes values for current clamp.

		If custom vector is played into amplitude, default or given is overriden
		Custom vector is input as 2-element list with first element being delay, and scond being stimulus
		"""
		if sect not in self.cell.nams: 
			print "The section you specified is not in the RGC"
			sys.exit(1)
		stim = neuron.h.IClamp(self.cell.parts[sect](region))

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
	
	def set_exstim(self,customVector=[],x_dist=0,y_dist=10,rho = 34.5,just_randflg = False):
		"""
		Initializes values for extracellular field.

		Default is square wave as specified in __init__. 
		Assumes linear resistance between point electrode and membrane.
		First section of soma is assumed to be at (x_dist,y_dist0) = (0,0)
		"""
		if not customVector: 
			[self.t_ext,self.i_ext] = make_square(self.delay,self.amp,self.dur,self.sim_time,self.dt)
		else:
			self.t_ext=neuron.h.Vector(customVector[0])
			self.i_ext=neuron.h.Vector(customVector[1])
#Calculate voltage at each compartment
		self.v_ext = []
		y_seg = 0
		x_seg = self.cell.seg_len/2.0 #take midpoint of first segment as its pos
		#add random extracellular vector 

		for i in self.cell.nams:
			for seg in self.cell.parts[i]:
				v_calculated = calc_v(rho,self.i_ext,x_dist,y_dist,x_seg,y_seg)
				if not self.cell.ex_rand: 
					self.v_ext.append(v_calculated)
				elif just_randflg:
					self.v_ext.append(neuron.h.Vector(self.cell.nchan_den[i]*np.random.normal(0,self.cell.ex_rand,self.sim_time/self.dt)))
				else:
					self.v_ext.append(neuron.h.Vector(np.array(v_calculated) + self.cell.nchan_den[i]*np.random.normal(0,self.cell.ex_rand,self.sim_time/self.dt)))

				x_seg+=self.cell.seg_len	
				self.v_ext[-1].play(seg._ref_e_extracellular, self.t_ext)


	def show(self,titl="Voltage vs. Time",showAx = 1, showStim=False):
		if self.go_already:
			plt.figure()
			x = np.array(self.rec_t)
			y = np.array(self.rec_v)
			#z = np.array(self.rec_vax)
			if showAx != 2:
				plt.plot(x, y)
			if showAx == 1 or showAx == 2: 
				z = np.array(self.rec_vax)
				plt.plot(x, z)
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
    
	def set_recording(self,spaceflg):
# Record Time
		self.rec_t = neuron.h.Vector()
		self.rec_t.record(neuron.h._ref_t)
# Record Voltage
		if self.cell.dflag: 
			if not spaceflg:
				self.rec_v = neuron.h.Vector()
				self.rec_v.record(self.cell.parts['soma'](0.5)._ref_v)
				self.rec_vax = neuron.h.Vector()
				self.rec_vax.record(self.cell.parts['distal'](0.6)._ref_v)
			else:
				self.rec_v = neuron.h.Vector()
				self.rec_v.record(self.cell.parts['soma'](0.5)._ref_v)
				self.rec_vax = [neuron.h.Vector() for i in range(self.cell.parts['distal'].nseg)] #array of vectors, one for each segment in section
				for iseg,seg in enumerate(self.cell.parts['distal']):
					self.rec_vax[iseg].record(seg._ref_v)
		else:

			self.rec_v = neuron.h.Vector()
			self.rec_v.record(self.cell.parts['soma'](0.5)._ref_v)

	def get_recording(self):
		time = np.array(self.rec_t)
		voltage = np.array(self.rec_v)
		if self.cell.dflag: 
			voltage_axon = np.array(self.rec_vax)
			return time, voltage, voltage_axon
		else: 
			return time, voltage

	def go(self, sim_time=None, spaceflg=False):
		self.set_recording(spaceflg)
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

#Random functions--------------------------------------------------

#Waveform making
def make_triphasic(delay,dur,maxAmp,total,dt):
	t = np.arange(total/dt)*dt
	#part=np.floor((dur/dt)/3.0)
	part=np.ceil((dur/dt)/3.0)
	rem = (np.ceil(dur/dt))%3.0
	unit=np.ones(part)
	tp = np.concatenate((unit*maxAmp*(2.0/3),unit*-maxAmp,unit*maxAmp*(1.0/3),np.zeros(rem))) #add remainder
	return t,np.concatenate((np.zeros(delay/dt),tp,np.zeros((total/dt-delay/dt-dur/dt))))

def make_dpp(delay,dur1,dur2,dpAmp,stimAmp,total,dt):
	t = np.arange(total/dt)*dt
	unit1 = np.ones(np.round(dur1/dt))
	unit2 = np.ones(dur2/dt)
	dpp = np.concatenate((unit1*dpAmp, unit2*stimAmp))
	return t, np.concatenate((np.zeros(delay/dt), dpp, np.zeros((total/dt-delay/dt-dur1/dt-dur2/dt))))

def make_dppbal(delay,dur1,dur2,dpAmp,stimAmp,total,dt):
	t = np.arange(total/dt)*dt
	unit1 = np.ones(np.round(dur1/dt))
	unit2 = np.ones(dur2/dt)
	dpp = np.concatenate((unit1*dpAmp, unit2*stimAmp, -unit2*((len(unit1)*dpAmp + len(unit2)*stimAmp)/len(unit2))))
	return t, np.concatenate((np.zeros(delay/dt), dpp, np.zeros((total/dt-delay/dt-dur1/dt-2*dur2/dt))))

def make_dppmod(delay,dur1,dur2,dpAmp,stimAmp,total,dt):
	t = np.arange(total/dt)*dt
	unit1 = np.ones(dur1/dt)
	unit2 = np.ones(dur2/dt)
	dpp = np.concatenate((unit1*dpAmp,unit2*0,unit2*stimAmp))
	return t, np.concatenate((np.zeros(delay/dt), dpp, np.zeros((total-delay-dur1-dur2-dur2)/dt)))

def make_dppmodbal(delay,dur1,dur2,dpAmp,stimAmp,total,dt):
	t = np.arange(total/dt)*dt
	unit1 = np.ones(dur1/dt)
	unit2 = np.ones(dur2/dt)
	dpp = np.concatenate((unit1*dpAmp,unit2*0,unit2*stimAmp,-unit2*((len(unit1)*dpAmp + len(unit2)*stimAmp)/len(unit2))))
	return t, np.concatenate((np.zeros(delay/dt), dpp, np.zeros((total-delay-dur1-dur2-dur2-dur2)/dt)))

def make_square(delay,dur,amp,total,dt):
	t = np.arange(total/dt)*dt
	i = np.concatenate((np.zeros(delay/dt),np.ones(np.round(dur/dt))*amp,np.zeros((total/dt-delay/dt-dur/dt))))
	return t,i

#Others
def read_param_file(param_file):
	"""
	Reads a CSV formatted parameters file and turns it into an array
	"""
	grot=[]
	with open(param_file,'rb') as f:
		wri=csv.reader(f,delimiter=',')
		for row in wri:
			grot.append(row)
	grot=np.transpose(np.array([i[1:] for i in grot[1:]]))
	params=[[] for i in range(grot.shape[0])] 
	for i in range(grot.shape[0]):
		for j in range(grot.shape[1]):
			try: 
				params[i].append(int(grot[i,j]))
			except ValueError:
				params[i].append(float(grot[i,j]))
	return params

def calc_v(rho,i,x,y,x0,y0):
	dist = np.sqrt((x - x0)**2 + (y - y0)**2)
	return neuron.h.Vector(np.array(i) * rho / (4 * np.pi * dist ) * .01) #in mV

