import pyneurlib as pdl

#globals
fdict = {'sq':pdl.make_square, 'tp': pdl.make_triphasic, 'dpp':pdl.make_dpp, 'dppmod':pdl.make_dppmod, 'dppbal':pdl.make_dppbal}

def get_max(cell, typmax, amp, x,y,stim_params, dppamp = 0): 
	if dppamp == 0: 
		[delay, dur, simtime, dt, myrho] = stim_params
		[t,i]=fdict[typmax](delay,dur,amp,simtime,dt)
	else:
		[delay, dur1, dur2, simtime, dt, myrho] = stim_params
		[t,i]=fdict[typmax](delay,dur1,dur2,dppamp,amp,simtime,dt)
	sim = pdl.Simulation(cell,dt,sim_time = simtime)
	sim.set_exstim([t,i],x_dist=x,y_dist = y,rho = myrho)
	sim.go()
	_,som_rec,ax_rec = sim.get_recording()
	return som_rec,ax_rec

def get_thresh(cell, typ, init_amp, init_step, x, y, stim_params, dpp_amp = 0, spike_thr = 0, amp_thr = .025e3, cnt_lim = 15):
	amp = init_amp if type(init_amp) != list else init_amp[1]
	step = init_step
	downflg = 0
	cnt = 0
	while step >= amp_thr: 
		if dpp_amp == 0:
			[msom_rec, max_rec] = get_max(cell, typ, amp, x,y, stim_params)
		else: 
			[msom_rec, max_rec] = get_max(cell, typ, amp, x,y, stim_params,dppamp = dpp_amp)
		if (max(msom_rec) > spike_thr or max(max_rec) > spike_thr): 
			amp -= step
			downflg = 1
		elif downflg:
			step *= 0.5
			amp += step
		else: amp +=step #don't reduce ampstep if the thing has never gone down
		cnt += 1
		if cnt > cnt_lim: break #make sure things don't go overboard
	return amp
