import pyneurlib_temp as pdl

#globals
fdict = {'sq':pdl.make_square, 'tp': pdl.make_triphasic, 'dpp':pdl.make_dpp, 'dppmod':pdl.make_dppmod, 'dppbal':pdl.make_dppbal, 'dppmodbal':pdl.make_dppmodbal}

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

def get_thresh(cell, typ, init_amp, init_step, x, y, stim_params, dpp_amp = 0, spike_thr = 0, amp_step_thr = .025e3, cnt_lim = 15, cathodic = 1):
	if cathodic: #default is that you want negative threshold
		amp = init_amp if type(init_amp) != list else init_amp[1]
		step = init_step
		#downflg = 0
		upflg = 0
		cnt = 0
		while step >= amp_step_thr: 
			[msom_rec, max_rec] = get_max(cell, typ, amp, x,y, stim_params,dppamp = dpp_amp)
			if (max(msom_rec) > spike_thr or max(max_rec) > spike_thr): 
				#amp -= step
				amp += step
				#downflg = 1
				upflg = 1
			#elif downflg:
			elif upflg:
				step *= 0.5
				amp -= step
			#else: amp +=step #don't reduce ampstep if the thing has never gone down
			else: amp -=step #don't reduce ampstep if the thing has never gone down
			cnt += 1
			if cnt > cnt_lim: #make sure things don't go overboard
				#now go the other direction in a simple way
				orig_amp = amp
				amp = init_amp
				step = init_step
				cnt = 0
				while step >= amp_step_thr: 
					#print 'hi'
					[msom_rec, max_rec] = get_max(cell, typ, amp, x,y, stim_params,dppamp = dpp_amp)
					if (max(msom_rec) > spike_thr or max(max_rec) > spike_thr): break
					#else: amp -= step
					else: amp += step
					if cnt > cnt_lim: return orig_amp #return the original super high one
					if amp > 0: return orig_amp
					cnt += 1
				break
			if amp > 0: break #can't have positive stimulus here
		return amp
	else:
		amp = init_amp if type(init_amp) != list else init_amp[1]
		step = init_step
		downflg = 0
		cnt = 0
		while step >= amp_step_thr: 
			[msom_rec, max_rec] = get_max(cell, typ, amp, x,y, stim_params,dppamp = dpp_amp)
			if (max(msom_rec) > spike_thr or max(max_rec) > spike_thr): 
				amp -= step
				downflg = 1
			elif downflg:
				step *= 0.5
				amp += step
			else: amp +=step #don't reduce ampstep if the thing has never gone down
			cnt += 1
			if cnt > cnt_lim: #make sure things don't go overboard
				#now go the other direction in a simple way
				orig_amp = amp
				amp = init_amp
				step = init_step
				cnt = 0
				while step >= amp_thr: 
					[msom_rec, max_rec] = get_max(cell, typ, amp, x,y, stim_params,dppamp = dpp_amp)
					if (max(msom_rec) > spike_thr or max(max_rec) > spike_thr): break
					else: amp -= step
					if cnt > cnt_lim: return orig_amp #return the original super high one
					if amp < 0: return orig_amp
					cnt += 1
				break
			if amp < 0: break #can't have negative stimulus here 
		return amp
