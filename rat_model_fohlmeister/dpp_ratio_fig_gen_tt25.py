import pickle
#reference for cmaps: http://matplotlib.org/examples/color/colormaps_reference.html

#definitions
def plot_3d(data, list1, list2, title, list1name, list2name, cbarlabel = 'Stimulation Threshold [uA]', myvmin = 0, myvmax = 1):
	fig, ax = plt.subplots()
	#im = ax.imshow(np.flipud(data), cmap=plt.cm.RdBu, vmin=myvmin, vmax=myvmax, extent=[min(list1),max(list1),min(list2),max(list2)], aspect = 'auto')
	im = ax.imshow(np.flipud(data), cmap='gist_rainbow', vmin=myvmin, vmax=myvmax, extent=[min(list1),max(list1),min(list2),max(list2)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	cb.set_label(cbarlabel, rotation = 270)
	plt.xlabel(list1name)
	plt.ylabel(list2name)
	plt.title(title)
	plt.show()

#load data files
with open('outsoma-tt25.p', 'rb') as f: outsoma = pickle.load(f)
with open('outaxon-tt25.p', 'rb') as f: outaxon = pickle.load(f)

#make lists


dppperclist = np.linspace(.5,1,20)
distlist = np.linspace(1,100,20)

cdict = {'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)), 'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)), 'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))}

#plots
#plot_3d(outsoma/1000, dppperclist, distlist, 'Threshold vs. DPP (Balanced) Electrode Distance and Amplitude (Soma)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]')
#plot_3d(outaxon/1000, dppperclist, distlist, 'Threshold vs. DPP (Balanced) Electrode Distance and Amplitude (Axon)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]')
plot_3d(abs(outaxon/outsoma), dppperclist, distlist, 'Ratio of Threshold vs. DPP (Axon vs. Soma)', 'DPP Percentage of Stimulus', 'Electrode Distance [microns]', 'Stimulation Threshold Ratio', myvmax = (outaxon/outsoma).max())

