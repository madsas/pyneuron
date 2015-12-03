import matplotlib.pyplot as plt
import numpy as np

def plot_3d(data, list1, list2, title, list1name, list2name, cbarlabel = 'Stimulation Threshold [uA]', myvmin = 0, myvmax = 1):
	fig, ax = plt.subplots()
	im = ax.imshow(np.flipud(data), cmap=plt.cm.RdBu, vmin=myvmin, vmax=myvmax, extent=[min(list1),max(list1),min(list2),max(list2)], aspect = 'auto')
	im.set_interpolation('none')
	cb = fig.colorbar(im, ax=ax)
	cb.set_label(cbarlabel, rotation = 270)
	plt.xlabel(list1name)
	plt.ylabel(list2name)
	plt.title(title)
	plt.show()

x = np.arange(10)
xl = 'hithere'
plot_3d(x,x,x,xl,xl,xl)
