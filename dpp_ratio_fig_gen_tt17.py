import pickle

#load data files
with open('outsoma-tt17.p', 'rb') as f: outsoma = pickle.load(f)
with open('outaxon-tt17.p', 'rb') as f: outaxon = pickle.load(f)

#make lists

dpplist = np.linspace(.01e3,3e3,20)
distlist = np.linspace(1,100,20)

#original plots
plot_3d(outsoma/1000, dpplist, distlist, 'Threshold vs. DPP (Balanced) Distance and Duration (Soma)', 'DPP Amplitude [nA]', 'Electrode Distance [microns]')
plot_3d(outaxon/1000, dpplist, distlist, 'Threshold vs. DPP (Balanced) Distance and Duration (Axon)', 'DPP Amplitude [nA]', 'Electrode Distance [microns]')
plot_3d(abs(outaxon/outsoma), dpplist, distlist, 'Ratio of Threshold vs. DPP (Axon vs. Soma)', 'DPP Amplitude [nA]', 'Electrode Distance [microns]', 'Stimulation Threshold Ratio', myvmax = 50)

#percentage plots
