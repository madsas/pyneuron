import pyneurlib as pdl
import matplotlib.pyplot as plt

plt.figure()
t,x = pdl.make_triphasic(100, 150, 1.5, 450, .1)
plt.plot(t,x-2)
t,x = pdl.make_dppbal(100, 200, 50, 100, 1, 1.5, 500, .1)
plt.plot(t,x+2)
plt.ylabel('Amplitude (nA)')
plt.xlabel('Time (us)')
plt.title('Triphasic and Depolarizing Pre-Pulse Waveforms')
plt.show()

savefig('waveform_2.png')
savefig('waveform_2.svg')
