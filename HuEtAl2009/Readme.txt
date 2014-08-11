This directory contains the Neuron source code for cortical Layer 5
pyramidal cell model and experiments employed in:

Distinct Contributions of Na(V)1.6 and Na(V)1.2 in Action Potential
Initiation and Backpropagation Wenqin Hu, Cuiping Tian, Tun Li, Mingpo
Yang, Han Hou & Yousheng Shu (2009) Nat Neurosci 12(8): 996-1002.


Part of model is based on:

Mainen, Z. F. and Sejnowski, T. J.  Nature 382: 363-6 (1996)
Yu, Y., Shu, Y., et al.  J Neurosci 28: 7260-72 (2008)
Shu, Y., Hasenstaub, A., et al.  Nature 441: 761-5. (2006)


===============================================

BRIEF OVERVIEW OF THE CONTENTS


Three different but related models are involved in this package:

1). A realistic model of Layer 5 pyramidal cell with sophisticatedly
described voltage-dependent sodium channels at the axon initial
segment. Either action potential's initiation site (figure not shown
in the aforementioned paper, see its main text) or backpropagation
failure threshold (Supplementary Fig.4 and Fig.8) can be tested
here. This section produced Fig.5d and 5e.

2). A model of "uniform axon" for addressing the relative contribution
of Nav1.6 and Nav1.2 in action potential's initiation. This produced
Fig.5b and 5c.

3). A single-compartment model for examining the activation and
inactivation properties of these two channel subtypes. This produced
Supplementary Fig.3.

All above are accessible through running "init.hoc". 

See "Overview.png" for detailed organization of files.

Note: The names "naaxon" and "nasoma" appeared in Online Methods have
been replaced by "na16" and "na12" here in this model for
clarity. They stand for low- and high-threshold Nav channels at the
axon initial segment, respectively, whereas those located in
somatodendritic region are termed "na", without any suffix.

===============================================

USAGE

1. If needed, generate your own nrnmech.dll with *.mod files in the
directory MECHANISM, and replace any old dll files with the new one.

2. "init.hoc" can be directly executed from Neuron.exe.  (double click
if Neuron has been installed in your computer)

3. For a better experience of editing and debugging, it is strongly
recommended to use the open source software PSPad editor to open the
project file "PSPad_Project.ppr".  (see
http://www.neuron.yale.edu/phpBB/viewtopic.php?f=5&t=1219)

20120312 a mosinit.hoc file was added for auto-launch and the solve
method was updated from euler to derivimplicit in cad.mod.  See
http://www.neuron.yale.edu/phpBB/viewtopic.php?f=28&t=592
20140119 init.hoc and experiment/UniformAxonIClamp_ChangeDens.hoc
updated xopens for unix case sensitive matches of filenames
experiment/UniformAxon_main.hoc and lib/U_dvdt.hoc respectively.
