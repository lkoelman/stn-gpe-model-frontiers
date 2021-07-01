Readme file

The present files complement the paper:

Balbi P, Martinoia S and Massobrio P (2015) Axon-somatic
back-propagation in detailed models of spinal alpha
motoneurons. Front. Comput. Neurosci. 9:15. doi:
10.3389/fncom.2015.00015

Pietro Balbi, december 2014

A morphologically detailed model of cat spinal alpha-motoneuron.  The
present implementation takes advantage of the somato-dendritic
morphologically detailed 3D reconstructions available on-line at
NeuroMorpho.org.  The present simulation also extensively adopts, with
few changes:

a) the channels mechanisms by "Powers RK, ElBasiouny SM, Rymer WZ,
Heckman CJ. Contribution of intrinsic properties and synaptic inputs
to motoneuron discharge patterns: a simulation study. J Neurophysiol
(2012) 107: 808-823";

b) a model of myelinated motor axon by "McIntyre CC, Richardson AG,
Grill WM. Modeling the excitability of Mammalian nerve fibers:
influence of afterpotentials on the recovery cycle. J Neurophysiol
(2002) 87:995-1006".  Both these models are available at ModelDB,
accession numbers 143671 and 3810, respectively.

'1_mosinit.hoc' initializes the simulation and displays the panel for
choosing amongst 14 detailed motoneurons. When a choice is made, the
corresponding 3D morphologically detailed reconstructions are
implemented into the model. The file also loads a previously saved
steady-state (300 ms without stimuli, to stabilize the resting
potential across the cell membrane). If a steady-state is not desired,
please set flag_svstate=1 from the interpreter.

'2_complete_cell.hoc' loads the 3_ to 7_ files in sequence.

'3_ins_ch.hoc' inserts the ionic channels at the appropriate locations
(see paper for details) onto soma and dendrites, then sets the
relative conductances.

'4_AHIS.hoc' creates an axon hillock and an axonal initial segment,
joins them to the soma, inserts the ionic channel and the relative
conductances onto them.

'5_axon.hoc' creates the myelinated axon with the corresponding
geometric and biophysical features, and attaches it to the distal end
of the initial segment.

'6_ax_term.hoc' creates the unmyelinated ending of the axon.

'7_morphometry' displays on the default output some morphometric
values of the model (in Italian).

The directory 'channels' contains the .mod files to instantiate the
mechanisms of passive and active channels.

The directory 'cat_spinal_mn' contains the 3D data imported from
NeuroMorho.org.

The directory 'States' contains the saved steady-states of each neuron.

Use of the model: load the mosinit.hoc file to start (in the channels
folder).

To replicate the simulations in the paper, after choosing the
preferred motoneuron a stimulating electrode should be inserted in a
distal node (node[15], for example), and the recording ones both in a
proximal node (only avoid node[0], because it is deleted during the
neuron 3D construction), the axonal initial segment and the soma, thus
detecting the voltage attenuation of the antidromic wave.  Changing
the geometric and sodium conductance parameters of the initial segment
affects the voltage attenuation at a point when the antidromic spike
invasion of the soma can happen.  It is convenient to change the
geometric and conductance parameters from the .hoc file (4_AHIS.hoc),
to be sure that the axon hillock parameters consistently change.  To
speed up the simulations, the VariableStepControl can be set on.  The
model also carries the extracellular mechanism.

The model only represents a first, coarse approximation of a spinal
alpha-motoneuron. In particular, the electrophysiological behaviour of
the dendritic part has not been extensively tested.

P.S.: apologize for multiple untranslated comments in Italian inside
the .hoc files.
