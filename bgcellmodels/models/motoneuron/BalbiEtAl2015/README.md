# Changes to source code

Changes with respect to original source code:

- file `kca2.mod`: change line `USEION caL READ icaL` to `USEION caL READ icaL VALENCE 2` as in `L_Ca.mod` to prevent NEURON 7.5 from throwing error


# Model structure

- Following sections areavailable on Hoc interpreter:
	+ `soma[N] <Section>` somatic sections
	+ `dend[M] <Section>` dendritic sections
	+ `AH[1] <Section>` axon hillock
	+ `IS[1] <Section>` axon initial segment
	+ `node[axonnodes]` <Section>
	+ `MYSA[paranodes1]` <Section>
	+ `FLUT[paranodes2]` <Section>
	+ `STIN[axoninter]` <Section>

- Following `SectionList` objects are made to group sections for inserting mechanisms and setting parameters
	+ 

# Files overview

`'1_mosinit.hoc'` initializes the simulation and displays the panel for
choosing amongst 14 detailed motoneurons. When a choice is made, the
corresponding 3D morphologically detailed reconstructions are
implemented into the model. The file also loads a previously saved
steady-state (300 ms without stimuli, to stabilize the resting
potential across the cell membrane). If a steady-state is not desired,
please set flag_svstate=1 from the interpreter.

`'2_complete_cell.hoc'` loads the 3_ to 7_ files in sequence.

`'3_ins_ch.hoc'` inserts the ionic channels at the appropriate locations
(see paper for details) onto soma and dendrites, then sets the
relative conductances.

`'4_AHIS.hoc' `creates an axon hillock and an axonal initial segment,
joins them to the soma, inserts the ionic channel and the relative
conductances onto them.

`'5_axon.hoc'` creates the myelinated axon with the corresponding
geometric and biophysical features, and attaches it to the distal end
of the initial segment.

`'6_ax_term.hoc'` creates the unmyelinated ending of the axon.

`'7_morphometry'` displays on the default output some morphometric
values of the model (in Italian).

`channels/` The directory 'channels' contains the .mod files to instantiate the
mechanisms of passive and active channels.

`cat_spinal_mn/` The directory 'cat_spinal_mn' contains the 3D data imported from
NeuroMorho.org.

`states/` The directory 'States' contains the saved steady-states of each neuron.