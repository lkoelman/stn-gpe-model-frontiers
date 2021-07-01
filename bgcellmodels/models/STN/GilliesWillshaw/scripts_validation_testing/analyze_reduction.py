"""
Analysis utilities for Marasco reduction

@author Lucas Koelman
@date	20-12-2016
"""

# Python modules
import math
import matplotlib.pyplot as plt

from neuron import h

# Own modules
from bgcellmodels.common import electrotonic
from reducemodel.redutils import ExtSecRef, getsecref

import gillies_model
from reducemodel import reduce_cell


def plot_tree_ppty(secref, allsecrefs, propname, secfilter, labelfunc, y_range=None, fig=None, ax=None):
	"""
	Descend the given dendritic tree start from the given root section
	and plot the given propertie.


	@param	secref		SectionRef to root section

	@param	allsecrefs	list(SectionRef) with references to all sections
						in tree

	@param	propname	segment property to plot

	@param	secfilter	filter function applied to each SectionRef
	"""
	if secref is None:
		return
	first_call = fig is None
	if first_call:
		fig = plt.figure()
	if y_range is None:
		y_range = [float('inf'), float('-inf')]

	# Plot current node
	if secfilter(secref):
		if not ax:
			nax = len(fig.axes)
			for i, ax in enumerate(fig.axes):
				ax.change_geometry(nax+1, 1, i+1) # change grid and position in grid
			ax = fig.add_subplot(nax+1, 1, nax+1)

		# get data to plot
		sec = secref.sec
		xg = [(seg.x, getattr(seg, propname)) for seg in sec]
		xvals, gvals = zip(*xg)
		if min(gvals) < y_range[0]:
			y_range[0] = min(gvals)
		if max(gvals) > y_range[1]:
			y_range[1] = max(gvals)

		# plot it
		ax.plot(xvals, gvals, 'o')
		ax.plot(xvals, gvals, '-')
		# ax.set_xlabel('x')
		ax.set_ylabel(labelfunc(secref))

		if y_range is not None:
			ax.set_ylim(y_range)

	# plot children
	for childsec in secref.sec.children():
		childref = getsecref(childsec, allsecrefs)
		plot_tree_ppty(childref, allsecrefs, propname, secfilter, labelfunc, y_range, fig)

	if first_call:
		plt.suptitle(propname)
		y_span = y_range[1]-y_range[0]
		for ax in fig.axes:
			ax.set_ylim((y_range[0]-0.1*y_span, y_range[1]+0.1*y_span))
		plt.show(block=False)
	return fig


def calc_subtree_Rin(root_secs):
	"""
	Compare model properties
	"""

	# Compare input resistance of large/left tree 
	Rin_DC = [redutils.inputresistance_tree(sec, 0., 'gpas_STh') for sec in root_secs]
	Rin_AC = [redutils.inputresistance_tree(sec, 100., 'gpas_STh') for sec in root_secs]
	
	for i, root_sec in enumerate(root_secs):
		print("\n=== INPUT RESISTANCE ===\
			\nSubtree of {}:\
			\nRin_AC={:.3f} \tRin_DC={:.3f}".format(root_sec, Rin_AC[i], Rin_DC[i]))



def inspect_passive_electrotonic_structure(reduced=True):
	"""
	Inspect passive electrotonic strucutre in NEURON GUI.

	After running this function, plot attenuation properties via Tools > Electrotonic Analysis
	"""
	if reduced:
		soma_refs, dend_refs = reduce_cell.fold_gillies_marasco(False)
		allsecrefs = soma_refs + dend_refs

	else:
		somaref, dendLrefs, dendRrefs = gillies_model.get_stn_refs()
		allsecrefs = [somaref] + dendLrefs + dendRrefs

	# Disable all active conductances (can also use 'uninsert' on all sections)
	gbar_active = [g for g in gillies_model.gillies_glist if (g != gillies_model.gleak_name)]

	sec_modified = 0
	seg_modified = 0

	# for sec in h.allsec():
	for secref in allsecrefs:
		sec = secref.sec
		secref.orig_range_props = redutils.get_range_props(secref, gbar_active)

		for seg in sec:
			for gbar in gbar_active:
				setattr(seg, gbar, 0.0)
			seg_modified += 1
		
		sec_modified += 1

	print("Set gbar to zero in {} sections and {} segments".format(sec_modified, seg_modified))
	print("Please plot electrotonic structure via Tools > Impedance ")
	from neuron import gui # opens GUI in another thread


def inspect_diameter(mech_param_names=None):
	"""
	Inspect diameter using NEURON GUI

	USAGE

	- after calling this function, open NEURON Main Menu > ModelView > plot mechanism 'd_dum'
	"""

	if mech_param_names is None:

		# THIS WORKS BUT MAKES MODELVIEW CRASH
		# Create a dummy mechanism for tracking diameter
		# h("""
		# begintemplate Dum
		# public d
		# d = 0
		# endtemplate Dum

		# make_mechanism(\"dum\", \"Dum\", \"d\")
		# """)

		mech_param_names = {'hh', 'gl'}

	mech = mech_param_names.keys()[0]
	param = mech_param_names[mech]
	pname = '{}_{}'.format(param, mech)

	# In each segment: assign diameter to dummy parameter
	for sec in h.allsec():
		sec.insert(mech)
		for seg in sec:
			setattr(seg, pname, seg.diam)

	# Plot dummy parameter in GUI
	print("Please plot parameter '{}' via Tools > ModelView ".format(pname))
	from neuron import gui

if __name__ == '__main__':
	inspect_passive_electrotonic_structure()
