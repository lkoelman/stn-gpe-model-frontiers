"""
Reduction of Gillies & Willshaw (2006) STN neuron model


@author Lucas Koelman
@date	03-11-2016
@note	must be run from script directory or .hoc files not found

"""


import numpy as np
import matplotlib.pyplot as plt
import math

import neuron
from neuron import h

# Load NEURON function libraries
h.load_file("stdlib.hoc") # Load the standard library
h.load_file("stdrun.hoc") # Load the standard run library

import sys
import os.path
scriptdir, scriptfile = os.path.split(__file__)
modulesbase = os.path.normpath(os.path.join(scriptdir, '..'))
sys.path.append(modulesbase)

# Load NEURON mechanisms
# add this line to nrn/lib/python/neuron/__init__.py/load_mechanisms()
# from sys import platform as osplatform
# if osplatform == 'win32':
# 	lib_path = os.path.join(path, 'nrnmech.dll')
NRN_MECH_PATH = os.path.normpath(os.path.join(scriptdir, 'nrn_mechs'))
neuron.load_mechanisms(NRN_MECH_PATH)

from interpolation import loadgstruct, gillies_gstructs

# Global variables
soma = None

class Tree:
	""" A binary tree node """
	def __init__(self, payload, left=None, right=None):
		self.payload = payload
		self.left  = left
		self.right = right

	def __str__(self):
		return str(self.payload)

	def setall(self, key, value):
		self.payload[key] = value
		if self.left is not None:
			self.left.setall(key, value)
		if self.right is not None:
			self.right.setall(key, value)

def gillies_treestruct():
	""" Return the two dendritic tree structures from the paper
		
		The 'index' refers to the branch/section indices in the
		.dat files that come with the paper. To get a reference to
		the corresponding section, use index-1 in the Vectors stored
		in the SThcell object.

		e.g. in tree1: branch5 = h.SThcell[0].dend1[4]
	"""
	# Tree structure (map structure to indices in tree1-nom.dat)
	## Fig. 1, right dendritic tree
	dend1tree = Tree({'index':1}, 
					Tree({'index':2}, 
						Tree({'index':4}, Tree({'index':6}), Tree({'index':7})), 
						Tree({'index':5})), 
					Tree({'index':3}, 
						Tree({'index':8}), 
						Tree({'index':9}, Tree({'index':10}), Tree({'index':11}))))

	## Fig. 1, left dendritic tree
	dend0upper = Tree({'index':2}, 
					Tree({'index':4}, 
						Tree({'index':6}, Tree({'index':8}), Tree({'index':9})), 
						Tree({'index':7})), 
					Tree({'index':5}, 
						Tree({'index':10}), 
						Tree({'index':11}, Tree({'index':12}), Tree({'index':13}))))

	dend0lower = Tree({'index':3}, 
					Tree({'index':14}, 
						Tree({'index':16}, Tree({'index':18}), Tree({'index':19})), 
						Tree({'index':17})), 
					Tree({'index':15}, 
						Tree({'index':20}), 
						Tree({'index':21}, Tree({'index':22}), Tree({'index':23}))))

	dend0tree = Tree({'index':1}, dend0upper, dend0lower)
	return dend0tree, dend1tree




def treechannelstruct():
	""" Return tree structure containing channel distributions.

		Each tree node's payload will contain an entry gname
		and xgname with the conductance values for each x value
	"""

	# Load conductance matrices
	gstructs = gillies_gstructs()

	# Initialize tree structure
	dend0tree, dend1tree = gillies_treestruct()

	def filltreegstruct(tree, dendidx):
		if tree is None: return

		# Fill this branch/section
		branchidx = tree.payload['index']
		for gname, gmat in gstructs:
			branchrows = (gmat['dendidx']==dendidx) & (gmat['branchidx']==branchidx-1)
			tree.payload['x'+gname] = gmat[branchrows]['x']
			tree.payload[gname] = gmat[branchrows]['g']

		# Fill left and right children
		filltreegstruct(tree.left, dendidx)
		filltreegstruct(tree.right, dendidx)

	# Fill the two dendritic trees
	filltreegstruct(dend0tree, 0)
	filltreegstruct(dend1tree, 1)
	return dend0tree, dend1tree


def getbranchproperties(branchidx, dendsecs):
	""" Return branch properties for the given branch index """
	secidx = branchidx-1

	# section properties
	L = dendsecs[secidx].L
	Ra = dendsecs[secidx].Ra

	# segment properties
	diam = dendsecs[secidx](0.5).diam
	cm = dendsecs[secidx](0.5).cm

	# Electrotonic length
	lambda100 = 1e5*np.sqrt(diam/(4*np.pi*100*Ra*cm))
	l_elec = L/lambda100

	return L, Ra, diam, cm, lambda100, l_elec

def printproperties(tree, dendsecs):
	""" Recursive depth first traversal of tree """
	if tree is None: return 0

	# Get branch properties
	branchidx = tree.payload['index']
	L, Ra, diam, cm, lambda100, l_elec = getbranchproperties(branchidx, dendsecs)

	# Print branch properties
	print("Branch index {0}: L={1}; Ra={2}; diam(0.5)={3}; cm(0.5)={4}; lambda100={5}; l_elec={6}".format(
		branchidx, L, Ra, diam, cm, lambda100, l_elec))
	printproperties(tree.left, dendsecs)
	printproperties(tree.right, dendsecs)

def rallreduce(tree, dendsecs):
	""" Reduce thre tree according to Rall algorithm

	@return		the equivalent electrotonic length of the tree 
				(i.e. that should be added to parent tree of this branch)
	@return		diameter^(3/2) of the given tree branch

	The following properties must be met:
		1. Rm and Ra must be the same in parent + children
		2. Same boundary conditions on all branches
		3. Child branches have same electrotonic length
		4. 3/2 diameter rule satisfied betweeen parents - children
	Properties 1 & 2 are satisfied according to model specification, 
	so only 3 & 4 must be checked.

	NOTE: if tree is the parent tree, the equivalent length is the 
	      returned value times its lambda100
	"""
	if tree is None: return 0, 0

	# Get own electrotonic length and diameter
	branchidx = tree.payload['index']
	L, Ra, diam, cm, lambda100, l_elec = getbranchproperties(branchidx, dendsecs)
	selfd32 = diam**(3./2.)

	# Get equivalent electrotonic length of children
	left_elec, leftd32 = rallreduce(tree.left, dendsecs)
	right_elec, rightd32 = rallreduce(tree.right, dendsecs)

	# Check if they are (almost) the same
	if (left_elec==0 and right_elec==0) or (left_elec/right_elec >= 0.95 and left_elec/right_elec <= 1.05):
		print('OK: electrotonic length of child branches is within 5%% of eachother')
	else:
		print('WARNING: electrotonic length of branches to be merged is not close enough!')
	d32ratio = (leftd32+rightd32)/selfd32
	if (leftd32==0 and rightd32==0) or (d32ratio >= 0.95 and d32ratio <= 1.05):
		print('OK: 3/2 rule satisfied within 5%% accuracy')
	else:
		print('WARNING: 3/2 rule not satisfied within 5%% accuracy" ratio={0}'.format(d32ratio))

	return ((left_elec+right_elec)/2.0 + l_elec), selfd32

def gillies_rallreduction():
	""" Measure electrical/geometrical ppties and reduce dendritic
		trees using Rall algorithm """

	# Create cell and three IClamp in soma
	h.xopen("gillies_cell_singleton.hoc")
	somasec = h.SThcell[0].soma
	dend0secs = h.SThcell[0].dend0
	dend1secs = h.SThcell[0].dend1
	dend0tree, dend1tree = gillies_treestruct()
	

	# Check from properties if reducable (checked by hand: YES)
	print('== Checking right dendritic tree (smaller one) ===')
	printproperties(dend1tree, dend1secs)

	# equivalent length
	print('== Reducing right dendritic tree (smaller one) ===')
	L, Ra, diam, cm, lambda_root, l_elec = getbranchproperties(1, dend1secs)
	l_equiv, _ = rallreduce(dend1tree, dend1secs)
	L_equiv = l_equiv * lambda_root
	print('Equivalent length is {0}'.format(L_equiv)) # 549.86

	# LEFT TREE
	print('== Checking left dendritic tree (larger one) ===')
	printproperties(dend0tree, dend0secs)

	print('== Reducing left dendritic tree (larger one) ===')
	L, Ra, diam, cm, lambda_root, l_elec = getbranchproperties(1, dend0secs)
	l_equiv, _ = rallreduce(dend0tree, dend0secs)
	L_equiv = l_equiv * lambda_root
	print('Equivalent length is {0}'.format(L_equiv)) # 703.34

def gillies_twosec_reduction():
	""" Collapse all dendritic trees into one equivalent section

	ALGORITHM:
		- use algorithm Sterrat p. 85-86 to compute equivalent imput resistance
		  of two trees in parallel
		- the equivalent Rin has two free parameters: d and L. Use d of the 
		  largest tree and calculate L from this and the equivalent resistance value.
	"""
	# use min. amount of sections to reproduce behaviour
	pass

def plotconductances(treetopo, treeidx, gstruct, includebranches=None):
	""" Recursively plot conductances in each branch by descending the tree """

	if treetopo is None: return

	# Get branch properties
	branchidx = treetopo.payload['index']

	# Plot own distribution
	if includebranches is None or branchidx in includebranches:

		branchrows = (gstruct['dendidx']==treeidx) & (gstruct['branchidx']==branchidx-1)
		x = gstruct[branchrows]['x']
		g = gstruct[branchrows]['g']

		plt.figure()
		plt.plot(x,g,'o')
		plt.ylim(0, 1e-2)
		plt.xlabel('x (normalized length)')
		plt.ylabel('gbar')
		plt.suptitle('Conductance distribution in branch {}'.format(branchidx))
		plt.show(block=False)

		# print it
		print('## branch {} conductances ##'.format(branchidx))
		print('\t'.join((str(l) for l in x)))
		print('\t'.join((str(gb) for gb in g)))

	plotconductances(treetopo.left, treeidx, gstruct, includebranches)
	plotconductances(treetopo.right, treeidx, gstruct, includebranches)

if __name__ == '__main__':
	plotconductances(gillies_treestruct()[1], 1, loadgstruct('gcaT_CaT'), includebranches=[1,2,5])
	# plotchanneldist(0, 'gcaL_HVA')
	# dend0tree, dend1tree = treechannelstruct()
	# gtstruct = loadgeotopostruct(0)