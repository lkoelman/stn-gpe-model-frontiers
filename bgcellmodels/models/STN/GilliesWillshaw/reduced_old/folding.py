"""
Utility functions for folding / collapsing operations.

@author Lucas Koelman

@date	28-08-2016

"""

# Own modules
import redutils

# logging of DEBUG/INFO/WARNING messages
import logging
logging.basicConfig(format='%(levelname)s:%(message)s @%(filename)s:%(lineno)s', level=logging.DEBUG)
logname = "folding" # __name__
logger = logging.getLogger(logname) # create logger for this module


def find_collapsable(allsecrefs, i_pass, Y_criterion='highest_level', zips_per_pass=1e9):
	"""
	Find branch points with child branches that can be collapsed.

	@return		list of SectionRef to the branchpoint-Sections with
				collapsable children
	"""

	if Y_criterion=='max_electrotonic':
		
		# Find Y with longest electrotonic length to first Y among its children.
		# This corresponds to  Y with child section that has longest L/lambda.
		# (Collapsing this Y will eliminate most compartments.)

		# Get minimum distance to next branch point (determines collapsable length)
		for secref in allsecrefs:
			
			child_secs = secref.sec.children()
			secref.collapsable_L_elec = 0.0
			
			if any(child_secs):
				min_child_L = min(sec.L for sec in child_secs)
			
			for chi_sec in child_secs:
				# Get segment that is that distance away from child
				furthest_seg = chi_sec(min_child_L/chi_sec.L)
				furthest_L_elec = redutils.seg_path_L_elec(furthest_seg, f_lambda, gleak_name)
				secref.collapsable_L_elec += (secref.pathLelec1 - furthest_L_elec)

		# Get maximum collapsable length
		max_L_collapsable = max(ref.collapsable_L_elec for ref in allsecrefs)

		# Find all branch points that have +/- same collapsable length
		low, high = 0.95*max_L_collapsable, 1.05*max_L_collapsable
		candidate_Y_secs = [ref for ref in allsecrefs if (any(ref.sec.children())) and (
															low<ref.collapsable_L_elec<high)]
		
		# Off these, pick the one with highest level, and select all branch points at this level
		target_level = max(ref.level for ref in candidate_Y_secs)
		target_Y_secs = [ref for ref in candidate_Y_secs if ref.level==target_level]

		# Report findings
		logger.debug("The maximal collapsable electrotonic length is %f", max_L_collapsable)
		logger.debug("Found %i sections with collapsable length within 5%% of this value", 
						len(candidate_Y_secs))
		logger.debug(("The highest level of a parent node with this value is %i, "
						"and the number of nodes at this level with the same value is %i"), 
						target_level, len(target_Y_secs))


	elif Y_criterion=='highest_level':

		# Simplify find all branch points at highest level
		max_level = max(ref.level for ref in allsecrefs)
		target_Y_secs = [ref for ref in allsecrefs if (ref.level==max_level-1) and (
						 i_pass+1 <= ref.max_passes) and (len(ref.sec.children()) >= 2)]

	else:
		raise Exception("Unknow Y-section selection criterion '{}'".format(Y_criterion))


	# Prune branchpoints: can't collapse more branchpoints than given maximum
	n_Y = len(target_Y_secs)
	target_Y_secs = [ref for i,ref in enumerate(target_Y_secs) if i<zips_per_pass]

	logger.debug("Found {0} Y-sections that meet selection criterion. Keeping {1}/{2}\n\n".format(
					n_Y, min(n_Y, zips_per_pass), n_Y))
	

	# Return branch points identified for collapsing
	return target_Y_secs