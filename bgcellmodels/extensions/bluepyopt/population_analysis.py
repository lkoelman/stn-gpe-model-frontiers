"""
Analysis of BluePyOpt optimisation results: population metrics etc.

@author	Lucas Koelman

@date	4/11/2017

@see	based on scripts:
			https://github.com/BlueBrain/BluePyOpt/blob/master/examples/l5pc/l5pc_analysis.py
			https://github.com/BlueBrain/BluePyOpt/blob/master/examples/l5pc/opt_l5pc.py
"""

import pickle

# Scipy
import numpy as np
from matplotlib import pyplot as plt


def plot_log(log, ymax='max'):
	"""
	Plot logbook

	@param	ymax	metric to use for axes y limits
	"""

	fig, axes = plt.subplots(facecolor='white')

	gen_numbers = log.select('gen')
	mean = np.array(log.select('avg'))
	std = np.array(log.select('std'))
	minimum = np.array(log.select('min'))
	maximum = np.array(log.select('max'))

	metrics = {
		'avg': mean,
		'std': std,
		'min': minimum,
		'max': maximum,
	}

	stdminus = mean - std
	stdplus = mean + std
	axes.plot(
		gen_numbers,
		mean,
		color='black',
		linewidth=2,
		label='population average')

	axes.fill_between(
		gen_numbers,
		stdminus,
		stdplus,
		color='lightgray',
		linewidth=2,
		label=r'population standard deviation')

	axes.plot(
		gen_numbers,
		minimum,
		color='green',
		linewidth=2,
		label='population minimum')

	axes.plot(
		gen_numbers,
		maximum,
		color='red',
		linewidth=2,
		label='population maximum')

	axes.set_xlim(min(gen_numbers) - 1, max(gen_numbers) + 1)
	axes.set_xlabel('Generation #')
	axes.set_ylabel('Sum of objectives')
	axes.set_ylim([0, max(metrics[ymax])])
	axes.legend()

	fig.tight_layout()

	return fig, axes


def plot_history(history):
	"""
	Plot the history (geneaology tree) of the individuals

	@param	history		deap.tools.History object: see
						https://github.com/BlueBrain/deap/blob/bbp-specific/deap/tools/support.py#L21
	"""

	import networkx
	import matplotlib.pyplot as plt

	plt.figure()

	graph = networkx.DiGraph(history.genealogy_tree)
	graph = graph.reverse()     # Make the grah top-down
	# colors = [\
	#        toolbox.evaluate(history.genealogy_history[i])[0] for i in graph]
	positions = networkx.graphviz_layout(graph, prog="dot")
	networkx.draw(graph, positions)


def plot_individual_fitness(objectives, box=None):
	"""
	Plot objectives of the cell model

	@param objectives	dict { (str) objective_name: (float) objective_score}

	USAGE:
		objectives = optimisation.evaluator.fitness_calculator.calculate_scores(responses)
		plot_fitness_scores(objectives)
	"""

	import collections
	objectives = collections.OrderedDict(sorted(objectives.items()))

	fig, axes = plt.subplots(facecolor='white')
	
	# left_margin = box['width'] * 0.4
	# right_margin = box['width'] * 0.05
	# top_margin = box['height'] * 0.05
	# bottom_margin = box['height'] * 0.1

	# axes = fig.add_axes(
	# 	(box['left'] + left_margin,
	# 	 box['bottom'] + bottom_margin,
	# 	 box['width'] - left_margin - right_margin,
	# 	 box['height'] - bottom_margin - top_margin))

	ytick_pos = [x + 0.5 for x in range(len(objectives.keys()))]

	axes.barh(ytick_pos,
			  objectives.values(),
			  height=0.5,
			  align='center',
			  color='#779ECB')
	axes.set_yticks(ytick_pos)
	axes.set_yticklabels(objectives.keys(), size='x-small')
	axes.set_ylim(-0.5, len(objectives.values()) + 0.5)
	axes.set_xlabel('Objective value (# std)')
	axes.set_ylabel('Objectives')

	return fig, axes


def plot_individual_params(
		opt,
		ax,
		params,
		marker,
		color,
		markersize=40,
		plot_bounds=False,
		fitness_cut_off=None): 
	'''
	plot the individual parameter values
	'''
	observations_count = len(params)
	param_count = len(params[0])

	results = np.zeros((observations_count, param_count))
	good_fitness = 0
	for i, param in enumerate(params):
		if fitness_cut_off < max(param.fitness.values):
			continue
		results[good_fitness] = param
		good_fitness += 1

	results = results

	for c in range(good_fitness):
		x = np.arange(param_count)
		y = results[c, :]
		ax.scatter(x=x, y=y, s=float(markersize), marker=marker, color=color)

	if plot_bounds:
		def plot_tick(column, y):
			col_width = 0.25
			x = [column - col_width,
				 column + col_width]
			y = [y, y]
			ax.plot(x, y, color='black')

		# plot min and max
		for i, parameter in enumerate(opt.evaluator.params):
			min_value = parameter.lower_bound
			max_value = parameter.upper_bound
			plot_tick(i, min_value)
			plot_tick(i, max_value)


def plot_diversity(opt, checkpoint_file, param_names):
	'''
	plot the whole history, the hall of fame, and the best individual
	from a unpickled checkpoint
	'''
	checkpoint = pickle.load(open(checkpoint_file, "r"))

	fig = plt.figure(facecolor='white')
	ax = fig.add_subplot(1, 1, 1)

	import copy
	best_individual = copy.deepcopy(checkpoint['halloffame'][0])

	# for index, param_name in enumerate(opt.evaluator.param_names):
	# 	best_individual[index] = release_params[param_name]
	
	plot_individual_params(
		opt,
		ax,
		checkpoint['history'].genealogy_history.values(),
		marker='.',
		color='grey',
		plot_bounds=True) 
	
	plot_individual_params(opt, ax, checkpoint['halloffame'],
						   marker='o', color='black')
	
	plot_individual_params(opt,
						   ax,
						   [checkpoint['halloffame'][0]],
						   markersize=150,
						   marker='x',
						   color='blue') 
	
	plot_individual_params(opt, ax, [best_individual], markersize=150,
						   marker='x', color='red')

	labels = [name.replace('.', '\n') for name in param_names]

	param_count = len(checkpoint['halloffame'][0])
	x = range(param_count)
	for xline in x:
		ax.axvline(xline, linewidth=1, color='grey', linestyle=':')

	plt.xticks(x, labels, rotation=80, ha='center', size='small')
	ax.set_xlabel('Parameter names')
	ax.set_ylabel('Parameter values')
	ax.set_yscale('log')
	ax.set_ylim(bottom=1e-7)

	plt.tight_layout()
	plt.plot()
	ax.set_autoscalex_on(True)

	return fig, ax