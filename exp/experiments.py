# parameters
ENABLE  = 1;
DISABLE = 0;

# solver
BB = 1;
BC = 2;

# node selection strategy
CPX_NODESEL_DFS         = 0 # Depth-first search
CPX_NODESEL_BESTBOUND   = 1 # Best-bound search; default
CPX_NODESEL_BESTEST     = 2 # Best-estimate search
CPX_NODESEL_BESTEST_ALT = 3 # Alternative best-estimate search 

# variable selection
CPX_VARSEL_MININFEAS 	 = -1 # Branch on variable with minimum infeasibility
CPX_VARSEL_DEFAULT   	 = 0  # Automatic: let CPLEX choose variable to branch on; default
CPX_VARSEL_MAXINFEAS 	 = 1  # Branch on variable with maximum infeasibility
CPX_VARSEL_PSEUDO 	   	 = 2  # Branch based on pseudo costs
CPX_VARSEL_STRONG 		 = 3  # Strong branching
CPX_VARSEL_PSEUDOREDUCED = 4  # Branch based on pseudo reduced costs

experiments = [
	{
		'enable': DISABLE,
		'type': 1,
		'description:': 'Different combinations of amount of nodes, partitions and densities.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [40,150,20],
		'density': 0,
		'density_range': [10,100,10],
		'partition_size': 0,
		'partitions_sizes': [20,70,10],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'clique_only': ENABLE,
		'load_limit': 30,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Graph density',
			'ylabel': 'Time (secs)',
			'title' : 'Branch & Bound vs. Branch & Cut ({vertex_size} nodes, {partition_size} partitions)',
			'xticks': ('10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', 'clique'),
			'filename': "../docs/img/bb_vs_bc_v{vertex_size}_p{partition_size}_i{iterations}_co{clique_only}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': DISABLE,
		'type': 2,
		'description:': 'Effect of breaking the symmetry.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [20,30,5],
		'density': 0,
		'density_range': [10,100,10],
		'partition_size': 0,
		'partitions_sizes': [5,15,5],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'clique_only': ENABLE,
		'load_limit': 30,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Bound Breaking Symmetry',
			'xlabel': 'Graph density',
			'ylabel': 'Time (secs)',
			'title' : 'Symmetry Breakers ({vertex_size} nodes, {partition_size} partitions)',
			'xticks': [('10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', 'clique')],
			'filename': "../docs/img/symmetry_v{vertex_size}_p{partition_size}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 3,
		'description:': 'Performance as the number of partitions increase.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 50,
		'vertex_sizes': [30,40,10],
		'density': 20,
		'density_range': [50,60,10],
		'partition_size': 0,
		'partitions_sizes': [20, 50, 5],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'clique_only': ENABLE,
		'load_limit': 30,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Number of partitions',
			'ylabel': 'Time (secs)',
			'title' : 'Effect of number of partitions on runtime ({vertex_size} nodes)',
			'xticks': range(20, 50, 5),
			'filename': "../docs/img/partitions_v{vertex_size}_d{density}_p{partition_size}_i{iterations}_co{clique_only}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': DISABLE,
		'type': 4,
		'description:': 'Effect of increasing the number of nodes.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [20,105,5],
		'density': 50,
		'density_range': [10,10,10],
		'partition_size': 10,
		'partitions_sizes': [10,15,5],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'clique_only': ENABLE,
		'load_limit': 30,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Nodes',
			'ylabel': 'Time (secs)',
			'title' : 'Effect of increasing the number of nodes ({density}% density, {partition_size} partitions)',
			'xticks': range(20,105,5),
			'filename': "../docs/img/nodes_d{density}_p{partition_size}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
]

import os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt

def findText(start, end, string):
	s = string.find(start)
	if s == -1:
		return "0"

	e = string[s+len(start):].find(end)
	if e == -1:
		return "0"

	return string[s+len(start):s+len(start)+e]

def generateGraph(experiment):
	# generate graph
	command = './generate {graph_type} {vertex_size} {density}'
	res = subprocess.check_output(command.format(**experiment), shell=True, universal_newlines=False)

def solveLP(experiment):

	# run graph coloring
	command = './coloring {graph_filename} {solver} {partition_size} {symmetry_breaker} {iterations} {clique_only} {load_limit} {traversal_strategy} {branching_strategy}'
	res = subprocess.check_output(command.format(**experiment), shell=True, universal_newlines=False)

	# find time taken to add cutting planes
	if experiment['solver'] == 2:
		time_taken_cp = float(findText("Time taken to add cutting planes: ","\n", res))
	else:
		time_taken_cp = 0

	# find time taken by LP
	time_taken_lp = float(findText("Time taken to solve final LP: ", "\n", res))

	total_time = time_taken_lp + time_taken_cp

	colors_used = int(findText("Colors used: ", "\n", res))

	if experiment['solver'] == 1:
		print "Result of Branch & Bound"
	else:
		print "Result of Branch & Cut"
		print "Time taken by CP: %f" % (time_taken_cp)
		print "Time taken by LP: %f" % (time_taken_lp)

		clique_res = findText("Loaded ", " unsatisfied", res)
		print "Clique restrictions loaded: %s" % (clique_res)
	
	print "Colors used: %d" % (colors_used)
	print "Total time: %f" % (total_time)

	return total_time

def createGraph(time1, time2, experiment):

	n_groups = len(time1)

	fig, ax = plt.subplots()

	index = np.arange(n_groups)
	bar_width = 0.35

	opacity = 0.4
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index, time1, bar_width,
	                 alpha=opacity,
	                 color='b',
	                 error_kw=error_config,
	                 label=experiment['settings']['label1'])

	if len(time2) > 0:
		rects2 = plt.bar(index + bar_width, time2, bar_width,
		                 alpha=opacity,
		                 color='r',
		                 error_kw=error_config,
		                 label=experiment['settings']['label2'])

	plt.xlabel(experiment['settings']['xlabel'])
	plt.ylabel(experiment['settings']['ylabel'])
	plt.title(experiment['settings']['title'].format(**experiment))
	plt.xticks(index + bar_width, experiment['settings']['xticks'])
	plt.legend(loc='upper left')

	plt.tight_layout()
	plt.savefig(experiment['settings']['filename'].format(**experiment))
	# plt.show(block=False)

def runExperiment1(experiment):

	for partition_size in range(*experiment['partitions_sizes']):

		experiment['partition_size'] = partition_size

		for vertex_size in range(*experiment['vertex_sizes']):

			if partition_size > vertex_size:
				continue;

			experiment['vertex_size'] = vertex_size

			if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
				print "Experiment already exists."
				continue

			time1 = []
			time2 = []

			print "\nvertex_size: {vertex_size}, partition_size: {partition_size}".format(**experiment)

			for density in range(*experiment['density_range']):

				print "Density %d" % density

				experiment['density'] = density

				generateGraph(experiment)

				experiment['solver'] = BB
				ellapsed_time = solveLP(experiment)
				time1.append(ellapsed_time)

				print "-"*20
				experiment['solver'] = BC
				ellapsed_time = solveLP(experiment)
				time2.append(ellapsed_time)
				print "\n"

			createGraph(time1, time2, experiment)

def runExperiment2(experiment):

	for partition_size in range(*experiment['partitions_sizes']):

		experiment['partition_size'] = partition_size

		for vertex_size in range(*experiment['vertex_sizes']):

			if partition_size > vertex_size:
				continue;

			experiment['vertex_size'] = vertex_size

			print vertex_size

			if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
				print "Experiment already exists."
				continue

			time1 = []
			time2 = []

			print "\nvertex_size: {vertex_size}, partition_size: {partition_size}".format(**experiment)

			for density in range(*experiment['density_range']):

				print "Density %d" % density

				experiment['density'] = density

				generateGraph(experiment)

				experiment['symmetry_breaker'] = DISABLE
				ellapsed_time = solveLP(experiment)
				time1.append(ellapsed_time)

				print "-"*20
				experiment['symmetry_breaker'] = ENABLE
				ellapsed_time = solveLP(experiment)
				time2.append(ellapsed_time)
				print "\n"

			createGraph(time1, time2, experiment)

def runExperiment3(experiment):

	generateGraph(experiment)

	time1 = []
	time2 = []

	for partition_size in range(*experiment['partitions_sizes']):

		experiment['partition_size'] = partition_size

		print "Partition size: %d" % partition_size

		experiment['solver'] = BB
		ellapsed_time = solveLP(experiment)
		time1.append(ellapsed_time)

		experiment['solver'] = BC
		ellapsed_time = solveLP(experiment)
		time2.append(ellapsed_time)
		print "\n"

	createGraph(time1, time2, experiment)

def runExperiment4(experiment):

	time1 = []
	time2 = []

	for vertex_size in range(*experiment['vertex_sizes']):

		experiment['vertex_size'] = vertex_size

		generateGraph(experiment)

		print "Vertex size: %d" % vertex_size

		experiment['solver'] = BB
		ellapsed_time = solveLP(experiment)
		time1.append(ellapsed_time)

		experiment['solver'] = BC
		ellapsed_time = solveLP(experiment)
		time2.append(ellapsed_time)
		print "\n"

	createGraph(time1, time2, experiment)

if __name__ == "__main__":

	for experiment in experiments:

		if experiment['enable'] == DISABLE:
			continue

		print "Experiment type %s" % experiment['type']

		if experiment['type'] == 1:
			runExperiment1(experiment)

		if experiment['type'] == 2:
			runExperiment2(experiment)

		if experiment['type'] == 3:
			runExperiment3(experiment)

		if experiment['type'] == 4:
			runExperiment4(experiment)