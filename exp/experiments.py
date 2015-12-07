# Experiments for LP configurations
# Advice: Next time, generate csv dataset first, then generate graphs from the dataset,
# the code ends up being much cleaner and it's easier to generate new graphs without
# having to re-generate the whole thing.

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
CPX_VARSEL_MININFEAS     = -1 # Branch on variable with minimum infeasibility
CPX_VARSEL_DEFAULT       = 0  # Automatic: let CPLEX choose variable to branch on; default
CPX_VARSEL_MAXINFEAS     = 1  # Branch on variable with maximum infeasibility
CPX_VARSEL_PSEUDO        = 2  # Branch based on pseudo costs
CPX_VARSEL_STRONG        = 3  # Strong branching
CPX_VARSEL_PSEUDOREDUCED = 4  # Branch based on pseudo reduced costs

# select cuts
CUTS_CLIQUE_ONLY  = 0
CUTS_ODDHOLE_ONLY = 1
CUTS_ALL          = 2

# select cplex config
CPLEX_DEFAULT = 0
CPLEX_CUSTOM  = 1

experiments = [
	{
		'enable': ENABLE,
		'type': 1,
		'description': 'Different combinations of amount of vertices, partitions and densities.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [20,40],
		'density': 0,
		'density_range': [30, 50],
		'partition_size': 0,
		'partitions_sizes': [10,20],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Graph density',
			'ylabel': 'Time (secs)',
			'title' : 'Branch & Bound vs. Branch & Cut Time ({vertex_size} vertices, {partition_size} partitions)',
			'xticks': ('30%', '50%', '70%', '90%'),
			'filename': "../docs/img/{type}-bb_vs_bc_v{vertex_size}_p{partition_size}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 2,
		'description': 'Effect of breaking the symmetry.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [20,30,5],
		'density': 0,
		'density_range': [10,90,10],
		'partition_size': 0,
		'partitions_sizes': [5,15,5],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Bound Breaking Symmetry',
			'xlabel': 'Graph density',
			'ylabel': 'Time (secs)',
			'title' : 'Symmetry Breakers ({vertex_size} vertices, {partition_size} partitions)',
			'xticks': ('10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', 'clique'),
			'filename': "../docs/img/{type}-symmetry_v{vertex_size}_p{partition_size}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 3,
		'description': 'Performance as the number of partitions increase.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 30,
		'density': 20,
		'density_range': [30, 50],
		'partition_size': 0,
		'partitions_sizes': [10, 35, 5],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Partitions',
			'ylabel': 'Time (secs)',
			'title' : 'Effect of number of partitions on runtime ({vertex_size} vertices, {density} density)',
			'xticks': range(20, 50, 5),
			'filename': "../docs/img/{type}-partitions_v{vertex_size}_d{density}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 4,
		'description': 'Effect of increasing the number of vertices.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [20,105,5],
		'density': 50,
		'partition_size': 10,
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Nodes',
			'ylabel': 'Time (secs)',
			'title' : 'Effect of increasing the number of vertices ({density}% density, {partition_size} partitions)',
			'xticks': range(20,105,5),
			'filename': "../docs/img/{type}-nodes_d{density}_p{partition_size}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 5,
		'description': 'Cutting Plane Strategies.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 40,
		'density': 0,
		'density_range': [10, 30, 50],
		'partition_size': 20,
		'solver': BC,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'C&B clique cuts',
			'label2': 'C&B oddhole cuts',
			'label3': 'C&B both',
			'label4': 'B&B',
			'xlabel': 'Graph density',
			'ylabel': 'Time (secs)',
			'title' : 'Cutting Plane Strategies ({vertex_size} vertices, {partition_size} partitions)',
			'xticks': ('10%', '30%', '50%', '70%'),
			'filename': "../docs/img/{type}-cuts_v{vertex_size}_p{partition_size}_i{iterations}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 6,
		'description': 'Cutting Plane Thresholds.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 40,
		'density': 50,
		'partition_size': 20,
		'partitions_sizes': [10,20],
		'solver': BC,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'load_limit_sizes': [10, 20, 30, 40, 50],
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'C&B clique cuts',
			'label2': 'C&B oddhole cuts',
			'label3': 'C&B both',
			'label4': 'B&B',
			'xlabel': 'Threshold',
			'ylabel': 'Time (secs)',
			'title' : 'Cutting Plane Thresholds ({vertex_size} vertices, {partition_size} partitions)',
			'xticks': [10, 20, 30, 40, 50],
			'filename': "../docs/img/{type}-thresholds_v{vertex_size}_p{partition_size}_i{iterations}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 7,
		'description': 'Cutting Plane Iterations.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 40,
		'density': 70,
		'partition_size': 10,
		'solver': BC,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'iteration_list': [1,2,3],
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'load_limit_sizes': [10, 20, 30, 40, 50],
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'C&B clique cuts',
			'label2': 'C&B oddhole cuts',
			'label3': 'C&B both',
			'label4': 'B&B',
			'xlabel': 'Iterations',
			'ylabel': 'Time (secs)',
			'title' : 'Iterations ({vertex_size} vertices, {density}% density, {partition_size} partitions, {load_limit} threshold)',
			'xticks': (1,2,3),
			'filename': "../docs/img/{type}-iterations_v{vertex_size}_p{partition_size}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	},
	{
		'enable': ENABLE,
		'type': 8,
		'description': 'B&B vs C&B vs Default CPLEX.',
		'graph_type': 1,
		'graph_filename': 'graph',
		'vertex_size': 0,
		'vertex_sizes': [20,40],
		'density': 0,
		'density_range': [30, 50, 70],
		'partition_size': 0,
		'partitions_sizes': [10,20],
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'load_limit': 40,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT,	
		'settings': {
			'label1': 'CPLEX default',
			'label2': 'Branch & Bound',
			'label3': 'C&B clique only',
			'xlabel': 'Graph density',
			'ylabel': 'Time (secs)',
			'title' : 'B&B vs C&B vs Default CPLEX ({vertex_size} vertices, {partition_size} partitions)',
			'xticks': ('30%', '50%', '70%', '90%'),
			'filename': "../docs/img/{type}-compare_v{vertex_size}_p{partition_size}_i{iterations}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}
	}
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
	command = './coloring {graph_filename} {solver} {partition_size} {symmetry_breaker} {iterations} {select_cuts} {load_limit} {custom_config} {traversal_strategy} {branching_strategy}'
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

	nodes_traversed = int(findText("Nodes traversed: ", "\n", res))

	if experiment['solver'] == 1:
		print "Result of Branch & Bound"
	else:
		print "Result of Branch & Cut"
		print "Time taken by CP: %f" % time_taken_cp
		print "Time taken by LP: %f" % time_taken_lp

		clique_res = findText("Loaded ", " unsatisfied clique restrictions!", res)
		print "Clique restrictions loaded: %s" % clique_res
	
		# oddhole_res = findText("Loaded ", " unsatisfied oddhole restrictions!", res) # fix substr
		# print "Oddhole restrictions loaded: %s" % oddhole_res


	print "Colors used: %d" % colors_used
	print "Nodes traversed: %d" % nodes_traversed
	print "Total time: %f" % total_time

	return {'total_time': total_time, 'nodes_traversed': nodes_traversed}

def createGraph(time1, time2, time3, time4, experiment):

	n_groups = len(time1)

	fig, ax = plt.subplots()

	index = np.arange(n_groups)
	bar_width = 0.2

	opacity = 0.4
	error_config = {'ecolor': '0.3'}


	if len(time4) > 0:
		s1 = -bar_width
		s2 = 0
		s3 = bar_width
		s4 = 2*bar_width
	elif len(time3) > 0:
		s1 = -bar_width*1/2
		s2 = bar_width/2
		s3 = bar_width*3/2
	else:
		s1 = 0
		s2 = bar_width

	rects1 = plt.bar(index + s1, time1, bar_width,
	                 alpha=opacity,
	                 color='b',
	                 error_kw=error_config,
	                 label=experiment['settings']['label1'])

	if len(time2) > 0:
		rects2 = plt.bar(index + s2, time2, bar_width,
		                 alpha=opacity,
		                 color='r',
		                 error_kw=error_config,
		                 label=experiment['settings']['label2'])

	if len(time3) > 0:
		rects3 = plt.bar(index + s3, time3, bar_width,
		                 alpha=opacity,
		                 color='g',
		                 error_kw=error_config,
		                 label=experiment['settings']['label3'])	

	if len(time4) > 0:
		rects3 = plt.bar(index + s4, time4, bar_width,
		                 alpha=opacity,
		                 color='y',
		                 error_kw=error_config,
		                 label=experiment['settings']['label4'])	

	plt.xlabel(experiment['settings']['xlabel'])
	plt.ylabel(experiment['settings']['ylabel'])
	plt.title(experiment['settings']['title'].format(**experiment))
	plt.xticks(index + bar_width, experiment['settings']['xticks'])
	plt.legend(loc='upper left')

	plt.tight_layout()
	plt.savefig(experiment['settings']['filename'].format(**experiment))
	# plt.show(block=False)

def runExperiment1(experiment):

	for partition_size in experiment['partitions_sizes']:

		experiment['partition_size'] = partition_size

		for vertex_size in experiment['vertex_sizes']:

			if partition_size > vertex_size:
				continue;

			experiment['vertex_size'] = vertex_size

			if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
				print "Experiment already exists."
				continue

			time1 = []
			time2 = []
			nodes1 = []
			nodes2 = []

			print "\nvertex_size: {vertex_size}, partition_size: {partition_size}".format(**experiment)

			for density in experiment['density_range']:

				print "Density %d" % density

				experiment['density'] = density

				generateGraph(experiment)

				experiment['solver'] = BB
				sol = solveLP(experiment)
				time1.append(sol['total_time'])
				nodes1.append(sol['nodes_traversed'])

				print "-"*20
				experiment['solver'] = BC
				sol = solveLP(experiment)
				time2.append(sol['total_time'])
				nodes2.append(sol['nodes_traversed'])
				print "\n"

			createGraph(time1, time2, [], [], experiment)

			experiment['settings'] = {
				'label1': 'Branch & Bound',
				'label2': 'Branch & Cut',
				'xlabel': 'Graph density',
				'ylabel': 'Nodes',
				'title' : 'Branch & Bound vs. Branch & Cut Nodes ({vertex_size} vertices, {partition_size} partitions)',
				'xticks': ('30%', '50%', '70%', '90%'),
				'filename': "../docs/img/{type}-bb_vs_bc_v{vertex_size}_p{partition_size}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}_nodes"
			}

			createGraph(nodes1, nodes2, [], [], experiment)

			experiment['settings'] = {
				'label1': 'Branch & Bound',
				'label2': 'Branch & Cut',
				'xlabel': 'Graph density',
				'ylabel': 'Time (secs)',
				'title' : 'Branch & Bound vs. Branch & Cut Time ({vertex_size} vertices, {partition_size} partitions)',
				'xticks': ('30%', '50%', '70%', '90%'),
				'filename': "../docs/img/{type}-bb_vs_bc_v{vertex_size}_p{partition_size}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
			}


def runExperiment2(experiment):

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

				experiment['symmetry_breaker'] = DISABLE
				sol = solveLP(experiment)
				time1.append(sol['total_time'])

				print "-"*20
				experiment['symmetry_breaker'] = ENABLE
				sol = solveLP(experiment)
				time2.append(sol['total_time'])
				print "\n"

			createGraph(time1, time2, [], [], experiment)

def runExperiment3(experiment):

	for density in experiment['density_range']:

		print "Density %d" % density

		experiment['density'] = density

		print experiment['settings']['filename'].format(**experiment)+".png"

		if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
			print "Experiment already exists."
			continue

		generateGraph(experiment)

		time1 = []
		time2 = []

		nodes1 = []
		nodes2 = []

		for partition_size in range(*experiment['partitions_sizes']):

			experiment['partition_size'] = partition_size

			print "Partition size: %d" % partition_size

			experiment['solver'] = BB
			sol = solveLP(experiment)
			time1.append(sol['total_time'])
			nodes1.append(sol['nodes_traversed'])

			experiment['solver'] = BC
			sol = solveLP(experiment)
			time2.append(sol['total_time'])
			nodes2.append(sol['nodes_traversed'])
			print "\n"

		createGraph(time1, time2, [], [], experiment)

		experiment['settings'] = {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Partitions',
			'ylabel': 'Nodes',
			'title' : 'Effect of number of partitions on nodes traversed ({vertex_size} vertices, {density} density)',
			'xticks': range(20, 50, 5),
			'filename': "../docs/img/{type}-partitions_v{vertex_size}_d{density}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}_nodes"
		}

		createGraph(nodes1, nodes2, [], [], experiment)

		experiment['settings'] = {
			'label1': 'Branch & Bound',
			'label2': 'Branch & Cut',
			'xlabel': 'Partitions',
			'ylabel': 'Time (secs)',
			'title' : 'Effect of number of partitions on runtime ({vertex_size} vertices, {density} density)',
			'xticks': range(20, 50, 5),
			'filename': "../docs/img/{type}-partitions_v{vertex_size}_d{density}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}"
		}

def runExperiment4(experiment):

	if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
		print "Experiment already exists."
		return

	time1 = []
	time2 = []

	for vertex_size in range(*experiment['vertex_sizes']):

		experiment['vertex_size'] = vertex_size

		generateGraph(experiment)

		print "Vertex size: %d" % vertex_size

		experiment['solver'] = BB
		sol = solveLP(experiment)
		time1.append(sol['total_time'])

		experiment['solver'] = BC
		sol = solveLP(experiment)
		time2.append(sol['total_time'])
		print "\n"

	createGraph(time1, time2, [], [], experiment)

def runExperiment5(experiment):

	print experiment['settings']['filename'].format(**experiment)+".png"

	if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
		print "Experiment already exists."
		return

	time1 = []
	time2 = []
	time3 = []
	time4 = []

	nodes1 = []
	nodes2 = []
	nodes3 = []
	nodes4 = []

	for density in experiment['density_range']:

		experiment['density'] = density

		generateGraph(experiment)

		print "Density: %d" % density

		experiment['solver'] = BC

		print "Clique only:"
		experiment['select_cuts'] = CUTS_CLIQUE_ONLY
		sol = solveLP(experiment)
		time1.append(sol['total_time'])
		nodes1.append(sol['nodes_traversed'])

		print "Oddhole only:"
		experiment['select_cuts'] = CUTS_ODDHOLE_ONLY
		sol = solveLP(experiment)
		time2.append(sol['total_time'])
		nodes2.append(sol['nodes_traversed'])

		print "All cuts:"
		experiment['select_cuts'] = CUTS_ALL
		sol = solveLP(experiment)
		time3.append(sol['total_time'])
		nodes3.append(sol['nodes_traversed'])

		print "Branch and Bound:"
		experiment['solver'] = BB
		sol = solveLP(experiment)
		time4.append(sol['total_time'])
		nodes4.append(sol['nodes_traversed'])

		print "\n"

	createGraph(time1, time2, time3, time4, experiment)

	experiment['settings'] = {
		'label1': 'C&B clique cuts',
		'label2': 'C&B oddhole cuts',
		'label3': 'C&B both',
		'label4': 'B&B',
		'xlabel': 'Graph density',
		'ylabel': 'Nodes',
		'title' : 'Cutting Plane Strategies ({vertex_size} vertices, {partition_size} partitions)',
		'xticks': ('10%', '30%', '50%', '70%'),
		'filename': "../docs/img/{type}-cuts_v{vertex_size}_p{partition_size}_i{iterations}_co{select_cuts}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}_nodes"
	}

	createGraph(nodes1, nodes2, nodes3, nodes4, experiment)

def runExperiment6(experiment):

	if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
		print "Experiment already exists."
		return

	generateGraph(experiment)

	time1 = []
	time2 = []
	time3 = []
	time4 = []

	nodes1 = []
	nodes2 = []
	nodes3 = []
	nodes4 = []

	for limit in experiment['load_limit_sizes']:

		experiment['load_limit'] = limit;

		print "Threshold: %d" % limit

		experiment['solver'] = BC

		print "Clique only:"
		experiment['select_cuts'] = CUTS_CLIQUE_ONLY
		sol = solveLP(experiment)
		time1.append(sol['total_time'])
		nodes1.append(sol['nodes_traversed'])

		print "Oddhole only:"
		experiment['select_cuts'] = CUTS_ODDHOLE_ONLY
		sol = solveLP(experiment)
		time2.append(sol['total_time'])
		nodes2.append(sol['nodes_traversed'])

		print "All cuts:"
		experiment['select_cuts'] = CUTS_ALL
		sol = solveLP(experiment)
		time3.append(sol['total_time'])
		nodes3.append(sol['nodes_traversed'])

		print "Branch and Bound:"
		experiment['solver'] = BB
		sol = solveLP(experiment)
		time4.append(sol['total_time'])
		nodes4.append(sol['nodes_traversed'])

		print "\n"

	createGraph(time1, time2, time3, time4, experiment)

	experiment['settings'] = {
		'label1': 'C&B clique cuts',
		'label2': 'C&B oddhole cuts',
		'label3': 'C&B both',
		'label4': 'B&B',
		'xlabel': 'Threshold',
		'ylabel': 'Nodes',
		'title' : 'Cutting Plane Thresholds ({vertex_size} vertices, {partition_size} partitions)',
		'xticks': [10, 20, 30, 40, 50],
		'filename': "../docs/img/{type}-thresholds_v{vertex_size}_p{partition_size}_i{iterations}_t{traversal_strategy}_b{branching_strategy}_nodes"
	}

	createGraph(nodes1, nodes2, nodes3, nodes4, experiment)

def runExperiment7(experiment):


	if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
		print "Experiment already exists."
		return

	generateGraph(experiment)

	time1 = []
	time2 = []
	time3 = []
	time4 = []

	nodes1 = []
	nodes2 = []
	nodes3 = []
	nodes4 = []

	for iterations in experiment['iteration_list']:

		experiment['iterations'] = iterations;

		print "Clique only:"
		experiment['select_cuts'] = CUTS_CLIQUE_ONLY
		sol = solveLP(experiment)
		time1.append(sol['total_time'])
		nodes1.append(sol['nodes_traversed'])

		print "Oddhole only:"
		experiment['select_cuts'] = CUTS_ODDHOLE_ONLY
		sol = solveLP(experiment)
		time2.append(sol['total_time'])
		nodes2.append(sol['nodes_traversed'])

		print "All cuts:"
		experiment['select_cuts'] = CUTS_ALL
		sol = solveLP(experiment)
		time3.append(sol['total_time'])
		nodes3.append(sol['nodes_traversed'])

		print "Branch and Bound:"
		experiment['solver'] = BB
		sol = solveLP(experiment)
		time4.append(sol['total_time'])
		nodes4.append(sol['nodes_traversed'])

		print "\n"

	createGraph(time1, time2, time3, time4, experiment)

	experiment['settings'] = {
		'label1': 'C&B clique cuts',
		'label2': 'C&B oddhole cuts',
		'label3': 'C&B both',
		'label4': 'B&B',
		'xlabel': 'Iterations',
		'ylabel': 'Nodes',
		'title' : 'Iterations ({vertex_size} vertices, {density}% density, {partition_size} partitions, {load_limit} threshold)',
		'xticks': (1,2,3),
		'filename': "../docs/img/{type}-iterations_v{vertex_size}_p{partition_size}_l{load_limit}_t{traversal_strategy}_b{branching_strategy}_nodes"
	}

	# createGraph(nodes1, nodes2, nodes3, nodes4, experiment)

def runExperiment8(experiment):

	for partition_size in experiment['partitions_sizes']:

		experiment['partition_size'] = partition_size

		for vertex_size in experiment['vertex_sizes']:

			if partition_size > vertex_size:
				continue;

			experiment['vertex_size'] = vertex_size

			if os.path.isfile(experiment['settings']['filename'].format(**experiment)+".png"):
				print "Experiment already exists."
				continue

			time1 = []
			time2 = []
			time3 = []

			print "\nvertex_size: {vertex_size}, partition_size: {partition_size}".format(**experiment)

			for density in experiment['density_range']:

				print "Density %d" % density

				experiment['density'] = density

				generateGraph(experiment)

				experiment['custom_config'] = CPLEX_DEFAULT

				print "-"*20
				experiment['solver'] = BB
				sol = solveLP(experiment)
				time1.append(sol['total_time'])

				experiment['custom_config'] = CPLEX_CUSTOM

				print "-"*20
				experiment['solver'] = BB
				sol = solveLP(experiment)
				time2.append(sol['total_time'])

				print "-"*20
				experiment['solver'] = BC
				sol = solveLP(experiment)
				time3.append(sol['total_time'])
				print "\n"

			createGraph(time1, time2, time3, [], experiment)

if __name__ == "__main__":

	for experiment in experiments:

		if experiment['enable'] == DISABLE:
			continue

		print "Experiment #%s" % experiment['type']
		print "Description: %s" % experiment['description']

		if experiment['type'] == 1:
			runExperiment1(experiment)

		if experiment['type'] == 2:
			runExperiment2(experiment)

		if experiment['type'] == 3:
			runExperiment3(experiment)

		if experiment['type'] == 4:
			runExperiment4(experiment)

		if experiment['type'] == 5:
			runExperiment5(experiment)

		if experiment['type'] == 6:
			runExperiment6(experiment)

		if experiment['type'] == 7:
			runExperiment7(experiment)

		if experiment['type'] == 8:
			runExperiment8(experiment)