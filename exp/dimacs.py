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

import os, fnmatch, csv, xlrd
import os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

def findText(start, end, string):
	s = string.find(start)
	if s == -1:
		return "0"

	e = string[s+len(start):].find(end)
	if e == -1:
		return "0"

	return string[s+len(start):s+len(start)+e]

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

	vertex_size = int(findText("vertex_size: ", ",", res))
	edge_size   = int(findText("edge_size: ", ",", res))

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

	return {'total_time': total_time, 'nodes_traversed': nodes_traversed, 'colors_used': colors_used, 'vertex_size': vertex_size, 'edge_size': edge_size}

if __name__ == "__main__":

	experiment = {
		'problem': '',
		'vertex_size': '',
		'edge_size': '',
		'total_time': '',
		'colors_used': '',
		'nodes_traversed': 'f',
		'partition_size': 20,
		'graph_filename': 'h',
		'solver': BB,
		'symmetry_breaker': ENABLE,
		'iterations': 1,
		'load_limit': 30,
		'select_cuts': CUTS_CLIQUE_ONLY,
		'custom_config': CPLEX_CUSTOM,
		'traversal_strategy': CPX_NODESEL_BESTBOUND,
		'branching_strategy': CPX_VARSEL_DEFAULT
	}

	csvFile = open('dimacs.csv', 'w+')
	csvFile.write(",".join(map(str,experiment.keys()))+"\n")

	for solver in [BB, BC]:

		experiment['solver'] = solver

		for filename in sorted(os.listdir('DIMACS')):
			#parse_csv(asset, filename, start, end)
			print "Processing %s" % filename

			experiment['problem'] = filename
			experiment['graph_filename'] = "DIMACS/" + filename
			res = solveLP(experiment)

			for key, value in res.iteritems():
				experiment[key] = value

			print ",".join(map(str,experiment.values()))+"\n"

			csvFile.write(",".join(map(str,experiment.values()))+"\n")