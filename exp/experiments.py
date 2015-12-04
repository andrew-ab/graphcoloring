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

def generateGraph(graph_type, vertex_size, density):
	# generate graph
	generate_graph = ['./generate', str(graph_type), str(vertex_size), str(density)]
	res = subprocess.check_output(generate_graph, universal_newlines=False)

def solveLP(solver, partitions, iterations):

	# run graph coloring
	color_graph    = ['./coloring', 'graph', str(solver), str(partitions), str(iterations)]
	res = subprocess.check_output(color_graph, universal_newlines=False)

	# find time taken to add cutting planes
	if solver == 2:
		time_taken_cp = float(findText("Time taken to add cutting planes: ","\n", res))
	else:
		time_taken_cp = 0

	# find time taken by LP
	time_taken_lp = float(findText("Time taken to solve final LP: ", "\n", res))

	total_time = time_taken_lp + time_taken_cp

	colors_used = int(findText("Colors used: ", "\n", res))

	if solver == 1:
		print "Result of Branch & Bound"
	else:
		print "Result of Branch & Cut"
		print "Time taken by CP: %f" % (time_taken_cp)
		print "Time taken by LP: %f" % (time_taken_lp)

		clique_res = int(findText("Found ", " unsatisfied clique", res))
		print "Clique restrictions found: %d" % (clique_res)
	
	print "Colors used: %d" % (colors_used)
	print "Total time: %f" % (total_time)

	return total_time

def createGraph(time_bb, time_bc, vertex_size, filename):

	n_groups = len(time_bc)

	fig, ax = plt.subplots()

	index = np.arange(n_groups)
	bar_width = 0.35

	opacity = 0.4
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index, time_bc, bar_width,
	                 alpha=opacity,
	                 color='b',
	                 error_kw=error_config,
	                 label='Branch & Cut')

	rects2 = plt.bar(index + bar_width, time_bb, bar_width,
	                 alpha=opacity,
	                 color='r',
	                 error_kw=error_config,
	                 label='Branch & Bound')

	plt.xlabel('Graph density')
	plt.ylabel('Time (ms)')
	plt.title('Branch & Bound vs. Branch & Cut ('+ str(vertex_size) +' nodes)')
	plt.xticks(index + bar_width, ('10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', 'clique'))
	plt.legend(loc='upper left')

	plt.tight_layout()
	plt.savefig("../docs/img/"+filename)
	plt.show()

if __name__ == "__main__":

	# 1st experiment
	graph_type = 1   # unused
	vertex_size = 50
	iterations = 3
	partitions = 10;

	time_bb = [];    # branch & bouund
	time_bc = [];    # branch & cut

	for density in range(10,100,10):

		print "Density: %d" % (density)
		print "-"*20

		generateGraph(graph_type, vertex_size, density)

		ellapsed_time = solveLP(1, partitions, iterations)
		time_bb.append(ellapsed_time)
		print "-"*20
		ellapsed_time = solveLP(2, partitions, iterations)
		time_bc.append(ellapsed_time)
		print "\n"

	createGraph(time_bb, time_bc, vertex_size, "all")

	# print "-"*20

	# 2nd experiment
	# graph_type = 1   # unused
	# vertex_size = 55

	# time_bb = [];    # branch & bouund
	# time_bc = [];    # branch & cut

	# time_bb = []; # branch & bouund
	# time_bc = []; # branch & cut

	# for density in range(30,40,10):
	# 	ellapsed_time = solveLP(2, graph_type, vertex_size, density)
	# 	time_bc.append(ellapsed_time)

	# for density in range(30,40,10):
	# 	ellapsed_time = solveLP(1, graph_type, vertex_size, density)
	# 	time_bb.append(ellapsed_time)

	# createGraph(time_bb, time_bc, vertex_size, "bc")