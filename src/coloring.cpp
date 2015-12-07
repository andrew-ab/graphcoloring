#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

#include <stdlib.h>
#include <cassert>

#include <algorithm>
#include <string>
#include <vector>
#include <set>

#define TOL 1e-05

ILOSTLBEGIN // macro to define namespace

// helper functions
int getVertexIndex(int id, int color, int partition_size);
inline int fromMatrixToVector(int from, int to, int vertex_size);
inline bool isAdyacent(int from, int to, int vertex_size, bool* adjacencyList);
bool adyacentToAll(int id, int vertex_size, bool* adjacencyList, const set<int>& clique);
bool cliqueNotContained(const set<int>& clique, int color, const vector<tuple<double, int, set<int> > >& clique_familly);

// load LP
int loadObjectiveFunction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size, char vtype);
int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList);
int loadSingleColorInPartitionRestriction(CPXENVptr& env, CPXLPptr& lp, vector<vector<int> >& partitions, int partition_size);
int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size);
int loadSymmetryBreaker(CPXENVptr& env, CPXLPptr& lp, int partition_size);

// cutting planes
int loadCuttingPlanes(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList, int iterations, int load_limit, int select_cuts);
int maximalCliqueFamillyHeuristic(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList, double* sol, int load_limit);
int loadUnsatisfiedCliqueRestriction(CPXENVptr& env, CPXLPptr& lp, int partition_size, const set<int>& clique, int color);

int oddholeFamillyHeuristic(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList, double* sol, int load_limit);
int loadUnsatisfiedOddholeRestriction(CPXENVptr& env, CPXLPptr& lp, int partition_size, const set<int>& path, int color);

// cplex functions
int solveLP(CPXENVptr& env, CPXLPptr& lp, int edge_size, int vertex_size, int partition_size);
int convertVariableType(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size, char vtype);
int setTraversalStrategy(CPXENVptr& env, int strategy);
int setBranchingVariableStrategy(CPXENVptr& env, int strategy);
int setBranchAndBoundConfig(CPXENVptr& env);

int checkStatus(CPXENVptr& env, int status);

// colors array!
const char* colors[] = {"Blue", "Red", "Green", "Yellow", "Grey", "Green", "Pink", "AliceBlue","AntiqueWhite","Aqua","Aquamarine","Azure","Beige",
"Bisque","Black","BlanchedAlmond","BlueViolet","Brown","BurlyWood","CadetBlue","Chartreuse","Chocolate","Coral","CornflowerBlue",
"Cornsilk","Crimson","Cyan","DarkBlue","DarkCyan","DarkGoldenRod","DarkGray","DarkGrey","DarkGreen","DarkKhaki","DarkMagenta","DarkOliveGreen",
"Darkorange","DarkOrchid","DarkRed","DarkSalmon","DarkSeaGreen","DarkSlateBlue","DarkSlateGray","DarkSlateGrey","DarkTurquoise",
"DarkViolet","DeepPink","DeepSkyBlue","DimGray","DimGrey","DodgerBlue","FireBrick","FloralWhite","ForestGreen","Fuchsia",
"Gainsboro","GhostWhite","Gold","GoldenRod","Gray","GreenYellow","HoneyDew","HotPink","IndianRed","Indigo",
"Ivory","Khaki","Lavender","LavenderBlush","LawnGreen","LemonChiffon","LightBlue","LightCoral","LightCyan","LightGoldenRodYellow",
"LightGray","LightGrey","LightGreen","LightPink","LightSalmon","LightSeaGreen","LightSkyBlue","LightSlateGray","LightSlateGrey",
"LightSteelBlue","LightYellow","Lime","LimeGreen","Linen","Magenta","Maroon","MediumAquaMarine","MediumBlue","MediumOrchid",
"MediumPurple","MediumSeaGreen","MediumSlateBlue","MediumSpringGreen","MediumTurquoise","MediumVioletRed","MidnightBlue",
"MintCream","MistyRose","Moccasin","NavajoWhite","Navy","OldLace","Olive","OliveDrab","Orange","OrangeRed","Orchid",
"PaleGoldenRod","PaleGreen","PaleTurquoise","PaleVioletRed","PapayaWhip","PeachPuff","Peru","Plum","PowderBlue",
"Purple","RosyBrown","RoyalBlue","SaddleBrown","Salmon","SandyBrown","SeaGreen","SeaShell","Sienna","Silver","SkyBlue",
"SlateBlue","SlateGray","SlateGrey","Snow","SpringGreen","SteelBlue","Tan","Teal","Thistle","Tomato","Turquoise","Violet",
"Wheat","White","WhiteSmoke","YellowGreen"};

int main(int argc, char **argv) {

	if (argc != 11) {
		printf("Usage: %s inputFile solver partitions symmetry_breaker iterations select_cuts load_limit custom_config traversal_strategy branching_strategy\n", argv[0]);
		exit(1);
	}

	int solver = atoi(argv[2]);
	int partition_size = atoi(argv[3]);
	bool symmetry_breaker = (atoi(argv[4]) == 1);
	int iterations = atoi(argv[5]);
	int select_cuts = atoi(argv[6]);              // 0: clique only, 1: oddhole only, 2: both
	int load_limit = atoi(argv[7]);
	int custom_config = atoi(argv[8]);           // 0: default, 1: custom
	int traversal_strategy = atoi(argv[9]);
	int branching_strategy = atoi(argv[10]);

	if (solver == 1) {
		printf("Solver: Branch & Bound\n");
	} else {
		printf("Solver: Cut & Branch\n");
	}

	/* read graph input file
	 * format: http://mat.gsia.cmu.edu/COLOR/instances.html
	 * graph representation chosen in order to load the LP easily.
	 * - vector of edges
	 * - vector of partitions
	 */
	FILE* fp = fopen(argv[1], "r");

	if (fp == NULL) {
		printf("Invalid input file.\n");
		exit(1);
	}

	char buf[100];
	int vertex_size, edge_size;

	set<pair<double,int> > edges; // sometimes we have to filter directed graphs

	while (fgets(buf, sizeof(buf), fp) != NULL) {
		if (buf[0] == 'c') continue;
		else if (buf[0] == 'p') {
			sscanf(&buf[7], "%d %d", &vertex_size, &edge_size);
		}
		else if (buf[0] == 'e') {
			int from, to;
			sscanf(&buf[2], "%d %d", &from, &to);
			if (from < to) {
				edges.insert(pair<double,int>(from, to));
			} else {
				edges.insert(pair<double,int>(to, from));
			}
		}
	}

	// build adyacency list
	edge_size = edges.size();
	int adyacency_size = vertex_size*vertex_size - ((vertex_size+1)*vertex_size/2);
	bool* adjacencyList = new bool[adyacency_size]; // can be optimized even more with a bitfield.
	fill_n(adjacencyList, adyacency_size, false);
	for (set<pair<double,int> >::iterator it = edges.begin(); it != edges.end(); ++it) {
		adjacencyList[fromMatrixToVector(it->first, it->second, vertex_size)] = true;
	}

	// set random seed
	// srand(time(NULL));

	// asign every vertex to a partition
	// int partition_size = rand() % vertex_size + 1;
	vector<vector<int> > partitions(partition_size, vector<int>());

	for (int i = 0; i < vertex_size; ++i) {
		partitions[i % partition_size].push_back(i+1);
	}

	// warning: this procedure doesn't guarantee every partition will have an element.
	// for (int i = 1; i <= vertex_size; ++i) {
	// 	int assign_partition = rand() % partition_size;
	// 	partitions[assign_partition].push_back(i);
	// }

	// // update partition_size
	// for (std::vector<vector<int> >::iterator it = partitions.begin(); it != partitions.end(); ++it) {
	// 	if (it->size() == 0) --partition_size;
	// }

	printf("Graph: vertex_size: %d, edge_size: %d, partition_size: %d\n", vertex_size, edge_size, partition_size);

	// start loading LP using CPLEX
	int status;
	CPXENVptr env; // pointer to enviroment
	CPXLPptr lp;   // pointer to the lp.

	env = CPXopenCPLEX(&status); // create enviroment
	checkStatus(env, status);

	// create LP
	lp = CPXcreateprob(env, &status, "Instance of partitioned graph coloring.");
	checkStatus(env, status);

	if (custom_config == 1) setBranchAndBoundConfig(env);
	setTraversalStrategy(env, traversal_strategy);
	setBranchingVariableStrategy(env, branching_strategy);

	if (solver == 1) { // pure branch & bound
		loadObjectiveFunction(env, lp, vertex_size, partition_size, CPX_BINARY);
	} else {
		loadObjectiveFunction(env, lp, vertex_size, partition_size, CPX_CONTINUOUS);
	}

	loadAdyacencyColorRestriction(env, lp, vertex_size, edge_size, partition_size, adjacencyList);
	loadSingleColorInPartitionRestriction(env, lp, partitions, partition_size);
	loadAdyacencyColorRestriction(env, lp, vertex_size, partition_size);

	if (symmetry_breaker) loadSymmetryBreaker(env, lp, partition_size);

	if (solver != 1) loadCuttingPlanes(env, lp, vertex_size, edge_size, partition_size, adjacencyList, iterations, load_limit, select_cuts);
		
	// write LP formulation to file, great to debug.
	status = CPXwriteprob(env, lp, "graph.lp", NULL);
	checkStatus(env, status);

	convertVariableType(env, lp, vertex_size, partition_size, CPX_BINARY); 
	
	solveLP(env, lp, edge_size, vertex_size, partition_size);

	delete[] adjacencyList;

	return 0;
}

int getVertexIndex(int id, int color, int partition_size) {
	return partition_size + ((id-1)*partition_size) + (color-1);
}

/* since the adyacency matrix is symmetric and the diagonal is not needed, we can simply
 * store the upper diagonal and get adyacency from a list. the math is quite simple, it
 * just uses the formula for the sum of integers. ids are numbered starting from 1.
 */
inline int fromMatrixToVector(int from, int to, int vertex_size) {

	// for speed, many parts of this code are commented, since by our usage we always
	// know from < to and are in range.

	// assert(from != to && from <= vertex_size && to <= vertex_size);

	// if (from < to)
		return from*vertex_size - (from+1)*from/2 - (vertex_size - to) - 1;
	// else
	// 	return to*vertex_size - (to+1)*to/2 - (vertex_size - from) - 1;
}

inline bool isAdyacent(int from, int to, int vertex_size, bool* adjacencyList) {
	return adjacencyList[fromMatrixToVector(from, to, vertex_size)];
}

bool adyacentToAll(int id, int vertex_size, bool* adjacencyList, const set<int>& clique) {
	for (set<int>::iterator it = clique.begin(); it != clique.end(); ++it) {
		if (!isAdyacent(*it, id, vertex_size, adjacencyList)) return false;
	}
	return true;
}

bool cliqueNotContained(const set<int>& clique, int color, const vector<tuple<double, int, set<int> > >& clique_familly) {
	for (vector<tuple<double, int, set<int> > >::const_iterator it = clique_familly.begin(); it != clique_familly.end(); ++it) {
		// by construction, sets are already ordered.
		if (get<1>(*it) == color && includes(get<2>(*it).begin(), get<2>(*it).end(), clique.begin(), clique.end())) return false;
	}
	return true;
}

int loadObjectiveFunction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size, char vtype) {

	// load objective function
	int n = partition_size + (vertex_size*partition_size);
	double *objfun	= new double[n];
	double *ub      = new double[n];
	char	 *ctype	= new char[n];
	char **colnames = new char*[n];

	for (int i = 0; i < partition_size; ++i) {
		objfun[i] = 1;
		ub[i] = 1;
		ctype[i]	= vtype;
		colnames[i] = new char[10];
		sprintf(colnames[i], "w_%d", (i+1));
	}

	for (int id = 1; id <= vertex_size; ++id) {
		for (int color = 1; color <= partition_size; ++color) {
			int index = getVertexIndex(id, color, partition_size);
			objfun[index]   = 0;
			ub[index] = 1;
			ctype[index]    = vtype;
			colnames[index] = new char[10];
			sprintf(colnames[index], "x%d_%d", id, color);
		}
	}

	// CPLEX bug? If you set ctype, it doesn't identify the problem as continous.
	int status = CPXnewcols(env, lp, n, objfun, NULL, ub, NULL, colnames);
	checkStatus(env, status);
	
	// free memory
	for (int i = 0; i < n; ++i) {
		delete[] colnames[i];
	}
	
	delete[] objfun;
	delete[] ub;
	delete[] ctype;
	delete[] colnames;

	return 0;
}

int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList) {

	// load first restriction
	int ccnt = 0;                          // new columns being added.
	int rcnt = edge_size * partition_size; // new rows being added.
	int nzcnt = rcnt*2;                    // nonzero constraint coefficients being added.

	double *rhs = new double[rcnt];        // independent term in restrictions.
	char *sense = new char[rcnt];          // sense of restriction inequality.

	int *matbeg = new int[rcnt];           // array position where each restriction starts in matind and matval.
	int *matind = new int[rcnt*2];         // index of variables != 0 in restriction (each var has an index defined above)
	double *matval  = new double[rcnt*2];  // value corresponding to index in restriction.
	char **rownames = new char*[rcnt];     // row labels.

	int i = 0;
	for (int from = 1; from <= vertex_size; ++from) {
		for (int to = from + 1; to <= vertex_size; ++to) {

			if (!isAdyacent(from, to, vertex_size, adjacencyList)) continue;

			for (int color = 1; color <= partition_size; ++color) {
				matbeg[i] = i*2;

				matind[i*2]   = getVertexIndex(from, color, partition_size);
				matind[i*2+1] = getVertexIndex(to  , color, partition_size);

				matval[i*2]   = 1;
				matval[i*2+1] = 1;

				rhs[i] = 1;
				sense[i] = 'L';
				rownames[i] = new char[40];
				sprintf(rownames[i], "%s", colors[color-1]);

				++i;
			}
		}
	}

	// debug flag
	// status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);

	// add restriction
	int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);
	checkStatus(env, status);

	// free memory
	for (int i = 0; i < rcnt; ++i) {
		delete[] rownames[i];
	}
	
	delete[] rhs;
	delete[] sense;
	delete[] matbeg;
	delete[] matind;
	delete[] matval;
	delete[] rownames;

	return 0;
}


int loadSingleColorInPartitionRestriction(CPXENVptr& env, CPXLPptr& lp, vector<vector<int> >& partitions, int partition_size) {

	// load second restriction
	int p = 1;
	for (std::vector<vector<int> >::iterator it = partitions.begin(); it != partitions.end(); ++it) {

		int size = it->size();                 // current partition size.
		if (size == 0) continue;               // skip empty partitions.

		int ccnt = 0;                          // new columns being added.
		int rcnt = 1;                          // new rows being added.
		int nzcnt = size*partition_size;       // nonzero constraint coefficients being added.

		double *rhs = new double[rcnt];        // independent term in restrictions.
		char *sense = new char[rcnt];          // sense of restriction inequality.

		int *matbeg = new int[rcnt];           // array position where each restriction starts in matind and matval.
		int *matind = new int[nzcnt];          // index of variables != 0 in restriction (each var has an index defined above)
		double *matval  = new double[nzcnt];   // value corresponding to index in restriction.
		char **rownames = new char*[rcnt];     // row labels.

		matbeg[0] = 0;
		sense[0]  = 'E';
		rhs[0]    = 1;
		rownames[0] = new char[40];
		sprintf(rownames[0], "partition_%d", p);

		int i = 0;
		for (std::vector<int>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
			for (int color = 1; color <= partition_size; ++color) {
				matind[i] = getVertexIndex(*it2, color, partition_size);
				matval[i] = 1;
				++i;
			}
		}

		// add restriction
		int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);
		checkStatus(env, status);

		// free memory
		delete[] rownames[0];
		delete[] rhs;
		delete[] sense;
		delete[] matbeg;
		delete[] matind;
		delete[] matval;
		delete[] rownames;

		++p;
	}

	return 0;
}

int loadSymmetryBreaker(CPXENVptr& env, CPXLPptr& lp, int partition_size) {

	int ccnt = 0;                          // new columns being added.
	int rcnt = partition_size - 1;         // new rows being added.
	int nzcnt = 2*rcnt;                    // nonzero constraint coefficients being added.

	double* rhs = new double[rcnt];        // independent term in restrictions.
	char *sense = new char[rcnt];          // sense of restriction inequality.

	int *matbeg = new int[rcnt];           // array position where each restriction starts in matind and matval.
	int *matind = new int[rcnt*2];         // index of variables != 0 in restriction (each var has an index defined above)
	double *matval  = new double[rcnt*2];  // value corresponding to index in restriction.
	char **rownames = new char*[rcnt];     // row labels.

	int i = 0;
	for (int color = 0; color < partition_size - 1; ++color) {
		matbeg[i] = i*2;
		matind[i*2]   = color;
		matind[i*2+1] = color + 1;
		matval[i*2]   = -1;
		matval[i*2+1] = 1;

		rhs[i] = 0;
		sense[i] = 'L';
		rownames[i] = new char[40];
		sprintf(rownames[i], "%s", "symmetry_breaker");

		++i;
	}


	// add restriction
	int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);
	checkStatus(env, status);

	// free memory
	for (int i = 0; i < rcnt; ++i) {
		delete[] rownames[i];
	}
	
	delete[] rhs;
	delete[] sense;
	delete[] matbeg;
	delete[] matind;
	delete[] matval;
	delete[] rownames;

	return 0;
}

int loadCuttingPlanes(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList, int iterations, int load_limit, int select_cuts) {

	printf("Finding Cutting Planes.\n");

	// calculate runtime
	double inittime, endtime;
	int status = CPXgettime(env, &inittime);

	int n = partition_size + (vertex_size*partition_size);

	double *sol = new double[n];
	int i = 1;
	int unsatisfied_restrictions = 0;
	while (i <= iterations) {

		printf("Iteration %d\n", i);

		// solve LP
		status = CPXlpopt(env, lp);
		checkStatus(env, status);

		status = CPXgetx(env, lp, sol, 0, n - 1);
		checkStatus(env, status);

		// print relaxation result
		// for (int id = 1; id <= vertex_size; ++id) {
		// 	for (int color = 1; color <= partition_size; ++color) {
		// 		int index = getVertexIndex(id, color, partition_size);
		// 		if (sol[index] == 0) continue;
		// 		cout << "x" << id << "_" << color << " = " << sol[index] << endl;
		// 	}
		// }

		if (select_cuts == 0 || select_cuts == 2) unsatisfied_restrictions += maximalCliqueFamillyHeuristic(env, lp, vertex_size, edge_size, partition_size, adjacencyList, sol, load_limit);
		if (select_cuts == 1 || select_cuts == 2) unsatisfied_restrictions += oddholeFamillyHeuristic(env, lp, vertex_size, edge_size, partition_size, adjacencyList, sol, load_limit);

		if (unsatisfied_restrictions == 0) break;

		unsatisfied_restrictions = 0;
		++i;
	}

	status = CPXgettime(env, &endtime);
	double elapsed_time = endtime-inittime;
	cout << "Time taken to add cutting planes: " << elapsed_time << endl;

	return 0;
}

int maximalCliqueFamillyHeuristic(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList, double* sol, int load_limit) {

	printf("Generating clique familly.\n");

	int loaded = 0;

	vector<tuple<double, int, set<int> > > clique_familly;

	for (int color = 1; color <= partition_size; ++color) {

		for (int id = 1; id <= vertex_size; id++) {

			if (sol[getVertexIndex(id, color, partition_size)] == 0) continue;

			double sum = sol[getVertexIndex(id, color, partition_size)];
			set<int> clique;
			clique.insert(id);
			for (int id2 = id + 1; id2 <= vertex_size; ++id2) {
				if (sol[getVertexIndex(id2, color, partition_size)] == 0) continue;

				if (adyacentToAll(id2, vertex_size, adjacencyList, clique)) {
					clique.insert(id2);
					sum += sol[getVertexIndex(id2, color, partition_size)];
				}
			}
			if (clique.size() > 2 && sum > sol[color-1] + TOL) {
				if (cliqueNotContained(clique, color, clique_familly)) {
					double score = sum - sol[color-1];
					clique_familly.push_back(tuple<double, int, set<int> >(score,color,clique));
				}
			}
		}
	}

	sort(clique_familly.begin(), clique_familly.end(), greater<tuple<double, int, set<int> > >());   

	//print the familly
	for (vector<tuple<double, int, set<int> > >::const_iterator it = clique_familly.begin(); 
		it != clique_familly.end() && loaded < load_limit; ++loaded, ++it) {

		loadUnsatisfiedCliqueRestriction(env, lp, partition_size, get<2>(*it), get<1>(*it));
		cout << "Score: " << get<0>(*it) << " - ";
		for (set<int>::iterator it2 = get<2>(*it).begin(); it2 != get<2>(*it).end(); ++it2) {
			cout << *it2 << " ";
		}
		cout << endl;
	}

	printf("Loaded %d/%d unsatisfied clique restrictions! (all colors)\n", loaded, (int) clique_familly.size());

	return loaded;
}

int loadUnsatisfiedCliqueRestriction(CPXENVptr& env, CPXLPptr& lp, int partition_size, const set<int>& clique, int color) {

	int ccnt = 0;
	int rcnt = 1;
	int nzcnt = clique.size() + 1;

	double rhs = 0;
	char sense = 'L';

	int matbeg = 0;
	int* matind    = new int[clique.size() + 1];
	double* matval = new double[clique.size() +1]; 
	char **rowname = new char*[rcnt];
	rowname[0] = new char[40];
	sprintf(rowname[0], "unsatisfied_clique");

	matind[0] = color - 1;
	matval[0] = -1;

	int i = 1;
	for (set<int>::iterator it = clique.begin(); it != clique.end(); ++it) {
		matind[i] = getVertexIndex(*it, color, partition_size);
		matval[i] = 1;
		++i;
	}

	// add restriction
	int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, &rhs, &sense, &matbeg, matind, matval, NULL, rowname);
	checkStatus(env, status);

	// free memory
	delete[] matind;
	delete[] matval;
	delete rowname[0];
	delete rowname;

	return 0;
}

int oddholeFamillyHeuristic(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adjacencyList, double* sol, int load_limit) {

	printf("Generating oddhole familly.\n");

	int loaded = 0;

	vector<tuple<double, int, set<int> > > path_familly; // dif, color, path

	for (int color = 1; color <= partition_size; ++color) {

		for (int id = 1; id <= vertex_size; id++) {

			if (sol[getVertexIndex(id, color, partition_size)] == 0) continue;

			double sum = 0;
			set<int> path;
			path.insert(id);
			for (int id2 = id + 1; id2 <= vertex_size; ++id2) {
				if (sol[getVertexIndex(id2, color, partition_size)] == 0) continue;

				if (isAdyacent(*(--path.end()), id2, vertex_size, adjacencyList)) {
					path.insert(id2);
				}

			}

			while (path.size() >= 3 && (path.size() % 2 == 0 || 
				!isAdyacent(*path.begin(), *(--path.end()), vertex_size, adjacencyList))) {
				path.erase(--path.end());
			}

			for (set<int>::iterator it = path.begin(); it != path.end(); ++it) {
				sum += sol[getVertexIndex(*it, color, partition_size)];
			}

			int k = (path.size() - 1) / 2;
			if (path.size() > 2 && sum > k*sol[color-1] + TOL) {
				double score = sum - k*sol[color-1];
				path_familly.push_back(tuple<double, int, set<int> >(score,color,path));
			}
		}
	}

	sort(path_familly.begin(), path_familly.end(), greater<tuple<double, int, set<int> > >());   

	//print the familly
	for (vector<tuple<double, int, set<int> > >::const_iterator it = path_familly.begin(); 
		it != path_familly.end() && loaded < load_limit; ++loaded, ++it) {
		loadUnsatisfiedOddholeRestriction(env, lp, partition_size, get<2>(*it), get<1>(*it));
		cout << "Score: " << get<0>(*it) << " - ";
		for (set<int>::iterator it2 = get<2>(*it).begin(); it2 != get<2>(*it).end(); ++it2) {
			cout << *it2 << " ";
		}
		cout << endl;
	}

	printf("Loaded %d/%d unsatisfied oddhole restrictions! (all colors)\n", loaded, (int) path_familly.size());

	return loaded;
}

int loadUnsatisfiedOddholeRestriction(CPXENVptr& env, CPXLPptr& lp, int partition_size, const set<int>& path, int color) {

	int ccnt = 0;
	int rcnt = 1;
	int nzcnt = path.size() + 1;

	double rhs = 0;
	char sense = 'L';

	int matbeg = 0;
	int* matind    = new int[path.size() + 1];
	double* matval = new double[path.size() +1]; 
	char **rowname = new char*[rcnt];
	rowname[0] = new char[40];
	sprintf(rowname[0], "unsatisfied_oddhole");

	int k = (path.size() - 1) / 2;

	matind[0] = color - 1;
	matval[0] = -k;

	int i = 1;
	for (set<int>::iterator it = path.begin(); it != path.end(); ++it) {
		matind[i] = getVertexIndex(*it, color, partition_size);
		matval[i] = 1;
		++i;
	}

	// add restriction
	int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, &rhs, &sense, &matbeg, matind, matval, NULL, rowname);
	checkStatus(env, status);

	// free memory
	delete[] matind;
	delete[] matval;
	delete rowname[0];
	delete rowname;

	return 0;
}

int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size)  {

	// load third restriction
	int ccnt = 0;                            // new columns being added.
	int rcnt = vertex_size * partition_size; // new rows being added.
	int nzcnt = rcnt*2;                      // nonzero constraint coefficients being added.

	double *rhs = new double[rcnt];          // independent term in restrictions.
	char *sense = new char[rcnt];            // sense of restriction inequality.

	int *matbeg = new int[rcnt];             // array position where each restriction starts in matind and matval.
	int *matind = new int[rcnt*2];           // index of variables != 0 in restriction (each var has an index defined above)
	double *matval  = new double[rcnt*2];    // value corresponding to index in restriction.
	char **rownames = new char*[rcnt];       // row labels.

	int i = 0;
	for (int v = 1; v <= vertex_size; ++v) {
		for (int color = 1; color <= partition_size; ++color) {
			matbeg[i] = i*2;

			matind[i*2]   = getVertexIndex(v, color, partition_size);
			matind[i*2+1] = color-1;

			matval[i*2]   = 1;
			matval[i*2+1] = -1;

			rhs[i] = 0;
			sense[i] = 'L';
			rownames[i] = new char[40];
			sprintf(rownames[i], "color_res");

			++i;
		}
	}

	// add restriction
	int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);
	checkStatus(env, status);

	// free memory
	for (int i = 0; i < rcnt; ++i) {
		delete[] rownames[i];
	}
	
	delete[] rhs;
	delete[] sense;
	delete[] matbeg;
	delete[] matind;
	delete[] matval;
	delete[] rownames;

	return 0;
}

int solveLP(CPXENVptr& env, CPXLPptr& lp, int edge_size, int vertex_size, int partition_size) {

	printf("\nSolving MIP.\n");

	int n = partition_size + (vertex_size*partition_size); // amount of total variables

	// calculate runtime
	double inittime, endtime;
	int status = CPXgettime(env, &inittime);
	checkStatus(env, status);

	// solve LP
	status = CPXmipopt(env, lp);
	checkStatus(env, status);

	status = CPXgettime(env, &endtime);
	checkStatus(env, status);

	// check solution state
	int solstat;
	char statstring[510];
	CPXCHARptr p;
	solstat = CPXgetstat(env, lp);
	p = CPXgetstatstring(env, solstat, statstring);
	string statstr(statstring);
	if (solstat != CPXMIP_OPTIMAL && solstat != CPXMIP_OPTIMAL_TOL &&
		solstat != CPXMIP_NODE_LIM_FEAS && solstat != CPXMIP_TIME_LIM_FEAS) {
		// printf("Optimization failed.\n");
		cout << "Optimization failed: " << solstat << endl;
		exit(1);
	}

	double objval;
	status = CPXgetobjval(env, lp, &objval);
	checkStatus(env, status);

	// get values of all solutions
	double *sol = new double[n];
	status = CPXgetx(env, lp, sol, 0, n - 1);
	checkStatus(env, status);

 	int nodes_traversed = CPXgetnodecnt(env, lp);

	// write solutions to current window
	cout << "Optimization result: " << statstring << endl;
	cout << "Time taken to solve final LP: " << (endtime - inittime) << endl;
	cout << "Colors used: " << objval << endl;
	cout << "Nodes traversed: " << nodes_traversed << endl;
	for (int color = 1; color <= partition_size; ++color) {
		if (sol[color-1] == 1) {
			cout << "w_" << color << " = " << sol[color-1] << " (" << colors[color-1] << ")" << endl;
		}
	}

	for (int id = 1; id <= vertex_size; ++id) {
		for (int color = 1; color <= partition_size; ++color) {
			int index = getVertexIndex(id, color, partition_size);
			if (sol[index] == 1) {
				cout << "x_" << id << " = " << colors[color-1] << endl;
			}
		}
	}

	delete[] sol;

	return 0;
}

int convertVariableType(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size, char vtype) {

	int n = partition_size + (vertex_size*partition_size);
	int* indices = new int[n];
	char* xctype = new char[n];

	for (int i = 0; i < n; i++) {
		indices[i] = i;
		xctype[i]  = vtype;
	}
	CPXchgctype(env, lp, n, indices, xctype);

	delete[] indices;
	delete[] xctype;

	return 0;
}

int setTraversalStrategy(CPXENVptr& env, int strategy) {

	// MIP node selection strategy
	// http://www-01.ibm.com/support/knowledgecenter/SSSA5P_12.3.0/ilog.odms.cplex.help/Content/Optimization/Documentation/Optimization_Studio/_pubskel/ps_refparameterscplex2299.html

	// 0 CPX_NODESEL_DFS          Depth-first search
	// 1 CPX_NODESEL_BESTBOUND    Best-bound search; default
	// 2 CPX_NODESEL_BESTEST      Best-estimate search
	// 3 CPX_NODESEL_BESTEST_ALT  Alternative best-estimate search 

	CPXsetintparam(env, CPX_PARAM_NODESEL, strategy);

	return 0;
}

int setBranchingVariableStrategy(CPXENVptr& env, int strategy) {

	// MIP variable selection strategy
	// http://www-01.ibm.com/support/knowledgecenter/SS9UKU_12.4.0/com.ibm.cplex.zos.help/Parameters/topics/VarSel.html

	// -1 	CPX_VARSEL_MININFEAS 	   Branch on variable with minimum infeasibility
	// 0 	CPX_VARSEL_DEFAULT   	   Automatic: let CPLEX choose variable to branch on; default
	// 1 	CPX_VARSEL_MAXINFEAS 	   Branch on variable with maximum infeasibility
	// 2 	CPX_VARSEL_PSEUDO 	   	   Branch based on pseudo costs
	// 3 	CPX_VARSEL_STRONG 		   Strong branching
	// 4 	CPX_VARSEL_PSEUDOREDUCED   Branch based on pseudo reduced costs

	CPXsetintparam(env, CPX_PARAM_VARSEL, strategy);

	return 0;
}

int setBranchAndBoundConfig(CPXENVptr& env) {

	// CPLEX config
	// http://www-01.ibm.com/support/knowledgecenter/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX916.html
	
	// deactivate pre-processing
	CPXsetintparam(env, CPX_PARAM_PRESLVND, -1);
	CPXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 0);
	CPXsetintparam(env, CPX_PARAM_RELAXPREIND, 0);
	CPXsetintparam(env, CPX_PARAM_REDUCE, 0);
	CPXsetintparam(env, CPX_PARAM_LANDPCUTS, -1);

	// maximize objective function
	// CPXchgobjsen(env, lp, CPX_MAX);

	// enable/disable screen output
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);

	// set excecution limit
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600);
	
	// disable presolve
	// CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);

	// enable traditional branch and bound
	CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);

	// use only one thread for experimentation
	// CPXsetintparam(env, CPX_PARAM_THREADS, 1);

	// do not add cutting planes
	CPXsetintparam(env, CPX_PARAM_EACHCUTLIM, CPX_OFF);

	// disable gomory fractional cuts
	CPXsetintparam(env, CPX_PARAM_FRACCUTS, -1);

	// measure time in CPU time
	// CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, CPX_ON);

	return 0;
}


int checkStatus(CPXENVptr& env, int status) {
	if (status) {
		char buffer[100];
		CPXgeterrorstring(env, status, buffer);
		printf("%s\n", buffer);
		exit(1);
	}
	return 0;
}