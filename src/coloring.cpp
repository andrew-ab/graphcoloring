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

struct edge {
	int from;
	int to;

	edge(int a, int b) {
		from = a;
		to = b;
	}
};

int getVertexIndex(int id, int color, int partition_size);
inline int fromMatrixToVector(int from, int to, int edge_size);
inline bool isAdyacent(int from, int to, int edge_size, bool* adyacencyList);
bool adyacentToAll(int id, int edge_size, bool* adyacencyList, const set<int>& clique);
bool cliqueNotContained(const set<int>& clique, const set<set<int> >& clique_set);

int loadObjectiveFunction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size, char vtype);
int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, vector<edge>& edges, int edge_size, int partition_size);
int loadSingleColorInPartitionRestriction(CPXENVptr& env, CPXLPptr& lp, vector<vector<int> >& partitions, int partition_size);
int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size);

int loadCuttingPlanes(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adyacencyList);
int maximalCliqueFamillyHeuristic(set<set<int> >& clique_familly, int vertex_size, int edge_size, bool* adjacencyList);
int findUnsatisfiedCliqueRestrictions(CPXENVptr& env, CPXLPptr& lp, set<set<int> >& clique_familly, int vertex_size, int partition_size, int n, double* sol);
int loadUnsatisfiedCliqueRestriction(CPXENVptr& env, CPXLPptr& lp, int partition_size, const set<int>& clique, int color, double* sol);

int solveLP(CPXENVptr& env, CPXLPptr& lp, int edge_size, int vertex_size, int partition_size);
int convertVariableType(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int partition_size, char vtype);
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

	if (argc != 3) {
		printf("Usage: type (1,2) %s inputFile\n", argv[0]);
		exit(1);
	}

	int solver = atoi(argv[1]);

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
	FILE* fp = fopen(argv[2], "r");

	if (fp == NULL) {
		printf("Invalid input file.\n");
		exit(1);
	}

	char buf[100];
	int vertex_size, edge_size;

	vector<edge> edges;
	bool* adyacencyList;

	while (fgets(buf, sizeof(buf), fp) != NULL) {
		if (buf[0] == 'c') continue;
		else if (buf[0] == 'p') {
			sscanf(&buf[7], "%d %d", &vertex_size, &edge_size);

			int adyacency_size = edge_size*edge_size - ((edge_size+1)*edge_size/2);
			adyacencyList = new bool[adyacency_size]; // can be optimized even more with a bitfield.
			fill_n(adyacencyList, adyacency_size, false);
			// printf("vertex_size: %d, edge_size: %d \n", vertex_size, edge_size);
			// printf("Adding edges! \n");
		}
		else if (buf[0] == 'e') {
			int from, to;
			sscanf(&buf[2], "%d %d", &from, &to);
			// printf("Edge: (%d,%d) \n", from, to);
			edges.push_back(edge(from, to));
			adyacencyList[fromMatrixToVector(from, to, edge_size)] = true;
		}
	}

	// set random seed
	// srand(time(NULL));

	// asign every vertex to a partition
	int partition_size = rand() % vertex_size + 1;
	vector<vector<int> > partitions(partition_size, vector<int>());

	// warning: this procedure doesn't guarantee every partition will have an element.
	for (int i = 1; i <= vertex_size; ++i) {
		int assign_partition = rand() % partition_size;
		partitions[assign_partition].push_back(i);
	}

	// update partition_size
	for (std::vector<vector<int> >::iterator it = partitions.begin(); it != partitions.end(); ++it) {
		if (it->size() == 0) --partition_size;
	}

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

	setBranchAndBoundConfig(env);

	if (solver == 1) { // pure branch & bound
		loadObjectiveFunction(env, lp, vertex_size, partition_size, CPX_BINARY);
	} else {
		loadObjectiveFunction(env, lp, vertex_size, partition_size, CPX_CONTINUOUS);
	}

	loadAdyacencyColorRestriction(env, lp, edges, edge_size, partition_size);
	loadSingleColorInPartitionRestriction(env, lp, partitions, partition_size);
	loadAdyacencyColorRestriction(env, lp, vertex_size, partition_size);
	//loadSymmetryBreaker

	if (solver != 1) loadCuttingPlanes(env, lp, vertex_size, edge_size, partition_size, adyacencyList);
		
	// write LP formulation to file, great to debug.
	status = CPXwriteprob(env, lp, "graph.lp", NULL);
	checkStatus(env, status);

	convertVariableType(env, lp, vertex_size, partition_size, CPX_BINARY); 
	
	solveLP(env, lp, edge_size, vertex_size, partition_size);

	delete[] adyacencyList;

	return 0;
}

int getVertexIndex(int id, int color, int partition_size) {
	return partition_size + ((id-1)*partition_size) + (color-1);
}

/* since the adyacency matrix is symmetric and the diagonal is not needed, we can simply
 * store the upper diagonal and get adyacency from a list. the math is quite simple, it
 * just uses the formula for the sum of integers. ids are numbered starting from 1.
 */
inline int fromMatrixToVector(int from, int to, int edge_size) {

	// for speed, many parts of this code are commented, since by our usage we always
	// know from < to and are in range.

	// assert(from != to && from <= edge_size && to <= edge_size);

	// if (from < to)
		return from*edge_size - (from+1)*from/2 - (edge_size - to) - 1;
	// else
	// 	return to*edge_size - (to+1)*to/2 - (edge_size - from) - 1;
}

inline bool isAdyacent(int from, int to, int edge_size, bool* adyacencyList) {
	return adyacencyList[fromMatrixToVector(from, to, edge_size)];
}

bool adyacentToAll(int id, int edge_size, bool* adyacencyList, const set<int>& clique) {
	for (set<int>::iterator it = clique.begin(); it != clique.end(); ++it) {
		if (!isAdyacent(id, *it, edge_size, adyacencyList)) return false;
	}
	return true;
}

bool cliqueNotContained(const set<int>& clique, const set<set<int> >& clique_set) {
	for (set<set<int> >::iterator it = clique_set.begin(); it != clique_set.end(); ++it) {
		// by construction, sets are already ordered.
		if (includes(it->begin(), it->end(), clique.begin(), clique.end())) return false;
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
			sprintf(colnames[index], "x_%d%d", id, color);
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

int loadAdyacencyColorRestriction(CPXENVptr& env, CPXLPptr& lp, vector<edge>& edges, int edge_size, int partition_size) {

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
	for (std::vector<edge>::iterator it = edges.begin(); it != edges.end(); ++it) {
		int from = it->from;
		int to   = it->to;
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


int loadCuttingPlanes(CPXENVptr& env, CPXLPptr& lp, int vertex_size, int edge_size, int partition_size, bool* adyacencyList) {

	printf("Finding Cutting Planes.\n");

	int n = partition_size + (vertex_size*partition_size);

	set<set<int> > clique_familly;
	maximalCliqueFamillyHeuristic(clique_familly, vertex_size, edge_size, adyacencyList);

	double *sol = new double[n];
	int cutting_plane_iterations = 1;
	while (cutting_plane_iterations > 0) {

		// solve LP
		int status = CPXlpopt(env, lp);
		checkStatus(env, status);

		status = CPXgetx(env, lp, sol, 0, n - 1);
		checkStatus(env, status);

		// for (int id = 1; id <= vertex_size; ++id) {
		// 	for (int color = 1; color <= partition_size; ++color) {
		// 		int index = getVertexIndex(id, color, partition_size);
		// 		if (sol[index] == 0) continue;
		// 		cout << "x_" << id << " " << color << " = " << sol[index] << endl;
		// 	}
		// }

		// check which elements in the familly do not satisfy the inequality
		findUnsatisfiedCliqueRestrictions(env, lp, clique_familly, vertex_size, partition_size, n, sol);

		cutting_plane_iterations--;
	}

	return 0;
}

int maximalCliqueFamillyHeuristic(set<set<int> >& clique_familly, int vertex_size, int edge_size, bool* adjacencyList) {

	printf("Maximal Clique Heuristic.\n");

	for (int id = 1; id <= vertex_size; id++) {
		set<int> clique;
		clique.insert(id);
		for (int id2 = id + 1; id2 <= vertex_size; id2++) {
			if (adyacentToAll(id2, edge_size, adjacencyList, clique)) {
				clique.insert(id2);
			}
		}
		if (clique.size() > 2) {
			if (cliqueNotContained(clique, clique_familly)) {
				clique_familly.insert(clique);
			}
		}
	}

	// print the familly
	// for (set<set<int> >::iterator it = clique_familly.begin(); it != clique_familly.end(); ++it) {
	// 	cout << "Clique: ";
	// 	for (set<int>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
	// 		cout << *it2 << " ";
	// 	}
	// 	cout << endl;
	// }

	printf("Clique Familly Generated.\n");

	return 0;
}

int findUnsatisfiedCliqueRestrictions(CPXENVptr& env, CPXLPptr& lp, set<set<int> >& clique_familly, int vertex_size, int partition_size, int n, double* sol) {
	for (set<set<int> >::iterator it = clique_familly.begin(); it != clique_familly.end(); ++it) {

		for (int color = 1; color <= partition_size; ++color) {
			double sum = 0;
			for (set<int>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				double coef = sol[getVertexIndex(*it2, color, partition_size)];
				sum += sol[getVertexIndex(*it2, color, partition_size)];
			}
			if (sum > sol[color-1]) {
				// printf("Unsatisfied Clique Restriction Found!\n");
				loadUnsatisfiedCliqueRestriction(env, lp, partition_size, *it, color, sol);
			}
		}
	}

	return 0;
}

int loadUnsatisfiedCliqueRestriction(CPXENVptr& env, CPXLPptr& lp, int partition_size, const set<int>& clique, int color, double* sol) {

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
	matval[0] = -sol[0];

	int i = 1;
	for (set<int>::iterator it = clique.begin(); it != clique.end(); ++it) {
		matind[i] = getVertexIndex(*it, color, partition_size);
		matval[i] = sol[getVertexIndex(*it, color, partition_size)];
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

	printf("Solving MIP.\n");

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

	// write solutions to current window
	cout << endl << "Optimization result: " << statstring << endl;
	cout << "Runtime: " << (endtime - inittime) << endl;
	cout << "Colors used: " << objval << endl;
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
	CPXsetintparam(env, CPX_PARAM_THREADS, 1);

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