#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

#include <stdlib.h>

#include <string>
#include <vector>

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

int main(int argc, char **argv) {

	if (argc != 2) {
		printf("Usage: %s inputFile\n", argv[0]);
		exit(1);
	}

	/* read graph input file
	 * format: http://mat.gsia.cmu.edu/COLOR/instances.html
	 * graph representation chosen in order to load the LP easily.
	 * - vector of edges
	 * - vector of partitions
	 */
	FILE* fp = fopen(argv[1], "r");

	if (fp == NULL) {
		printf("Invalid input file. \n");
		exit(1);
	}

	char buf[100];
	int vertex_size, edge_size;

	vector<edge> edges;

	while (fgets(buf, sizeof(buf), fp) != NULL) {
		if (buf[0] == 'c') continue;
		else if (buf[0] == 'p') {
			sscanf(&buf[7], "%d %d", &vertex_size, &edge_size);
			printf("vertex_size: %d, edge_size: %d \n", vertex_size, edge_size);
			// printf("Adding edges! \n");
		}
		else if (buf[0] == 'e') {
			int from, to;
			sscanf(&buf[2], "%d %d", &from, &to);
			// printf("Edge: (%d,%d) \n", from, to);
			edges.push_back(edge(from, to));
		}
	}

	// set random seed
 	srand(time(NULL));

	// asign every vertex to a partition
	int partition_size = rand() % vertex_size;

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
	
	// start loading LP using CPLEX
	int status;
	CPXENVptr env; // pointer to enviroment
	CPXLPptr lp;	 // pointer to the lp.

	env = CPXopenCPLEX(&status); // create enviroment

	if (env == NULL) {
		printf("Error creating enviroment.\n");
		exit(1);
	}

	// create LP
	lp = CPXcreateprob(env, &status, "Instance of partitioned graph coloring.");

	if (lp == NULL) {
		printf("Error creating the LP.\n");
		exit(1);
	}

	// load objective function
	int n = partition_size + (vertex_size*partition_size);
	double *objfun	= new double[n];
	char	 *ctype	 = new char[n];
	char **colnames = new char*[n];

	for (int i = 0; i < partition_size; ++i) {
		objfun[i] = 1;
		ctype[i]	= 'B';
		colnames[i] = new char[10];
		sprintf(colnames[i], "w_%d", (i+1));
	}

	for (int id = 1; id <= vertex_size; ++id) {
		for (int color = 1; color <= partition_size; ++color) {
			int index = getVertexIndex(id, color, partition_size);
			objfun[index]   = 0;
			ctype[index]    = 'B';
			colnames[index] = new char[10];
			sprintf(colnames[index], "x_%d%d", id, color);
			// cout << colnames[index] << " (" << index << ")" << endl;
		}
	}

	status = CPXnewcols(env, lp, n, objfun, NULL, NULL, ctype, colnames);

	if (status) {
		printf("Problem adding variables with CPXnewcols.\n");
		exit(1);
	}
	
	// free memory
	for (int i = 0; i < n; ++i) {
		delete[] colnames[i];
	}
	
	delete[] objfun;
	delete[] ctype;
	delete[] colnames;

	// colors array!
	const char* colors[] = {"Blue", "Red", "AliceBlue","AntiqueWhite","Aqua","Aquamarine","Azure","Beige","Bisque","Black","BlanchedAlmond",
	"BlueViolet","Brown","BurlyWood","CadetBlue","Chartreuse","Chocolate","Coral","CornflowerBlue","Cornsilk","Crimson",
	"Cyan","DarkBlue","DarkCyan","DarkGoldenRod","DarkGray","DarkGrey","DarkGreen","DarkKhaki","DarkMagenta","DarkOliveGreen",
	"Darkorange","DarkOrchid","DarkRed","DarkSalmon","DarkSeaGreen","DarkSlateBlue","DarkSlateGray","DarkSlateGrey","DarkTurquoise",
	"DarkViolet","DeepPink","DeepSkyBlue","DimGray","DimGrey","DodgerBlue","FireBrick","FloralWhite","ForestGreen","Fuchsia",
	"Gainsboro","GhostWhite","Gold","GoldenRod","Gray","Grey","Green","GreenYellow","HoneyDew","HotPink","IndianRed","Indigo",
	"Ivory","Khaki","Lavender","LavenderBlush","LawnGreen","LemonChiffon","LightBlue","LightCoral","LightCyan","LightGoldenRodYellow",
	"LightGray","LightGrey","LightGreen","LightPink","LightSalmon","LightSeaGreen","LightSkyBlue","LightSlateGray","LightSlateGrey",
	"LightSteelBlue","LightYellow","Lime","LimeGreen","Linen","Magenta","Maroon","MediumAquaMarine","MediumBlue","MediumOrchid",
	"MediumPurple","MediumSeaGreen","MediumSlateBlue","MediumSpringGreen","MediumTurquoise","MediumVioletRed","MidnightBlue",
	"MintCream","MistyRose","Moccasin","NavajoWhite","Navy","OldLace","Olive","OliveDrab","Orange","OrangeRed","Orchid",
	"PaleGoldenRod","PaleGreen","PaleTurquoise","PaleVioletRed","PapayaWhip","PeachPuff","Peru","Pink","Plum","PowderBlue",
	"Purple","RosyBrown","RoyalBlue","SaddleBrown","Salmon","SandyBrown","SeaGreen","SeaShell","Sienna","Silver","SkyBlue",
	"SlateBlue","SlateGray","SlateGrey","Snow","SpringGreen","SteelBlue","Tan","Teal","Thistle","Tomato","Turquoise","Violet",
	"Wheat","White","WhiteSmoke","Yellow","YellowGreen"};

	// load first restriction
	{
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
	status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);

	if (status) {
		printf("Problem adding restriction with CPXaddrows.\n");
		exit(1);
	}

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
	}

	// load second restriction
	{
	printf("partition_size: %d\n", partition_size);
	int p = 1;
	for (std::vector<vector<int> >::iterator it = partitions.begin(); it != partitions.end(); ++it) {

		int size = it->size();                 // current partition size.
		// if (size == 0) continue;

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
		sense[0]  = 'L';
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
		status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);

		if (status) {
			printf("Problem adding restriction with CPXaddrows.\n");
			exit(1);
		}

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
	}

	// load third restriction
	{
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
	status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, matbeg, matind, matval, NULL, rownames);

	if (status) {
		printf("Problem adding restriction with CPXaddrows.\n");
		exit(1);
	}

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
	}

	// CPLEX by default minimizes the objective function. Just in case you want to maximize.
	// CPXchgobjsen(env, lp, CPX_MAX);

	// enable screen output
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
		
	if (status) {
		printf("Problem setting CPX_PARAM_SCRIND\n");
		exit(1);
	}
		
	// set excecution limit
	status = CPXsetdblparam(env, CPX_PARAM_TILIM, 3600);
	
	if (status) {
		cerr << "Problema seteando el tiempo limite" << endl;
		exit(1);
	}
 
	// write LP formulation to file
	status = CPXwriteprob(env, lp, "graph.lp", NULL);
		
	if (status) {
		printf("Problem writing LP problem to file.");
		exit(1);
	}
		
	// calculate runtime
	double inittime, endtime;
	status = CPXgettime(env, &inittime);

	// solve LP
	status = CPXlpopt(env, lp);

	status = CPXgettime(env, &endtime);

	if (status) {
		printf("Optimization problem.\n");
		exit(1);
	}

	// check solution state
	int solstat;
	char statstring[510];
	CPXCHARptr p;
	solstat = CPXgetstat(env, lp);
	p = CPXgetstatstring(env, solstat, statstring);
	string statstr(statstring);
	cout << endl << "Optimization result: " << statstring << endl;
	if (solstat != CPX_STAT_OPTIMAL){
		 exit(1);
	}	
		
	double objval;
	status = CPXgetobjval(env, lp, &objval);
		
	if (status) {
		printf("Problem obtaining optimal solution.\n");
		exit(1);
	}
		
	cout << "Solution: " << "\t" << objval << "\t" << (endtime - inittime) << endl; 

	// write solutions to file
	std::string outputfile = "graph.sol";
	ofstream solfile(outputfile.c_str());

	// get values of all solutions
	double *sol = new double[n];
	status = CPXgetx(env, lp, sol, 0, n - 1);

	if (status) {
		printf("Problem obtaining the solution of the LP.\n");
		exit(1);
	}

	// Write solutions different than 0 (tolerance: 1e-05).
	solfile << "Solution status: " << statstr << endl;
	for (int i = 0; i < n; i++) {
		if (sol[i] > TOL) {
			solfile << "x_" << i << " = " << sol[i] << endl;
		}
	}

	delete[] sol;
	solfile.close();

	return 0;
}

int getVertexIndex(int id, int color, int partition_size) {
	return partition_size + ((id-1)*partition_size) + (color-1);
}