#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>

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
			printf("Adding edges! \n");
		}
		else if (buf[0] == 'e') {
			int from, to;
			sscanf(&buf[2], "%d %d", &from, &to);
			printf("Edge: (%d,%d) \n", from, to);
			edges.push_back(edge(from, to));
		}
	}

	// asign every vertex to a partition
	int partition_size = 1 + (rand() % vertex_size);
	vector<vector<int> > partitions(partition_size, vector<int>(vertex_size/partition_size));

	// warning: this procedure doesn't guarantee every partition will have an element.
	for (int i = 1; i <= vertex_size; ++i) {
		int assign_partition = rand() % partition_size;
		partitions[assign_partition].push_back(i);
	}

	// start loading LP using CPLEX
	int status;
	CPXENVptr env; // pointer to enviroment
	CPXLPptr lp;   // pointer to the lp.

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
	double *objfun  = new double[partition_size];
	char   *ctype   = new char[partition_size];
	char **colnames = new char*[partition_size];

	for (int i = 0; i < partition_size; i++) {
		objfun[i] = 1;
		ctype[i]  = 'B';
		colnames[i] = new char[10];
		sprintf(colnames[i], "w_%d", (i+1));
	}

	status = CPXnewcols(env, lp, partition_size, objfun, NULL, NULL, ctype, colnames);

	if (status) {
		printf("Problem adding variables with CPXnewcols.\n");
		exit(1);
	}
	
	// free memory
	for (int i = 0; i < partition_size; ++i) {
		delete[] colnames[i];
	}
	
	delete[] objfun;
	delete[] ctype;
	delete[] colnames;


	// CPLEX by default minimizes the objective function. Just in case you want to maximize.
	// CPXchgobjsen(env, lp, CPX_MAX);

	// write LP formulation to file
	status = CPXwriteprob(env, lp, "graph.lp", NULL);
		
	if (status) {
		printf("Problem writing LP problem to file.");
		exit(1);
	}

	return 0;
}
