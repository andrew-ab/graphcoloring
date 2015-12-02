#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

unsigned choose(unsigned n, unsigned k);

int main(int argc, char** argv) {

	if (argc != 4) {
		cout << "Usage: " << argv[0] << " graph_type vertex_size density" << endl;
		return 0;
	}

	int type = atoi(argv[1]);
	int vertex_size = atoi(argv[2]);
	int density = atoi(argv[3]);
	int clique_size = choose(vertex_size, 2);
	int edge_size = clique_size * (density/100.0);

	/* I was looking for a data structure with efficient random access and removal.
	 * So far it seems the best option is a vector. However, I'm sure some sort of
	 * bijection can be found to improve space complexity and that there must exist
	 * a more efficient data structure to sample edges to generate a random graph.
	 * One good option was to build an array of bools, shuffle it and then use the
	 * biyeccion to generate the edges. Remember the problem is when generating graphs
	 * that are not sparse.
	 */
	vector<pair<int,int> > possible_edges;

	cout << "Type: " << type << " Vertex Size: " << vertex_size << " Density: " << density << " Edges: " << edge_size << endl;

	for (int i = 1; i <= vertex_size; ++i) {
		for (int j = i + 1; j <= vertex_size; ++j) {
			possible_edges.push_back(pair<int,int>(i,j));
		}
	}

	char filename[] = "graph";
	ofstream myfile;
	myfile.open(filename);
	myfile << "c FILE: " << filename << endl;
	myfile << "c Randomly generated graph" << endl;
	myfile << "c Density: " << density << "%" << endl;
	myfile << "p edge " << edge_size << " " << vertex_size << endl;

	for (int i = 0; i < edge_size; ++i) {
		int index = rand() % possible_edges.size();
		pair<int, int>& p = possible_edges[index];
		myfile << "e " << p.first << " " << p.second << endl;
		swap(p, possible_edges.back());
		possible_edges.pop_back();
	}

	myfile.close();

	return 0;
}

unsigned choose(unsigned n, unsigned k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}