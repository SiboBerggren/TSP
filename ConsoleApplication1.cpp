#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <iostream>     
#include <algorithm>    
#include <vector>  
#include <map>

// these can namespace usages can be ignored if reader is unfamiliar with C++
using namespace std;
using namespace std::chrono;

// constants used by the 3-opt
int THREE_OPT_SIT_TWO_OPT = 1;
int THREE_OPT_SIT_A = 2;
int THREE_OPT_SIT_B = 3;

// we store these as static variables so they don't have to passed around by all functions
int n = 0; // number of nodes in graph G
double x[1000]; // x-coords of nodes
double y[1000]; // y-coords of nodees
int dist[1000][1000]; // matrix that will contain distance between each pair of nodes
vector<int> all_nodes; // utility vector will contain {0,1,...,n-1}
int num_edges; // utility int that holds the number of edges in G, i.e., n(n-1)/2.

int D[65536][16]; // matrix used by dynamic programming construction algorithm
int D_opt[65536][16]; // matrix used by dynamic programming optimisation

int tour[1000] = { 0 }; // array for storing a tour in G
int min_tour[1000]; // array for storing a the shortest tour in G we find when running algorithms multiple times
int min_sum = -1; // the length (or sum of distances) of min_tour

vector<pair<int,int>> edges; // used by some optimisation algorithms, will contain the edges of G

vector<vector<pair<int,int>>> closest_neighbours; // this will be used by optimisation algorithms to look up closest neighbours of a node. for 
// sorting reasons this is a vector<vector<pair<int,int>>> and not a vector<vector<int>>. each vector<pair<int,int>> are the closest neighbours
// of the first int in the pair, the second int in the pair is a neighbour.
vector<vector<pair<int, int>>> closest_neighbours_two_opt; // same thing as above, but here we only consider neighbours of higher node index. this is used by 2-opt 
vector<vector<pair<int, int>>> closest_neighbours_quadrants; // in this one we store the 10 closest neighbours from each quadrant around a node. this
// is used in an optimised version of 3-opt adapted especially for a 2D-plane, see report for more information.

int hub; // used by Clarke-Wright function, must be a static int since it is also used by sorting functions.

// entry in a (single) linked list
class Entry {
public:
	int node = -1;
	Entry *next = NULL;
};

// entry in a (double) linked list
class DoubleEntry {
public:
	int node = -1;
	DoubleEntry *prev = NULL;
	DoubleEntry *next = NULL;
};

// read G from input stream as specified on Kattis
void readInput() {
	cin >> n;
	for (int i = 0; i < n; i++)
		cin >> x[i] >> y[i];
}

// used for debugging, place a file "in.txt" with input as specified on Kattis in same folder as code
void readInputFromFile() {
	ifstream fin("in.txt");
	fin >> n;
	for (int i = 0; i < n; i++)
		fin >> x[i] >> y[i];
}

// populate dist-matrix
void computeDistances() {
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			dist[i][j] = round(sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2)));
			dist[j][i] = dist[i][j];
		}
	}
}

// updates min_tour[] to tour[] if tour has shorter length (sum of edges) than min_tour
void UpdateMinTour() {
	int sum = 0;

	for (int i = 1; i < n; i++) {
		sum += dist[tour[i - 1]][tour[i]];
	}
	sum += dist[tour[0]][tour[n - 1]];

	if (min_sum == -1 || sum < min_sum) {
		min_sum = sum;
		for (int i = 0; i < n; i++) {
			min_tour[i] = tour[i];
		}
	}
}

// write the best tour we found min_tour to output stream
void OutputMinTour() {
	for (int i = 0; i < n; i++)
		cout << min_tour[i] << '\n';
}

// populate utility vector all_nodes
void addAllNodes() {
	all_nodes.reserve(n);
	for (int i = 0; i < n; i++)
		all_nodes.push_back(i);
}

// return true if and only if edge e1 is shorter than edge e2
bool CompareEdges(pair<int,int> &e1, pair<int, int> &e2) {
	return dist[e1.first][e1.second] < dist[e2.first][e2.second];
}

// populate and sort edges vector
void sortEdges() {
	edges.reserve(num_edges);

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			edges.push_back(make_pair(i,j));
		}
	}

	sort(begin(edges), begin(edges) + num_edges, CompareEdges);
}

// sort the neighbours of each node in increasing order by distance
void computeClosestNeighbours() { 
	for (int i = 0; i < n; i++) {
		vector<pair<int, int>> closest_neighbours_of_i;
		for (int j = 0; j < n; j++) { // this can be made faster with array copy, but does it matter for runtime overall?
			if (j != i) {
				closest_neighbours_of_i.push_back(make_pair(i, j));
			}
		}
		sort(closest_neighbours_of_i.begin(), closest_neighbours_of_i.end(), CompareEdges);
		closest_neighbours.push_back(closest_neighbours_of_i);
	}
}

// sort the neighbours of each node in increasing order by distance
// only neighbours of higher node index are considered, this is used specifically for the 2-opt
void computeClosestNeighboursTwoOpt() { 
	for (int i = 0; i < n; i++) {
		vector<pair<int, int>> closest_neighbours_of_i;
		for (int j = i + 1; j < n; j++) {
			closest_neighbours_of_i.push_back(make_pair(i, j));
		}
		sort(closest_neighbours_of_i.begin(), closest_neighbours_of_i.end(), CompareEdges);
		closest_neighbours_two_opt.push_back(closest_neighbours_of_i);
	}
}

// this function tells us which half plane (0 or 1) a node with coordinate coord lies in
inline int coord_to_half_plane(int coord) { // returns 1 if coord > 0, and 0 otherwise
	return coord > 0 ? 1 : 0;
}

// sort the neighbours of each node in increasing order by distance
// in this computeClosestNeighbours variant: for every node i we pick some number of neighbours from each quadrant around i
void computeClosestNeighboursQuadrants() {
	for (int i = 0; i < n; i++) {
		vector<pair<int, int>> closest_neighbours_of_i;
		for (int j = 0; j < n; j++) {
			if (j != i) {
				closest_neighbours_of_i.push_back(make_pair(i, j));
			}
		}
		sort(closest_neighbours_of_i.begin(), closest_neighbours_of_i.end(), CompareEdges);

		vector<pair<int, int>> closest_neighbours_of_i_from_quads;
		int i_x = x[i]; // x-coord of node i
		int i_y = y[i]; // y-coord of node i
		int nodes_left_from_quad[2][2]; // pick (at most) 8 close neighbours from each quadrant AROUND i, these are the quadrants when we think of i's coordinates as the origin.
		nodes_left_from_quad[0][0] = 8;
		nodes_left_from_quad[0][1] = 8;
		nodes_left_from_quad[1][0] = 8;
		nodes_left_from_quad[1][1] = 8;
		for (int j = 0; j < closest_neighbours_of_i.size(); j++) {
			int neighbour = closest_neighbours_of_i[j].second;
			int x_half_plane = coord_to_half_plane(i_x - x[neighbour]); // x_half_plane and y_half_plane combined tells us which quadrant around i that neighbour lies in.
			int y_half_plane = coord_to_half_plane(i_y - y[neighbour]);

			if (nodes_left_from_quad[x_half_plane][y_half_plane] > 0) {
				closest_neighbours_of_i_from_quads.push_back(make_pair(i, neighbour)); // (for code consistency we create a pair as we do in other computeClosestNeighbours-functions, if this comment is confusing: ignore it.)

				nodes_left_from_quad[x_half_plane][y_half_plane]--;
				if (nodes_left_from_quad[x_half_plane][y_half_plane] == 0 && // do we have all the neighbours we want?
					nodes_left_from_quad[0][0] == 0 &&
					nodes_left_from_quad[0][1] == 0 &&
					nodes_left_from_quad[1][0] == 0 &&
					nodes_left_from_quad[1][1] == 0) {
					break;
				}
			}
		}

		closest_neighbours_quadrants.push_back(closest_neighbours_of_i_from_quads);
	}
}

// *** CONSTRUCTION ALGOS

// standard nearest neighbour from kattis description (called greedy there).
void NearestNeighbour(int start_node) {
	vector<bool> used(n, false); // keeps track of which nodes we have used so far in tour, set all values to false initially.

	tour[0] = start_node;
	used[start_node] = true;
	for (int i = 1; i < n; i++) {
		int closest_neighbour = -1; // this will be the closest neighbour (that is not already used) of previous node (i.e. tour[i-1]) on tour
		for (int j = 0; j < n; j++) {
			if (used[j])
				continue;
			if (closest_neighbour == -1 || dist[tour[i - 1]][j] < dist[tour[i - 1]][closest_neighbour])
				closest_neighbour = j;
		}

		tour[i] = closest_neighbour;
		used[closest_neighbour] = true;
	}
}

// Backtracking optimal solution from minimum, in standard dynamic programming fashion
void DynprogBacktrack(int subset, int tour_idx) { 
	if (subset == 0) // are we done?
		return;

	int last = tour[tour_idx + 1]; // note: the first time this function is called, last will be n-1.
	// we want to find a minimum path from node n-1 that ends in last. all nodes in subset should be on the path.
	// what we do below is find the node (min_node) that is before last on a minimum path.
	int min = -1;
	int min_node = -1;
	int c_flag = 1; // as in DynProg(), c_flag encodes c.
	for (int c = 0; c_flag <= subset; c++, c_flag = c_flag << 1) { // note that c can only be in subset if c_flag <= subset
		if ((subset & c_flag) == 0) // continue if c is not in subset
			continue;

		if (min == -1 || D[subset][c] + dist[c][last] < min) {
			min = D[subset][c] + dist[c][last];
			min_node = c;
		}
	}

	tour[tour_idx] = min_node;
	DynprogBacktrack(subset ^ (1 << min_node), --tour_idx); // find minimum path that uses the rest of the nodes in subset and ends in min_node...
}

// This function finds a minimum (i.e. optimal) tour (that starts in node n-1)
void DynProg() {
	// the int set_of_all_nodes, whose binary encoding is 111...1 (n-1 number of 1's), encodes the set of all nodes {0,...,n-2}. element c is encoded by 2^c.
	int set_of_all_nodes = pow(2, n - 1) - 1; 

	// loop over all subsets encoded by i
	for (int subset = 1; subset <= set_of_all_nodes; subset++) {
		// try removing each element c from subset
		int c_flag = 1; // c_flag encodes the node c, i.e. c_flag = 2^c.
		for (int c = 0; c_flag <= subset; c++, c_flag = c_flag << 1) { // note that c can only be in subset if c_flag <= subset
			if ((subset & c_flag) == 0) // continue if c is not in subset
				continue;

			if (subset == c_flag) { // base case: subset = {c}
				D[subset][c] = dist[n-1][c];
				break;
			}

			// remove c from subset
			int subset_minus_c = subset ^ c_flag;
			// find min over all nodes x in subset_minus_c
			int min = -1;
			int x_flag = 1; // x_flag encodes the node x, i.e. x_flag = 2^x.
			for (int x = 0; x_flag <= subset_minus_c; x++, x_flag = x_flag << 1) { // note that x can only be in subset_minus_c if x_flag <= subset_minus_c
				if ((x_flag & subset_minus_c) == 0) // continue if x is not in subset
					continue;

				if (min == -1 || D[subset_minus_c][x] + dist[x][c] < min)
					min = D[subset_minus_c][x] + dist[x][c];
			}
			D[subset][c] = min;
		}
	}

	// find optimal tour by backtracking
	tour[n-1] = n-1;
	DynprogBacktrack(set_of_all_nodes, n - 2);

	UpdateMinTour();
}

// greedy algorithm that starts with a tour consisting of 0 and node. Then inserts each node (into the tour) where the node adds the least to 
// the cycle's total length.
void GreedyTour(int node) {
	// allocate n number of linked list entries. we will use a linked list to represent the tour. entries(i) is entry for node i.
	vector<Entry> entries(n);

	// place 0 and node on tour to start with
	Entry* zero_entry = &entries[0];
	Entry* node_entry = &entries[node];
	zero_entry->node = 0;
	node_entry->node = node;
	zero_entry->next = node_entry;
	node_entry->next = zero_entry;

	// loop over nodes and insert them where they add the least to the length of the current tour
	for (int k = 1; k < n; k++) {
		int i = all_nodes[k];
		if (i == node_entry->node) // make sure we do not add node twice
			continue;

		// find best place on tour to insert node i at
		int min = -2; // this will contain the minimum cost for which we can insert i. note that cost of inserting i at some place can actually be -1 due to broken triangle inequality.
		Entry *min_entry = NULL; // min_entry is the entry on the tour that we should insert node i after at min cost
		Entry *entry = zero_entry;
		do {
			// try fitting node i between entry and entry->next
			int cost = dist[entry->node][i] + dist[entry->next->node][i] - dist[entry->node][entry->next->node]; // cost of inserting i between entry and entry->next
			if (min == -2 || cost < min) {
				min = cost;
				min_entry = entry;
			}

			entry = entry->next;
		} while (entry != zero_entry);

		// insert node i after min_entry on tour
		entries[i].node = i;
		entries[i].next = min_entry->next;
		min_entry->next = &entries[i];
	}

	// copy linked list to tour-array
	Entry *entry = zero_entry;
	int tour_idx = 0;
	do {
		tour[tour_idx++] = entry->node;
		entry = entry->next;
	} while (entry != zero_entry);
}

// greedy algorithm that starts with a tour consisting of 0 and a random node. Then inserts each node (into the tour) where the node adds the least to 
// the cycle's total length. The nodes, considered for insertion, are considered in a random order.
void GreedyTourRandom() {
	// allocate n number of linked list entries. we will use a linked list to represent the tour. entries(i) is entry for node i.
	vector<Entry> entries(n);

	// place 0 and a random node on tour to start with
	Entry* zero_entry = &entries[0];
	int rand_node = (rand() % (n - 1)) + 1; // random node in {1,...,n-1}
	Entry* rand_node_entry = &entries[rand_node];
	zero_entry->node = 0;
	rand_node_entry->node = rand_node;
	zero_entry->next = rand_node_entry;
	rand_node_entry->next = zero_entry;

	// create a vector of all nodes (except 0) in shuffled order
	vector<int> all_nodes_copy(all_nodes.begin()+1, all_nodes.end());
	random_shuffle(all_nodes_copy.begin(), all_nodes_copy.end());

	// loop over nodes and insert them where they add the least to the length of the current tour
	for (int k = 0; k < n-1; k++) {
		int i = all_nodes_copy[k];
		if (i == rand_node) // make sure we do not add rand_node twice
			continue;

		// find best place on tour to insert node i at
		int min = -2; // this will contain the minimum cost for which we can insert i. note that cost of inserting i at some place can actually be -1 due to broken triangle inequality.
		Entry *min_entry = NULL; // min_entry is the entry on the tour that we should insert node i after at min cost
		Entry *entry = zero_entry;
		do {
			// try fitting node i between entry and entry->next
			int cost = dist[entry->node][i] + dist[entry->next->node][i] - dist[entry->node][entry->next->node]; // cost of inserting i between entry and entry->next
			if (min == -2 || cost < min) {
				min = cost;
				min_entry = entry;
			}

			entry = entry->next;
		} while (entry != zero_entry);

		// insert node i after min_entry on tour
		entries[i].node = i;
		entries[i].next = min_entry->next;
		min_entry->next = &entries[i];
	}

	// copy linked list to tour-array
	Entry *entry = zero_entry;
	int tour_idx = 0;
	do {
		tour[tour_idx++] = entry->node;
		entry = entry->next;
	} while (entry != zero_entry);
}

void Greedy() {
	// we will greedily add all edges of the tour. we will save these edges as a collection of paths since any node will have at most two neighbours at any point of this process.
	vector<vector<int>> neighbours(n); // holds the neighbours of each node.
	int num_nodes_of_deg_2 = 0; // when all nodes except two have deg 2 we can end for loop over edges.
	vector<int> endpoint_to_endpoint(n); // lookup table from an endpoint of a path to the other endpoint. If a node is the endpoint of a path, then we can use this table to find the other endpoint.
	for (int i = 0; i < n; i++)
		endpoint_to_endpoint[i] = -1; // -1 represents that a node is not the endpoint of a path
	// loop over all sorted edges {v,u}
	for (int i = 0; i < edges.size(); i++) {
		int v = edges[i].first;
		int u = edges[i].second;
		int v_degree = neighbours[v].size();
		int u_degree = neighbours[u].size();

		if (v_degree == 2 || u_degree == 2)
			continue;

		if (v_degree == 0 && u_degree == 0) {
			// add new path {v,u}
			neighbours[v].push_back(u);
			neighbours[u].push_back(v);
			endpoint_to_endpoint[v] = u;
			endpoint_to_endpoint[u] = v;
			continue;
		}

		// at least one of v and u have exactly degree 1. for what follows, make sure v has degree 1.
		if (v_degree == 0) {
			int tmp = v;
			v = u;
			u = tmp;
			u_degree = 0;
			v_degree = 1;
		}
		if (u_degree == 1) { // both u and v have degree 1
			// make sure u and v are not endpoints of the same path as that would create a proper subcycle 
			// (actually it will not create a proper subcycle if all nodes are on a single big path, but we add this last edge afterwards instead)
			if (endpoint_to_endpoint[v] == u)
				continue;
			// connect paths that have u and v as endpoints
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = endpoint_to_endpoint[u];
			endpoint_to_endpoint[endpoint_to_endpoint[u]] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
			endpoint_to_endpoint[u] = -1;
			// u will have deg 2 after this
			num_nodes_of_deg_2++;
		}
		else { // only v has degree 1
			// u is new endpoint of path that v was formerly endpoint of
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = u;
			endpoint_to_endpoint[u] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
		}
		neighbours[v].push_back(u);
		neighbours[u].push_back(v);
		// v will have deg 2 after this
		num_nodes_of_deg_2++;
		if (num_nodes_of_deg_2 == n - 2)
			break;
	}

	// save the tour
	// 1. find an endpoint of the single path left
	int endpoint;
	for (int i = 0; i < n; i++)
		if (endpoint_to_endpoint[i] != -1) {
			endpoint = i;
			break;
		}
	// 2. walk through tour and save it to tour-array
	tour[0] = endpoint;
	int node = neighbours[endpoint][0]; // walk to only neighbour of endpoint
	int previous_node = endpoint; // this is used to make sure we don't go backwards.
	for (int i = 1; i < n - 1; i++) {
		tour[i] = node;
		int next_node = neighbours[node][0] == previous_node ? neighbours[node][1] : neighbours[node][0]; // step forward
		previous_node = node;
		node = next_node;
	}
	tour[n - 1] = node; // node is now other endpoint
}

void GreedyRandomShuffle() {
	// randomisation: we do a (slight) shuffle of the sorted edges. loop over the edges and with probability 1/3 swap edge i with edge i+1.
	// since edges is used over multiple runs of this algorithm we first need to copy edges to edges_copy.
	vector<pair<int, int>> edges_copy(edges);
	for (int i = 0; i < edges_copy.size() - 1; i++)
		if (i < edges_copy.size() - 2 && rand() % 3 == 2)
			swap(edges_copy[i], edges_copy[i + 1]);

	// we will now greedily add all edges of the tour. we will save these edges as a collection of paths since any node will have at most two neighbours at any point of this process.
	vector<vector<int>> neighbours(n); // holds the neighbours of each node.
	int num_nodes_of_deg_2 = 0; // when all nodes except two have deg 2 we can end for loop over edges.
	vector<int> endpoint_to_endpoint(n); // lookup table from an endpoint of a path to the other endpoint. If a node is the endpoint of a path, then we can use this table to find the other endpoint.
	for (int i = 0; i < n; i++)
		endpoint_to_endpoint[i] = -1; // -1 represents that a node is not the endpoint of a path
	// loop over all sorted edges {v,u}
	for (int i = 0; i < edges_copy.size(); i++) {
		int v = edges_copy[i].first;
		int u = edges_copy[i].second;
		int v_degree = neighbours[v].size();
		int u_degree = neighbours[u].size();

		if (v_degree == 2 || u_degree == 2)
			continue;

		if (v_degree == 0 && u_degree == 0) {
			// add new path {v,u}
			neighbours[v].push_back(u);
			neighbours[u].push_back(v);
			endpoint_to_endpoint[v] = u;
			endpoint_to_endpoint[u] = v;
			continue;
		}

		// at least one of v and u have exactly degree 1. for what follows, make sure v has degree 1.
		if (v_degree == 0) {
			int tmp = v;
			v = u;
			u = tmp;
			u_degree = 0;
			v_degree = 1;
		}
		if (u_degree == 1) { // both u and v have degree 1
			// make sure u and v are not endpoints of the same path as that would create a proper subcycle
			// (actually it will not create a proper subcycle if all nodes are on a single big path, but we add this last edge afterwards instead)
			if (endpoint_to_endpoint[v] == u)
				continue;
			// connect paths that have u and v as endpoints
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = endpoint_to_endpoint[u];
			endpoint_to_endpoint[endpoint_to_endpoint[u]] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
			endpoint_to_endpoint[u] = -1;
			// u will have deg 2 after this
			num_nodes_of_deg_2++;
		}
		else { // only v has degree 1
			// u is new endpoint of path that v was formerly endpoint of
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = u;
			endpoint_to_endpoint[u] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
		}
		neighbours[v].push_back(u);
		neighbours[u].push_back(v);
		// v will have deg 2 after this
		num_nodes_of_deg_2++;
		if (num_nodes_of_deg_2 == n - 2)
			break;
	}

	// save the tour
	// 1. find an endpoint of the single path left
	int endpoint;
	for (int i = 0; i < n; i++)
		if (endpoint_to_endpoint[i] != -1) {
			endpoint = i;
			break;
		}
	// 2. walk through tour and save it to tour-array
	tour[0] = endpoint;
	int node = neighbours[endpoint][0]; // walk to only neighbour of endpoint
	int previous_node = endpoint; // this is used to make sure we don't go backwards.
	for (int i = 1; i < n - 1; i++) {
		tour[i] = node;
		int next_node = neighbours[node][0] == previous_node ? neighbours[node][1] : neighbours[node][0]; // step forward
		previous_node = node;
		node = next_node;
	}
	tour[n - 1] = node; // node is now other endpoint
}

void GreedyRandomNode(int rand_node) {
	// we will greedily add all edges of the tour. we will save these edges as a collection of paths since any node will have at most two neighbours at any point of this process.
	vector<vector<int>> neighbours(n); // holds the neighbours of each node.
	int num_nodes_of_deg_2 = 0; // when all nodes except two have deg 2 we can end for loop over edges.
	vector<int> endpoint_to_endpoint(n); // lookup table from an endpoint of a path to the other endpoint. If a node is the endpoint of a path, then we can use this table to find the other endpoint.
	for (int i = 0; i < n; i++)
		endpoint_to_endpoint[i] = -1; // -1 represents that a node is not the endpoint of a path

	// first, add the shortest edge for node. this is what enables this algorithm to be randomized.
	int closest_neighbour = closest_neighbours[rand_node][0].second;
	neighbours[rand_node].push_back(closest_neighbour);
	neighbours[closest_neighbour].push_back(rand_node);
	endpoint_to_endpoint[rand_node] = closest_neighbour;
	endpoint_to_endpoint[closest_neighbour] = rand_node;

	// now add the rest of the edges of the tour.
	// loop over all sorted edges {v,u}
	for (int i = 0; i < edges.size(); i++) {
		int v = edges[i].first;
		int u = edges[i].second;
		int v_degree = neighbours[v].size();
		int u_degree = neighbours[u].size();

		if (v_degree == 2 || u_degree == 2)
			continue;

		if (v_degree == 0 && u_degree == 0) {
			// add new path {v,u}
			neighbours[v].push_back(u);
			neighbours[u].push_back(v);
			endpoint_to_endpoint[v] = u;
			endpoint_to_endpoint[u] = v;
			continue;
		}

		// at least one of v and u have exactly degree 1. for what follows, make sure v has degree 1.
		if (v_degree == 0) {
			int tmp = v;
			v = u;
			u = tmp;
			u_degree = 0;
			v_degree = 1;
		}
		if (u_degree == 1) { // both u and v have degree 1
			// make sure u and v are not endpoints of the same path as that would create a proper subcycle
			// (actually it will not create a proper subcycle if all nodes are on a single big path, but we add this last edge afterwards instead)
			if (endpoint_to_endpoint[v] == u)
				continue;
			// connect paths that have u and v as endpoints
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = endpoint_to_endpoint[u];
			endpoint_to_endpoint[endpoint_to_endpoint[u]] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
			endpoint_to_endpoint[u] = -1;
			// u will have deg 2 after this
			num_nodes_of_deg_2++;
		}
		else { // only v has degree 1
			   // u is new endpoint of path that v was formerly endpoint of
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = u;
			endpoint_to_endpoint[u] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
		}
		neighbours[v].push_back(u);
		neighbours[u].push_back(v);
		// v will have deg 2 after this
		num_nodes_of_deg_2++;
		if (num_nodes_of_deg_2 == n - 2)
			break;
	}

	// save the tour
	// 1. find an endpoint of the single path left
	int endpoint;
	for (int i = 0; i < n; i++)
		if (endpoint_to_endpoint[i] != -1) {
			endpoint = i;
			break;
		}
	// 2. walk through tour and save it to tour-array
	tour[0] = endpoint;
	int node = neighbours[endpoint][0]; // walk to only neighbour of endpoint
	int previous_node = endpoint; // this is used to make sure we don't go backwards.
	for (int i = 1; i < n - 1; i++) {
		tour[i] = node;
		int next_node = neighbours[node][0] == previous_node ? neighbours[node][1] : neighbours[node][0]; // step forward
		previous_node = node;
		node = next_node;
	}
	tour[n - 1] = node; // node is now other endpoint
}

bool CompareSavings(pair<int, int> &e1, pair<int, int> &e2) {
	return dist[hub][e1.first] + dist[hub][e1.second] - dist[e1.first][e1.second] > 
		dist[hub][e2.first] + dist[hub][e2.second] - dist[e2.first][e2.second];
}

void ClarkeWright(int hub_idx) {
	// set the hub, this is a static variable as it is used also by CompareSavings
	hub = hub_idx;

	// we will build a vector non_hub_pairs consisting of all pairs {v,u} where neither v nor u is the hub
	int num_non_hub_pairs = num_edges - (n - 1); // number of such pairs {v,u}
	vector<pair<int, int>> non_hub_pairs;
	non_hub_pairs.reserve(num_non_hub_pairs);
	for (int i = 0; i < n; i++) {
		if (i == hub)
			continue;
		for (int j = i + 1; j < n; j++) {
			if (j == hub)
				continue;
			non_hub_pairs.push_back(make_pair(i, j));
		}
	}
	// sort (in decreasing order) non_hub_pairs {v,u} by how much we save if we replace the edges {v,hub} and {hub,u} with {v,u}
	sort(begin(non_hub_pairs), begin(non_hub_pairs) + num_non_hub_pairs, CompareSavings);

	// we will now add all non-hub edges of the tour. we will save these edges as a collection of paths since any node will have at most two neighbours at any point of this process.
	vector<vector<int>> neighbours(n); // holds the neighbours of each node. hub will have zero neighbours since we only care about non-hub edges.
	int num_nodes_of_deg_2 = 0; // when all non-hub nodes except two have deg 2 we can end for loop over edges.
	int endpoint_to_endpoint[1000]; // lookup table from an endpoint of a path to the other endpoint. If a node is the endpoint of a path, then we can use this table to find the other endpoint.
	for (int i = 0; i < n; i++)
		endpoint_to_endpoint[i] = -1; // -1 represents that a node is not the endpoint of a path
	// loop over all sorted non_hub_pairs
	for (int i = 0; i < num_non_hub_pairs; i++) {
		pair<int, int> non_hub_pair = non_hub_pairs[i];
		int v = non_hub_pair.first;
		int u = non_hub_pair.second;
		int v_degree = neighbours[v].size();
		int u_degree = neighbours[u].size();

		if (v_degree == 2 || u_degree == 2)
			continue;

		if (v_degree == 0 && u_degree == 0) { 
			// add new path {v,u}
			neighbours[v].push_back(u);
			neighbours[u].push_back(v);
			endpoint_to_endpoint[v] = u;
			endpoint_to_endpoint[u] = v;
			continue;
		}
		
		// at least one of v and u have exactly degree 1. for what follows, make sure v has degree 1.
		if (v_degree == 0) {
			int tmp = v;
			v = u;
			u = tmp;
			u_degree = 0;
			v_degree = 1;
		}
		if (u_degree == 1) { // both u and v have degree 1
			// make sure u and v are not endpoints of the same path as that would create a proper subcycle
			if (endpoint_to_endpoint[v] == u)
				continue;
			// connect paths that have u and v as endpoints
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = endpoint_to_endpoint[u];
			endpoint_to_endpoint[endpoint_to_endpoint[u]] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
			endpoint_to_endpoint[u] = -1;
			// u will have deg 2 after this
			num_nodes_of_deg_2++;
		}
		else { // only v has degree 1
			// u is new endpoint of path that v was formerly endpoint of
			endpoint_to_endpoint[endpoint_to_endpoint[v]] = u;
			endpoint_to_endpoint[u] = endpoint_to_endpoint[v];
			endpoint_to_endpoint[v] = -1;
		}
		neighbours[v].push_back(u);
		neighbours[u].push_back(v);
		// v will have deg 2 after this
		num_nodes_of_deg_2++;
		if (num_nodes_of_deg_2 == n - 3)
			break;
	}

	// save the tour
	// 1. find an endpoint, this is a node that is still connected to the hub. there are exactly two such nodes
	int endpoint;
	for (int i = 0; i < n; i++)
		if (endpoint_to_endpoint[i] != -1) {
			endpoint = i;
			break;
		}
	// 2. walk through tour and save it to tour-array
	tour[0] = endpoint;
	int node = neighbours[endpoint][0]; // walk to only neighbour of endpoint
	int previous_node = endpoint; // this is used to make sure we don't go backwards.
	for (int i = 1; i < n-2; i++) {
		tour[i] = node;
		int next_node = neighbours[node][0] == previous_node ? neighbours[node][1] : neighbours[node][0]; // step forward
		previous_node = node;
		node = next_node;
	}
	tour[n - 2] = node; // node is now other endpoint
	tour[n - 1] = hub;
}

// ***

// *** OPTMIZATION ALGOS

void TwoOpt() {
	// allocate n double linked list entries. we will work with the tour as a linked list.
	vector<DoubleEntry> entries(n);
	// build linked list of tour
	for (int i = 0; i < n; i++) {
		int node = tour[i];
		entries[node].node = node;
		entries[node].next = &entries[tour[i < n - 1 ? i + 1 : 0]];
		entries[node].prev = &entries[tour[i > 0 ? i - 1 : n - 1]];
	}

	bool made_improvement = false; // we do 2-opt until we longer make an improvement for any node i below
	do {
		made_improvement = false;
		for (int i = 0; i < n; i++) {
			DoubleEntry *t1 = &entries[i];
			DoubleEntry *t2 = t1->next;

			int max_saving = 0; // this will keep track of the best saving a 3-opt move involving t1 and t2 can give us
			DoubleEntry *max_saving_t3 = NULL; // these are needed when we are going to perform the move that gives the best saving (if we find a move that gives saving > 0 at all that is)
			DoubleEntry *max_saving_t4 = NULL; //

			// loop over all t4-candidates (these are close neighbours of t1 (node i))
			int num_t4_to_check = min(30, (int) closest_neighbours_two_opt[t1->node].size()); // optimisation: only look at 30 closest neighbours
			for (int j = 0; j < num_t4_to_check; j++) { 
				int t4_node = closest_neighbours_two_opt[t1->node][j].second; // keep in mind closestNeighbours[t1->node] is a vector of pairs where t1 is the first in each pair
				DoubleEntry *t4 = &entries[t4_node];
				DoubleEntry *t3 = t4->next;

				int saving = dist[t1->node][t2->node] + dist[t3->node][t4->node] - (dist[t1->node][t4->node] + dist[t2->node][t3->node]); // old cost - new cost
				if (saving > max_saving) {
					max_saving = saving;
					max_saving_t3 = t3;
					max_saving_t4 = t4;
				}
			}

			// make move if improvement found
			if (max_saving > 0) {
				DoubleEntry *t3 = max_saving_t3;
				DoubleEntry *t4 = max_saving_t4;
				
				// note that t1 != t3 and t2 != t4 since otherwise 2-opt would not be profitable (saving would be zero)
				t1->next = t4;
				DoubleEntry *entry = t4->prev;
				t4->prev = t1;
				t4->next = entry;
				while (entry != t2) {
					DoubleEntry *tmp = entry->next;
					entry->next = entry->prev;
					entry->prev = tmp;

					entry = entry->next;
				}
				t2->prev = t2->next;
				t2->next = t3;
				t3->prev = t2;

				made_improvement = true;
			}
		}
	} while (made_improvement);

	// copy tour to tour-array
	DoubleEntry *start = &entries[0];
	DoubleEntry *entry = start;
	int tour_idx = 0;
	do {
		tour[tour_idx++] = entry->node;
		entry = entry->next;
	} while (entry != start);

	UpdateMinTour();
}

inline void CheckDegenerateThreeOptMove(
	DoubleEntry* &t1, 
	DoubleEntry* &t2,
	DoubleEntry* &t3,
	DoubleEntry* &t4,
	int &max_saving,
	DoubleEntry* &max_saving_t1,
	DoubleEntry* &max_saving_t2,
	DoubleEntry* &max_saving_t3, 
	DoubleEntry* &max_saving_t4,
	int &situation
) {
	int saving = dist[t1->node][t2->node] + dist[t3->node][t4->node] - (dist[t1->node][t4->node] + dist[t2->node][t3->node]); // old cost - new cost
	if (saving > max_saving) {
		max_saving = saving;
		max_saving_t1 = t1;
		max_saving_t2 = t2;
		max_saving_t3 = t3;
		max_saving_t4 = t4;
		situation = THREE_OPT_SIT_TWO_OPT;
	}
}

inline void CheckThreeOptMoveByT5(
	DoubleEntry* &t1,
	DoubleEntry* &t2,
	DoubleEntry* &t3,
	DoubleEntry* &t4,
	DoubleEntry* &t5,
	vector<int> &tour_pos,
	int &max_saving,
	DoubleEntry* &max_saving_t1,
	DoubleEntry* &max_saving_t2,
	DoubleEntry* &max_saving_t3,
	DoubleEntry* &max_saving_t4,
	DoubleEntry* &max_saving_t5,
	DoubleEntry* &max_saving_t6,
	int &situation
) {
	// we have t5 and we need t6. whether t6 is t5->next or t5->prev depends on how nodes are ordered on tour.
	DoubleEntry* t6 = NULL;
	// keep in mind that we know t4 is before t1 on tour and that t5 != t4 and t5 != t3.
	if (tour_pos[t5->node] <= tour_pos[t1->node] && // check if t5 is between t3 and t1 on tour
		tour_pos[t5->node] > tour_pos[t3->node]) {
		t6 = t5->prev; // this is situation B from report
	}
	else { // t5 is not between t3 and t1 on tour
		t6 = t5->next; // this is situation A
	}
	// check how much we save with this move
	int saving = dist[t1->node][t2->node] + dist[t3->node][t4->node] + dist[t5->node][t6->node] -
		(dist[t4->node][t1->node] + dist[t3->node][t5->node] + dist[t2->node][t6->node]); // save = old cost - new cost
	if (saving > max_saving) {
		max_saving = saving;
		max_saving_t1 = t1;
		max_saving_t2 = t2;
		max_saving_t3 = t3;
		max_saving_t4 = t4;
		max_saving_t5 = t5;
		max_saving_t6 = t6;
		situation = t6 == t5->next ? THREE_OPT_SIT_A : THREE_OPT_SIT_B;
	}
}

inline void CheckThreeOptMoveByT6(
	DoubleEntry* &t1,
	DoubleEntry* &t2,
	DoubleEntry* &t3,
	DoubleEntry* &t4,
	DoubleEntry* &t6,
	vector<int> &tour_pos,
	int &max_saving,
	DoubleEntry* &max_saving_t1,
	DoubleEntry* &max_saving_t2,
	DoubleEntry* &max_saving_t3,
	DoubleEntry* &max_saving_t4,
	DoubleEntry* &max_saving_t5,
	DoubleEntry* &max_saving_t6,
	int &situation
) {
	// we have t6 and we need t5. whether t5 is t6->next or t6->prev depends on how nodes are ordered on tour.
	DoubleEntry* t5 = NULL;
	// keep in mind that we know t4 is before t1 on tour and that t6 != t1 (by earlier check) and t6 != t2 (since t6 is close neighbour of t2).
	if (tour_pos[t6->node] < tour_pos[t1->node] && // check if t6 is between t3 and t1 on tour
		tour_pos[t6->node] >= tour_pos[t3->node]) {
		t5 = t6->next; // this is situation B from report
	}
	else { // t6 is not between t3 and t1 on tour
		t5 = t6->prev; // this is situation A
	}
	// check how much we save with this move
	int saving = dist[t1->node][t2->node] + dist[t3->node][t4->node] + dist[t5->node][t6->node] -
		(dist[t4->node][t1->node] + dist[t3->node][t5->node] + dist[t2->node][t6->node]); // save = old cost - new cost
	if (saving > max_saving) {
		max_saving = saving;
		max_saving_t1 = t1;
		max_saving_t2 = t2;
		max_saving_t3 = t3;
		max_saving_t4 = t4;
		max_saving_t5 = t5;
		max_saving_t6 = t6;
		situation = t6 == t5->next ? THREE_OPT_SIT_A : THREE_OPT_SIT_B;
	}
}

void ThreeOpt() {
	// allocate n double linked list entries. we will work with the tour as a linked list.
	vector<DoubleEntry> entries(n);

	// build linked list of tour
	for (int i = 0; i < n; i++) {
		int node = tour[i];
		entries[node].node = node;
		entries[node].next = &entries[tour[i < n - 1 ? i + 1 : 0]];
		entries[node].prev = &entries[tour[i > 0 ? i - 1 : n - 1]];
	}

	bool made_improvement = false; // we do 3-opt until we longer make an improvement for any node i below
	int iteration_counter = 0;
	do {
		iteration_counter++;
		made_improvement = false;
		vector<int> tour_pos(n); // tour_pos is used when we need to know how nodes are ordered on tour. it is a lookup table such 
		// that tour_pos[i] gives the position of node i on the tour. we let node 0 have position 0.
		bool tour_pos_outdated = true; // we need to update tour_pos initially, after that tour_pos will only be outdated after we made an improvement (move) for some i.
		for (int i = 0; i < n; i++) {

			// do we need to update tour_pos?
			if (tour_pos_outdated) {
				DoubleEntry *start = &entries[0];
				DoubleEntry *entry = start;
				int counter = 0;
				do {
					tour_pos[entry->node] = counter++;
					entry = entry->next;
				} while (entry != start);
				tour_pos_outdated = false;
			}

			int max_saving = 0; // this will keep track of the best saving a 3-opt move involving t1 and t2 can give us
			DoubleEntry *max_saving_t1 = NULL; // these are needed when we are going to perform the move that gives the best saving (if we find a move that gives saving > 0 at all that is)
			DoubleEntry *max_saving_t2 = NULL; //
			DoubleEntry *max_saving_t3 = NULL; //
			DoubleEntry *max_saving_t4 = NULL; //
			DoubleEntry *max_saving_t5 = NULL; //
			DoubleEntry *max_saving_t6 = NULL; //
			int situation = 0; // tells us which one of the three 3-opt situations we are in. see report and constants THREE_OPT_SIT at beginning of file
			
			// loop over all t4-candidates (these are close neighbours of t1 (node i))
			int num_t4_to_check = min(30, n-1); // optimisation: only look at 30 closest neighbours
			for (int j = 0; j < num_t4_to_check; j++) { 
				int t4_node = closest_neighbours[i][j].second; // keep in mind closestNeighbours[i] is a vector of pairs where t1 (i) is the first in each pair
				DoubleEntry *t1 = NULL;
				DoubleEntry *t4 = NULL;
				if (tour_pos[t4_node] < tour_pos[i]) { // for case handeling we always want t1 after t4 on tour, if this is not the case we switch the variables.
					t1 = &entries[i];
					t4 = &entries[t4_node];
				}
				else {
					t1 = &entries[t4_node];
					t4 = &entries[i];
				}
				DoubleEntry *t2 = t1->next;
				DoubleEntry *t3 = t4->next;

				// check for a two-opt move
				CheckDegenerateThreeOptMove(t1, t2, t3, t4, max_saving, max_saving_t1, max_saving_t2, max_saving_t3, max_saving_t4, situation);

				// now check for 3-opt-moves. first, check for a neighbour of t3
				int num_t5_to_check = min(20, n - 1);
				for (int k = 0; k < num_t5_to_check; k++) {
					int t5_node = closest_neighbours[t3->node][k].second; // keep in mind closestNeighbours[t3->node] is a vector of pairs where t3 is the first in each pair
					if (t5_node == t4->node) // this special case will result in a two-opt, which we already covered above
						continue;
					DoubleEntry *t5 = &entries[t5_node];
					CheckThreeOptMoveByT5(t1, t2, t3, t4, t5, tour_pos, max_saving, max_saving_t1, max_saving_t2, max_saving_t3, max_saving_t4, max_saving_t5, max_saving_t6, situation);
				}
				// second, check for neighbour of t2
				for (int k = 0; k < num_t5_to_check; k++) {
					int t6_node = closest_neighbours[t2->node][k].second; // keep in mind closestNeighbours[t2->node] is a vector of pairs where t2 is the first in each pair
					if (t6_node == t1->node) // this special case will result in a two-opt, which we already covered above
						continue;
					DoubleEntry *t6 = &entries[t6_node];
					CheckThreeOptMoveByT6(t1, t2, t3, t4, t6, tour_pos, max_saving, max_saving_t1, max_saving_t2, max_saving_t3, max_saving_t4, max_saving_t5, max_saving_t6, situation);
				}
			}
				
			// make move if improvement found
			if (max_saving > 0) {
				DoubleEntry *t1 = max_saving_t1;
				DoubleEntry *t2 = max_saving_t2;
				DoubleEntry *t3 = max_saving_t3;
				DoubleEntry *t4 = max_saving_t4;
				DoubleEntry *t5 = max_saving_t5;
				DoubleEntry *t6 = max_saving_t6;
				if (situation == THREE_OPT_SIT_TWO_OPT) {
					// note that t1 != t3 and t2 != t4 since otherwise 2-opt would not be profitable (saving would be zero)
					t1->next = t4;
					DoubleEntry *entry = t4->prev;
					t4->prev = t1;
					t4->next = entry;
					while (entry != t2) {
						DoubleEntry *tmp = entry->next;
						entry->next = entry->prev;
						entry->prev = tmp;

						entry = entry->next;
					}
					t2->prev = t2->next;
					t2->next = t3;
					t3->prev = t2;
				}
				else if (situation == THREE_OPT_SIT_A) { // reconnect tour for situation A
					t1->next = t4;
					if (t4 == t6) {
						t4->prev = t1;
					}
					else {
						DoubleEntry *entry = t4->prev;
						t4->prev = t1;
						t4->next = entry;
						while (entry != t6) {
							DoubleEntry *tmp = entry->next;
							entry->next = entry->prev;
							entry->prev = tmp;

							entry = entry->next;
						}
						t6->prev = t6->next;
					}
					t6->next = t2;
					t2->prev = t6;
					t5->next = t3;
					t3->prev = t5;
				}
				else { // reconnect tour for situation B
					t6->next = t2;
					t2->prev = t6;
					t4->next = t1;
					if (t1 == t5) {
						t1->prev = t4;
					}
					else {
						DoubleEntry *entry = t1->prev;
						t1->prev = t4;
						t1->next = entry;
						while (entry != t5) {
							DoubleEntry *tmp = entry->next;
							entry->next = entry->prev;
							entry->prev = tmp;

							entry = entry->next;
						}
						t5->prev = t5->next;
					}
					t5->next = t3;
					t3->prev = t5;
				}

				made_improvement = true;
				tour_pos_outdated = true; // need to update tour pos next run
			}
		}
	} while (made_improvement);

	// copy tour to tour-array
	DoubleEntry *start = &entries[0];
	DoubleEntry *entry = start;
	int tour_idx = 0;
	do {
		tour[tour_idx++] = entry->node;
		entry = entry->next;
	} while (entry != start);

	UpdateMinTour();
}

void ThreeOptQuadrants() { // this 3-opt is the same as the previous one except that it uses closest_neighbours_quadrants instead
	// of the regular closest_neighbours. There is really no need to read this code if the reader has already read ThreeOpt().
	// the reason there are two so similar functions has to do with time limitations.

	// allocate n double linked list entries. we will work with the tour as a linked list.
	vector<DoubleEntry> entries(n);

	// build linked list of tour
	for (int i = 0; i < n; i++) {
		int node = tour[i];
		entries[node].node = node;
		entries[node].next = &entries[tour[i < n - 1 ? i + 1 : 0]];
		entries[node].prev = &entries[tour[i > 0 ? i - 1 : n - 1]];
	}

	bool made_improvement = false; // we do 3-opt until we longer make an improvement for any node i below
	int iteration_counter = 0;
	do {
		iteration_counter++;
		made_improvement = false;
		vector<int> tour_pos(n); // tour_pos is used when we need to know how nodes are ordered on tour. it is a lookup table such 
								 // that tour_pos[i] gives the position of node i on the tour. we let node 0 have position 0.
		bool tour_pos_outdated = true; // we need to update tour_pos initially, after that tour_pos will only be outdated after we made an improvement (move) for some i.
		for (int i = 0; i < n; i++) {

			// do we need to update tour_pos?
			if (tour_pos_outdated) {
				DoubleEntry *start = &entries[0];
				DoubleEntry *entry = start;
				int counter = 0;
				do {
					tour_pos[entry->node] = counter++;
					entry = entry->next;
				} while (entry != start);
				tour_pos_outdated = false;
			}

			int max_saving = 0; // this will keep track of the best saving a 3-opt move involving t1 and t2 can give us
			DoubleEntry *max_saving_t1 = NULL; // these are needed when we are going to perform the move that gives the best saving (if we find a move that gives saving > 0 at all that is)
			DoubleEntry *max_saving_t2 = NULL; //
			DoubleEntry *max_saving_t3 = NULL; //
			DoubleEntry *max_saving_t4 = NULL; //
			DoubleEntry *max_saving_t5 = NULL; //
			DoubleEntry *max_saving_t6 = NULL; //
			int situation = 0; // tells us which one of the three 3-opt situations we are in. see report and constants THREE_OPT_SIT at beginning of file

							   // loop over all t4-candidates (these are close neighbours of t1 (node i))
			int num_t4_to_check = min(30, (int) closest_neighbours_quadrants[i].size()); // optimisation: only look at 30 closest neighbours
			for (int j = 0; j < num_t4_to_check; j++) {
				int t4_node = closest_neighbours_quadrants[i][j].second; // keep in mind closestNeighbours[i] is a vector of pairs where t1 (i) is the first in each pair
				DoubleEntry *t1 = NULL;
				DoubleEntry *t4 = NULL;
				if (tour_pos[t4_node] < tour_pos[i]) { // for case handeling we always want t1 after t4 on tour, if this is not the case we switch the variables.
					t1 = &entries[i];
					t4 = &entries[t4_node];
				}
				else {
					t1 = &entries[t4_node];
					t4 = &entries[i];
				}
				DoubleEntry *t2 = t1->next;
				DoubleEntry *t3 = t4->next;

				// check for a two-opt move
				CheckDegenerateThreeOptMove(t1, t2, t3, t4, max_saving, max_saving_t1, max_saving_t2, max_saving_t3, max_saving_t4, situation);

				// now check for 3-opt-moves. first, check for a neighbour of t3
				int num_t5_to_check = min(20, (int) closest_neighbours_quadrants[t3->node].size());
				for (int k = 0; k < num_t5_to_check; k++) {
					int t5_node = closest_neighbours_quadrants[t3->node][k].second; // keep in mind closestNeighbours[t3->node] is a vector of pairs where t3 is the first in each pair
					if (t5_node == t4->node) // this special case will result in a two-opt, which we already covered above
						continue;
					DoubleEntry *t5 = &entries[t5_node];
					CheckThreeOptMoveByT5(t1, t2, t3, t4, t5, tour_pos, max_saving, max_saving_t1, max_saving_t2, max_saving_t3, max_saving_t4, max_saving_t5, max_saving_t6, situation);
				}
				// second, check for neighbour of t2
				int num_t6_to_check = min(20, (int) closest_neighbours_quadrants[t2->node].size());
				for (int k = 0; k < num_t6_to_check; k++) {
					int t6_node = closest_neighbours_quadrants[t2->node][k].second; // keep in mind closestNeighbours[t2->node] is a vector of pairs where t2 is the first in each pair
					if (t6_node == t1->node) // this special case will result in a two-opt, which we already covered above
						continue;
					DoubleEntry *t6 = &entries[t6_node];
					CheckThreeOptMoveByT6(t1, t2, t3, t4, t6, tour_pos, max_saving, max_saving_t1, max_saving_t2, max_saving_t3, max_saving_t4, max_saving_t5, max_saving_t6, situation);
				}
			}

			// make move if improvement found
			if (max_saving > 0) {
				DoubleEntry *t1 = max_saving_t1;
				DoubleEntry *t2 = max_saving_t2;
				DoubleEntry *t3 = max_saving_t3;
				DoubleEntry *t4 = max_saving_t4;
				DoubleEntry *t5 = max_saving_t5;
				DoubleEntry *t6 = max_saving_t6;
				if (situation == THREE_OPT_SIT_TWO_OPT) {
					// note that t1 != t3 and t2 != t4 since otherwise 2-opt would not be profitable (saving would be zero)
					t1->next = t4;
					DoubleEntry *entry = t4->prev;
					t4->prev = t1;
					t4->next = entry;
					while (entry != t2) {
						DoubleEntry *tmp = entry->next;
						entry->next = entry->prev;
						entry->prev = tmp;

						entry = entry->next;
					}
					t2->prev = t2->next;
					t2->next = t3;
					t3->prev = t2;
				}
				else if (situation == THREE_OPT_SIT_A) { // reconnect tour for situation A
					t1->next = t4;
					if (t4 == t6) {
						t4->prev = t1;
					}
					else {
						DoubleEntry *entry = t4->prev;
						t4->prev = t1;
						t4->next = entry;
						while (entry != t6) {
							DoubleEntry *tmp = entry->next;
							entry->next = entry->prev;
							entry->prev = tmp;

							entry = entry->next;
						}
						t6->prev = t6->next;
					}
					t6->next = t2;
					t2->prev = t6;
					t5->next = t3;
					t3->prev = t5;
				}
				else { // reconnect tour for situation B
					t6->next = t2;
					t2->prev = t6;
					t4->next = t1;
					if (t1 == t5) {
						t1->prev = t4;
					}
					else {
						DoubleEntry *entry = t1->prev;
						t1->prev = t4;
						t1->next = entry;
						while (entry != t5) {
							DoubleEntry *tmp = entry->next;
							entry->next = entry->prev;
							entry->prev = tmp;

							entry = entry->next;
						}
						t5->prev = t5->next;
					}
					t5->next = t3;
					t3->prev = t5;
				}

				made_improvement = true;
				tour_pos_outdated = true; // need to update tour pos next run
			}
		}
	} while (made_improvement);

	// copy tour to tour-array
	DoubleEntry *start = &entries[0];
	DoubleEntry *entry = start;
	int tour_idx = 0;
	do {
		tour[tour_idx++] = entry->node;
		entry = entry->next;
	} while (entry != start);

	UpdateMinTour();
}

void DynprogOptBacktrack(int subset, int &s, vector<int> &optimal_path, vector<int> &path) {
	if (subset == 0)
		return;

	int last = optimal_path[optimal_path.size() - 1];

	int min = -1;
	int min_idx = -1;
	int min_flag = -1;
	// flag encodes the element c below
	int flag = 1;
	for (int c = 0; flag <= subset; c++) {
		if ((subset & flag) == 0) { // continue if c is not in subset
			flag = flag << 1;
			continue;
		}

		if (min == -1 || D_opt[subset][c] + dist[path[c]][path[last]] < min) {
			min = D_opt[subset][c] + dist[path[c]][path[last]];
			min_idx = c;
			min_flag = flag;
		}

		flag = flag << 1;
	}

	optimal_path.push_back(min_idx);
	DynprogOptBacktrack(subset ^ min_flag, s, optimal_path, path);
}

// ******

//C8
int main()
{
	high_resolution_clock::time_point time = high_resolution_clock::now();

	readInputFromFile();

	if (n <= 3) // handle trivial corner cases
	{
		for (int i = 0; i < n; i++)
			cout << i << '\n';
	}
	else {
		computeDistances();
		num_edges = (n*(n - 1)) / 2;

		ClarkeWright(2);
		UpdateMinTour();

		OutputMinTour();
	}

	return 0;
}
