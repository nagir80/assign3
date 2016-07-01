// This C++ program creates a graph from an input file,
// and if it is connected, calculates the
// minimum spanning tree (MST) using Prim's algorithm

#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <queue>
using namespace std;

// Set this switch to true to print more details
const bool verbose = true;

// Set the input filename
const string input_filename("test_data.txt");

// Define infinity
const double inf = 1.0/0.0;
const int fill_value = -1;

// Define the edge structure
struct Edge {
  int source_id;
  int target_id;
  double cost;
};

struct CompareEdgeFunction {
  inline bool operator()(const Edge& p, const Edge& q) {
    if (p.source_id != q.source_id) {
      return p.source_id < q.source_id;
    }
    return p.target_id < q.target_id;
  }
};

// Compare the cost of edge p against the cost of edge q.
struct CompareEdgeCostFunction {
  inline bool operator() (const Edge& p, const Edge& q) {
    return (p.cost < q.cost);
  }
};

// Find the index of the edge in the vector of edges.
// Return fill_value if not found.
int find_edge_index(const vector<Edge>& edges, const Edge& edge) {
  int num_edges = static_cast<int>(edges.size());
  for (int i = 0; i < num_edges; ++i) {
    Edge this_edge = edges[i];
    if (this_edge.source_id == edge.source_id &&
      this_edge.target_id == edge.target_id &&
      this_edge.cost == edge.cost) {
      return i;
    }
  }
  return fill_value;
}

// Overload the << operator to be able to print the edge to screen
ostream& operator<< (ostream& out, const Edge& e) {
  out << "(" << e.source_id << ", " << e.target_id
      << ", " << e.cost << ")";
  return out;
}

// We use the list representation of a graph.
// We make use of vectors to store:
// 1) the indices of the connected vertices, and
// 2) the weights (or costs) of each connection.
class Graph {
  public:
    Graph(int num_vertices_in_graph); // Default constructor for an empty graph
    Graph(const string& filename); // Constructor for graph read from file
    bool is_connected();
    int get_num_edges();
    Graph get_minimal_spanning_tree();
    double get_total_cost();
    void show();
  private:
    int num_vertices;
    vector<vector<int> > vertices;
    vector<vector<double> > edge_costs;
    void add_edge(int v1_id, int v2_id, double cost);
    vector<Edge> get_all_edges(void);
    vector<Edge> get_unique_edges(void);
};

// Add the edge to the undirected graph
void Graph::add_edge(int v1_id, int v2_id, double cost) {
  // If v2_id is not in list for vertex 1, add it
  if (find(vertices[v1_id].begin(), vertices[v1_id].end(), v2_id) == vertices[v1_id].end()) {
    vertices[v1_id].push_back(v2_id);
    edge_costs[v1_id].push_back(cost);
  }
  // If v1_id is not in list for vertex 2, add it
  if (find(vertices[v2_id].begin(), vertices[v2_id].end(), v1_id) == vertices[v2_id].end()) {
    vertices[v2_id].push_back(v1_id);
    edge_costs[v2_id].push_back(cost);
  }
}

// Show the graph. For each vertex, we show a list of adjacent vertices.
// For example, 3: 4[9] 5[8] 6[2] means
// that vertex 3 is connected to vertex 4, 5, and 6, and
// the corresponding weights (between brackets) are 9, 8, and 2, respectively.
void Graph::show(void) {
  for (int i = 0; i < num_vertices; ++i) {
    cout << i << ": ";
    int num_adjacent_vertices = static_cast<int>(vertices[i].size());
    for (int j = 0; j < num_adjacent_vertices; ++j) {
        cout << vertices[i][j] << "[" << edge_costs[i][j] << "]" << " ";
    }
    cout << endl;
  }
}

// Split the string into separate words.
inline vector<string> split(const string& s, char delimiter) {
  vector<string> words;
  stringstream ss(s);
  string word;
  while (getline(ss, word, delimiter)) {
    words.push_back(word);
  }
  return words;
}

// Create an empty graph with 'num_vertices_in_graph' vertices.
Graph::Graph(int num_vertices_in_graph) {
  // Allow roomt to store one empty vector for each vertex
  num_vertices = num_vertices_in_graph;
  vertices.resize(num_vertices_in_graph);
  edge_costs.resize(num_vertices_in_graph);
}

// Constructor to create a graph from a file.
// The file format is: an initial integer
// that is the node size of the graph and the
// further values are integer triples: (i, j, cost).
Graph::Graph(const string& input_filename) {
  // Open the input file
  ifstream input_file(input_filename.c_str());

  // Read header
  string line;
  getline(input_file, line);
  int num_nodes = atoi(line.c_str());

  // Set the number of vertices
  num_vertices = num_nodes;

  // Store one empty vector for each vertex
  vertices.resize(num_vertices);
  edge_costs.resize(num_vertices);

  while (!input_file.eof())
  {
    getline(input_file, line);
    vector<string> words = split(line, ' ');
    if (words.size() == 3) {
      // We have a valid edge.
      // Extract the edge information
      int v1_id = atoi(words[0].c_str());
      int v2_id = atoi(words[1].c_str());
      double cost = static_cast<double>(atoi(words[2].c_str()));

      // Add the edge to the graph
      add_edge(v1_id, v2_id, cost);
    }
  }

  input_file.close();
}

// Return the edge with the smallest cost.
Edge get_best_edge(const vector<Edge>& edges) {

  // Copy the vector
  vector<Edge> edges_work = edges;

  // Sort the vector
  sort(edges_work.begin(), edges_work.end(), CompareEdgeCostFunction());

  // Assume the first element has the smallest cost.
  Edge best_edge = edges_work[0];
  return best_edge;
}

// Given an input vector of edge and the vertex id,
// return the vector of neighbouring edges
vector<Edge> get_adjacent_edges(const vector<Edge>& edges, int vertex_id) {
  vector<Edge> adjacent_edges;
  int num_edges = static_cast<int>(edges.size());
  for (int i = 0; i < num_edges; ++i) {
    Edge e = edges[i];
    if (e.source_id == vertex_id) {
      adjacent_edges.push_back(e);
    }
  }
  return adjacent_edges;
}

// Check if graph is connected,
// i.e. all nodes can be reached from node 0
// using bread first search.
bool Graph::is_connected(void) {
  // Get the list of unique edges
  vector<Edge> edges = get_all_edges();

  // Keep track of nodes that have been visited.
  // Start with all nodes are unvisited.
  bool visited[num_vertices];
  for (int i = 0; i < num_vertices; ++i) {
    visited[i] = false;
  }

  // Start with node 0:
  // add edge with  source_id = 0 to the queue.
  queue<Edge> q;
  q.push(edges[0]);
  visited[0] = true;

  while (!q.empty()) {
    // Grab the first edge and remove it from the queue
    Edge e = q.front();
    q.pop();

    // Get the neigbouring vertices
    int source_id = e.source_id;
    vector<Edge> adjacent_edges = get_adjacent_edges(edges, source_id);

    // Loop over the neighbours
    int num_adjacent_edges = static_cast<int>(adjacent_edges.size());
    for (int j = 0; j < num_adjacent_edges; ++j) {
      int target_id = adjacent_edges[j].target_id;

      // If the node has not been visitited before add it to the queue
      // and mark it as visited.
      if (!visited[target_id]) {
        visited[target_id] = true;

        // Loop over the neighbours
        vector<Edge> next_adjacent_edges = get_adjacent_edges(edges, target_id);
        int num_next_adjacent_edges = static_cast<int>(next_adjacent_edges.size());
        for (int k = 0; k < num_next_adjacent_edges; ++k) {
          q.push(next_adjacent_edges[k]);
        }
      }
    }
  }

  // Count the number of visited nodes
  int num_visited = 0;
  for (int i = 0; i < num_vertices; ++i) {
    if (visited[i]) {
      num_visited++;
    }
  }

  if (num_visited == num_vertices) {
    // Yes, the graph is connected
    return true;
  } else {
    // If the number of visited nodes is less than
    // the number of vertices in the graph, the graph is unconnected
    return false;
  }
}

// Use Prim's algorithm to derive the minimum spanning tree
Graph Graph::get_minimal_spanning_tree(void) {
  // Check if the graph is connected.
  // If not, return an empty graph: a minimal spanning tree is not available.
  if (!is_connected()) {
    Graph m(0);
    return m;
  }

  // The minimum spanning tree has the same number of
  // vertices as the input graph
  int num_vertices_in_graph = num_vertices;
  Graph m(num_vertices_in_graph);

  // Associate with each vertex v of the graph a number C[v]
  // (the cheapest cost of a connection to v) and an edge E[v]
  // (the edge providing that cheapest connection). To initialize these values,
  // set all values of C[v] to +inf (or to any number larger than the maximum edge weight)
  // and set each E[v] to a special flag value indicating that
  // there is no edge connecting v to earlier vertices.
  vector<double> C(num_vertices_in_graph);
  vector<Edge> E(num_vertices_in_graph);
  // Store the finished edges in F.
  vector<Edge> F;
  for (int i = 0; i < num_vertices_in_graph; ++i) {
    C[i] = inf;
  }

  // Get the list edges
  vector<Edge> Q = get_all_edges();

  // Start with vertex 0
  C[0] = 0.0;
  E[0].source_id = 0;

  // Repeat the following steps until Q is empty
  while (static_cast<int>(F.size()) < num_vertices_in_graph-1) {
    // Find and remove a vertex v from Q having the minimum possible value of C[v]
    Edge best_edge = get_best_edge(Q);

    // Remove the edge from Q in both directions
    int ip = find_edge_index(Q, best_edge);
    if (ip != fill_value) {
      Q.erase(Q.begin()+ip);
    }
    Edge best_edge_flipped = {best_edge.target_id, best_edge.source_id, best_edge.cost};
    int ip_flipped = find_edge_index(Q, best_edge_flipped);
    if (ip_flipped != fill_value) {
      Q.erase(Q.begin()+ip_flipped);
    }

    // Add this edge to the minimal spanning tree
    m.add_edge(best_edge.source_id, best_edge.target_id, best_edge.cost);
    m.add_edge(best_edge.target_id, best_edge.source_id, best_edge.cost);
    F.push_back(best_edge);

    // Get adjacent vertices
    int source_id = best_edge.source_id;
    vector<Edge> adjacent_edges = get_adjacent_edges(Q, source_id);
    int num_adjacent_edges = static_cast<int>(adjacent_edges.size());
    for (int j = 0; j < num_adjacent_edges; ++j) {
      int target_id = adjacent_edges[j].target_id;
      double cost = adjacent_edges[j].cost;

      // Pick the edge with the smallest cost
      if (adjacent_edges[j].cost < C[target_id]) {
        C[target_id] = cost;
        E[target_id].source_id = source_id;
        E[target_id].target_id = target_id;
        E[target_id].cost = cost;
      }
    }
  }

  // Return the minimum spanning tree
  return m;
}

// Convert the graph data to a vector of unique edges in one direction.
vector<Edge> Graph::get_unique_edges(void) {
  // Create a set of unique edges
  set<Edge, CompareEdgeFunction> es;
  int i;
  for (i = 0; i < num_vertices; ++i) {
    int num_adjacent_vertices = static_cast<int>(vertices[i].size());
    for (int j = 0; j < num_adjacent_vertices; ++j) {
      int source_id = i;
      int target_id = vertices[i][j];
      double cost = edge_costs[i][j];
      Edge edge = {source_id, target_id, cost};
      es.insert(edge);
      Edge edge_flipped = {target_id, source_id, cost};
      es.insert(edge_flipped);
    }
  }

  // Remove the edges in the backward direction
  for (set<Edge>::iterator it = es.begin(); it != es.end(); ++it) {
    Edge edge = *it;
    Edge edge_flipped = {edge.target_id, edge.source_id, edge.cost};
    es.erase(edge_flipped);
  }

  // Convert the set to a vector
  vector<Edge> edges(es.begin(), es.end());
  return edges;
}

// Convert the graph data to a vector of unique edges in both directions
vector<Edge> Graph::get_all_edges(void) {
  // Create a set of unique edges
  set<Edge, CompareEdgeFunction> es;
  int i;
  for (i = 0; i < num_vertices; ++i) {
    int num_adjacent_vertices = static_cast<int>(vertices[i].size());
    for (int j = 0; j < num_adjacent_vertices; ++j) {
      int source_id = i;
      int target_id = vertices[i][j];
      double cost = edge_costs[i][j];
      Edge edge = {source_id, target_id, cost};
      es.insert(edge);
      Edge edge_flipped = {target_id, source_id, cost};
      es.insert(edge_flipped);
    }
  }

  // Keep edges in both directions.
  // Convert the set to a vector.
  vector<Edge> edges(es.begin(), es.end());
  return edges;
}

// Obtain the number of edges in this undirected graph
int Graph::get_num_edges(void) {
  // Obtain the edges as a vector
  // (i.e. edge in one direction only)
  vector<Edge> edges = get_unique_edges();
  return static_cast<int>(edges.size());
}

// Sum up the costs of all edges
double Graph::get_total_cost(void) {
  // Convert the graph data to a list of unique edges.
  // Make sure we have edges in one direction (i.e.
  // do not count them double).
  vector<Edge> edges = get_unique_edges();

  // Sum the costs of all edges in the array
  double total_cost = 0.0;
  int num_edges = static_cast<int>(edges.size());
  for (int i = 0; i < num_edges; ++i) {
     total_cost += edges[i].cost;
  }
  return total_cost;
}

int main() {
  // Create the graph
  cout << "\n";
  cout << "---------------------------------" << endl;
  cout << "Taking Input Graph from File..."<< endl;
  cout << "---------------------------------" << endl;
  cout << "\n";

  Graph g(input_filename);
  if (verbose) {
    g.show();
  }

  // First check if the graph is connected
  bool is_connected = g.is_connected();
  if (is_connected) {
    cout << "\n";
    cout << "\n";
    cout << "------------------------------" << endl;
    cout << "Graph is Connected..!!" << endl;
    // Obtain the minimum spanning tree (which is a new graph).
    cout << "------------------------------" << endl;
    cout << "\n";
    cout << "Calculating the MST for the Input Graph..!!"<< endl;
    Graph m = g.get_minimal_spanning_tree();
    if (verbose) {
      cout << "\n";
      cout << "Minimum Spanning Tree: " << endl;
      cout << "\n";
      m.show();
    }
    // Determine the cost of the minimum spanning tree
      cout << "Cost: " << m.get_total_cost() << endl;
  } else {
    cout << "Graph is unconnected." << endl;
  }

  return 0;
}
