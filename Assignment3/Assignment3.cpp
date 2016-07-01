// Assignment3.cpp 
// solution of Assignment based on Prim's alghorithm;


#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <list>
#include <time.h>
#include <queue>
#include <algorithm>


//class for data of separate Weighted Edge
class WeightedEdge
{   int v;
    int w;
	double weight;
public:
//constructor
	WeightedEdge(int v, int w, double weight):v(v), w(w), weight(weight){};
//Getting weight of edge
	double get_weight(void) const { return weight;}
// Getting points of edge
	int from(void) const {return v;}
	int to(void) const {return w;}
// friend functions for output of edge and checking equality
	friend std::ostream& operator<< (std::ostream& out, const WeightedEdge& we); 
	friend bool operator==(WeightedEdge const& lhs, WeightedEdge const& rhs);
	~WeightedEdge(void);
};

std::ostream& operator<< (std::ostream& out, const WeightedEdge& we){
	out << "(" << we.from() << "," << we.to() << "," << we.get_weight() << ")";
	return out;
}

bool operator==(WeightedEdge const& lhs, WeightedEdge const& rhs){
	if(((lhs.from() == rhs.from())&&(lhs.to() == rhs.to())&&(lhs.get_weight() == rhs.get_weight()))||
	((lhs.from() == rhs.to())&&(lhs.to() == rhs.from())&&(lhs.get_weight() == rhs.get_weight())))
		return true;
	else
		return false;
}

WeightedEdge::~WeightedEdge(void)
{
}

//Base class for creating Random Graph. 
class Graph{
private:
	bool ** table;

// number of vertices
	int size;
public:
	double **weights;
// constructor for type Graph
	Graph(unsigned int size = 10, unsigned int p = 20);
// overloaded operator<< for printing table of graph
	friend std::ostream& operator<<(std::ostream& out,  Graph const & gr);

//return number of vertices
	int v(void) const {return size;}
	bool** get_table(void) const{return table;}
//return number of edges
	int e(void);
// check if edge is between x and y 
	bool adjacent(int x, int y) const {return table[x][y];}
};

//constructor for random undirected graph
Graph::Graph(unsigned int size, unsigned int p):size(size){
	table = new bool*[size];
	weights = new double*[size];
	for(unsigned int i = 0 ; i < size ; i++){
		table[i] = new bool[size];
		weights[i] = new double[size];
		table [i][i] = false;
		weights[i][i] = 0;
	}

	for (int i = 0; i < size ; i++){
		for (int j = 0 ; j < i ; j++){
				table [i][j] = table[j][i] = (rand() % 100 < p );
				weights[i][j] = weights[j][i] = static_cast<double>((rand()%100))/100;
		}
	}
}

int Graph::e(void){
	int edges = 0;
	for (int i = 0; i < size ; i++){
		for (int j = 0 ; j < size ; j++){
			if (table[i][j] != 0)
				edges++;
		}
	}
	return edges;
}


std::ostream& operator<<(std::ostream& out,  Graph const & gr){
	for (int i = 0 ; i < gr.size ; i++){
		for (int j = 0 ; j < gr.size ; j++){
			out << gr.table[i][j] << " ";
		}
		out << "\r\n";
	}
	return out;
}



//List of lists of Weighted edges, which based on Random Graph (previous class)
// or read from file
class WeightedEdgeGraph{
	std::list<std::list<WeightedEdge>> *listEdges;
	std::list<double> *vert_values;
public:
//constructor for creating Weighted Graph based on random graph
	WeightedEdgeGraph(Graph *graph);
	//constructor for creating Weighted Graph based on file
	WeightedEdgeGraph(std::string pathfile);
	friend std::ostream& operator<<(std::ostream& out,  WeightedEdgeGraph const & gr);
	std::list<std::list<WeightedEdge>>& get(void) const {return *listEdges;}
//getting edges from inputted vertex
	std::list<WeightedEdge> edges(int num_vertex);
//getting quantity if vertex
	int get_num_vertex(void) {return vert_values->size();}
};



WeightedEdgeGraph::WeightedEdgeGraph(Graph *graph){
	listEdges =  new std::list<std::list<WeightedEdge>>();
	vert_values = new std::list<double>();
	for(int i=0 ; i < graph->v() ; i ++){

		std::list<WeightedEdge> *hlp_list = new std::list<WeightedEdge>();
		if (i == 0) vert_values->push_back(0);
		else vert_values->push_back(100000000.0);

		for (int j=0 ; j <  graph->v() ; j++){
			if (graph->adjacent(i,j)){
				
				WeightedEdge *wEdge = new WeightedEdge(i, j, graph->weights[i][j]);

				hlp_list->push_back(*wEdge);


				delete wEdge;
			}
		}
		listEdges->push_back(*hlp_list);
	}
}

WeightedEdgeGraph::WeightedEdgeGraph(std::string pathfile){
	listEdges = new std::list<std::list<WeightedEdge>>();
	auto hlp_list = new std::list<WeightedEdge>();
	std::ifstream graph_file(pathfile);
	std::istream_iterator<int> start(graph_file), end;
	std::vector<int> words(start, end);

	int num_list = 0;
	for (auto it = words.begin()+1; it != words.end(); it+=3){
		int from;
		int to;
		double wgt;
		from = *it;
		to = *(it + 1);
		wgt = static_cast<double>(*(it + 2));

		if (num_list != from){
			listEdges->push_back(*hlp_list);
			hlp_list->clear();
			num_list = from;
		}

		hlp_list->push_back(*new WeightedEdge(from, to, wgt));
	}
	listEdges->push_back(*hlp_list);
	delete hlp_list;
}

std::list<WeightedEdge> WeightedEdgeGraph::edges(int num_vertex){
	std::list<std::list<WeightedEdge>>::iterator it = listEdges->begin();
	int i = num_vertex;
	while (i >0){
		it++;
		i--;
	}
	
	return *it;

}

std::ostream& operator<<(std::ostream& out,  WeightedEdgeGraph const & grl){
	for (std::list<std::list<WeightedEdge>>::iterator it = grl.get().begin() ; it != grl.get().end() ; ++it){
		static int num = 0;
		 out << num <<":";
		for (std::list<WeightedEdge>::iterator it_hlp = it->begin() ; it_hlp != it->end() ; ++it_hlp){
			out << "(" << it_hlp->from() << " " << it_hlp->to() << " " << it_hlp->get_weight() << ")";
		}
		out << "\r\n";
		num++;
	}
	return out;
}

//This PriorityQueue is based on std::map<> from STL

class PriorityQueue : public std::vector<WeightedEdge>//public std::map<int, double>
{
public:
// This function insert pair of number of the point and distanñe to this point from the start
	void insert(WeightedEdge& we);
// getting point with the shortest path
	WeightedEdge& top(void);
// deleting point with the shortest path
	void minPrioirty(void);
//size of the std::map, containing points
	int size(void);
	PriorityQueue(void);
	~PriorityQueue(void);
};

PriorityQueue::PriorityQueue(void)
{
}


PriorityQueue::~PriorityQueue(void)
{
}


void PriorityQueue::minPrioirty(void){
	//std::priority_queue<double, std::vector<double>, std::greater<double>>::pop();
	double min = 100000000.0;
	int it_pos = -1;
	int pos = 0;
	auto it_mem = this->begin();
	for (auto it = this->begin() ; it != this->end() ; ++it){
		if (it->get_weight() < min){
			min = it->get_weight();
			it_mem = it;
		}
		pos++;
	}
		this->erase(it_mem);
}


void PriorityQueue::insert(WeightedEdge& we){
	//std::priority_queue<double, std::vector<double>, std::greater<double>>::push(i);
	std::vector<WeightedEdge>::push_back(we);
}


int PriorityQueue::size(void){
	//return std::priority_queue<double, std::vector<double>, std::greater<double>>::size();
	return std::vector<WeightedEdge>::size();
}


WeightedEdge& PriorityQueue::top(void){
	//return std::priority_queue<double, std::vector<double>, std::greater<double>>::top();
	double min = 100000000.0;
	int it_pos = -1;
	auto it_mem = this->begin();
	for (auto it = this->begin() ; it != this->end() ; ++it){
		if (it->get_weight() < min){
			min = it->get_weight();
			it_mem = it;
		}
	}
		
	return *it_mem;
}

// class for implementation of Prim's algorithm
class Prim_MST{
	PriorityQueue* pq;
public:
	std::list<int> out_vertex;
	std::list<WeightedEdge> out_edges;
	double cost_tree;
// constructor gets objects of WeightedEdgeGraph, process it by algorithm and fill lists out_vertex and out_edges
// which include vertexes and edges of resulting MST
	Prim_MST(WeightedEdgeGraph *in_graph):cost_tree(0.0){
		pq = new PriorityQueue();
// adding of edges from start vertex to priority queue
		for (auto we : in_graph->edges(0)){
			pq->insert(we);
		}
//adding the start vertex in our MST
		out_vertex.push_back(0);
		// variable for counting edges in minimum spanning tree

		while (!pq->empty()){
//getting edge with the minimal weight
			WeightedEdge we_h = pq->top();
			pq->minPrioirty();
//checking if both vertexes of this edge doesn't include in our MST
			int n = 0;
			for (auto i : out_vertex){
				if (i == we_h.from())
					n++;
			}
			for (auto i : out_vertex){
				if (i == we_h.to())
					n++;
			}
// if not (n<2) then process more...
			if (n<2){
// adding verexes and edges to our MST
				out_edges.push_back(we_h);
				out_vertex.push_back(we_h.to());
// checking if our edge haven't been already added earlier
				for (auto we : in_graph->edges(we_h.to())){

					int compare = 0;
					for (auto i : out_edges){
						if (i == we){
							compare = 1;
							break;
						}
					}
					for (auto i : out_vertex){
						//add to our MST
						if (i == we.to()){
							compare = 1;
							break;
						}
					}
					if (!compare){
						pq->insert(we);

					}
				}
			}

		}
	}
// output of result of algoritm working
	void out_results(void){
		std::cout << "Vertexes: " << "\r\n";
		for (auto i : out_vertex){
			std::cout << i << " ";
		}

		std::cout << "\r\n";
		std::cout << "Edges: " << "\r\n";
		for (auto i : out_edges){
			cost_tree += i.get_weight();
			std::cout << i << " ";
		}
		std::cout << "\r\n";
		std::cout << "Cost of MST: " << cost_tree << "\r\n";

		std::cout << "\r\n";
	}



};

int _tmain(int argc, _TCHAR* argv[])
{
	srand(time(0));
//create the random graph
	Graph* graph = new Graph(40,50);
//create WeightedEdgeGraph based on random graph
	WeightedEdgeGraph *weg = new WeightedEdgeGraph(graph);
//create WeightedEdgeGraph based on file from site
	WeightedEdgeGraph *inp_weg = new WeightedEdgeGraph("graph.txt");
// Apply Prim's algorithm to random graph
	Prim_MST random_graph_mst(weg);
// Apply Prim's algorithm to graph from site
	Prim_MST inp_graph_mst(inp_weg);

//output result for both graphs
	std::cout << *inp_weg << "\r\n";
	random_graph_mst.out_results();
	std::cout << "\r\n";
	inp_graph_mst.out_results();
	std::cout << "\r\n";


	int i;
	std::cin >> i;
	return 0;
}

