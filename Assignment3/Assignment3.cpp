// Assignment3.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <list>
#include <time.h>


#include <queue>
#include <functional>   // std::greater
#include <algorithm>
#include <map>

class WeightedEdge
{   int v;
    int w;
	double weight;
public:
	WeightedEdge(int v, int w, double weight):v(v), w(w), weight(weight){};
	double get_weight(void) const { return weight;}
	int from(void) const {return v;}
	int to(void) const {return w;}
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



class WeightedEdgeGraph{
	std::list<std::list<WeightedEdge>> *listEdges;
	std::list<double> *vert_values;
public:
	WeightedEdgeGraph(Graph *graph);
	friend std::ostream& operator<<(std::ostream& out,  WeightedEdgeGraph const & gr);
	std::list<std::list<WeightedEdge>>& get(void) const {return *listEdges;}
	std::list<WeightedEdge> edges(int num_vertex);
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
// This function insert pair of number of the point and distanсe to this point from the start
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






int _tmain(int argc, _TCHAR* argv[])
{
	srand(time(0));
	Graph* graph = new Graph(4,50);
	WeightedEdgeGraph *weg = new WeightedEdgeGraph(graph);
	PriorityQueue* pq = new PriorityQueue();

		std::cout << *graph << "\r\n";
	std::cout << *weg << "\r\n";
//points were included in MST
	std::list<int> spanTree;
	std::list<WeightedEdge> spanTree2;
//length of MST
	double length_tree = 0.0;
	for (auto we: weg->edges(0)){
		pq->insert(we);
	}
	spanTree.push_back(0);
// variable for counting edges in minimum spanning tree
	int eBegin = -1;
	while(!pq->empty()){
		WeightedEdge we_h = pq->top();
		pq->minPrioirty();
		int n = 0;
		for(auto i: spanTree){
			if (i == we_h.from())
				n++;
		}
		for(auto i: spanTree){
			if (i == we_h.to())
				n++;
		}
		if (n<2){
			spanTree2.push_back(we_h);
			spanTree.push_back(we_h.to());
			//	length_tree += 
			for(auto we: weg->edges(we_h.to())){

					int compare = 0;
					for (auto i: spanTree2){
	//add to our MST
						if (i == we){
							compare = 1;
							break;
						}
					}
					for (auto i: spanTree){
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

	for (auto i: spanTree){
		std::cout << i << " ";
	}

		std::cout << "\r\n";

	for (auto i: spanTree2){
		length_tree += i.get_weight();
		std::cout << i << " ";
	}
	std::cout << "\r\n";
	std::cout << "Length of MST: " << length_tree << "\r\n";

	std::cout << "\r\n";



	int i;
	std::cin >> i;
	return 0;
}

