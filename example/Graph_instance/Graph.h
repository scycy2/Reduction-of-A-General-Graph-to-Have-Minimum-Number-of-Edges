#ifndef INSTANCE_GRAPH_H_
#define INSTANCE_GRAPH_H_

#include "Vertex.h"

class Graph {

private:
	vector<Vertex *> vertices;
	int num_of_vertices;

public:
    int num_of_edges;
	vector<Vertex *> getVertexList();
	void setVertexList(vector<Vertex *> list_of_vertex);
	int getNumberOfVertices();
};

#endif /* INSTANCE_GRAPH_H_ */


