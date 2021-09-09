#include "Graph.h"
#include <iostream>
vector<Vertex *> Graph::getVertexList() {
	return vertices;
}

void Graph::setVertexList(vector<Vertex *> list_of_vertex) {

	vertices = list_of_vertex;
	num_of_vertices = vertices.size();

}

int Graph::getNumberOfVertices() {
	return num_of_vertices;
}


