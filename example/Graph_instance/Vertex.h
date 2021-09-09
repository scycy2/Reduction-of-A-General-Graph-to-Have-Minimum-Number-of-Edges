#ifndef INSTANCE_VERTEX_H_
#define INSTANCE_VERTEX_H_

#include <string>
#include <vector>
#include "Rect.h"
using namespace std;

class Vertex {

private:
	string name;
	vector<int> neighbour;
    
public:
    double x;
    double y;
    int vid;
    int num_of_neighbour;
    Rect rect;
    
	Vertex(string n);
	string getName();
	vector<int> getNeighbour();
	void addNeighbour(int new_neighbour);
	int getNumberOfNeighbour();
    void assignRect();
    int getNameLength();
};

#endif /* INSTANCE_VERTEX_H_ */
