#include "Vertex.h"

Vertex::Vertex(string n)
:name(n)
{
    num_of_neighbour = 0;
}

string Vertex::getName() {

	return name;
}

vector<int> Vertex::getNeighbour() {

	return neighbour;
}

void Vertex::addNeighbour(int new_neighbour) {
	neighbour.push_back(new_neighbour);
	num_of_neighbour = neighbour.size();
}

int Vertex::getNumberOfNeighbour() {

	return num_of_neighbour;
}

void Vertex::assignRect() {
    int l = name.length();
    rect = Rect(x-5, y-5, x+5*l, y+5);
}

int Vertex::getNameLength() {
    int n = 0;
    for (int i = 0; i < name.length(); i++) {
        if (name.at(i) != ' ') {
            n++;
        }
    }
    return n;
}
