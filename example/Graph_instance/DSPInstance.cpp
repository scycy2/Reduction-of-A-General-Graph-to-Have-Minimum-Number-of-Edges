#include "DSPInstance.h"
#include <random>
#include <iostream>
#include <fstream>

void DSPInstance::setgraph(Graph * thegraph){
	graph = thegraph;
}

Graph* DSPInstance::getgraph(){
	return graph;
}
