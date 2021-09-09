
#ifndef INSTANCE_DSPINSTANCE_H_
#define INSTANCE_DSPINSTANCE_H_

#include "Graph.h"
#include <random>
//#include "../DSPObjectiveFunction.h"

class DSPInstance {

private:
	Graph * graph;
	int *value;

public:
	void setgraph(Graph * graph);
	Graph* getgraph();
};

#endif /* INSTANCE_DSPINSTANCE_H_ */
