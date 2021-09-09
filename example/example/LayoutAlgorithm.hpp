#ifndef LayoutAlgorithm_hpp
#define LayoutAlgorithm_hpp

#include <stdio.h>
#include <vector>
#include "geometricFunctions.hpp"

class LayoutAlgorithm {
    
public:
    LayoutAlgorithm(Graph *g);
    void frOperator_(vector<Point*> positions, int s);
    void FR(vector<Point*> positions, int iters_count);

private:
    Graph *g;
    double k;
    double k_squared;
    double temp_;
    vector<Vector> disp;
};

#endif /* LayoutAlgorithm_hpp */
