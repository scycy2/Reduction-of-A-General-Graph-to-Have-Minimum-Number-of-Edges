#ifndef kamada_kawai_hpp
#define kamada_kawai_hpp

#include <stdio.h>
#include <vector>
#include "DSPInstanceReader.h"
#include "geometricFunctions.hpp"
//#include "main.hpp"

class KamadaKawai {
public:
    struct Spring {
        double length;
        double strength;
    };
    
    Graph *g;
    double energy_threshold;
    vector<vector<Spring>> springs_;
    
    KamadaKawai();
    KamadaKawai(Graph *g, double k, double energy_threshold);
    
    static vector<vector<int>> floyd_warshall(Graph *g);
    void kk(vector<Point*> positions);
    double find_max_vertex_energy(vector<Point*> positions, int *max_energy_v_id);
    double compute_vertex_energy(int v_id, vector<Point*> positions);
    void compute_next_vertex_position(int v_id, vector<Point*> positions);
};


#endif /* kamada_kawai_hpp */
