#ifndef INSTANCE_READER_DSPINSTANCEREADER_H_
#define INSTANCE_READER_DSPINSTANCEREADER_H_

#include <stdio.h>
#include <string>
#include "../DSPInstance.h"
#include "../Graph.h"
#include "../Vertex.h"
#include <map>

using namespace std;

class DSPInstanceReader {
public:
	 int number_vertex,number_edge;
	 string vertex_name,edgefrom,edgeto,separator;
	 DSPInstance * instance;
	 Graph * graph;
	 vector <Vertex *>  vertexs;
	 Vertex * v;
	 map<string , int > map;
    
    int sign_pos; // if it contains the position information, sign_pos = 1; else 0.

    DSPInstance* readDSPInstance(string filename);
    void outputDSPInstance(Graph* g, string filename, double angle, double d_p2e, double d_p2p, double f, double time_interval) ;
};


#endif /* INSTANCE_READER_DSPINSTANCEREADER_H_ */
