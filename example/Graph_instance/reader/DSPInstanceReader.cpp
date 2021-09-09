
#include "DSPInstanceReader.h"
#include <iostream>
#include <fstream>
#include "../Vertex.h"
#include <vector>
#include<stdlib.h>


DSPInstance* DSPInstanceReader::readDSPInstance(string filename){
    ifstream file;
    file.open(filename);
    if (!file) {
        cout << "Failed to open file" << endl;
        exit (0);
    }
    file >> separator;
//    cout << separator << endl;
    file >> separator;
//    cout << separator << endl;
    file >> separator;
//    cout << separator << endl;
    file >> number_vertex;
    file >> separator;
    file >> number_edge;
    file >> separator;

    for( int i = 0; i < number_vertex; i++ ){
        file >> vertex_name;
        for (int j = 0; j < vertex_name.length(); j++) {
            if (vertex_name.at(j) == '_') {
                vertex_name.at(j) = ' ';
            }
        }
        v = new Vertex(vertex_name);
        v->vid = i;
        v->x = 0;
        v->y = 0; // set it to 0
        vertexs.push_back(v);
        map.insert(pair<string, int>(vertex_name,i));
    }
    
//    for (int i = 0; i < vertexs.size(); i++) {
//        cout << vertexs.at(i)->getName() << endl;
//    }
    
    
    
    
    while(separator != "Edge" && separator != "Minimal_Angle:") {
        file >> separator;
    }

    // contains the position info
    if(separator == "Minimal_Angle:") {
        sign_pos = 1;
        file >> separator;
        cout << "Minimal angle:" << separator << endl;
        file >> separator;
        file >> separator;
        cout << "Minimal distance from p2e:" << separator << endl;
        file >> separator;
        file >> separator;
        cout << "Minimal distance from p2p:" << separator << endl;
        file >> separator;
        file >> separator;
        cout << "Value of evaluation:" << separator << endl;
        file >> separator;
        file >> separator;
        cout << "Time:" << separator << endl;
        
        file >> separator;
        cout << separator;
        if(separator == "Pos:") {
//            string a, b;
            double a,b;
            for( int i = 0; i < number_vertex; i++ ){
                file >> a;
                file >> b;
                v = vertexs[i];
                v->x = a;
                v->y = b;
//                printf("%f %f\n", a, b);
            }
        }
        file >> separator;
    }
    else {
        sign_pos = 0;
    }
    
//    cout << sign_pos << endl;
        
    // read the edges
//    file >> separator;
    if(separator != "Edge" ) {
        cout << separator;
        cout << "Error occurs!" << endl;
        exit(0);
    }

    for( int i = 0; i < number_edge; i++ ){
        file >> edgefrom;
        
        for (int j = 0; j < edgefrom.length(); ++j) {
            if (edgefrom.at(j) == '_') {
                edgefrom.at(j) = ' ';
            }
        }
        file >> edgeto;
        for (int j = 0; j < edgeto.length(); ++j) {
            if (edgeto.at(j) == '_') {
                edgeto.at(j) = ' ';
            }
        }
        cout << i << ": " << edgefrom << " to " << edgeto << endl;
        vertexs[map.find(edgefrom)->second]->addNeighbour(map.find(edgeto)->second);
        vertexs[map.find(edgeto)->second]->addNeighbour(map.find(edgefrom)->second);
    }

    graph = new Graph();
    graph->setVertexList(vertexs);
    graph->num_of_edges = number_edge;
    
    for (int i = 0; i < graph->getVertexList().size(); i++) {
        cout << graph->getVertexList().at(i)->getName() << endl;
    }

    instance = new DSPInstance();
    instance->setgraph(graph);
    file.close();
    return instance;
}


void DSPInstanceReader::outputDSPInstance(Graph* g, string filename, double angle, double d_p2e, double d_p2p, double f, double time_interval) {
    ofstream fileout;
    fileout.open(filename);
    if (!fileout) {
        cout << "Failed to create file" << endl;
        exit (0);
    }
    
    fileout << "filename: " << filename << "\n";
    fileout << "number_of_vertex: " << g->getNumberOfVertices() << "\n";
    fileout << "number_of_vertex: " << g->num_of_edges << "\n";

    // vertex list
    fileout << "Vertex:\n";
    for( int i = 0; i < g->getNumberOfVertices(); i++ ){
        fileout << "n" << i << " ";
    }

    fileout << "\nMinimal_Angle:\n" << angle << "\nMinimal_Distance_p2e:\n" << d_p2e << "\nMinimal_Distance_p2p:\n" << d_p2p <<  "\nValue_of_evaluation:\n" << f <<  "\nTime:\n" << time_interval << "\nPos:\n";

    // positions
    vector<Vertex *> vlist = g->getVertexList();
    for( int i = 0; i < g->getNumberOfVertices(); i++ ){
        fileout << vlist[i]->x << " " << vlist[i]->y << " ";
    }
    
    const int size = g->getNumberOfVertices();
    int sign[size][size];
    for( int i = 0; i < size; i++ ){
        for(int j = 0; j < size; j++) {
            sign[i][j] = 0; // initalize
        }
    }
     
    fileout << "\nEdge\n";
    for( int i = 0; i < g->getNumberOfVertices(); i++ ){
        for(int j = 0; j < vlist[i]->getNumberOfNeighbour(); j++) {
            int j_id = vlist[i]->getNeighbour()[j];
            if(sign[i][j_id] == 0 && sign[j_id][i] == 0 ) {
                fileout << "n" << i << " n" << vlist[i]->getNeighbour()[j] << "\n";
                sign[i][j_id] = 1;
                sign[j_id][i] = 1;
            }
        }
    }
    fileout.close();
}
