#ifndef main_hpp
#define main_hpp

#define GLUT_DISABLE_ATEXIT_HACK
#include <stdio.h>
#include "DSPInstanceReader.h"
#include "kamada_kawai.hpp"
#include "Newkk.hpp"
#include "LayoutAlgorithm.hpp"
#include "GenerateGraph.hpp"
#include <iostream>
#include <fstream>
#include <GL/glut.h>

#endif /* main_hpp */


using namespace std;

GLfloat Xplace;
GLfloat Yplace;
GLfloat keyboard_move;
int scale;
int psize;
Graph *g;
vector<Point*> drawing;
vector<Point*> init_drawing;
vector<Point*> *curt_drawing;

void drawGraph(vector<Point*> positions, Graph *g);
void drawPoint(GLfloat x, GLfloat y, GLfloat size);
void drawRect(GLfloat x, GLfloat y);
void draw_graph_algorithm();

void keyboard(unsigned char key, int x, int y);
void reshape(int w, int h);
void display(void);
void copyDrawing(vector<Point*> p, vector<Point*> des);
void printAllpoints(vector<Point*> positions);
