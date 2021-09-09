#include "main.hpp"
#include "geometricFunctions.hpp"
#include <GL/glut.h>

//vector<string> split(const string& str, const string& pattern)
//{
//    vector<string> ret;
//    if(pattern.empty()) return ret;
//    size_t start=0,index=str.find_first_of(pattern,0);
//    while(index!=str.npos)
//    {
//        if(start!=index)
//            ret.push_back(str.substr(start,index-start));
//        start=index+1;
//        index=str.find_first_of(pattern,start);
//    }
//    if(!str.substr(start).empty())
//        ret.push_back(str.substr(start));
//    return ret;
//}

void init()
{
    // load graph instance
    DSPInstanceReader a; // /result/test10_graphviz
//    a.readDSPInstance("/Users/apple/Desktop/example/Graph_instance/My_dataset/Group1/test0.txt");
    a.readDSPInstance("/Users/apple/Desktop/test1.txt");
    g = a.graph;

    // allocate memory
    drawing = initializeDrawing(g, 1);
    
//    ouputDotFile(g, "/Users/wangdanyun/Desktop/3.1g.dot");
//    if(a.sign_pos == 1) {
//        int vsize = g->getNumberOfVertices();
//        vector<Vertex*> vlist = g->getVertexList();
//
//        // copy the positions
//        for(int i=0; i<vsize; ++i) {
//            drawing[i]->x = vlist[i]->x;
//            drawing[i]->y = vlist[i]->y;
//        }
//    }
//
//    drawing = readpPlainFile(g, "/Users/wangdanyun/Desktop/DrawSphere-master/example的副本/Graph_instance/My_dataset/Group1/test0_graphvaiz.plain");
//    g = generate_Kn_graph(6);
//    int k = 4;
//    int h = 4;
//    g = generate_k_ary_tree(k, h);
//    experiement1();
//    experiement2();
    
    curt_drawing = &init_drawing;
    init_drawing = initializeDrawing(g, 1);
    copyDrawing(init_drawing, drawing);

    // set showing parameters
    scale = g->getNumberOfVertices()/100 + 2;
    psize = 5;
    Xplace = 1;
    Yplace = 1;
    keyboard_move = 0.5;
    
    glClearColor(1.f, 1.f, 1.f, 1.f); // set white as background
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_DITHER);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


void draw_graph_algorithm() {
//    not use
//    KamadaKawai *ka = new KamadaKawai(g, 500.0, 0.0001);
//    ka->kk(drawing);
//    Newkk *ka = new Newkk(g, 500.0, 20);
//    ka->kk(drawing);
    
    
    // calculate positions
    clock_t start_time=clock();
    
    LayoutAlgorithm fr(g);
    fr.FR(drawing, 100);

    clock_t end_time=clock();
    double runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC; // the run time for the algorithm
    cout << "time: "<< runtime <<endl;
    
//    printAllpoints(drawing);
    printf("evaluate_drawing:%f\n", evaluate_drawing(drawing, g));
}

int main(int argc,char* argv[])
{
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(600,600);
    glutCreateWindow("Drawing graph");
       init();
    glutDisplayFunc(display);
       draw_graph_algorithm();
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    
    /* ------------ */
    glPushMatrix();
    glTranslatef(-Xplace, -Yplace, 0);
    glScalef(scale, scale, 1);
    

    // draw
    drawGraph(*curt_drawing, g);
    
    glPopMatrix();
    glutSwapBuffers();
}

void drawPoint(GLfloat x, GLfloat y, GLfloat size){
    glPointSize(size);
    glBegin(GL_POINTS);
    glColor3f(0.f, 0.f, 0.f);
    glVertex2f(x, y);
    glEnd();
}

void drawRect(GLfloat x, GLfloat y, GLfloat l, GLfloat w) {
    glBegin(GL_LINES);
    glColor3f(0.f, 0.f, 0.f);
    glVertex2f(x, y);
    glVertex2f(x+l, y);
    glVertex2f(x+l, y);
    glVertex2f(x+l, y+w);
    glVertex2f(x+l, y+w);
    glVertex2f(x, y+w);
    glVertex2f(x, y+w);
    glVertex2f(x, y);
    glEnd();
    
//    glClear(GL_COLOR_BUFFER_BIT);
//    glRectf(x, y, x+0.025f, y+0.025f);
    glColor4f(0.f, 0.75f, 0.75f, 0.5);
    glBegin(GL_POLYGON);
    glVertex2f(x, y);
    glVertex2f(x+l, y);
    glVertex2f(x+l, y+w);
    glVertex2f(x, y+w);
    glEnd();
//    glFlush();
}

void drawLabel(GLfloat x, GLfloat y, string name) {
    glColor3f(0, 0, 0);
    glRasterPos2f(x + 0.01f, y + 0.01f);
    char nameArray[] = {0};
    strcpy(nameArray, name.c_str());
    for (int i = 0; i < name.length(); ++i) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, nameArray[i]);
//        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, '\n');
    }
}

void drawGraph(vector<Point*> position, Graph *g){
    int vsize = g->getNumberOfVertices();
    vector<Vertex *> vertices = g->getVertexList();
    
    for(int i=0;i<vsize;i++){
        drawLabel(position[i]->x, position[i]->y, vertices[i]->getName());
        Vertex *v = g->getVertexList().at(i);
        drawRect(position[i]->x, position[i]->y, v->getNameLength() * 0.04f, 0.04f);
    }
    
    // draw edges
    for(int i=0;i<vsize;i++) {
//        drawPoint(position[i]->x, position[i]->y, 0.1);
        Vertex *v = g->getVertexList().at(i);
        long nsize = v->getNeighbour().size();
        
        for(int j=0;j<nsize;j++) {
            int neibor = v->getNeighbour().at(j);
            
            glBegin(GL_LINES);
            glColor3f(0.f, 0.f, 0.f);
            glVertex2f(position[i]->x + v->getNameLength() * 0.04f/2, position[i]->y);
            glVertex2f(position[neibor]->x + g->getVertexList().at(neibor)->getNameLength() * 0.04f/2, position[neibor]->y);
            glEnd();
        }
//        drawRect(position[i]->x, position[i]->y, v->getNameLength() * 0.04f, 0.04f);
    }
}

void reshape(int w, int h) {
    GLfloat aspect = (GLfloat)w / (GLfloat)h;
    GLfloat nRange = 2.0f;
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (w<=h)
    {   glOrtho(-nRange, nRange, -nRange * aspect, nRange * aspect, -nRange, nRange);  }
    else
    {  glOrtho(-nRange, nRange, -nRange / aspect, nRange / aspect, -nRange, nRange);  }
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'w':
            Yplace += keyboard_move;
            glutPostRedisplay();
            break;
        case 's':
            Yplace -= keyboard_move;
            glutPostRedisplay();
            break;
        case 'a':
            Xplace -= keyboard_move;
            glutPostRedisplay();
            break;
        case 'd':
            Xplace += keyboard_move;
            glutPostRedisplay();
            break;
        case 'j':
            scale += 1;
            psize += 1;
            glutPostRedisplay();
            break;
        case 'k':
            scale -= 1;
            psize -= 1;
            glutPostRedisplay();
            break;
        case 'n':
            if(curt_drawing == &drawing)
                curt_drawing = &init_drawing;
            else
                curt_drawing = &drawing;
            glutPostRedisplay();
            break;
        case 27:
            exit(0);
        default:
            break;
    }
}

void copyDrawing(vector<Point*> p, vector<Point*> des) {
    for(int i=0; i<p.size();i++) {
        p[i]->x = des[i]->x;
        p[i]->y = des[i]->y;
    }
}

