#include "Newkk.hpp"
#include <iostream>

#define MAX_VERTEX_ITERS_COUNT 50
float k = 0.5;
int times = 0;


Newkk::Newkk(Graph *g, double k, double energy_threshold) {
    vector<vector<int>> distances = floyd_warshall(g);
    this->g = g;
    this->energy_threshold =  energy_threshold;
    int size = g->getNumberOfVertices();
    
    // find biggest distance
    unsigned int biggest_distance = 0;
    for (int v_id = 0; v_id < size; v_id++) {
        for (int other_id = 0; other_id < size; other_id++) {
            if (distances[v_id][other_id] > biggest_distance) {
                biggest_distance = distances[v_id][other_id];
            }
        }
    }

    // Ideal length for all edges. we don't really care, the layout is going to be scaled.
    // Let's chose 1.0 as the initial positions will be on a 1.0 radius circle, so we're
    // on the same order of magnitude
    double length = 1.4 / biggest_distance;
    printf("ideal: %f\n", length);

    // initialize springs lengths and strengths matrices
    for (int v_id = 0; v_id < size; v_id++) {
        vector<Spring> v_springs;

        for (int other_id = 0; other_id < size; other_id++) {
            Spring spring;
            if (v_id == other_id) {
                spring.length = 0.0;
                spring.strength = 0.0;
            } else {
                unsigned int distance = distances[v_id][other_id];
                spring.length = distance * length;
                spring.strength = k / (distance * distance);
            }

            v_springs.push_back(spring);
        }
        springs_.push_back(v_springs);
    }
}

vector<vector<int>> Newkk::floyd_warshall(Graph *g) {
    // build adjacency matrix (infinity = no edge, 1 = edge)
    int infinity = std::numeric_limits<unsigned int>::max()/2;
    int size = g->getNumberOfVertices();
    vector<vector<int>> distances(size, vector<int>(size, infinity)); // allocate a 2d vector
    
    for (int v_id = 0; v_id < size; v_id++) {
        distances[v_id][v_id] = 0;
        Vertex *v = g->getVertexList().at(v_id);
        for (int i = 0; i < v->num_of_neighbour; i++) {
            int adj_id = v->getNeighbour().at(i);
            if (adj_id > v_id) {
                distances[v_id][adj_id] = 1;
                distances[adj_id][v_id] = 1;
            }
        }
    }

    // floyd warshall itself, find length of shortest path for each pair of vertices
    for (int k = 0; k < size; k++) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int d;
                if(distances[i][k] == infinity || distances[k][j] == infinity)
                    d = infinity;
                else
                    d = distances[i][k] + distances[k][j];
                
                if(distances[i][j] > d) {
                    distances[i][j] = d;
                }
            }
        }
    }
//    for(int i=0;i<size;i++){
//        for(int j=0;j<size;j++){
//            printf("%d ", distances[i][j]);
//        }
//        printf("\n");
//    }
    return distances;
}

/* Reduce the energy of the next vertex with most energy
** until all the vertices have a energy below energy_threshold
*/
void Newkk::kk(vector<Point*> positions) {
    int v_id;
    while (find_max_vertex_energy(positions, &v_id) > energy_threshold) {
        // move vertex step by step until its energy goes below threshold
        // (apparently this is equivalent to the newton raphson method)
        unsigned int count = 0;
//        printf("vid:%d\n",v_id);
        do {
            compute_next_vertex_position(v_id, positions);
            count++;
        } while (compute_vertex_energy(v_id, positions) > energy_threshold && count < MAX_VERTEX_ITERS_COUNT);
        
        times++;
        if(times == 500)
            break;
    }
}

/* Find max_energy_v_id with the most potential energy and return its energy
 */
double Newkk::find_max_vertex_energy(vector<Point*> positions, int *max_energy_v_id) {
    double max_energy = -1.0;
    for (int v_id = 0; v_id < this->g->getNumberOfVertices(); v_id++) {
        double energy = compute_vertex_energy(v_id, positions);
        
        if(times >= 50) {
            k = 0.01;
            energy +=  compute_vertex_energy(v_id, positions);
        }
        
        if (energy > max_energy) {
            *max_energy_v_id = v_id;
            max_energy = energy;
        }
    }
    printf("max_energy:%f\n", max_energy);
    return max_energy;
}

/* return the potential energies of springs between v_id and all other vertices */
double Newkk::compute_vertex_energy(int v_id, vector<Point*> positions) {
    double x_energy = 0.0;
    double y_energy = 0.0;

    for (int other_id = 0; other_id < g->getNumberOfVertices(); other_id++) {
        if(v_id == other_id) {
            continue;
        }

        Vector delta;
        delta.dx = positions[v_id]->x - positions[other_id]->x;
        delta.dy = positions[v_id]->y - positions[other_id]->y;
        double distance = delta.norm();

        // delta * k * (1 - l / distance)
        Spring spring = springs_[v_id][other_id];
        x_energy += delta.dx * spring.strength * (1.0 - spring.length / distance);
        y_energy += delta.dy * spring.strength * (1.0 - spring.length / distance);
    }

    return sqrt(x_energy * x_energy + y_energy * y_energy);
}

/* returns next position for v_id reducing its potential energy
** ie the energy in the whole graph caused by its position.
** The position's delta depends on K (TODO bigger K = faster?).
*/
void Newkk::compute_next_vertex_position(int v_id, vector<Point*> positions){
    double xx_energy = 0.0, xy_energy = 0.0, yx_energy = 0.0, yy_energy = 0.0;
    double x_energy = 0.0, y_energy = 0.0;

    for (int other_id = 0; other_id <  g->getNumberOfVertices(); other_id++) {
        if (v_id == other_id) {
            continue;
        }

        Vector delta;
        delta.dx = positions[v_id]->x - positions[other_id]->x;
        delta.dy = positions[v_id]->y - positions[other_id]->y;
        
        double distance = delta.norm();
        double cubed_distance = distance * distance * distance;

        Spring spring = springs_[v_id][other_id];
        
        x_energy += delta.dx * spring.strength * (1.0 - spring.length / distance);
        y_energy += delta.dy * spring.strength * (1.0 - spring.length / distance);
        xy_energy += spring.strength * spring.length * delta.dx * delta.dy / cubed_distance;
        xx_energy += spring.strength * (1.0 - spring.length * delta.dy * delta.dy / cubed_distance);
        yy_energy += spring.strength * (1.0 - spring.length * delta.dx * delta.dx / cubed_distance);
    }
    yx_energy = xy_energy;

    double denom = xx_energy * yy_energy - xy_energy * yx_energy;
    double dx = (xy_energy * y_energy - yy_energy * x_energy) / denom;
    double dy = (xy_energy * x_energy - xx_energy * y_energy) / denom;
    
    if(times >= 50) {
        k = 0.2;
    }
    else {
        k = 0.5;
    }
    double mag = sqrt(dx*dx + dy*dy);
    Vector *edgeRepulsive = compute_edgeRepulsion(v_id, positions);
    dx = dx + edgeRepulsive->dx;
    dy = dy + edgeRepulsive->dy;
    positions[v_id]->x += dx;
    positions[v_id]->y += dy;
    
//    printf("%lf %lf %lf denom:%f\n", dx, dy, sqrt(dx*dx + dy*dy), denom);
//    printf("d:%lf %lf\n", edgeRepulsive->dx, edgeRepulsive->dy);
    
//    positions[v_id]->x += (xy_energy * y_energy - yy_energy * x_energy) / denom;
////    printf("%f\n", positions[v_id].x);
//
//    positions[v_id]->y += (xy_energy * x_energy - xx_energy * y_energy) / denom;
//    return position;
}

//double calStrightDis(double x0, double y0, double x1, double y1, double x2, double y2) {
//    double c = (x0 - x1) * (y2 - y1) + (x1 - x2) * (y0 - y1);
//    c *= c;
//    c /= (pow(y2 - y1, 2) + pow(x1 - x2, 2));
//    return sqrt(c);
//}

//// check if the orthogonal line of the point to the edge intersect the edge
//// if not, the crossing point is outside the edge
//double isPointInEdgeRange(Point *p, Point *a, Point *b, double range) {
//    double vx = a->x - b->x;
//    double vy = a->y - b->y;
//    double norm = sqrt(vx*vx + vy*vy);
////    cout <<"\n"<< vx << " ?"<< vy  << " ?"<<norm <<endl;
////    // vector ab is (vx, vy), so the perpendicular vectors: (vy, -vx) (-vy, vx)
//
//    Point A, B, C, D;
//    A.x = a->x + vy/norm*range; A.y = a->y - vx/norm*range;
//    B.x = a->x - vy/norm*range; B.y = a->y + vx/norm*range;
//    C.x = b->x - vy/norm*range; C.y = b->y + vx/norm*range;
//    D.x = b->x + vy/norm*range; D.y = b->y - vx/norm*range;
////    cout << A.x << B.x << C.x << D.x <<endl;
////    cout << calStrightDis(A.x, A.y, a->x, a->y, b->x, b->y) <<endl;
//
//    float x = p->x;
//    float y = p->y;
//    float a1 = (B.x - A.x)*(y - A.y) - (B.y - A.y)*(x - A.x);
//    float b1 = (C.x - B.x)*(y - B.y) - (C.y - B.y)*(x - B.x);
//    float c1 = (D.x - C.x)*(y - C.y) - (D.y - C.y)*(x - C.x);
//    float d1 = (A.x - D.x)*(y - D.y) - (A.y - D.y)*(x - D.x);
//    if((a1 >= 0 && b1 >= 0 && c1 >= 0 && d1 >= 0) || (a1 <= 0 && b1 <= 0 && c1 <= 0 && d1 <= 0)) {
//        return calStrightDis(p->x, p->y, a->x, a->y, b->x, b->y);
//    }
//    return -1;
//}




Vector* Newkk::compute_edgeRepulsion(int vid, vector<Point*> positions) {
    double l = 0.30;
    
    Vector *vec = new Vector;
    vec->dx = 0;
    vec->dy = 0;
    
    Vertex *v = g->getVertexList().at(vid);
    
    // find all edges
    for(int i=0;i<g->getVertexList().size();i++) {
        Vertex *v1 = g->getVertexList().at(i);
        if(i == vid) continue;
        
        for (int j = 0; j < v1->num_of_neighbour; j++) {
            
            int v2id = v1->getNeighbour().at(j);
            Vertex *v2 = g->getVertexList().at(v2id);
             if(v2id == vid) continue;
            
            double straightdis = isPointInEdgeRange(positions[vid], positions[i], positions[v2id], l);
            if(straightdis == -1) {
                continue;
            }
            else {
                //  GetFootOfPerpendicular
                Point retVal;
                double dx = positions[i]->x - positions[v2id]->x;
                double dy = positions[i]->y - positions[v2id]->y;
                if(abs(dx) < 0.00001 && abs(dy) < 0.0001 ) {
                    retVal.x = positions[i]->x; retVal.y = positions[i]->y;
                }
                double u = (positions[vid]->x - positions[i]->x)*(positions[i]->x - positions[v2id]->x) +
                    (positions[vid]->y - positions[i]->y)*(positions[i]->y - positions[v2id]->y);
                u = u/((dx*dx)+(dy*dy));
                retVal.x = positions[i]->x + u*dx;
                retVal.y = positions[i]->y + u*dy;
             
                // distance == 0
                if(abs(straightdis) <= 0.0000001) {
                    double norm = sqrt(dx * dx + dy * dy);
                    dx = -dy/norm;
                    dy = dx/norm;
                }
                else {
                    // the vector from foot point to the vertex
                    dx = (positions[i]->x - retVal.x)/straightdis;
                    dy = (positions[i]->y - retVal.y)/straightdis;
                }
                
                float c = 0.01;
                double repulsion = c /(straightdis*straightdis + k);
//                printf("%f %f\n", positions[i]->x, retVal.x);
//                printf("%f %f\n", positions[i]->y, retVal.y);
//                printf("straightdis:%f %f %f %f %f\n",straightdis, repulsion, positions[i]->x - retVal.x, positions[i]->y - retVal.y, sqrt(dx*dx + dy*dy) );

                vec->dx += dx * repulsion;
                vec->dy += dx * repulsion;
//
//                printf("d:%f %f %f\n", repulsion, vec->dx, vec->dy);
                
//                // test
//                Point *p = new Point[3];
//                p[0].x = 2.8;
//                p[0].y = 0;
//                p[1].x = 0;
//                p[1].y = 0;
//                p[2].x = 1;
//                p[2].y = 1;
//                double x = isPointInEdgeRange(&p[0], &p[1], &p[2], 1);
//                cout <<"!! "<< x << endl;
            }
        }
    }
    return vec;
}

