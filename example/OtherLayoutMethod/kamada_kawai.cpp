#include "kamada_kawai.hpp"

#define MAX_VERTEX_ITERS_COUNT 50


KamadaKawai::KamadaKawai(Graph *g, double k, double energy_threshold) {
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

vector<vector<int>> KamadaKawai::floyd_warshall(Graph *g) {
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
void KamadaKawai::kk(vector<Point*> positions) {
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
    }
}

/* Find max_energy_v_id with the most potential energy and return its energy
 */
double KamadaKawai::find_max_vertex_energy(vector<Point*> positions, int *max_energy_v_id) {
    double max_energy = -1.0;
    for (int v_id = 0; v_id < this->g->getNumberOfVertices(); v_id++) {
        double energy = compute_vertex_energy(v_id, positions);
        if (energy > max_energy) {
            *max_energy_v_id = v_id;
            max_energy = energy;
        }
    }
    return max_energy;
}

/* return the potential energies of springs between v_id and all other vertices */
double KamadaKawai::compute_vertex_energy(int v_id, vector<Point*> positions) {
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
void KamadaKawai::compute_next_vertex_position(int v_id, vector<Point*> positions){
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
//    printf("%f\n", positions[v_id].x);
    
    positions[v_id]->x += (xy_energy * y_energy - yy_energy * x_energy) / denom;
//    printf("%f\n", positions[v_id].x);
    
    positions[v_id]->y += (xy_energy * x_energy - xx_energy * y_energy) / denom;
//    return position;
}
