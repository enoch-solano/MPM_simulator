#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>



template<class T, int dim>
class GridNode {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    T  m;   // node mass
    TV vn;  // node transfered velocity
    TV v;   // node updated velocity
    TV f;   // node force

    GridNode() : m((T) 0), v(TV::Zero()), vn(TV::Zero()), f(TV::Zero()) {}

    GridNode(const GridNode &g) : m(g.m), v(g.v), vn(g.vn), f(g.f) {}

    void print() {
        printf("m: %f\n", m);
        printf("vn.x: %f, vn.y: %f, vn.z: %f\n", vn(0), vn(1), vn(2));
        printf("v.x: %f, v.y: %f, v.z: %f\n",     v(0),  v(1),  v(2));
        printf("f.x: %f, f.y: %f, f.z: %f\n\n",   f(0),  f(1),  f(2));
    }
};

template<class T, int dim>
class Grid {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using IV = Eigen::Matrix<int,dim,1>;

    TV min_corner, max_corner;
    T dx;   // width between nodes

    IV res; // resolution in each direction
    int num_nodes;

    std::vector<GridNode<T,dim>> grid;
    std::vector<int> active_nodes;      // indices of active nodes (i.e. non-zero mass)

    Grid(TV min, TV max, T dx)
        : min_corner(min), max_corner(max), dx(dx), num_nodes(1)
    {
        TV tmp_res = (max_corner - min_corner) / dx + TV::Ones();

        for (int i = 0; i < dim; i++) {
            res(i) = static_cast<int>(tmp_res(i));
            num_nodes *= res(i);
        }

        grid = std::vector<GridNode<T,dim>>(num_nodes, GridNode<T,dim>());
    }

    // gets mass of grid node (i,j,k)
    T m(int i, int j, int k) const {
        return grid[to1D(i,j,k)].m;
    }

    T& m(int i, int j, int k) {
        return grid[to1D(i,j,k)].m;
    }

    T m(int i) const {
        return grid[i].m;
    }

    T& m(int i) {
        return grid[i].m;
    }

    // gets updated velocity of grid node (i,j,k)
    TV v(int i, int j, int k) const {
        return grid[to1D(i,j,k)].v;
    }

    TV& v(int i, int j, int k) {
        return grid[to1D(i,j,k)].v;
    }

    TV v(int i) const {
        return grid[i].v;
    }

    TV& v(int i) {
        return grid[i].v;
    }

    // gets transfered velocity of grid node (i,j,k)
    TV vn(int i, int j, int k) const {
        return grid[to1D(i,j,k)].vn;
    }

    TV& vn(int i, int j, int k) {
        return grid[to1D(i,j,k)].vn;
    }

    TV vn(int i) const {
        return grid[i].vn;
    }

    TV& vn(int i) {
        return grid[i].vn;
    }

    // gets force of grid node (i,j,k)
    TV f(int i, int j, int k) const {
        return grid[to1D(i,j,k)].f;
    }

    TV& f(int i, int j, int k) {
        return grid[to1D(i,j,k)].f;
    }

    TV f(int i) const {
        return grid[i].f;
    }

    TV& f(int i) {
        return grid[i].f;
    }

    // helper function to convert 3D index to 1D
    int to1D(int i, int j, int k) {
        return i + j * res(0) + k * res(0) * res(1);
    }

    // clears entire grid
    void clear() {
        grid = std::vector<GridNode<T,dim>>(num_nodes, GridNode<T,dim>());
        active_nodes.clear();
    }

    void print() {
        for (auto g : grid) {
            g.print();
        }

        print_active_nodes();
    }

    void print_active_nodes() {
        for (auto i : active_nodes) {
            grid[i].print();
        }
        printf("number of active nodes: %d\n", (int) active_nodes.size());
    }
};

template<class T, int dim>
class Particle {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    T  m;   // particle mass
    TV x;   // particle position
    TV v;   // particle velocity
    T vol;  // particle volume

    std::array<int, dim> base_node;         // index of base node
    std::array<T, 3*dim> weights;             // weights of nodes (in 1d along all dim)
    std::array<T, 3*dim> weight_gradients;    // weight gradients of nodes (in 1d along all dim)


    std::array<T, dim*dim*dim> weights_ijk;             // weights of nodes in 3D
    std::array<T, dim*dim*dim> weight_gradients_ijk;    // weight gradients of nodes in 3D

    Particle() : m((T) 0), x(TV::Zero()), v(TV::Zero()), vol((T) 0) {
        for (int i = 0; i < dim; i++) {
            base_node[i] = 0;
        }

        for (int i = 0; i < 3*dim; i++) {
            weights[i] = 0;
            weight_gradients[i] = 0;
        }

        for (int i = 0; i < dim*dim*dim; i++) {
            weights_ijk[i] = 0;
            weight_gradients_ijk[i] = 0;
        }
    }

    Particle(TV x, T m, T vol)
        : m(m), x(x), v(TV::Zero()), vol(vol)
        {
            for (int i = 0; i < dim; i++) {
                base_node[i] = 0;
            }

            for (int i = 0; i < 3*dim; i++) {
                weights[i] = 0;
                weight_gradients[i] = 0;
            }

            for (int i = 0; i < dim*dim*dim; i++) {
                weights_ijk[i] = 0;
                weight_gradients_ijk[i] = 0;
            }
        }

    Particle(const Particle &p)
        : m(p.m), x(p.x), v(p.v), vol(p.vol)
        {
            for (int i = 0; i < dim; i++) {
                base_node[i] = p.base_node[i];
            }

            for (int i = 0; i < 3*dim; i++) {
                weights[i] = p.weights[i];
                weight_gradients[i] = p.weight_gradients[i];
            }

            for (int i = 0; i < dim*dim*dim; i++) {
                weights_ijk[i] = p.weights_ijk[i];
                weight_gradients_ijk[i] = p.weight_gradients_ijk[i];
            }
        }


    // return base index in dimension 'dim'
    int b(int d) const {
        return base_node[d];
    }

    int& b(int d) {
        return base_node[d];
    }

    // returns weight of node 'n' along dimension 'd'
    T w(int d, int n) const {
        return weights[d*dim + n];
    }

    T& w(int d, int n) {
        return weights[d*dim + n];
    }

    // returns weight gradient of node 'n' along dimension 'd'
    T dw(int d, int n) const {
        return weight_gradients[d*dim + n];
    }

    T& dw(int d, int n) {
        return weight_gradients[d*dim + n];
    }

    // compute 1D quadratic B spline weights
    // x_idx is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
    void compute_weights_1D(int d, T x_idx) {
        T base_node_i = std::floor(x_idx - 0.5) + 1;
        b(d) = static_cast<int>(base_node_i) - 1;

        T d0 = x_idx - base_node_i + 1;
        T z = 1.5 - d0;
        T z2 = z * z;
        w(d, 0) = 0.5 * z2;

        T d1 = d0 - 1;
        w(d, 1) = 0.75 - d1 * d1;

        T d2 = 1 - d1;
        T zz = 1.5 - d2;
        T zz2 = zz * zz;
        w(d, 2) = 0.5 * zz2;
    }

    T w_(int i, int j, int k) const {
        return weights_ijk[to1D(i,j,k)];
    }

    T& w_(int i, int j, int k) {
        return weights_ijk[to1D(i,j,k)];
    }

    // helper function to convert 3D index to 1D
    int to1D(int i, int j, int k) {
        return i + j * 3 + k * 3 * 3;
    }
};

template<class T, int dim>
class MPMSystem {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    T youngs_modulus;
    T damping_coeff;

    Grid<T, dim> grid;
    std::vector<Particle<T, dim>> particles;

    MPMSystem(TV min, TV max, T dx)
        : grid(min, max, dx) {}

    void populate_cube(TV min, TV max) {
        T rho = 1000;   // density
        T vol = (grid.dx * grid.dx) / 8;
        T mass = vol * rho;

        for (int i = 1; i < grid.res(0)-1; i++) {
            for (int j = 1; j < grid.res(1)-1; j++) {
                for (int k = 1; k < grid.res(2)-1; k++) {
                    // node position
                    TV n_x = grid.min_corner + TV(i,j,k) * grid.dx;

                    // stratified sampling
                    for (int n = 0; n < 8; n++) {
                        TV p_x = n_x;

                        int is_within_box = 1;

                        for (int d = 0; d < dim; d++) {
                            T r = (((T) rand()) / RAND_MAX) * grid.dx;
                            r -= grid.dx / 2;
                            p_x(d) += r;

                            is_within_box = is_within_box && (p_x(d) > min(d) && p_x(d) < max(d));
                        }

                        if (is_within_box) {
                            Particle<T,dim> p(p_x, mass, vol);
                            particles.push_back(p);
                        }
                    }
                }
            }
        }
    }



    void compute_weights_3D() {
        for (auto &p : particles) {
            TV x_idx = p.x / grid.dx;

            for (int d = 0; d < dim; d++) {
                p.compute_weights_1D(d, x_idx(d));
            }

            for (int i = 0; i < 3; i++) {
                T w_i = p.w(0, i);

                for (int j = 0; j < 3; j++) {
                    T w_ij = w_i * p.w(1, j);

                    for (int k = 0; k < 3; k++) {
                        T w_ijk = w_ij * p.w(2, k);

                        p.w_(i,j,k) = w_ijk;
                    }
                }
            }
        }

        for (auto &p : particles) {
            TV pos = TV::Zero();
            TV x_idx = p.x / grid.dx;

            for (int i = 0; i < 3; i++) {
                int n_i = p.b(0) + i;
                for (int j = 0; j < 3; j++) {
                    int n_j = p.b(1) + j;
                    for (int k = 0; k < 3; k++) {
                        int n_k = p.b(2) + k;

                        pos(0) = pos(0) + n_i*p.w_(i,j,k);
                        pos(1) = pos(1) + n_j*p.w_(i,j,k);
                        pos(2) = pos(2) + n_k*p.w_(i,j,k);
                    }
                }
            }

            // printf("position: (%f, %f, %f)\n", x_idx(0), x_idx(1), x_idx(2));
            // printf("interpolated position: (%f, %f, %f)\n\n", pos(0), pos(1), pos(2));
        }
    }


    void transfer_P2G() {
        int count = 0;
        for (auto &p : particles) {
            for (int i = 0; i < 3; i++) {
                int n_i = p.b(0) + i;

                for (int j = 0; j < 3; j++) {
                    int n_j = p.b(1) + j;

                    for (int k = 0; k < 3; k++) {
                        int n_k = p.b(2) + k;

                        // interpolate (splat) mass
                        grid.m(n_i, n_j, n_k)  = grid.m(n_i, n_j, n_k)  + (p.w_(i,j,k) * p.m);

                        // interpolate (splat) momentum
                        grid.vn(n_i, n_j, n_k) = grid.vn(n_i, n_j, n_k) + (p.w_(i,j,k) * p.m * p.v);
                    }
                }
            }
        }

        for (int i = 0; i < grid.num_nodes; i++) {
            if (grid.m(i) == 0) {
                grid.vn(i) = TV::Zero();
            } else {
                grid.active_nodes.push_back(i);
                grid.vn(i) = grid.vn(i) / grid.m(i);
            }
        }
    }

    void addGravity(TV gravity) {
        for (int i : grid.active_nodes) {
            // printf("f before gravity: (%.3f, %.3f, %.3f)\n", grid.f(i)[0], grid.f(i)[1], grid.f(i)[2]);
            grid.f(i) = grid.f(i) + (grid.m(i) * gravity);
            // printf("f after gravity:  (%.3f, %.3f, %.3f)\n\n", grid.f(i)[0], grid.f(i)[1], grid.f(i)[2]);
        }
    }

    void addElasticity() {
        // TODO: implement this
    }

    void updateGridVelocity(T dt) {
        for (int i : grid.active_nodes) {
            grid.v(i) = grid.vn(i) + (dt * grid.f(i) / grid.m(i));
        }
    }

    void setBoundaryVelocities(int thickness) {
        for (int i = 0; i < grid.res(0); i++) {
            for (int j = 0; j < grid.res(1); j++) {
                for (int k = 0; k < thickness; k++) {
                    grid.v(i,j,k) = TV::Zero();
                }
            }
        }

        for (int i = 0; i < grid.res(0); i++) {
            for (int j = 0; j < thickness; j++) {
                for (int k = 0; k < grid.res(2); k++) {
                    grid.v(i,j,k) = TV::Zero();
                }
            }
        }

        for (int i = 0; i < thickness; i++) {
            for (int j = 0; j < grid.res(1); j++) {
                for (int k = 0; k < grid.res(2); k++) {
                    grid.v(i,j,k) = TV::Zero();
                }
            }
        }

        for (int i = grid.res(0)-thickness; i < grid.res(0); i++) {
            for (int j = 0; j < grid.res(1); j++) {
                for (int k = 0; k < grid.res(2); k++) {
                    grid.v(i,j,k) = TV::Zero();
                }
            }
        }

        for (int i = 0; i < grid.res(0); i++) {
            for (int j = grid.res(1)-thickness; j < grid.res(1); j++) {
                for (int k = 0; k < grid.res(2); k++) {
                    grid.v(i,j,k) = TV::Zero();
                }
            }
        }

        for (int i = 0; i < grid.res(0); i++) {
            for (int j = 0; j < grid.res(1); j++) {
                for (int k = grid.res(2)-thickness; k < grid.res(2); k++) {
                    grid.v(i,j,k) = TV::Zero();
                }
            }
        }
    }

    void transfer_G2P(T flip, T dt) {
        for (auto &p : particles) {
            TV vpic = TV::Zero();
            TV vflip = p.v;

            for (int i = 0; i < 3; i++) {
                int n_i = p.b(0) + i;

                for (int j = 0; j < 3; j++) {
                    int n_j = p.b(1) + j;

                    for (int k = 0; k < 3; k++) {
                        int n_k = p.b(2) + k;

                        // interpolates (splat) grid velocity
                        vpic  = vpic  + p.w_(i,j,k) * grid.v(n_i, n_j, n_k);
                        vflip = vflip + p.w_(i,j,k) * (grid.v(n_i, n_j, n_k) - grid.vn(n_i, n_j, n_k));
                    }
                }
            }

            p.v = (1 - flip)*vpic + flip*vflip;
            p.x = p.x + dt * vpic;
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto p : particles) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++) {
                fs << " " << p.x(i);
            }

            if (dim == 2) {
                fs << " 0";
            }

            fs << "\n";
        }

        fs << "POLYS\n";
        fs << "END\n";

        fs.close();
    }

};
