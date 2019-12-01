#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "mesh_query0.1/mesh_query.h"

#include "SimulationDriver.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    int num_vert = 2;
    const double vert[] = {1, 2, 3};
    int num_tri = 5;
    const int tri[] = {1, 2, 3, 1, 3, 2};
    MeshObject* m = construct_mesh_object(num_vert, vert, num_tri, tri);
    int as = 2;
    as--;

    // grid parameters
    TV min_corner = TV::Zero();
    TV max_corner = TV::Ones();
    T dx = 0.05;

    T E = 1e4;
    T nu = 0.3;

    SimulationDriver<T,dim> driver(min_corner, max_corner, dx, E, nu);

    // simulate
    driver.run(240);

    return 0;
}
