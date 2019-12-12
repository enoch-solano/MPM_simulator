#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>
#include <boost/tokenizer.hpp>

#include "mesh_query0.1/mesh_query.h"

#include "SimulationDriver.h"

int import_mesh(std::vector<MeshObject*> &meshes, char *filename) {
    std::ifstream in(filename);
    if(!in) {
        std::cout << "Cannot open input file.\n";
        return 0;
    }

    boost::char_separator<char> space_sep(" ");

    std::string line;
    while (std::getline(in, line)) {
        boost::tokenizer<boost::char_separator<char>> tokens(line, space_sep);
        if (tokens.begin() == tokens.end()) {
            printf("skipping: %s", line.c_str());
            continue;
        }

        for (const auto& t : tokens) {
            std::cout << t << "\n";
        }

    }

    in.close();

    return 1;
}

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    std::vector<MeshObject*> meshes;

    if (argc > 1) {
        printf("importing meshes...\n");

        for (int i = 1; i < argc; i++) {
            printf("importing mesh with the file path: %s\n", argv[i]);
            import_mesh(meshes, argv[i]);
        }

        printf("finished importing meshes.\n");
    }

    int num_vert = 4;
    std::vector<double> vert = { 0, 0, 0,
                                 0, 1, 0,
                                 1, 1, 0,
                                 1, 0, 0 };
    int num_tri = 2;
    std::vector<int> tri = { 1, 2, 3, 1, 3, 4 };
    MeshObject* m = construct_mesh_object(num_vert, vert.data(), num_tri, tri.data());

    // grid parameters
    TV min_corner = TV::Zero();
    TV max_corner = TV::Ones();
    T dx = 0.05;

    T E = 1e4;
    T nu = 0.3;

    SimulationDriver<T,dim> driver(min_corner, max_corner, dx, E, nu);

    // simulate
    //driver.run(240);

    destroy_mesh_object(m);
    return 0;
}
