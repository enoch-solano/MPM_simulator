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

void triangulate(std::vector<int> &idx, std::vector<int> &face_idx) {
    for (int i = 0; i < face_idx.size()-2; i++) {
        idx.push_back(face_idx[0]);
        idx.push_back(face_idx[i+1]);
        idx.push_back(face_idx[i+2]);
    }
}

int import_mesh(std::vector<MeshObject*> &meshes, char *filename) {
    std::ifstream in(filename);
    if(!in) {
        std::cout << "Cannot open input file.\n";
        return 0;
    }

    boost::char_separator<char> space_sep(" ");

    std::vector<int> idx;
    std::vector<double> verts;

    std::string line;
    while (std::getline(in, line)) {
        boost::tokenizer<boost::char_separator<char>> tokens(line, space_sep);

        // continues if there are no tokens in this line
        if (tokens.begin() == tokens.end()) {
            continue;
        }

        // checks to see if line corresponds to a face line
        if ((*tokens.begin()).compare("f") == 0) {
            std::vector<int> face_idx;
            for (auto tok = ++tokens.begin(); tok != tokens.end(); ++tok) {
                boost::char_separator<char> sep("/");
                boost::tokenizer<boost::char_separator<char>> tmp_tokens(*tok, sep);

                // gets index of vertex position (minus 1 because objs are 1-indexed)
                face_idx.push_back(std::stoi(*tmp_tokens.begin()) - 1);
            }

            // triangulates face vertices as needed
            triangulate(idx, face_idx);
        }

        // checks to see if token corresponds to a vertex line
        if ((*tokens.begin()).compare("v") == 0) {
            int vert_count = 0;
            for (auto tok = ++tokens.begin(); tok != tokens.end(); ++tok) {
                // doesn't allow more than three values per vertex
                if (vert_count == 3) {
                    continue;
                }

                verts.push_back(std::stod(*tok));
                vert_count++;
            }
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
