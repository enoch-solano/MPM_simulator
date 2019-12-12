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

using T = float;
constexpr int dim = 3;
using TV = Eigen::Matrix<T,dim,1>;

void triangulate(std::vector<int> &idx, std::vector<int> &face_idx) {
    for (int i = 0; i < face_idx.size()-2; i++) {
        idx.push_back(face_idx[0]);
        idx.push_back(face_idx[i+1]);
        idx.push_back(face_idx[i+2]);
    }
}

void compute_bounds(TV &min, TV &max, std::vector<double> &positions) {
    for (int v = 0; v < positions.size(); v += 3) {
        for (int i = 0; i < 3; i++) {
            if (min(i) > positions[v+i]) {
                min(i) = positions[v+i];
            }

            if (max(i) < positions[v+i]) {
                max(i) = positions[v+i];
            }
        }
    }
}

int parse_mesh(char *filename, std::vector<double> &positions, std::vector<int> &idx) {
    std::ifstream in(filename);
    if(!in) {
        std::cout << "Cannot open input file.\n";
        return 0;
    }

    boost::char_separator<char> space_sep(" ");

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
            int pos_count = 0;
            for (auto tok = ++tokens.begin(); tok != tokens.end(); ++tok) {
                // doesn't allow more than three values per vertex
                if (pos_count == 3) {
                    continue;
                }

                positions.push_back(std::stod(*tok));
                pos_count++;
            }
        }
    }

    in.close();
    return 1;
}

int main(int argc, char* argv[])
{
    // grid parameters
    TV min_corner = TV::Zero();
    TV max_corner = TV::Zero();

    for (int i = 0; i < dim; i++) {
        min_corner(i) = std::numeric_limits<float>::max();
        max_corner(i) = std::numeric_limits<float>::min();
    }

    std::vector<MeshObject*> meshes;


    if (argc > 1) {
        printf("importing meshes...\n");

        for (int i = 1; i < argc; i++) {
            printf("importing mesh with the file path: %s\n", argv[i]);

            std::vector<double> positions;
            std::vector<int> idx;
            parse_mesh(argv[i], positions, idx);

            int num_vertices = positions.size() / 3;
            int num_tri = idx.size() / 3;
            MeshObject *m = construct_mesh_object(num_vertices, positions.data(),
                                                  num_tri, idx.data());
            meshes.push_back(m);

            compute_bounds(min_corner, max_corner, positions);
        }

        printf("finished importing meshes.\n");
    }

    std::cout << "min_corner" << std::endl;
    std::cout << min_corner << std::endl;
    std::cout << "max_corner" << std::endl;
    std::cout << max_corner << std::endl;

    T dx = 0.05;

    T E = 1e4;
    T nu = 0.3;

    SimulationDriver<T,dim> driver(min_corner, max_corner, dx, E, nu);

    // simulate
    //driver.run(240);
    return 0;
}
