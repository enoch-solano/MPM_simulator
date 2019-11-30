#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MassSpringSystem.h"

template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MPMSystem<T,dim> mpm;
    T dt;
    TV gravity;

    TV sphere_center;
    T sphere_radius;
    T ground;

    SimulationDriver(TV min, TV max, T dx)
      : mpm(min, max, dx), dt((T)1e-3)   // TODO: choose a better/worse dt?
    {
        gravity.setZero();
        gravity(1) = -9.8 * 0.15;

        sphere_center = TV::Ones()*0.5;
        sphere_radius = 0.2;
        ground = 0.1;

        TV min_box = TV::Ones()*0.3;
        TV max_box = TV::Ones()*0.7;

        printf("populating particles...\n");
        mpm.populate_cube(min_box, max_box);
        printf("finished populating particles.\n");
    }

    void run(const int max_frame)
    {
        for(int frame = 1; frame < max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int) (((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                advanceMPMSim();
            }
            mkdir("output/", 0777);
            std::string filename = "output/" + std::to_string(frame) + ".poly";
            mpm.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

    void advanceMPMSim() {
        // init grid data to zero
        mpm.grid.clear();

        // for (int i = 0; i < mpm.grid.grid.size(); i++) {
        //     if (mpm.grid.f(i)[0] > 0 || mpm.grid.f(i)[1] > 0 || mpm.grid.f(i)[2] > 0) {
        //         printf("f: (%.3f, %.3f, %.3f)\n", mpm.grid.f(i)[0], mpm.grid.f(i)[1], mpm.grid.f(i)[2]);
        //     }
        // }

        // TODO: (sanity checks)
        // compute total particle momentum before p2g
        // compute total grid momentum after p2g

        mpm.compute_weights_3D();
        // mpm.print_particles();



        // transfer particle momentum to grid
        mpm.transfer_P2G();
        // mpm.grid.print_active_nodes();


        // compute force on grid
        mpm.addGravity(gravity);
        mpm.addElasticity();

        // update velocity on grid
        mpm.updateGridVelocity(dt);

        // boundary conditions
        mpm.setBoundaryVelocities(3);

        // TODO: (sanity checks)
        // compute total grid momentum after g2p
        // compute total particle momentum before g2p

        // evolve force

        // transfer grid moment to particles
        T flip = 0.95;
        mpm.transfer_G2P(flip, dt);
    }
};
