
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include "FEMSolver.h"

int main(int argc, char *argv[]) {

    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set OpenMP threads based on available cores
    int num_threads = 0;

    #pragma omp parallel
    {
        #pragma omp master
        {
            num_threads = omp_get_num_threads();
        }
    }

    std::cout << "MPI Process Rank: " << rank << " out of " << size << " processes." << std::endl;
    if (rank == 0) {
        std::cout << "Running with " << num_threads << " OpenMP threads per process." << std::endl;
    }

    // Set parameters for FEM problem (grid size, domain size, force)
    int Nx = 10;  // number of elements in x-direction
    int Ny = 10;  // number of elements in y-direction

    double Lx = 1.0;  // length in x-direction
    double Ly = 1.0;  // length in y-direction

    double f = 1.0;  // force applied

    // Instantiate and run the FEM solver
    FEMSolver femSolver(Nx, Ny, Lx, Ly, f, rank, size);
    femSolver.run();

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
