
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <petscksp.h>
#include "FEMSolver.h"


#define DEBUG_PRINT(rank, msg) std::cout << "Rank " << rank << ": " << msg << std::endl;


void test() {
    int N = 4;
    CSRMatrix A(N, N);

    int nnz = 12;
    A.reserve(nnz);

    A.insert(0, 0, 10);
    A.insert(1, 0, 3);
    A.insert(3, 0, 3);
    A.insert(1, 1, 9);
    A.insert(2, 1, 7);
    A.insert(2, 2, 8);
    A.insert(2, 3, 7);
    A.insert(3, 2, 8);

    // Finalize CSR structure
    A.finalize();

    // Printing matrix in CSR format
    A.print();
}


int main(int argc, char *argv[]) {

    // Initialize MPI and PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set the number of threads to 4 programmatically
    omp_set_num_threads(4);

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

    DEBUG_PRINT(rank, "Entering main function");
 

    // Set parameters for FEM problem (grid size, domain size, force)
    int Nx = 10;  // number of elements in x-direction
    int Ny = 10;  // number of elements in y-direction

    double Lx = 1.0;  // length in x-direction
    double Ly = 1.0;  // length in y-direction

    double f = 100.0;  // force applied

    // Instantiate and run the FEM solver

    FEMSolver femSolver(Nx, Ny, Lx, Ly, f, rank, size);

    for (int k = 0; k< 1; k++)
        femSolver.run();

    femSolver.printSolution();  // Call the new function
    DEBUG_PRINT(rank, "Exiting main function");

    // Finalize PETSc and MPI
    PetscFinalize();

    return 0;
}
