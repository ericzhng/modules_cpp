
#include <iostream>
#include <vector>

#include "pardiso.h"


/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso PARDISO
#else
#define PARDISO pardiso
#endif

#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif


int main(void) {
    // Define the sparse matrix in CSR format
    MKL_INT n = 4;
    std::vector<MKL_INT> ia = {1, 5, 8, 10, 13};
    std::vector<MKL_INT> ja = {1, 2, 3, 4, 1, 3, 4, 2, 3, 1, 2, 4};
    std::vector<double> a = {10, -1, 2, 0, -1, 9, 3, -1, 7, 3, -2, 6};

    // Define the right-hand side vector
    std::vector<double> b = {6, 25, -11, 15};

    // Define the solution vector
    std::vector<double> x(n, 0.0);

    // PARDISO control parameters
    MKL_INT iparm[64];
    for (int i = 0; i < 64; ++i) {
        iparm[i] = 0;
    }
    iparm[0] = 1;   // No default values for solver
    iparm[1] = 2;   // Fill-in reducing ordering (MKL_DSS_DEFAULTS)
    iparm[2] = 3;   // Numbers of processors, value of MKL_DSS_DEFAULTS is environment dependent.
    iparm[3] = 0;   // No iterative-refinement step
    iparm[4] = 0;   // No check for singularity
    iparm[5] = 0;   // Zero-based indexing for matrix (MKL_DSS_ZERO_BASED_INDEXING)
    iparm[6] = 0;   // Not in use
    iparm[7] = 2;   // Maximum number of iterative refinement steps
    iparm[9] = 13;  // Perturb the pivot elements with 1E-13
    iparm[10] = 1;  // Use nonsymmetric permutation
    iparm[12] = 0;  // Output: number of perturbed pivots
    iparm[34] = 1;  // Activate the output of the amount of memory consumed

    MKL_INT maxfct = 1; // Maximum number of factors
    MKL_INT mnum = 1;   // Which factorization to use
    MKL_INT mtype = 11;  // Real and nonsymmetric matrix
    void* pt[64];       // Pointer to the solver data structure
    MKL_INT nrhs = 1;   // Number of right-hand sides
    MKL_INT error = 0;  // Error indicator
    
    // Initialize PARDISO
    pardisoinit(pt, &mtype, iparm);

    // Phase 1: Analysis
    MKL_INT phase = 11;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], nullptr, &nrhs, iparm, &error, nullptr, nullptr);
    if (error != 0) {
        std::cerr << "Error during analysis phase: " << error << std::endl;
        return 1;
    }

    // Phase 2: Factorization
    phase = 22;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], nullptr, &nrhs, iparm, &error, nullptr, nullptr);
    if (error != 0) {
        std::cerr << "Error during factorization phase: " << error << std::endl;
        return 1;
    }

    // Phase 3: Solve
    phase = 33;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], nullptr, &nrhs, iparm, &error, &b[0], &x[0]);
    if (error != 0) {
        std::cerr << "Error during solve phase: " << error << std::endl;
        return 1;
    }

    // Print the solution
    std::cout << "Solution:" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    // Phase -1: Release internal memory
    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], nullptr, &nrhs, iparm, &error, nullptr, nullptr);
    if (error != 0) {
        std::cerr << "Error during memory release phase: " << error << std::endl;
        return 1;
    }

    return 0;

}
