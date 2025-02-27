
#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <algorithm>

#include <mkl.h>
#include <mkl_pardiso.h>

#include "FEMSolver.h"


// ======== FEMSolver Methods ========

Eigen::MatrixXd FEMSolver::unit_stiffness(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &p3)
{
    // Compute element area
    double A = 0.5 * abs(p1.x() * (p2.y() - p3.y())
        + p2.x() * (p3.y() - p1.y())
        + p3.x() * (p1.y() - p2.y())
    );

    // Shape function derivatives
    double b1 = p2.y() - p3.y(), b2 = p3.y() - p1.y(), b3 = p1.y() - p2.y();
    double c1 = p3.x() - p2.x(), c2 = p1.x() - p3.x(), c3 = p2.x() - p1.x();

    // Strain-displacement matrix B (3x6)
    Eigen::MatrixXd B(3, 6);
    B << b1, 0, b2, 0, b3, 0,
        0, c1, 0, c2, 0, c3,
        c1, b1, c2, b2, c3, b3;
    B *= (1.0 / (2.0 * A));

    // Plane stress constitutive matrix D (3x3)
    Eigen::MatrixXd D(3, 3);
    D << 1, nu, 0,
        nu, 1, 0,
        0, 0, (1 - nu) / 2;
    D *= (E / (1 - nu * nu));

    // Compute element stiffness matrix Ke = t * A * B^T * D * B
    Eigen::MatrixXd Ke = thick * A * (B.transpose() * D * B);

    return Ke;
}


void FEMSolver::assemble()
{
    auto nodes = mesh.nodes;
    auto elements = mesh.elements;

    int numNodes = nodes.size();
    F.setZero(numNodes * 2);
    U.setZero(numNodes * 2);

    // Parallelize element matrix computation with OpenMP
    #pragma omp parallel for num_threads(num_threads)
    for (int e = 0; e < elements.size(); ++e) {
		auto elem = elements[e];

        Eigen::MatrixXd Ke = unit_stiffness( nodes[elem(0)], nodes[elem(1)], nodes[elem(2)] );

        // Print the matrix
        std::cout << "Ke:\n" << Ke << std::endl;

        #pragma omp critical
        for (int i = 0; i < 3; ++i)
		{
            int row = elem(i);
			for (int j = 0; j < 3; ++j)
			{
				int col = elem(j);
				K.insert(2 * row, 2 * col, Ke(2 * i, 2 * j));
				K.insert(2 * row + 1, 2 * col, Ke(2 * i + 1, 2 * j));
				K.insert(2 * row, 2 * col + 1, Ke(2 * i, 2 * j + 1));
				K.insert(2 * row + 1, 2 * col + 1, Ke(2 * i + 1, 2 * j + 1));
			}
		}
	}

    // Apply force
    for (int i = 0; i < numNodes; ++i) {

        std::cout << nodes[i].y() << "\n";

        if (abs(nodes[i].y()) < 1E-2) {
            F(2 * i + 1) = f;
        }
    }

	// Apply boundary conditions
	for (int i = 0; i < numNodes; ++i) {

		if (abs(nodes[i].y() - Ly) < 1E-2) {

            // Apply BC for x-direction (row and column)
            K.modifyRow(i);

			F(2 * i) = 0;
			F(2 * i + 1) = 0;
		}
	}

	K.print_rows();

    // Finalize CSR structure
    K.finalize();

    K.print();

    std::cout << F << std::endl;

    // MPI synchronization for shared nodes
    MPI_Allreduce(MPI_IN_PLACE, F.data(), numNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void FEMSolver::solve()
{
    void* pt[64] = { nullptr };
    MKL_INT maxfct = 1, mnum = 1;
    MKL_INT mtype = 11; // Real and nonsymmetric matrix

    MKL_INT phase;
    MKL_INT nsize = K.get_rows();

    MKL_INT nrhs = 1, error = 0;

    // Solver-specific parameters that control various aspects of the Pardiso solver
    MKL_INT iparm[64] = { 0 };
    iparm[0] = 0;   // Solver default parameters overridden with provided by iparm
    /*
    iparm[1] = 0;   // You can control the parallel execution of the solver, The minimum degree algorithm
    iparm[5] = 0;  // Write solution into x
    iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS
    iparm[34] = 1;  // Zero-based indexing: columns and rows indexing in arrays ia, ja, and perm starts from 0 (C-style indexing).
    */
    std::vector<MKL_INT> ia_long(K.row_ptr.size());
    std::vector<MKL_INT> ja_long(K.col_ind.size());

    // Use std::transform to cast from int to long long
    std::transform(K.row_ptr.begin(), K.row_ptr.end(), ia_long.begin(),
        [](int val) { return static_cast<long long>(val); });

    std::transform(K.col_ind.begin(), K.col_ind.end(), ja_long.begin(),
        [](int val) { return static_cast<long long>(val); });

    phase = 13; // Analysis, numerical factorization, solve
    std::vector<double> tt(K.val);

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nsize, &tt[0], &ia_long[0], &ja_long[0],
        nullptr, &nrhs, iparm, nullptr, F.data(), U.data(), &error);

    if (error != 0) {
        std::cerr << "PARDISO error: " << error << std::endl;
        MPI_Abort(MPI_COMM_WORLD, error);
    }

    phase = -1; // Release all internal memory for all matrices
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nsize, nullptr, nullptr, nullptr, nullptr, &nrhs, iparm, nullptr, nullptr, nullptr, &error);

    if (error != 0) {
        std::cerr << "PARDISO error: " << error << std::endl;
        MPI_Abort(MPI_COMM_WORLD, error);
    }
}
