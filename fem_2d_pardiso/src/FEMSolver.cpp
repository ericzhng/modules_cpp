
#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <algorithm>

#include <mkl.h>
#include <mkl_pardiso.h>

#include "FEMSolver.h"


// ======== CSRMatrix Methods ========
CSRMatrix::CSRMatrix(int n) : rows(n), cols(n) {
    ia.resize(n + 1, 0);
}

void CSRMatrix::reserve(int nnz) {
    ja.reserve(nnz);
    a.reserve(nnz);
}

void CSRMatrix::addValue(int row, int col, double value) {
    ja.push_back(col);
    a.push_back(value);
    ia[row + 1]++;
}

void CSRMatrix::finalize() {
    for (size_t i = 1; i < ia.size(); ++i) {
        ia[i] += ia[i - 1];
    }
}

// ======== Mesh Class ========
Mesh::Mesh(int nx, int ny, double lx, double ly, int rank, int size) {
    int numNodesX = nx + 1, numNodesY = ny + 1;
    int totalRows = numNodesY;
    int rowsPerProc = totalRows / size;
    startRow = rank * rowsPerProc;
    endRow = (rank == size - 1) ? totalRows - 1 : (rank + 1) * rowsPerProc - 1;

    for (int j = startRow; j <= endRow; ++j)
        for (int i = 0; i < numNodesX; ++i)
            nodes.emplace_back(i * lx / nx, j * ly / ny);

    for (int j = startRow; j < endRow; ++j)
        for (int i = 0; i < nx; ++i) {
            int n1 = j * numNodesX + i, n2 = n1 + 1;
            int n3 = n1 + numNodesX, n4 = n3 + 1;
            elements.emplace_back(n1, n2, n3);
            elements.emplace_back(n2, n4, n3);
        }
}

// ======== FEMSolver Methods ========
FEMSolver::FEMSolver(int Nx, int Ny, double Lx, double Ly, double f, int rank, int size)
    : Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly), f(f), rank(rank), size(size),
      mesh(Nx, Ny, Lx, Ly, rank, size), K(mesh.nodes.size()), F(mesh.nodes.size()), U(mesh.nodes.size()) {
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
}

Eigen::Matrix3d FEMSolver::computeElementMatrix(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &p3) {
    Eigen::Matrix<double, 2, 3> B;
    Eigen::Matrix2d D = Eigen::Matrix2d::Identity();
    Eigen::Matrix3d M;
    M << 1, p1.x(), p1.y(),
         1, p2.x(), p2.y(),
         1, p3.x(), p3.y();
    double area = 0.5 * std::abs(M.determinant());

    B << (p2.y() - p3.y()), (p3.y() - p1.y()), (p1.y() - p2.y()),
         (p3.x() - p2.x()), (p1.x() - p3.x()), (p2.x() - p1.x());
    B /= (2.0 * area);

    return (area * (B.transpose() * D * B));
}

void FEMSolver::assembleSystem() {
    int numNodes = mesh.nodes.size();
    F.setZero(numNodes);

    // Parallelize element matrix computation with OpenMP
    #pragma omp parallel for num_threads(num_threads)
    for (int e = 0; e < mesh.elements.size(); ++e) {
		auto elem = mesh.elements[e];
        Eigen::Matrix3d Ke = computeElementMatrix(
            mesh.nodes[elem(0)], mesh.nodes[elem(1)], mesh.nodes[elem(2)]
        );

        #pragma omp critical
		for (int i = 0; i < 3; ++i) {
			int row = mesh.elements[e](i);
			for (int j = 0; j < 3; ++j) {
				int col = mesh.elements[e](j);
				K.addValue(row, col, Ke(i, j));
			}
			F(row) += f / 3.0;
		}
    }

    // Finalize CSR structure
    K.finalize();

    // MPI synchronization for shared nodes
    MPI_Allreduce(MPI_IN_PLACE, F.data(), numNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void FEMSolver::solveWithPardiso() {

    long long n = K.rows, nrhs = 1, mtype = -2, phase, error = 0;

    void* pt[64] = { nullptr };
    long long iparm[64] = { 0 }, maxfct = 1, mnum = 1;

    iparm[0] = 1;
    iparm[1] = 2;
    iparm[10] = 13;
    iparm[39] = num_threads; // Set OpenMP threads

    std::vector<long long> ia_long(K.ia.begin(), K.ia.end());
    std::vector<long long> ja_long(K.ja.begin(), K.ja.end());

    phase = 13;
    pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, K.a.data(), ia_long.data(), ja_long.data(),
        nullptr, &nrhs, iparm, nullptr, F.data(), U.data(), &error);

    if (error != 0) {
        std::cerr << "PARDISO error on rank " << rank << ": " << error << std::endl;
        MPI_Abort(MPI_COMM_WORLD, error);
    }

    phase = -1;
    pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, nullptr, nullptr, nullptr, nullptr, &nrhs, iparm, nullptr, nullptr, nullptr, &error);
}


void FEMSolver::run() {

    assembleSystem();
    solveWithPardiso();

    if (rank == 0) {
        std::cout << "Solution U:\n" << U.transpose() << std::endl;
    }
}
