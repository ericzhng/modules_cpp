
#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <omp.h>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


class CSRMatrix {
public:
    std::vector<long long> ia, ja;
    std::vector<double> a;
    long long rows, cols;

    CSRMatrix(int n);
    void reserve(int nnz);
    void addValue(int row, int col, double value);
    void finalize();
};

class Mesh {
public:
    std::vector<Eigen::Vector2d> nodes;
    std::vector<Eigen::Vector3i> elements;
    int startRow, endRow;

    Mesh(int nx, int ny, double lx, double ly, int rank, int size);
};


class FEMSolver {
private:
    int rank, size, num_threads;
    int Nx, Ny;
    double Lx, Ly, f;
    Mesh mesh;
    CSRMatrix K;
    Eigen::VectorXd F, U;

    Eigen::Matrix3d computeElementMatrix(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &p3);
    void assembleSystem();
    void solveWithPardiso();

public:
    FEMSolver(int Nx, int Ny, double Lx, double Ly, double f, int rank, int size);
    void run();
};

#endif // FEM_SOLVER_H
