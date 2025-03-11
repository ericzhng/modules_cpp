
#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <omp.h>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


class CSRMatrix
{
public:
    std::vector<double> val{};       // Nonzero values
    std::vector<int> col_ind{};      // Column indices
    std::vector<int> row_ptr{};      // Row pointer

private:
    int rows = 0, cols = 0;

    // Temporary storage for unordered insertions (row-wise)
    std::vector<std::vector<std::pair<int, double>>> row_vec{};

public:
    CSRMatrix(int n) : rows(n), cols(n) {
        row_ptr.resize(rows + 1, 0);
        row_vec.resize(rows);
    }

    CSRMatrix(int rows, int cols) : rows(rows), cols(cols) {
        row_ptr.resize(rows + 1, 0);
        row_vec.resize(rows);
    }

    int get_rows() {
        return rows;
    }

    std::vector<int> get_row_ptr() {
        return row_ptr;
    }

    std::vector<int> get_col_ind() {
        return col_ind;
    }

    std::vector<double> get_val() {
        return val;
    }

    void modifyRow(int elemId)
    {
        // Set all values in the row to 0 except for the diagonal
        row_vec[elemId * 2].clear();
        row_vec[elemId * 2].emplace_back(elemId * 2, 1e10);

        row_vec[elemId * 2 + 1].clear();
        row_vec[elemId * 2 + 1].emplace_back(elemId * 2 + 1, 1e10);
    }


    // Insert a nonzero element (typically used in matrix assembly)
    void insert(int row, int col, double value) {

        // Iterate through existing values in the row to check for the column index
        bool found = false;
        for (auto& entry : row_vec[row]) {
            if (entry.first == col) {
                // If the column exists, add the value to the existing entry
                entry.second += value;
                found = true;
                break;
            }
        }

        // If no entry for this column, add a new entry
        if (!found) {
            row_vec[row].emplace_back(col, value);
        }
    }

    void CSRMatrix::reserve(int nnz) {
        col_ind.reserve(nnz);
        val.reserve(nnz);
    }

    // Finalize the CSR matrix: Sort columns and compute row_ptr offsets
    void finalize() {
        // Step 1: Compute `row_ptr` offsets
        for (int i = 0; i < rows; i++) {
            row_ptr[i + 1] = row_ptr[i] + row_vec[i].size();
        }

        // Reserve space for `val` and `col_ind`
        val.resize(row_ptr[rows]);
        col_ind.resize(row_ptr[rows]);

        // Step 2: Sort and copy values into CSR format
        for (int i = 0; i < rows; i++) {
            // Sort row data by column index
            std::sort(row_vec[i].begin(), row_vec[i].end());

            // Fill `val` and `col_ind`
            int index = row_ptr[i];
            for (const auto& entry : row_vec[i]) {
                col_ind[index] = entry.first;
                val[index] = entry.second;
                index++;
            }
        }

        // Clear temporary storage
        row_vec.clear();
        row_vec.reserve(0);
    }

    // Print CSR format
    void print() {
        std::cout << "val: ";
        for (double v : val) std::cout << v << " ";
        std::cout << "\ncol_ind: ";
        for (int c : col_ind) std::cout << c << " ";
        std::cout << "\nrow_ptr: ";
        for (int r : row_ptr) std::cout << r << " ";
        std::cout << "\n";
    }

    void print_rows() {
        for (size_t i = 0; i < row_vec.size(); ++i) {
            std::cout << "Row " << i << ": ";
            for (const auto& elem : row_vec[i]) {
                std::cout << "(" << elem.first << ", " << elem.second << ") ";
            }
            std::cout << std::endl;
        }
    }
};


class Mesh {
public:
    std::vector<Eigen::Vector2d> nodes{};
    std::vector<Eigen::Vector3i> elements{};

    Mesh(int nx, int ny, double lx, double ly, int rank, int size)
    {
        int numNodesX = nx + 1, numNodesY = ny + 1;

        int totalRows = numNodesY;
        int rowsPerProc = totalRows / size;

        int startRow = rank * rowsPerProc;
        int endRow = (rank == size - 1) ? totalRows - 1 : (rank + 1) * rowsPerProc - 1;

        for (int j = startRow; j <= endRow; ++j)
            for (int i = 0; i < numNodesX; ++i)
                nodes.emplace_back(i * lx / nx, j * ly / ny);

        // triangular mesh
        for (int j = startRow; j < endRow; ++j) {
            for (int i = 0; i < nx; ++i) {
                int n1 = j * numNodesX + i, n2 = n1 + 1;
                int n3 = n1 + numNodesX, n4 = n3 + 1;
                elements.emplace_back(n1, n2, n3);
                elements.emplace_back(n2, n4, n3);
            }
        }
    }
};


class FEMSolver {
public:
    double thick = 10E-3;
    double E = 70E9;
    double nu = 0.33;

private:
    int num_threads = 0;
    int Nx = 0, Ny = 0;
    double Lx = 0.0, Ly = 0.0, f = 0.0;

    Mesh mesh;
    CSRMatrix K;
    Eigen::VectorXd F{}, U{};

    Eigen::MatrixXd unit_stiffness(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &p3);

    void assemble();

    void solve();

public:
    FEMSolver(int Nx, int Ny, double Lx, double Ly, double f, int rank, int size)
        : Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly), f(f), mesh(Nx, Ny, Lx, Ly, rank, size), 
        K(mesh.nodes.size() * 2), F(mesh.nodes.size() * 2), U(mesh.nodes.size() * 2) {

        #pragma omp parallel
        {
            #pragma omp master
            {
                num_threads = omp_get_num_threads();
            }
        }
        int nn = mesh.nodes.size();
        K.reserve(nn*20);
    }

    void run() {
        assemble();
        solve();
        // std::cout << "Solution U:\n" << U << std::endl;
    }

};

#endif // FEM_SOLVER_H
