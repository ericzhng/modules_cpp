
#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <algorithm>
#include <petscksp.h>

#include "FEMSolver.h"

// Ensure PETSc works with C++
PetscErrorCode ierr;


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

    // Compute element stiffness matrix
    // Eigen::MatrixXd Ke = thick * A * (B.transpose() * D * B);
    Eigen::MatrixXd Ke = A * (B.transpose() * D * B);

    return Ke;
}


void FEMSolver::assemble()
{
    auto nodes = mesh.nodes;
    auto elements = mesh.elements;

    int numNodes = nodes.size();
    F.setZero(numNodes * 2);
    U.setZero(numNodes * 2);

    #pragma omp parallel
    {
        CSRMatrix K_local(numNodes * 2);
        K_local.reserve(elements.size() * 18);

        #pragma omp for schedule(dynamic)
        for (int e = 0; e < elements.size(); ++e) {
            auto elem = elements[e];
            Eigen::MatrixXd Ke = unit_stiffness(nodes[elem(0)], nodes[elem(1)], nodes[elem(2)]);
            for (int i = 0; i < 3; ++i) {
                int row = elem(i);
                for (int j = 0; j < 3; ++j) {
                    int col = elem(j);
                    K_local.insert(2 * row, 2 * col, Ke(2 * i, 2 * j));
                    K_local.insert(2 * row + 1, 2 * col, Ke(2 * i + 1, 2 * j));
                    K_local.insert(2 * row, 2 * col + 1, Ke(2 * i, 2 * j + 1));
                    K_local.insert(2 * row + 1, 2 * col + 1, Ke(2 * i + 1, 2 * j + 1));
                }
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < numNodes * 2; ++i) {
                for (const auto& entry : K_local.row_vec[i]) {
                    K.insert(i, entry.first, entry.second);
                }
            }
        }
    }

    // Apply force
    for (int i = 0; i < numNodes; ++i) {

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

    // Communicate ghost node data for K and F
    std::vector<MPI_Request> requests;
    for (const auto& [neighbor, ghost_indices] : mesh.ghost_map) {
        std::vector<double> send_vals, send_F;
        std::vector<int> send_cols;
        for (int i : ghost_indices) {
            int row_start = K.row_ptr[2 * i];
            int row_end = K.row_ptr[2 * i + 1];
            for (int k = row_start; k < row_end; ++k) {
                send_cols.push_back(K.col_ind[k]);
                send_vals.push_back(K.val[k]);
            }
            row_start = K.row_ptr[2 * i + 1];
            row_end = K.row_ptr[2 * i + 2];
            for (int k = row_start; k < row_end; ++k) {
                send_cols.push_back(K.col_ind[k]);
                send_vals.push_back(K.val[k]);
            }
            send_F.push_back(F(2 * i));
            send_F.push_back(F(2 * i + 1));
        }

        MPI_Request req;
        MPI_Isend(send_cols.data(), send_cols.size(), MPI_INT, neighbor, 0, MPI_COMM_WORLD, &req);
        requests.push_back(req);
        MPI_Isend(send_vals.data(), send_vals.size(), MPI_DOUBLE, neighbor, 1, MPI_COMM_WORLD, &req);
        requests.push_back(req);
        MPI_Isend(send_F.data(), send_F.size(), MPI_DOUBLE, neighbor, 2, MPI_COMM_WORLD, &req);
        requests.push_back(req);
    }

    for (const auto& [neighbor, ghost_indices] : mesh.ghost_map) {
        std::vector<int> recv_cols(ghost_indices.size() * 20); // Approximate
        std::vector<double> recv_vals(ghost_indices.size() * 20), recv_F(ghost_indices.size() * 2);
        MPI_Status status;
        int count;
        MPI_Recv(recv_cols.data(), recv_cols.size(), MPI_INT, neighbor, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &count);
        recv_cols.resize(count);
        MPI_Recv(recv_vals.data(), recv_vals.size(), MPI_DOUBLE, neighbor, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &count);
        recv_vals.resize(count);
        MPI_Recv(recv_F.data(), recv_F.size(), MPI_DOUBLE, neighbor, 2, MPI_COMM_WORLD, &status);

        int idx = 0, f_idx = 0;
        for (int i : ghost_indices) {
            for (int k = 0; k < 20 && idx < recv_cols.size(); ++k) { // Approximate entries per row
                K.insert(2 * i, recv_cols[idx], recv_vals[idx]);
                K.insert(2 * i + 1, recv_cols[idx], recv_vals[idx]);
                idx++;
            }
            F(2 * i) += recv_F[f_idx++];
            F(2 * i + 1) += recv_F[f_idx++];
        }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // Finalize CSR structure
    K.finalize();

    // MPI synchronization for shared nodes
    // MPI_Allreduce(MPI_IN_PLACE, F.data(), numNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


PetscErrorCode FEMSolver::solve()
{
    // PETSc variables
    Mat A;           // Stiffness matrix
    Vec b, u;        // RHS and solution vectors
    KSP ksp;         // Linear solver context
    PC pc;           // Preconditioner context
    PetscInt n = K.get_rows(); // Size of the system

    // Initialize PETSc matrix and vectors
    ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
    ierr = MatSetType(A, MATMPIAIJ); CHKERRQ(ierr); // Sparse parallel matrix
    ierr = MatMPIAIJSetPreallocation(A, 20, NULL, 20, NULL); CHKERRQ(ierr); // Preallocate 20 nonzeros per row

    // Fill the PETSc matrix from CSRMatrix
    for (PetscInt i = 0; i < n; i++) {
        PetscInt row_start = K.row_ptr[i];
        PetscInt row_end = K.row_ptr[i + 1];
        for (PetscInt k = row_start; k < row_end; k++) {
            PetscInt col = K.col_ind[k];
            PetscScalar val = K.val[k];
            ierr = MatSetValue(A, i, col, val, INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    // Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetType(b, VECMPI); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &u); CHKERRQ(ierr);
    ierr = VecSetSizes(u, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetType(u, VECMPI); CHKERRQ(ierr);

    // Fill the RHS vector
    for (PetscInt i = 0; i < n; i++) {
        ierr = VecSetValue(b, i, F(i), INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

    // Create and configure the linear solver
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr); // Example preconditioner
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);  // Conjugate Gradient solver
    ierr = KSPSetTolerances(ksp, 1e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr); // Allow runtime options

    // Solve the system
    ierr = KSPSolve(ksp, b, u); CHKERRQ(ierr);

    // Extract the solution into Eigen vector U
    const PetscScalar *u_array;
    ierr = VecGetArrayRead(u, &u_array); CHKERRQ(ierr);
    for (PetscInt i = 0; i < n; i++) {
        U(i) = u_array[i];
    }
    ierr = VecRestoreArrayRead(u, &u_array); CHKERRQ(ierr);

    // Synchronize ghost node displacements
    std::vector<MPI_Request> requests;
    for (const auto& [neighbor, ghost_indices] : mesh.ghost_map) {
        std::vector<double> send_U;
        for (int i : ghost_indices) {
            send_U.push_back(U(2 * i));
            send_U.push_back(U(2 * i + 1));
        }
        MPI_Request req;
        MPI_Isend(send_U.data(), send_U.size(), MPI_DOUBLE, neighbor, 3, MPI_COMM_WORLD, &req);
        requests.push_back(req);
    }

    for (const auto& [neighbor, ghost_indices] : mesh.ghost_map) {
        std::vector<double> recv_U(ghost_indices.size() * 2);
        MPI_Recv(recv_U.data(), recv_U.size(), MPI_DOUBLE, neighbor, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int idx = 0;
        for (int i : ghost_indices) {
            U(2 * i) = recv_U[idx++];
            U(2 * i + 1) = recv_U[idx++];
        }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    // Clean up PETSc objects
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);

    return 0; // Success
}


void FEMSolver::printSolution() {
    int numNodes = mesh.nodes.size();
    std::cout << "Rank " << rank << " Solution (Node, X, Y, Ux, Uy):\n";
    for (int i = 0; i < numNodes; ++i) {
        double x = mesh.nodes[i].x();
        double y = mesh.nodes[i].y();
        double ux = U(2 * i);     // x-displacement
        double uy = U(2 * i + 1); // y-displacement
        std::cout << "Node " << i << ": (" << x << ", " << y << ") -> ("
                  << ux << ", " << uy << ")\n";
    }
    std::cout << std::endl;
}