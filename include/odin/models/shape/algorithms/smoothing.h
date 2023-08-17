#ifndef ODIN_SMOOTHING_H
#define ODIN_SMOOTHING_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

namespace odin::algorithms {

using namespace Eigen;

template<typename T>
using sparse_matrix_type = SparseMatrix<T>;
typedef Triplet<double> T;

template<typename T>
sparse_matrix_type<T> build_variable_indices(const MatrixX<T> &band) {
    int      num_variables = band.count();
    MatrixXi variable_indices
            = MatrixXi::Constant(band.rows(), band.cols(), -1);
    int count = 0;
    for (int i = 0; i < band.rows(); ++i) {
        for (int j = 0; j < band.cols(); ++j) {
            if (band(i, j) != 0) {
                variable_indices(i, j) = count;
                count++;
            }
        }
    }
    return variable_indices;
}

template<typename T>
sparse_matrix_type<T> buildq2d(const MatrixXi &variable_indices) {
    int                   num_variables = variable_indices.maxCoeff() + 1;
    sparse_matrix_type<T> filterq(3 * num_variables, num_variables);

    MatrixXi padded_variable_indices = MatrixXi::Constant(
            variable_indices.rows() + 1, variable_indices.cols() + 1, -1);
    padded_variable_indices.block(
            0, 0, variable_indices.rows(), variable_indices.cols())
            = variable_indices;

    std::vector<T> triplet_list;
    triplet_list.reserve(3 * num_variables);

    for (int i = 0; i < variable_indices.rows(); ++i) {
        for (int j = 0; j < variable_indices.cols(); ++j) {
            int count = variable_indices(i, j);
            if (count >= 0) {
                triplet_list.emplace_back(2 * count, count, -2);
                int neighbor = padded_variable_indices(i, j + 1);
                if (neighbor >= 0) {
                    triplet_list.emplace_back(2 * count, neighbor, 1);
                } else {
                    triplet_list.emplace_back(2 * count, count, 1);
                }

                neighbor = padded_variable_indices(i + 1, j);
                if (neighbor >= 0) {
                    triplet_list.emplace_back(2 * count, neighbor, 1);
                } else {
                    triplet_list.emplace_back(2 * count, count, 1);
                }

                triplet_list.emplace_back(2 * count + 1, count, -2);
                neighbor = padded_variable_indices(i, j - 1);
                if (neighbor >= 0) {
                    triplet_list.emplace_back(2 * count + 1, neighbor, 1);
                } else {
                    triplet_list.emplace_back(2 * count + 1, count, 1);
                }

                neighbor = padded_variable_indices(i - 1, j);
                if (neighbor >= 0) {
                    triplet_list.emplace_back(2 * count + 1, neighbor, 1);
                } else {
                    triplet_list.emplace_back(2 * count + 1, count, 1);
                }
            }
        }
    }

    filterq.setFromTriplets(triplet_list.begin(), triplet_list.end());
    return filterq.transpose() * filterq;
}

template<typename T>
VectorXd jacobi(const sparse_matrix_type<T> &filterq,
        const VectorXd                      &x0,
        const VectorXd                      &lower_bound,
        const VectorXd                      &upper_bound,
        int                                  max_iters = 10,
        double                               rel_tol   = 1e-6,
        double                               weight    = 0.5) {
    sparse_matrix_type<T> jacobi_r = filterq;
    jacobi_r.diagonal().setZero();
    VectorXd jacobi_d = filterq.diagonal().cwiseInverse();
    VectorXd x        = x0;

    double cum_rel_tol = 1 - std::pow(1 - rel_tol, 10);

    double energy_now = 0.5 * x.dot(filterq * x);
    std::cout << "Energy at iter 0: " << energy_now << std::endl;

    for (int i = 0; i < max_iters; ++i) {
        VectorXd x_1 = -jacobi_d.cwiseProduct(jacobi_r * x);
        x            = weight * x_1 + (1 - weight) * x;

        x = x.cwiseMax(lower_bound);
        x = x.cwiseMin(upper_bound);

        if ((i + 1) % 10 == 0) {
            double energy_before = energy_now;
            energy_now           = 0.5 * x.dot(filterq * x);
            std::cout << "Energy at iter " << i + 1 << ": " << energy_now
                      << std::endl;

            double cum_rel_improvement
                    = (energy_before - energy_now) / energy_before;
            if (cum_rel_improvement < cum_rel_tol) { break; }
        }
    }

    return x;
}

//SparseMatrix<float> buildq3d(const MatrixXf& variable_indices) {
//    int                 num_variables = variable_indices.maxCoeff() + 1;
//    SparseMatrix<float> filterq(3 * num_variables, num_variables);
//    MatrixXf            padded_variable_indices = MatrixXf::Constant(
//            variable_indices.rows() + 2, variable_indices.cols() + 2, -1);
//    padded_variable_indices.block(
//            1, 1, variable_indices.rows(), variable_indices.cols())
//            = variable_indices;
//
//    std::vector<Triplet<float>> triplets;
//    triplets.reserve(7 * num_variables);
//
//    for (int i = 1; i < padded_variable_indices.rows() - 1; ++i) {
//        for (int j = 1; j < padded_variable_indices.cols() - 1; ++j) {
//            int count = padded_variable_indices(i, j);
//            if (count >= 0) {
//                triplets.emplace_back(3 * count, count, -2);
//                float neighbor;
//
//                neighbor = padded_variable_indices(i - 1, j);
//                if (neighbor >= 0) {
//                    triplets.emplace_back(3 * count, neighbor, 1);
//                } else {
//                    triplets.emplace_back(3 * count, count, 1);
//                }
//
//                neighbor = padded_variable_indices(i + 1, j);
//                if (neighbor >= 0) {
//                    triplets.emplace_back(3 * count, neighbor, 1);
//                } else {
//                    triplets.emplace_back(3 * count, count, 1);
//                }
//
//                // ... similar blocks for j-1, j+1, j, k-1 and j, k+1
//            }
//        }
//    }
//
//    filterq.setFromTriplets(triplets.begin(), triplets.end());
//    filterq.makeCompressed();
//
//    return filterq.transpose() * filterq;
//}
//
//VectorXf jacobi(SparseMatrix<float>& filterq,
//        VectorXf&                    x0,
//        VectorXf&                    lower_bound,
//        VectorXf&                    upper_bound,
//        int                          max_iters = 10,
//        float                        rel_tol   = 1e-6,
//        float                        weight    = 0.5) {
//
//    SparseMatrix<float> jacobi_r = filterq;
//    VectorXf            jacobi_d = filterq.diagonal().cwiseInverse();
//    for (int k = 0; k < jacobi_r.outerSize(); ++k)
//        for (SparseMatrix<float>::InnerIterator it(jacobi_r, k); it; ++it)
//            if (it.row() == it.col()) it.valueRef() = 0;
//
//    VectorXf x          = x0;
//    float    energy_now = 0.5 * x.transpose() * filterq * x;
//
//    for (int i = 0; i < max_iters; ++i) {
//        VectorXf x_1 = -jacobi_d.cwiseProduct(jacobi_r * x);
//        x            = weight * x_1 + (1 - weight) * x;
//
//        // Apply constraints
//        for (int j = 0; j < x.size(); ++j) {
//            x[j] = std::max(x[j], lower_bound[j]);
//            x[j] = std::min(x[j], upper_bound[j]);
//        }
//
//        // Stopping criterion
//        if ((i + 1) % 10 == 0) {
//            float energy_before = energy_now;
//            energy_now          = 0.5 * x.transpose() * filterq * x;
//            float cum_rel_improvement
//                    = (energy_before - energy_now) / energy_before;
//
//            if (cum_rel_improvement < rel_tol) break;
//        }
//    }
//
//    return x;
//}

}// namespace odin::algorithms
#endif// ODIN_SMOOTHING_H
