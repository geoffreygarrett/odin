/** @file
 * @brief Benchmarking tests for matrix inversion policies
 *
 * This file contains benchmarking tests for various matrix inversion
 * policies using the Google Benchmark library. These tests are designed to
 * measure both the execution time and the accuracy of the inversion policies.
 *
 * Each benchmark creates a random square matrix of a predefined size and
 * inverts it using a particular inversion policy. The execution time for the
 * inversion operation is measured and recorded by the benchmark.
 *
 * The accuracy of the inversion is measured by calculating the relative error
 * between the product of the original matrix and its inverse, and the identity
 * matrix. The relative error is computed as the norm of the difference between
 * the product and the identity matrix, divided by the norm of the identity matrix.
 * This value is calculated for each iteration and the average relative
 * error over all iterations is recorded in the benchmark's counters.
 *
 * The predefined size of the matrices used in the benchmarks can be changed by
 * modifying the value of the MATRIX_SIZE macro.
 *
 * The benchmarks are implemented using two macros named RUN_BENCHMARK_DYNAMIC and
 * RUN_BENCHMARK_STATIC, each taking the name of an inversion policy as a parameter.
 * To add a benchmark for a new inversion policy, simply add a call to RUN_BENCHMARK_DYNAMIC
 * and RUN_BENCHMARK_STATIC with the name of the policy at the end of the file. These macros generate
 * benchmarks for dynamically and statically sized matrices respectively. Comparing the results
 * of these two types of benchmarks can provide insights into the performance characteristics
 * of the inversion policies with respect to the size specification method of the matrices.
 *
 * REFERENCE
 * - https://nvlpubs.nist.gov/nistpubs/jres/78B/jresv78Bn2p65_A1b.pdf
 */

#include <Eigen/Dense>
#include <benchmark/benchmark.h>
#include <odin/eigen/inversion_policy.hpp>

/// The size of the matrix for the benchmark
#define MATRIX_SIZE 100
#define SCALAR_TYPE double


using MatrixXd = Eigen::MatrixXd;


/// Macro for running a benchmark test for dynamically sized matrices.
///
/// This macro generates a benchmark test for a matrix inversion policy using dynamically sized matrices.
/// The test creates a random square matrix of a predefined size and inverts it
/// using the given policy. Then it calculates and records the relative error
/// between the result and the expected identity matrix. The performance results
/// from this benchmark can provide insights into how well the policy performs
/// with dynamically sized matrices.
///
/// @param Policy The inversion policy to benchmark.
#define RUN_BENCHMARK_DYNAMIC(Policy)                                          \
    static void BM_##Policy##_dynamic(benchmark::State &state) {               \
        MatrixXd matrix = MatrixXd::Random(MATRIX_SIZE, MATRIX_SIZE);          \
        MatrixXd identity = MatrixXd::Identity(MATRIX_SIZE, MATRIX_SIZE);      \
        SCALAR_TYPE totalRelError = 0;                                              \
        for (auto _: state) {                                                  \
            auto invertedMatrix = policy::Policy::invert(matrix); \
            benchmark::DoNotOptimize(invertedMatrix);                          \
            MatrixXd product = matrix * invertedMatrix;                        \
            SCALAR_TYPE relError = (product - identity).norm() / identity.norm();   \
            totalRelError += relError;                                         \
        }                                                                      \
        state.counters["avg_rel_error"] = totalRelError / state.iterations();  \
    }                                                                          \
    BENCHMARK(BM_##Policy##_dynamic);

/// Macro for running a benchmark test for statically sized matrices.
///
/// This macro generates a benchmark test for a matrix inversion policy using statically sized matrices.
/// The test creates a random square matrix of a predefined size and inverts it
/// using the given policy. Then it calculates and records the relative error
/// between the result and the expected identity matrix. The performance results
/// from this benchmark can provide insights into how well the policy performs
/// with statically sized matrices.
///
/// @param Policy The inversion policy to benchmark.
#define RUN_BENCHMARK_STATIC(Policy)                                                                                            \
    static void BM_##Policy##_static(benchmark::State &state) {                                                                 \
        Eigen::Matrix<SCALAR_TYPE, MATRIX_SIZE, MATRIX_SIZE> matrix = Eigen::Matrix<SCALAR_TYPE, MATRIX_SIZE, MATRIX_SIZE>::Random();     \
        Eigen::Matrix<SCALAR_TYPE, MATRIX_SIZE, MATRIX_SIZE> identity = Eigen::Matrix<SCALAR_TYPE, MATRIX_SIZE, MATRIX_SIZE>::Identity(); \
        SCALAR_TYPE totalRelError = 0;                                                                                               \
        for (auto _: state) {                                                                                                   \
            auto invertedMatrix = policy::Policy::invert(matrix);                                                  \
            benchmark::DoNotOptimize(invertedMatrix);                                                                           \
            Eigen::Matrix<SCALAR_TYPE, MATRIX_SIZE, MATRIX_SIZE> product = matrix * invertedMatrix;                                  \
            SCALAR_TYPE relError = (product - identity).norm() / identity.norm();                                                    \
            totalRelError += relError;                                                                                          \
        }                                                                                                                       \
        state.counters["avg_rel_error"] = totalRelError / state.iterations();                                                   \
    }                                                                                                                           \
    BENCHMARK(BM_##Policy##_static);

/// Macro for running a benchmark test for both static and dynamic matrices.
///
/// This macro generates benchmark tests for a matrix inversion policy.
/// The tests create a random square matrix of a predefined size and inverts it
/// using the given policy. Then it calculates and records the relative error
/// between the result and the expected identity matrix.
///
/// The relative error is calculated as the norm of the difference between
/// the product of the matrix and its inverse and the identity matrix, divided
/// by the norm of the identity matrix. This value is recorded for each
/// iteration and then the average relative error over all iterations is
/// calculated and stored in the benchmark's counters.
///
/// @param Policy The inversion policy to benchmark.
#define RUN_BENCHMARK(Policy)     \
    RUN_BENCHMARK_DYNAMIC(Policy) \
    RUN_BENCHMARK_STATIC(Policy)

RUN_BENCHMARK(direct)
RUN_BENCHMARK(pinv)
RUN_BENCHMARK(direct_regularized)
RUN_BENCHMARK(partial_piv_lu)
RUN_BENCHMARK(full_piv_lu)
RUN_BENCHMARK(householder_qr)
RUN_BENCHMARK(col_piv_householder_qr)
RUN_BENCHMARK(full_piv_householder_qr)
RUN_BENCHMARK(complete_orthogonal_decomposition)
RUN_BENCHMARK(llt)
RUN_BENCHMARK(ldlt)
RUN_BENCHMARK(bdcs_svd)
RUN_BENCHMARK(jacobi_svd)

BENCHMARK_MAIN();
