/**
 * @file matrix_checker.hpp
 * @brief This file contains matrix dimensionality check functions
 */

#include <Eigen/Dense>
#include <stdexcept>
#include <type_traits>

/**
 * @enum MatrixOperation
 * @brief A set of possible operations for matrices.
 */
enum class MatrixOperation {
    MULTIPLY,
    ADD,
    SUBTRACT,
    ELEMENT_WISE_MULTIPLY,
    ELEMENT_WISE_DIVIDE
};

/**
 * @brief Checks matrix dimensions for a particular operation and throws an exception if they are inconsistent.
 * @tparam MatrixType1 Type of the first matrix
 * @tparam MatrixType2 Type of the second matrix
 * @param matrix1 First matrix
 * @param matrix2 Second matrix
 * @param operation The operation to be performed
 * @param error_msg Error message
 */
template<typename MatrixType1, typename MatrixType2>
void check_matrix_dimensions(const MatrixType1 &matrix1, const MatrixType2 &matrix2, MatrixOperation operation,
                             const std::string &error_msg = "") {
    std::string default_error_msg = "Inconsistent dimensions for matrix operation."
                                    + std::string("\n MatrixType1: ") + typeid(MatrixType1).name()
                                    + std::string("\n MatrixType2: ") + typeid(MatrixType2).name();
    switch (operation) {
        case MatrixOperation::MULTIPLY:
            if (matrix1.cols() != matrix2.rows()) {
                throw std::runtime_error(error_msg.empty() ? default_error_msg : error_msg);
            }
            break;
        case MatrixOperation::ADD:
        case MatrixOperation::SUBTRACT:
        case MatrixOperation::ELEMENT_WISE_MULTIPLY:
        case MatrixOperation::ELEMENT_WISE_DIVIDE:
            if (matrix1.rows() != matrix2.rows() || matrix1.cols() != matrix2.cols()) {
                throw std::runtime_error(error_msg.empty() ? default_error_msg : error_msg);
            }
            break;
    }
}

/**
 * @brief Helper structure to separate compile-time and runtime checks
 * @tparam MatrixType1 Type of the first matrix
 * @tparam MatrixType2 Type of the second matrix
 * @tparam Check Boolean value to check if dimensions are known at compile time
 */
template<typename MatrixType1, typename MatrixType2, bool Check = MatrixType1::ColsAtCompileTime != Eigen::Dynamic &&
                                                                  MatrixType2::RowsAtCompileTime != Eigen::Dynamic>
struct check_dimensions_helper;

template<typename MatrixType1, typename MatrixType2>
struct check_dimensions_helper<MatrixType1, MatrixType2, true> {
    static void check(const MatrixType1 &matrix1, const MatrixType2 &matrix2, MatrixOperation operation,
                      const std::string &error_msg = "") {
        switch (operation) {
            case MatrixOperation::MULTIPLY:
                static_assert(static_cast<int>(MatrixType1::ColsAtCompileTime) ==
                              static_cast<int>(MatrixType2::RowsAtCompileTime),
                              "Inconsistent dimensions for matrix multiplication");
                break;
            case MatrixOperation::ADD:
            case MatrixOperation::SUBTRACT:
            case MatrixOperation::ELEMENT_WISE_MULTIPLY:
            case MatrixOperation::ELEMENT_WISE_DIVIDE:
                static_assert(static_cast<int>(MatrixType1::RowsAtCompileTime) ==
                              static_cast<int>(MatrixType2::RowsAtCompileTime) &&
                              static_cast<int>(MatrixType1::ColsAtCompileTime) ==
                              static_cast<int>(MatrixType2::ColsAtCompileTime),
                              "Inconsistent dimensions for matrix operation");
                break;
        }

    }
};

template<typename MatrixType1, typename MatrixType2>
struct check_dimensions_helper<MatrixType1, MatrixType2, false> {
    static void check(const MatrixType1 &matrix1, const MatrixType2 &matrix2, MatrixOperation operation,
                      const std::string &error_msg = "") {
        check_matrix_dimensions(matrix1, matrix2, operation, error_msg);
    }
};

/**
 * @brief Checks the dimensions of two matrices and throws an error if they are inconsistent for the given operation
 * @tparam MatrixType1 Type of the first matrix
 * @tparam MatrixType2 Type of the second matrix
 * @param matrix1 First matrix
 * @param matrix2 Second matrix
 * @param operation The operation to be performed
 * @param error_msg Error message
 */
template<typename MatrixType1, typename MatrixType2>
void check_dimensions(const MatrixType1 &matrix1, const MatrixType2 &matrix2, MatrixOperation operation,
                      const std::string &error_msg = "") {
    check_dimensions_helper<MatrixType1, MatrixType2>::check(matrix1, matrix2, operation, error_msg);
}
