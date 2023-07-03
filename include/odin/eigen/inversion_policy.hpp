
#ifndef MATRIX_INVERSION_INVERSION_POLICIES_H
#define MATRIX_INVERSION_INVERSION_POLICIES_H
/**
 * @file
 *
 * | Decomposition Method           | Policy Name                       | Requirements on the Matrix       | Speed (small-to-medium)| Speed (large) | Accuracy |
 * |--------------------------------|-----------------------------------|----------------------------------|------------------------|---------------|----------|
 * | Inversion                      | direct                            | Square matrix                    | ++                     | ++            | ++       |
 * | Pseudoinverse                  | pinv                              | None                             | +                      | -             | +++      |
 * | Regularized Inversion          | regularized                       | Square matrix                    | +                      | +             | ++       |
 * | PartialPivLU                   | partial_piv_lu                    | Invertible                       | ++                     | ++            | +        |
 * | FullPivLU                      | full_piv_lu                       | None                             | -                      | --            | +++      |
 * | HouseholderQR                  | householder_qr                    | None                             | ++                     | ++            | +        |
 * | ColPivHouseholderQR            | col_piv_householder_qr            | None                             | +                      | -             | +++      |
 * | FullPivHouseholderQR           | full_piv_householder_qr           | None                             | -                      | --            | +++      |
 * | CompleteOrthogonalDecomposition| complete_orthogonal_decomposition | None                             | +                      | -             | +++      |
 * | LLT                            | llt                               | Positive definite                | +++                    | +++           | +        |
 * | LDLT                           | ldlt                              | Positive or negative semidefinite| +++                    | +             | ++       |
 * | BDCSVD                         | bdcs_svd                          | None                             | -                      | -             | +++      |
 * | JacobiSVD                      | jacobi_svd                        | None                             | -                      | ---           | +++      |
 */

namespace policy {

    /** Base policy for inversion. This must be overridden by a specific implementation. */
    struct inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(always_false<M>::value, "inversion_policy::invert() not implemented");
            return matrix;// this will never be executed
        }

    private:
        template<typename T>
        struct always_false : std::false_type {
        };
    };


    struct direct {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(M::RowsAtCompileTime == M::ColsAtCompileTime,
                          "Matrix must be square to compute the inverse directly");
            return matrix.inverse();
        }
    };

//    template<typename MatrixType, typename VectorType>
//    MatrixType pinv_helper(const MatrixType &matrix) {
//        using Scalar = typename MatrixType::Scalar;
//        constexpr Scalar epsilon = std::numeric_limits<Scalar>::epsilon();
//        Eigen::JacobiSVD<MatrixType> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
//
//        VectorType singularValues = svd.singularValues();
//        for (long i = 0; i < singularValues.rows(); i++) {
//            if (std::abs(singularValues(i)) > epsilon)
//                singularValues(i) = 1.0 / singularValues(i);
//            else
//                singularValues(i) = 0;
//        }
//
//        MatrixType singularValuesInverse = singularValues.asDiagonal();
//        return svd.matrixV() * singularValuesInverse * svd.matrixU().adjoint();
//    }
    template<typename MatrixType, typename VectorType>
    MatrixType pinv_helper(const MatrixType &matrix) {
        using Scalar = typename MatrixType::Scalar;
        constexpr Scalar epsilon = std::numeric_limits<Scalar>::epsilon();

        // Use a dynamic matrix for the SVD computation
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> dynamicMatrix = matrix;

        Eigen::JacobiSVD<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> svd(dynamicMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

        VectorType singularValues = svd.singularValues();
        for (long i = 0; i < singularValues.rows(); i++) {
            if (std::abs(singularValues(i)) > epsilon)
                singularValues(i) = 1.0 / singularValues(i);
            else
                singularValues(i) = 0;
        }

        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> singularValuesInverse = singularValues.asDiagonal();
        return svd.matrixV() * singularValuesInverse * svd.matrixU().adjoint();
    }

    struct pinv {
        template<typename MatrixType>
        static MatrixType invert(const MatrixType &matrix) {
            using Scalar = typename MatrixType::Scalar;
            using VectorType = Eigen::Vector<Scalar, Eigen::Dynamic>;
            return pinv_helper<MatrixType, VectorType>(matrix);
        }
    };

//    template<typename Scalar, int Rows, int Cols>
//    struct pinv<Eigen::Matrix<Scalar, Rows, Cols>> {
//        using MatrixType = Eigen::Matrix<Scalar, Rows, Cols>;
//        using DiagonalMatrixType = Eigen::DiagonalMatrix<Scalar, Rows, Cols>;
//        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
//
//        static MatrixType invert(const MatrixType &matrix) {
//            return pinv_helper<MatrixType, VectorType>(matrix);
//        }
//    };

    struct direct_regularized {
        template<typename M, int PowerOf10 = -3>
        static M invert(const M &matrix) {
            static_assert(M::RowsAtCompileTime == M::ColsAtCompileTime,
                          "Matrix must be square to compute the direct regularized inverse");

            using Scalar = typename M::Scalar;
            Scalar lambda = std::pow(Scalar(10), Scalar(PowerOf10));
            M identity = M::Identity(matrix.rows(), matrix.cols());
            return (matrix.transpose() * matrix + lambda * identity).ldlt().solve(matrix.transpose());
        }
    };


    struct ldlt {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(M::RowsAtCompileTime == M::ColsAtCompileTime,
                          "Matrix must be square to compute the LDLT decomposition");
            return matrix.ldlt().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };


    struct llt {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(M::RowsAtCompileTime == M::ColsAtCompileTime,
                          "Matrix must be square to compute the LLT decomposition");
            return matrix.llt().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };


    /** Policy for HouseholderQR decomposition-based inversion. Can be applied to any matrix. */
    struct householder_qr : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(
                    M::RowsAtCompileTime == M::ColsAtCompileTime,
                    "Matrix must be square to compute the Householder QR decomposition");
            return matrix.householderQr().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

    /** Policy for CompleteOrthogonalDecomposition-based inversion. Can be applied to any matrix. */
    struct complete_orthogonal_decomposition : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(
                    M::RowsAtCompileTime == M::ColsAtCompileTime,
                    "Matrix must be square to compute the CompleteOrthogonalDecomposition");
            return matrix.completeOrthogonalDecomposition().pseudoInverse();
        }
    };

    /** Policy for PartialPivLU decomposition-based inversion. Requires an invertible matrix. */
    struct partial_piv_lu : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(
                    M::RowsAtCompileTime == M::ColsAtCompileTime,
                    "Matrix must be invertible to compute the PartialPivLU decomposition");
            return matrix.partialPivLu().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

    /** Policy for FullPivLU decomposition-based inversion. Can be applied to any matrix. */
    struct full_piv_lu : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(
                    M::RowsAtCompileTime == M::ColsAtCompileTime,
                    "Matrix must be square to compute the FullPivLU decomposition");
            return matrix.fullPivLu().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

    /** Policy for ColPivHouseholderQR decomposition-based inversion. Can be applied to any matrix. */
    struct col_piv_householder_qr : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(
                    M::RowsAtCompileTime == M::ColsAtCompileTime,
                    "Matrix must be square to compute the ColPivHouseholderQR decomposition");
            return matrix.colPivHouseholderQr().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

    /** Policy for FullPivHouseholderQR decomposition-based inversion. Can be applied to any matrix. */
    struct full_piv_householder_qr : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            static_assert(
                    M::RowsAtCompileTime == M::ColsAtCompileTime,
                    "Matrix must be square to compute the FullPivHouseholderQR decomposition");
            return matrix.fullPivHouseholderQr().solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

    /** Policy for BDCSVD decomposition-based inversion. Can be applied to any matrix. */
    struct bdcs_svd : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            Eigen::BDCSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            return svd.solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

    /** Policy for JacobiSVD decomposition-based inversion. Can be applied to any matrix. */
    struct jacobi_svd : public inversion_policy {
        template<typename M>
        static M invert(const M &matrix) {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            return svd.solve(M::Identity(matrix.rows(), matrix.cols()));
        }
    };

}// namespace policy


#endif//MATRIX_INVERSION_INVERSION_POLICIES_H