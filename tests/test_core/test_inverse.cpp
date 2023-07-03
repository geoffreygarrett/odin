#include <Eigen/Dense>
#include <include/gtest/gtest.h>

#include <eigen/inversion_policy.hpp> // assume this is where you defined your policy namespace

namespace {

    using Eigen::MatrixXd;

    TEST(InversionPolicyTest, DirectTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::direct<MatrixXd>::invert(matrix));
    }

    TEST(InversionPolicyTest, PinvTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::pinv<MatrixXd>::invert(matrix));
    }

    TEST(InversionPolicyTest, RegularizedTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::direct_regularized<MatrixXd>::invert(matrix));
    }

    TEST(InversionPolicyTest, PartialPivLuTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::partial_piv_lu::invert(matrix));
    }

    TEST(InversionPolicyTest, FullPivLuTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::full_piv_lu::invert(matrix));
    }

    TEST(InversionPolicyTest, HouseholderQrTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::householder_qr::invert(matrix));
    }

    TEST(InversionPolicyTest, ColPivHouseholderQrTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::col_piv_householder_qr::invert(matrix));
    }

    TEST(InversionPolicyTest, FullPivHouseholderQrTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::full_piv_householder_qr::invert(matrix));
    }

    TEST(InversionPolicyTest, CompleteOrthogonalDecompositionTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::complete_orthogonal_decomposition::invert(matrix));
    }

    TEST(InversionPolicyTest, LltTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        matrix = matrix * matrix.transpose(); // Ensure positive definiteness
        EXPECT_NO_THROW(policy::llt<MatrixXd>::invert(matrix));
    }

    TEST(InversionPolicyTest, LdltTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        matrix = matrix * matrix.transpose(); // Ensure semidefiniteness
        EXPECT_NO_THROW(policy::ldlt<MatrixXd>::invert(matrix));
    }

    TEST(InversionPolicyTest, BdcsSvdTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::bdcs_svd::invert(matrix));
    }

    TEST(InversionPolicyTest, JacobiSvdTest) {
        MatrixXd matrix = MatrixXd::Random(2, 2);
        EXPECT_NO_THROW(policy::jacobi_svd::invert(matrix));
    }

    TEST(InversionPolicyTest, DirectNonSquareTest) {
        MatrixXd matrix = MatrixXd::Random(2, 3);
        EXPECT_DEATH(policy::direct<MatrixXd>::invert(matrix),
                     "Assertion `rows\\(\\) == cols\\(\\)' failed");
    }

    TEST(InversionPolicyTest, RegularizedNonSquareTest) {
        MatrixXd matrix = MatrixXd::Random(2, 3);
        EXPECT_DEATH(policy::direct_regularized<MatrixXd>::invert(matrix),
                     "Assertion `aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)' failed");
    }

//    TEST(InversionPolicyTest, PartialPivLuNonInvertibleTest) {
//        MatrixXd matrix = MatrixXd::Zero(2, 2);
//        EXPECT_DEATH(policy::partial_piv_lu::invert(matrix),
//                     "Assertion `aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)' failed");
//    }

//    TEST(InversionPolicyTest, LltNonPositiveDefiniteTest) {
//        MatrixXd matrix = MatrixXd::Identity(2, 2);
//        matrix(0, 0) = -1; // make it not positive definite
//        EXPECT_ANY_THROW(policy::llt::invert(matrix));
//    }
//
//    TEST(InversionPolicyTest, LdltNonSemiDefiniteTest) {
//        MatrixXd matrix = MatrixXd::Identity(2, 2);
//        matrix(0, 0) = -2; // make it not semi-definite
//        EXPECT_ANY_THROW(policy::ldlt::invert(matrix));
//    }

}  // namespace
