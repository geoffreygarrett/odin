#include <Eigen/Dense>
#include <eigen.hpp>
#include <include/gtest/gtest.h>

// Note that eigen assertions can be changes with this macro
// #define eigen_assert(X) do { if(!(X)) throw std::runtime_error(#X); } while(false);
// https://stackoverflow.com/questions/43108121/what-is-the-best-way-to-care-about-assertion-fail-in-eigen-c

// Test check_dimensions for fixed matrices
TEST(CheckDimensions, StaticMatrices) {
    Eigen::Matrix3d m1; // 3x3 matrix
    Eigen::Matrix3d m2; // 3x3 matrix
    Eigen::Matrix4d m3; // 4x4 matrix

    // Test for compatible operations
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));

    // The following lines should cause compile-time errors
    // EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::MULTIPLY));
    // EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::ADD));
    // EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::SUBTRACT));
    // EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    // EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::ELEMENT_WISE_DIVIDE));
}

// Test check_dimensions for dynamic matrices
TEST(CheckDimensions, DynamicMatrices) {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(5, 4);
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(4, 3);
    Eigen::MatrixXd m3 = Eigen::MatrixXd::Random(5, 4);
    Eigen::MatrixXd m4 = Eigen::MatrixXd::Random(4, 5);

    // Test multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::MULTIPLY));

    // Test addition
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::ADD));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    EXPECT_ANY_THROW(check_dimensions(m1, m4, MatrixOperation::ADD));

    // Test subtraction
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::SUBTRACT));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));

    // Test element-wise multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));

    // Test element-wise division
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::ELEMENT_WISE_DIVIDE));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
}

// Test check_dimensions with different matrix types
TEST(CheckDimensions, DifferentMatrixTypes) {
    Eigen::MatrixXf m1 = Eigen::MatrixXf::Random(5, 4);
    Eigen::MatrixXf m2 = Eigen::MatrixXf::Random(4, 3);
    Eigen::MatrixXf m3 = Eigen::MatrixXf::Random(5, 4);

    // Test multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m1, m3, MatrixOperation::MULTIPLY));

    // Test addition
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::ADD));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));

    // Test subtraction
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::SUBTRACT));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));

    // Test element-wise multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));

    // Test element-wise division
    EXPECT_NO_THROW(check_dimensions(m1, m3, MatrixOperation::ELEMENT_WISE_DIVIDE));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
}

// Test operations for mixed types that should fail
TEST(CheckDimensions, MixedTypeFailures) {
    Eigen::Matrix3d m_static; // 3x3 static matrix
    Eigen::MatrixXd m_dynamic = Eigen::MatrixXd::Random(4, 4); // 4x4 dynamic matrix

    // The following operations should throw a runtime exception
    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ADD));
    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::SUBTRACT));
    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ELEMENT_WISE_DIVIDE));
}

// Test operations for dynamic matrices that should fail
TEST(CheckDimensions, DynamicMatrixFailures) {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(3, 3); // 3x3 dynamic matrix
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(4, 4); // 4x4 dynamic matrix

    // The following operations should throw a runtime exception
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
}

// Test operations for static matrices that should fail
TEST(CheckDimensions, StaticMatrixFailures) {
    Eigen::Matrix3d m1; // 3x3 static matrix
    Eigen::Matrix4d m2; // 4x4 static matrix

    // The following lines should cause compile-time errors
    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
}

// Test operations for mixed types
TEST(CheckDimensions, MixedTypeOperations) {
    Eigen::Matrix3d m_static; // 3x3 static matrix
    Eigen::MatrixXd m_dynamic = Eigen::MatrixXd::Random(3, 3); // 3x3 dynamic matrix

    // Test multiplication
    EXPECT_NO_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::MULTIPLY));
    EXPECT_NO_THROW(m_static * m_dynamic);

    // Test addition
    EXPECT_NO_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ADD));
    EXPECT_NO_THROW(m_static + m_dynamic);

    // Test subtraction
    EXPECT_NO_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::SUBTRACT));
    EXPECT_NO_THROW(m_static - m_dynamic);

    // Test element-wise multiplication
    EXPECT_NO_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_NO_THROW(m_static.array() * m_dynamic.array());

    // Test element-wise division
    EXPECT_NO_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ELEMENT_WISE_DIVIDE));
    EXPECT_NO_THROW(m_static.array() / m_dynamic.array());
}

// Test operations for static matrices
TEST(CheckDimensions, StaticMatrixOperations) {
    Eigen::Matrix3d m1; // 3x3 static matrix
    Eigen::Matrix3d m2; // 3x3 static matrix

    // Test multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    EXPECT_NO_THROW(m1 * m2);

    // Test addition
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    EXPECT_NO_THROW(m1 + m2);

    // Test subtraction
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    EXPECT_NO_THROW(m1 - m2);

    // Test element-wise multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_NO_THROW(m1.array() * m2.array());

    // Test element-wise division
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
    EXPECT_NO_THROW(m1.array() / m2.array());
}

// Test operations for dynamic matrices
TEST(CheckDimensions, DynamicMatrixOperations) {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(3, 3); // 3x3 dynamic matrix
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3); // 3x3 dynamic matrix

    // Test multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    EXPECT_NO_THROW(m1 * m2);

    // Test addition
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    EXPECT_NO_THROW(m1 + m2);

    // Test subtraction
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    EXPECT_NO_THROW(m1 - m2);

    // Test element-wise multiplication
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    EXPECT_NO_THROW(m1.array() * m2.array());

    // Test element-wise division
    EXPECT_NO_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
    EXPECT_NO_THROW(m1.array() / m2.array());
}

// Test operations for mixed types that should fail
TEST(CheckDimensions, MixedTypeOperationsFailures) {
    Eigen::Matrix3d m_static; // 3x3 static matrix
    Eigen::MatrixXd m_dynamic = Eigen::MatrixXd::Random(4, 4); // 4x4 dynamic matrix

    // The following operations should throw a runtime exception in both the check and the operation
    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::MULTIPLY));
    ASSERT_DEATH(m_static * m_dynamic, "lhs.cols\\(\\) == rhs.rows\\(\\).*invalid matrix product");

    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ADD));
    ASSERT_DEATH(m_static + m_dynamic, "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::SUBTRACT));
    ASSERT_DEATH(m_static - m_dynamic, "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    ASSERT_DEATH(m_static.array() * m_dynamic.array(), "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

    EXPECT_ANY_THROW(check_dimensions(m_static, m_dynamic, MatrixOperation::ELEMENT_WISE_DIVIDE));
    ASSERT_DEATH(m_static.array() / m_dynamic.array(),  "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

}

// Test operations for static matrices that should fail
TEST(CheckDimensions, StaticMatrixOperationsFailures) {
    Eigen::Matrix3d m1; // 3x3 static matrix
    Eigen::Matrix4d m2; // 4x4 static matrix

    // The following lines should cause compile-time errors in both the check and the operation
    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    // EXPECT_ANY_THROW(m1 * m2);

    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    // EXPECT_ANY_THROW(m1 + m2);

    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    // EXPECT_ANY_THROW(m1 - m2);

    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    // EXPECT_ANY_THROW(m1.array() * m2.array());

    // EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
    // EXPECT_ANY_THROW(m1.array() / m2.array());
}

// Test operations for dynamic matrices that should fail
TEST(CheckDimensions, DynamicMatrixOperationsFailures) {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(3, 3); // 3x3 dynamic matrix
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(4, 4); // 4x4 dynamic matrix

    // The following operations should throw a runtime exception in both the check and the operation
    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::MULTIPLY));
    ASSERT_DEATH(m1 * m2, "lhs.cols\\(\\) == rhs.rows\\(\\).*invalid matrix product");

    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ADD));
    ASSERT_DEATH(m1 + m2, "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::SUBTRACT));
    ASSERT_DEATH(m1 - m2, "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_MULTIPLY));
    ASSERT_DEATH(m1.array() * m2.array(), "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");

    EXPECT_ANY_THROW(check_dimensions(m1, m2, MatrixOperation::ELEMENT_WISE_DIVIDE));
    ASSERT_DEATH(m1.array() / m2.array(),  "aLhs.rows\\(\\) == aRhs.rows\\(\\) && aLhs.cols\\(\\) == aRhs.cols\\(\\)");
}


