#include <include/gtest/gtest.h>

// autodiff include
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual.hpp>

using namespace autodiff;

// The vector function for which the Jacobian is needed
VectorXreal f1(const VectorXreal &x) {
    return x * x.sum();
}

// The vector function with parameters for which the Jacobian is needed
VectorXreal f2(const VectorXreal &x, const VectorXreal &p, const real &q) {
    return x * p.sum() * exp(q);
}

// The vector function for which the Jacobian is needed
VectorXreal f3(const VectorXreal &x) {
    return x * x.sum();
}


// Test the first Jacobian computation
TEST(AutoDiffTest, JacobianTest1) {
    using Eigen::MatrixXd;

    VectorXreal x(5);                           // the input vector x with 5 variables
    x << 1, 2, 3, 4, 5;                         // x = [1, 2, 3, 4, 5]

    VectorXreal F;                              // the output vector F = f(x) evaluated together with Jacobian matrix below

    MatrixXd J = jacobian(f1, wrt(x), at(x), F); // evaluate the output vector F and the Jacobian matrix dF/dx

    // Print the output vector and Jacobian matrix
//    std::cout << "F = \n" << F << std::endl;    // print the evaluated output vector F
//    std::cout << "J = \n" << J << std::endl;    // print the evaluated Jacobian matrix dF/dx

    // Assert the size of the output vector and Jacobian matrix
    ASSERT_EQ(F.size(), 5);
    ASSERT_EQ(J.rows(), 5);
    ASSERT_EQ(J.cols(), 5);
}

//// Test the second Jacobian computation
TEST(AutoDiffTest, JacobianTest2) {
    using Eigen::MatrixXd;

    VectorXreal x(5);    // the input vector x with 5 variables
    x << 1, 2, 3, 4, 5;  // x = [1, 2, 3, 4, 5]

    VectorXreal p(3);    // the input parameter vector p with 3 variables
    p << 1, 2, 3;        // p = [1, 2, 3]

    real q = -2;         // the input parameter q as a single variable

    VectorXreal F;       // the output vector F = f(x, p, q) evaluated together with Jacobian below

    MatrixXd Jx = jacobian(f2, wrt(x), at(x, p, q),
                           F);       // evaluate the function and the Jacobian matrix Jx = dF/dx
    MatrixXd Jp = jacobian(f2, wrt(p), at(x, p, q),
                           F);       // evaluate the function and the Jacobian matrix Jp = dF/dp
    MatrixXd Jq = jacobian(f2, wrt(q), at(x, p, q),
                           F);       // evaluate the function and the Jacobian matrix Jq = dF/dq
    MatrixXd Jqpx = jacobian(f2, wrt(q, p, x), at(x, p, q),
                             F); // evaluate the function and the Jacobian matrix Jqpx = [dF/dq, dF/dp, dF/dx]

    // Print the output vector and Jacobian matrices
//    std::cout << "F = \n" << F << std::endl;     // print the evaluated output vector F
//    std::cout << "Jx = \n" << Jx << std::endl;   // print the evaluated Jacobian matrix dF/dx
//    std::cout << "Jp = \n" << Jp << std::endl;   // print the evaluated Jacobian matrix dF/dp
//    std::cout << "Jq = \n" << Jq << std::endl;   // print the evaluated Jacobian matrix dF/dq
//    std::cout << "Jqpx = \n" << Jqpx << std::endl; // print the evaluated Jacobian matrix [dF/dq, dF/dp, dF/dx]

    // Assert the size of the output vector and Jacobian matrices
    ASSERT_EQ(F.size(), 5);
    ASSERT_EQ(Jx.rows(), 5);
    ASSERT_EQ(Jx.cols(), 5);
    ASSERT_EQ(Jp.rows(), 5);
    ASSERT_EQ(Jp.cols(), 3);
    ASSERT_EQ(Jq.rows(), 5);
    ASSERT_EQ(Jq.cols(), 1);
    ASSERT_EQ(Jqpx.rows(), 5);
    ASSERT_EQ(Jqpx.cols(), 9);
}

// Test the third Jacobian computation
TEST(AutoDiffTest, JacobianTest3) {
    using Eigen::Map;
    using Eigen::MatrixXd;

    VectorXreal x(5);                           // the input vector x with 5 variables
    x << 1, 2, 3, 4, 5;                         // x = [1, 2, 3, 4, 5]
    double y[25];                               // the output Jacobian as a flat array

    VectorXreal F;                              // the output vector F = f(x) evaluated together with Jacobian matrix below
    Map<MatrixXd> J(y, 5, 5);                   // the output Jacobian dF/dx mapped onto the flat array

    jacobian(f3, wrt(x), at(x), F, J);           // evaluate the output vector F and the Jacobian matrix dF/dx

    // Print the output vector, Jacobian matrix, and flat array
//    std::cout << "F = \n" << F << std::endl;    // print the evaluated output vector F
//    std::cout << "J = \n" << J << std::endl;    // print the evaluated Jacobian matrix dF/dx
//    std::cout << "y = ";                        // print the flat array
//    for (int i = 0; i < 25; ++i)
//        std::cout << y[i] << " ";
//    std::cout << std::endl;

    // Assert the size of the output vector, Jacobian matrix, and flat array
    ASSERT_EQ(F.size(), 5);
    ASSERT_EQ(J.rows(), 5);
    ASSERT_EQ(J.cols(), 5);
    ASSERT_EQ(sizeof(y) / sizeof(y[0]), 25);
}