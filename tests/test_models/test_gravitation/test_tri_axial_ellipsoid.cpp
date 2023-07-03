#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <odin/logging.hpp>
#include <odin/models/gravitational/ellipsoidal.hpp>


// Define a type-parameterized test case
template<class T>
class EllipsoidTest : public ::testing::Test {
protected:
    void SetUp() override {
        ellipsoid = new Ellipsoid<T>(2.0, 3.0, 4.0, 5.0);
        position  = Eigen::Vector3<T>(1.0, 1.0, 1.0);
    }

    void TearDown() override {
        delete ellipsoid;
    }

    Ellipsoid<T>     *ellipsoid;
    Eigen::Vector3<T> position;
};

// Types to test
typedef ::testing::Types<float, double> MyTypes;
TYPED_TEST_SUITE(EllipsoidTest, MyTypes);

TYPED_TEST(EllipsoidTest, PotentialTest) {
    auto potential = this->ellipsoid->potential(this->position);
    // Here you need to define the expected value based on your domain knowledge
    // or it might be some hand calculated values
    TypeParam expected_value = static_cast<TypeParam>(10.0);
    EXPECT_NEAR(potential, expected_value, 1e-5);
}

TYPED_TEST(EllipsoidTest, AccelerationTest) {
    auto acceleration = this->ellipsoid->acceleration(this->position);
    // Here you need to define the expected value based on your domain knowledge
    // or it might be some hand calculated values
    Eigen::Vector3<TypeParam> expected_value(2.0, 2.0, 2.0);
    EXPECT_TRUE(acceleration.isApprox(expected_value, 1e-5));
}

int main(int argc, char **argv) {
    INIT_ODIN_LOGGING("test_tri_axial_ellipsoid", "./log/test_tri_axial_ellipsoid.log");
    ODIN_SET_VLOG_LEVEL(100);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}