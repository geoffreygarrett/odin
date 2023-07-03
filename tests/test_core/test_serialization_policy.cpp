#include <gtest/gtest.h>
#include <odin/core/policy/serialization/serialization_policy.hpp>
//#include <core/policy/serialization/eigen_serialization_policy.hpp>
//#include <core/policy/serialization/basic_serialization_policy.hpp>
#include <Eigen/Core>
#include <list>
#include <odin/core/policy/container/stl_associative_container_policy.hpp>


struct Person {
    std::string name;
    int age{};

    template<class Archive>
    void serialize(Archive &ar) {
        ar(CEREAL_NVP(name), CEREAL_NVP(age));
    }
};

template<typename ArchiveType>
void runNestedStructureSerializationTest() {
    SerializationPolicy<std::vector<Person>, ArchiveType> sp;
    std::vector<Person> people;

    people.push_back({"Alice", 25});
    people.push_back({"Bob", 30});

    std::string serialized = sp.serialize(people);
    ASSERT_FALSE(serialized.empty());

    std::vector<Person> deserialized = sp.deserialize(serialized);
    ASSERT_EQ(people.size(), deserialized.size());

    for (size_t i = 0; i < people.size(); i++) {
        ASSERT_EQ(people[i].name, deserialized[i].name);
        ASSERT_EQ(people[i].age, deserialized[i].age);
    }
}

TEST(SerializationPolicyTest, CanSerializeNestedStructure_JSON) {
    runNestedStructureSerializationTest<archive::json>();
}

TEST(SerializationPolicyTest, CanSerializeNestedStructure_XML) {
    runNestedStructureSerializationTest<archive::xml>();
}

TEST(SerializationPolicyTest, CanSerializeNestedStructure_Binary) {
    runNestedStructureSerializationTest<archive::binary>();
}

TEST(SerializationPolicyTest, CanSerializeNestedStructure_PortableBinary) {
    runNestedStructureSerializationTest<archive::portable_binary>();
}


TEST(SerializationPolicyTest, CanSerializeBinaryString) {
    SerializationPolicy<std::string, archive::binary> sp;
    std::string test_string = "Hello, world!";
    std::string serialized = sp.serialize(test_string);
    ASSERT_FALSE(serialized.empty());
    std::string deserialized = sp.deserialize(serialized);
    ASSERT_EQ(test_string, deserialized);
}


TEST(SerializationPolicyTest, CanSerializeString) {
    SerializationPolicy<std::string, archive::json> sp;
    std::string test_string = "Hello, world!";
    std::string serialized = sp.serialize(test_string);
    ASSERT_FALSE(serialized.empty());
    std::string deserialized = sp.deserialize(serialized);
    ASSERT_EQ(test_string, deserialized);
}


TEST(SerializationPolicyTest, CanSerializeVector) {
    SerializationPolicy<std::vector<int>> sp;
    std::vector<int> test_vector = {1, 2, 3, 4, 5};
    std::string serialized = sp.serialize(test_vector);
    ASSERT_FALSE(serialized.empty());
    std::vector<int> deserialized = sp.deserialize(serialized);
    ASSERT_EQ(test_vector, deserialized);
}

TEST(SerializationPolicyTest, CanSerializeMap) {
    SerializationPolicy<std::map<std::string, int>> sp;
    std::map<std::string, int> test_map = {{"one", 1},
                                           {"two", 2},
                                           {"three", 3}};
    std::string serialized = sp.serialize(test_map);
    ASSERT_FALSE(serialized.empty());
    std::map<std::string, int> deserialized = sp.deserialize(serialized);
    ASSERT_EQ(test_map, deserialized);
}

TEST(SerializationPolicyTest, CanSerializeComplex) {
    SerializationPolicy<std::complex<float>> sp;
    std::complex<float> test_complex = {1.0f, 2.0f};
    std::string serialized = sp.serialize(test_complex);
    ASSERT_FALSE(serialized.empty());
    std::complex<float> deserialized = sp.deserialize(serialized);
    ASSERT_EQ(test_complex, deserialized);
}


#define MATRIX_DIM 3// Matrix dimensions for the test

template<typename MatrixType, std::enable_if_t<
                                      MatrixType::RowsAtCompileTime == Eigen::Dynamic && MatrixType::ColsAtCompileTime == Eigen::Dynamic> * = nullptr>
void resizeMatrix(MatrixType &matrix, int size) {
    matrix.resize(size, size);
}

template<typename MatrixType, std::enable_if_t<
                                      MatrixType::RowsAtCompileTime == Eigen::Dynamic && MatrixType::ColsAtCompileTime != Eigen::Dynamic> * = nullptr>
void resizeMatrix(MatrixType &matrix, int size) {
    matrix.resize(size, MatrixType::ColsAtCompileTime);
}

template<typename MatrixType, std::enable_if_t<
                                      MatrixType::RowsAtCompileTime != Eigen::Dynamic && MatrixType::ColsAtCompileTime == Eigen::Dynamic> * = nullptr>
void resizeMatrix(MatrixType &matrix, int size) {
    matrix.resize(MatrixType::RowsAtCompileTime, size);
}

template<typename MatrixType, std::enable_if_t<
                                      MatrixType::RowsAtCompileTime != Eigen::Dynamic && MatrixType::ColsAtCompileTime != Eigen::Dynamic> * = nullptr>
void resizeMatrix(MatrixType &matrix, int size) {
    // Fixed size, no need to resize. We don't need the 'size' argument in this case.
    // Maybe we could assert here that RowsAtCompileTime and ColsAtCompileTime are equal to the size, just to be safe?
}


template<typename TypeParam>
class EigenMatrixSerializationPolicyTest : public ::testing::Test {
protected:
    using ArchiveType = typename TypeParam::first_type;
    using MatrixType = typename TypeParam::second_type;
    SerializationPolicy<MatrixType, ArchiveType> sp;

    MatrixType matrix;
    static constexpr int RowsAtCompileTime = MatrixType::RowsAtCompileTime;
    static constexpr int ColsAtCompileTime = MatrixType::ColsAtCompileTime;

    void SetUp() override {
        resizeMatrix(matrix, MATRIX_DIM);
        // Fill the matrix with random values
        matrix.setRandom();
    }
};

// Archive types and Matrix types
using ConfigTypes = ::testing::Types<
        // Double precision types
        std::pair<archive::json, Eigen::MatrixXd>,
        std::pair<archive::xml, Eigen::MatrixXd>,
        std::pair<archive::binary, Eigen::MatrixXd>,
        std::pair<archive::portable_binary, Eigen::MatrixXd>,
        std::pair<archive::json, Eigen::VectorXd>,
        std::pair<archive::xml, Eigen::VectorXd>,
        std::pair<archive::binary, Eigen::VectorXd>,
        std::pair<archive::portable_binary, Eigen::VectorXd>,
        std::pair<archive::json, Eigen::RowVectorXd>,
        std::pair<archive::xml, Eigen::RowVectorXd>,
        std::pair<archive::binary, Eigen::RowVectorXd>,
        std::pair<archive::portable_binary, Eigen::RowVectorXd>,
        std::pair<archive::json, Eigen::Matrix3d>,
        std::pair<archive::xml, Eigen::Matrix3d>,
        std::pair<archive::binary, Eigen::Matrix3d>,
        std::pair<archive::portable_binary, Eigen::Matrix3d>,
        std::pair<archive::json, Eigen::Vector3d>,
        std::pair<archive::xml, Eigen::Vector3d>,
        std::pair<archive::binary, Eigen::Vector3d>,
        std::pair<archive::portable_binary, Eigen::Vector3d>,
        std::pair<archive::json, Eigen::RowVector3d>,
        std::pair<archive::xml, Eigen::RowVector3d>,
        std::pair<archive::binary, Eigen::RowVector3d>,
        std::pair<archive::portable_binary, Eigen::RowVector3d>,
        // Single precision types
        std::pair<archive::json, Eigen::MatrixXf>,
        std::pair<archive::xml, Eigen::MatrixXf>,
        std::pair<archive::binary, Eigen::MatrixXf>,
        std::pair<archive::portable_binary, Eigen::MatrixXf>,
        std::pair<archive::json, Eigen::VectorXf>,
        std::pair<archive::xml, Eigen::VectorXf>,
        std::pair<archive::binary, Eigen::VectorXf>,
        std::pair<archive::portable_binary, Eigen::VectorXf>,
        std::pair<archive::json, Eigen::RowVectorXf>,
        std::pair<archive::xml, Eigen::RowVectorXf>,
        std::pair<archive::binary, Eigen::RowVectorXf>,
        std::pair<archive::portable_binary, Eigen::RowVectorXf>,
        std::pair<archive::json, Eigen::Matrix3f>,
        std::pair<archive::xml, Eigen::Matrix3f>,
        std::pair<archive::binary, Eigen::Matrix3f>,
        std::pair<archive::portable_binary, Eigen::Matrix3f>,
        std::pair<archive::json, Eigen::Vector3f>,
        std::pair<archive::xml, Eigen::Vector3f>,
        std::pair<archive::binary, Eigen::Vector3f>,
        std::pair<archive::portable_binary, Eigen::Vector3f>,
        std::pair<archive::json, Eigen::RowVector3f>,
        std::pair<archive::xml, Eigen::RowVector3f>,
        std::pair<archive::binary, Eigen::RowVector3f>,
        std::pair<archive::portable_binary, Eigen::RowVector3f>,
        // Extra
        std::pair<archive::json, Eigen::MatrixXi>,
        std::pair<archive::xml, Eigen::MatrixXi>,
        std::pair<archive::binary, Eigen::MatrixXi>,
        std::pair<archive::portable_binary, Eigen::MatrixXi>,
        std::pair<archive::json, Eigen::MatrixXcd>,
        std::pair<archive::xml, Eigen::MatrixXcd>,
        std::pair<archive::binary, Eigen::MatrixXcd>,
        std::pair<archive::portable_binary, Eigen::MatrixXcd>
        //        std::pair<archive::json, Eigen::ArrayXd>
        //        std::pair<archive::xml, Eigen::ArrayXd>,
        //        std::pair<archive::binary, Eigen::ArrayXd>,
        //        std::pair<archive::portable_binary, Eigen::ArrayXd>
        //        std::pair<archive::json, Eigen::ArrayXXd>,
        //        std::pair<archive::xml, Eigen::ArrayXXd>,
        //        std::pair<archive::binary, Eigen::ArrayXXd>,
        //        std::pair<archive::portable_binary, Eigen::ArrayXXd>,
        //        std::pair<archive::json, Eigen::Array3d>,
        //        std::pair<archive::xml, Eigen::Array3d>,
        //        std::pair<archive::binary, Eigen::Array3d>
        //        std::pair<archive::portable_binary, Eigen::Array3d>
        >;

TYPED_TEST_SUITE(EigenMatrixSerializationPolicyTest, ConfigTypes);

TYPED_TEST(EigenMatrixSerializationPolicyTest, CanSerializeAndDeserialize) {
    // Get information about the current test
    const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();

    // Set the golden file name based on the currently running test
    std::string golden_file_name =
            std::string(test_info->test_suite_name()) + "_" + std::string(test_info->name()) + ".golden";

    // Print the name of the currently running test suite and test
    //    printf("We are in test %s of test suite %s.\n", test_info->name(), test_info->test_suite_name());

    // Serialize the matrix
    std::string serialized = this->sp.serialize(this->matrix);
    ASSERT_FALSE(serialized.empty());

    // Deserialize the matrix
    typename TestFixture::MatrixType deserialized = this->sp.deserialize(serialized);

    // Compare the matrices
    ASSERT_TRUE(this->matrix.isApprox(deserialized));
}


// Additional tests could include:
// - Serializing empty containers
// - Deserializing empty data
// - Handling serialization/deserialization of unsupported types
// - Error handling and exceptions
