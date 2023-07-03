#include <include/gtest/gtest.h>

#include <list>
#include <odin/core/policy/container/container_policy.hpp>
#include <odin/core/policy/container/stl_associative_container_policy.hpp>
#include <odin/core/policy/serialization/serialization_policy.hpp>

template<typename T>
T create_from_int(int val);

template<>
std::string create_from_int<std::string>(int val) {
    return std::to_string(val);
}

template<>
double create_from_int<double>(int val) {
    return static_cast<double>(val);
}

template<>
int create_from_int<int>(int val) {
    return val;
}

// Test Fixture for Map Container Policy
template<typename T>
class MapContainerPolicyTest : public ::testing::Test {
protected:
    using MapType = T;
    MapType container;
    ContainerPolicy<MapType> policy;
};

// Specify the types to be used for the tests
using MapTypes = ::testing::Types<std::map<int, std::string>, std::map<std::string, double>>;
TYPED_TEST_SUITE(MapContainerPolicyTest, MapTypes);

//TYPED_TEST(MapContainerPolicyTest, Serialization) {
//    using Key = typename TypeParam::key_type;
//    using Value = typename TypeParam::mapped_type;
//
//    for (int i = 0; i < 10; ++i) {
//        Key key = create_from_int<Key>(i);
//        Value value = create_from_int<Value>(i);
//        this->policy.add(this->container, std::make_pair(key, value));
//    }
//
//    std::stringstream ss;
//    {
//        cereal::JSONOutputArchive oarchive(ss);
//        this->policy.save_impl(oarchive, this->container);
//    }
//
////    std::cout<<ss.str()<<std::endl;
//
//    TypeParam new_container;
//    {
//        cereal::JSONInputArchive iarchive(ss);
//        this->policy.load_impl(iarchive, new_container);
//    }
//
//    ASSERT_EQ(this->container, new_container);
//}


TYPED_TEST(MapContainerPolicyTest, AddAndClear) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    const Key key1 = Key();// Make this const
    Value value1 = Value();
    this->policy.add(this->container, std::make_pair(key1, value1));

    ASSERT_EQ(this->container.size(), 1);
    ASSERT_EQ(this->container.begin()->first, key1);
    ASSERT_EQ(this->container.begin()->second, value1);

    this->policy.clear(this->container);
    ASSERT_EQ(this->container.size(), 0);
}

TYPED_TEST(MapContainerPolicyTest, MultipleAddAndClear) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    for (int i = 0; i < 10; ++i) {
        Key key = create_from_int<Key>(i);
        Value value = create_from_int<Value>(i);
        this->policy.add(this->container, std::make_pair(key, value));
    }

    ASSERT_EQ(this->container.size(), 10);
    this->policy.clear(this->container);
    ASSERT_EQ(this->container.size(), 0);
}

TYPED_TEST(MapContainerPolicyTest, SizeInBytes) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    Key key1 = Key();
    Value value1 = Value();
    this->policy.add(this->container, std::make_pair(key1, value1));

    // You might need to adjust the expected size depending on the implementation details and the types used.
    size_t expected_size = sizeof(std::pair<const Key, Value>) + 3 * sizeof(void *) + sizeof(bool);
    ASSERT_EQ(this->policy.size_in_bytes(this->container), expected_size);
}

TYPED_TEST(MapContainerPolicyTest, EmptyMap) {
    ASSERT_EQ(this->container.size(), 0);
    ASSERT_EQ(this->policy.size_in_bytes(this->container), 0);
}


TYPED_TEST(MapContainerPolicyTest, AddUniqueItems) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    for (int i = 0; i < 10; ++i) {
        Key key = create_from_int<Key>(i);
        Value value = create_from_int<Value>(i);
        this->policy.add(this->container, key, value);
    }

    ASSERT_EQ(this->container.size(), 10);
}

TYPED_TEST(MapContainerPolicyTest, AddDuplicateItems) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    Key key = create_from_int<Key>(1);
    Value value = create_from_int<Value>(1);

    for (int i = 0; i < 10; ++i) {
        this->policy.add(this->container, key, value);
    }

    ASSERT_EQ(this->container.size(), 1);
}

TYPED_TEST(MapContainerPolicyTest, SizeInBytesAfterMultipleAdd) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    for (int i = 0; i < 10; ++i) {
        Key key = create_from_int<Key>(i);
        Value value = create_from_int<Value>(i);
        this->policy.add(this->container, key, value);
    }

    size_t expected_size =
            this->container.size() * (sizeof(std::pair<const Key, Value>) + 3 * sizeof(void *) + sizeof(bool));
    ASSERT_EQ(this->policy.size_in_bytes(this->container), expected_size);
}

TYPED_TEST(MapContainerPolicyTest, ClearAndCheckMemory) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;

    for (int i = 0; i < 10; ++i) {
        Key key = create_from_int<Key>(i);
        Value value = create_from_int<Value>(i);
        this->policy.add(this->container, key, value);
    }

    this->policy.clear(this->container);
    ASSERT_EQ(this->policy.size_in_bytes(this->container), 0);
}

TYPED_TEST(MapContainerPolicyTest, SizeAtEachStep) {
    using Key = typename TypeParam::key_type;
    using Value = typename TypeParam::mapped_type;
    size_t expected_size = 0;

    for (int i = 0; i < 10; ++i) {
        Key key = create_from_int<Key>(i);
        Value value = create_from_int<Value>(i);
        this->policy.add(this->container, key, value);

        expected_size += sizeof(std::pair<const Key, Value>) + 3 * sizeof(void *) + sizeof(bool);
        ASSERT_EQ(this->container.size(), i + 1);
        ASSERT_EQ(this->policy.size_in_bytes(this->container), expected_size);
    }

    this->policy.clear(this->container);
    ASSERT_EQ(this->container.size(), 0);
    ASSERT_EQ(this->policy.size_in_bytes(this->container), 0);
}

#include <odin/core/policy/container/stl_sequence_container_policy.hpp>

template<typename T>
class SequenceContainerPolicyTest : public ::testing::Test {
protected:
    T container;
    ContainerPolicy<T> policy;
};

typedef ::testing::Types<
        std::vector<int>,
        std::list<int>,
        std::deque<int>,
        std::vector<std::string>,
        std::list<std::string>,
        std::deque<std::string>>
        SequenceContainerTypes;

TYPED_TEST_SUITE(SequenceContainerPolicyTest, SequenceContainerTypes);

TYPED_TEST(SequenceContainerPolicyTest, AddAndClear) {
    using ElementType = typename TypeParam::value_type;

    for (int i = 0; i < 10; ++i) {
        this->policy.add(this->container, create_from_int<ElementType>(i));
    }

    ASSERT_EQ(this->container.size(), 10);
    this->policy.clear(this->container);
    ASSERT_EQ(this->container.size(), 0);
}

TYPED_TEST(SequenceContainerPolicyTest, SizeInBytes) {
    using ElementType = typename TypeParam::value_type;

    for (int i = 0; i < 10; ++i) {
        this->policy.add(this->container, create_from_int<ElementType>(i));
    }

    // You might need to adjust the expected size depending on the implementation details and the types used.
    size_t expected_size = this->container.size() * sizeof(ElementType);
    ASSERT_EQ(this->policy.size_in_bytes(this->container), expected_size);
}

TYPED_TEST(SequenceContainerPolicyTest, EmptyContainer) {
    ASSERT_EQ(this->container.size(), 0);
    ASSERT_EQ(this->policy.size_in_bytes(this->container), 0);
}


//// Test Fixture for Set Container Policy
//template<typename SetType>
//class SetContainerPolicyTest : public ::testing::Test {
//protected:
//    SetType container;
//    ContainerPolicy<SetType> policy;
//};


//
//// Specify the types to be used for the tests
//using SetTypes = ::testing::Types<std::set<int>, std::set<std::string>>;
//TYPED_TEST_SUITE(SetContainerPolicyTest, SetTypes);
//
//TYPED_TEST(SetContainerPolicyTest, AddAndClear) {
//    using Key = typename TypeParam::key_type;
//
//    Key key1 = Key();
//    this->policy.add(this->container, key1);
//
//    ASSERT_EQ(this->container.size(), 1);
//    ASSERT_EQ(*this->container.begin(), key1);
//
//    this->policy.clear(this->container);
//    ASSERT_EQ(this->container.size(), 0);
//}