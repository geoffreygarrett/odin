#include <gtest/gtest.h>
#include <memory>
#include <odin/io.hpp>
#include <odin/tree/base.hpp>
#include <odin/tree/serialization.hpp>

using namespace odin;

TEST(NodeTest, BasicTest) {
    auto raw_node = std::make_unique<RawNode<int>>(10);
    EXPECT_EQ(raw_node->get_data(), 10);
    EXPECT_EQ(raw_node->get_parent(), nullptr);
    EXPECT_EQ(raw_node->get_children().size(), 0);

    auto safe_node = std::make_shared<SafeNode<int>>(20);
    EXPECT_EQ(safe_node->get_data(), 20);
    EXPECT_EQ(safe_node->get_parent().lock(), nullptr);
    EXPECT_EQ(safe_node->get_children().size(), 0);
}

TEST(NodeTest, ConversionTest) {
    auto raw_node  = std::make_unique<RawNode<int>>(10);
    auto safe_node = raw_node->to_safe();
    EXPECT_EQ(safe_node->get_data(), 10);
    EXPECT_EQ(safe_node->get_parent().lock(), nullptr);
    EXPECT_EQ(safe_node->get_children().size(), 0);

    auto raw_node_back = safe_node->to_raw();
    EXPECT_EQ(raw_node_back->get_data(), 10);
    EXPECT_EQ(raw_node_back->get_parent(), nullptr);
    EXPECT_EQ(raw_node_back->get_children().size(), 0);
}

TEST(NodeTest, ConversionWithChildrenTest) {
    auto raw_node  = std::make_unique<RawNode<int>>(10);
    auto raw_child = std::make_unique<RawNode<int>>(5);
    raw_node->add_child(std::move(raw_child));// Using std::move here

    auto safe_node = raw_node->to_safe();
    EXPECT_EQ(safe_node->get_data(), 10);
    EXPECT_EQ(safe_node->get_parent().lock(), nullptr);
    EXPECT_EQ(safe_node->get_children().size(), 1);
    EXPECT_EQ(safe_node->get_children()[0]->get_data(), 5);

    auto safe_child = safe_node->get_children()[0];
    EXPECT_EQ(safe_child->get_data(), 5);
    EXPECT_EQ(safe_child->get_parent().lock(), safe_node);
    EXPECT_EQ(safe_child->get_children().size(), 0);

    auto raw_node_back = safe_node->to_raw();
    EXPECT_EQ(raw_node_back->get_data(), 10);
    EXPECT_EQ(raw_node_back->get_parent(), nullptr);
    EXPECT_EQ(raw_node_back->get_children().size(), 1);
    EXPECT_EQ(raw_node_back->get_children()[0]->get_data(), 5);

    // Access the child directly from the vector without moving it
    EXPECT_EQ(raw_node_back->get_children()[0]->get_data(), 5);
    EXPECT_EQ(raw_node_back->get_children()[0]->get_parent(), raw_node_back.get());
    EXPECT_EQ(raw_node_back->get_children()[0]->get_children().size(), 0);
}

TEST(NodeTest, ParentChildRelationshipTest) {
    auto parent_node = std::make_shared<SafeNode<int>>(10);
    auto child_node  = std::make_shared<SafeNode<int>>(5);

    parent_node->add_child(child_node);

    EXPECT_EQ(parent_node->get_children().size(), 1);
    EXPECT_EQ(parent_node->get_children()[0], child_node);
    EXPECT_EQ(child_node->get_parent().lock(), parent_node);
}

TEST(NodeTest, ToRawConversionWithParentChildRelationshipTest) {
    auto parent_node = std::make_shared<SafeNode<int>>(10);
    auto child_node  = std::make_shared<SafeNode<int>>(5);
    parent_node->add_child(child_node);

    auto raw_parent = parent_node->to_raw();
    EXPECT_EQ(raw_parent->get_children().size(), 1);
    EXPECT_EQ(raw_parent->get_children()[0]->get_data(), 5);
    EXPECT_EQ(raw_parent->get_children()[0]->get_parent(), raw_parent.get());
}

TEST(NodeTest, SafeNodeCopyTest) {
    auto safe_node      = std::make_shared<SafeNode<int>>(10);
    auto safe_node_copy = safe_node->to_safe();

    EXPECT_NE(safe_node, safe_node_copy);
    EXPECT_EQ(safe_node->get_data(), safe_node_copy->get_data());
    EXPECT_EQ(safe_node->get_children().size(), safe_node_copy->get_children().size());
}

TEST(NodeTest, RawNodeCopyTest) {
    auto raw_node      = std::make_unique<RawNode<int>>(10);
    auto raw_node_copy = raw_node->to_raw();

    EXPECT_NE(raw_node.get(), raw_node_copy.get());
    EXPECT_EQ(raw_node->get_data(), raw_node_copy->get_data());
    EXPECT_EQ(raw_node->get_children().size(), raw_node_copy->get_children().size());
}

TEST(NodeTest, MultipleChildrenTest) {
    auto parent_node = std::make_shared<SafeNode<int>>(10);
    auto child_node1 = std::make_shared<SafeNode<int>>(5);
    auto child_node2 = std::make_shared<SafeNode<int>>(15);

    parent_node->add_child(child_node1);
    parent_node->add_child(child_node2);

    EXPECT_EQ(parent_node->get_children().size(), 2);
    EXPECT_EQ(parent_node->get_children()[0], child_node1);
    EXPECT_EQ(parent_node->get_children()[1], child_node2);
    EXPECT_EQ(child_node1->get_parent().lock(), parent_node);
    EXPECT_EQ(child_node2->get_parent().lock(), parent_node);
}

TEST(NodeTest, MultilevelTreeConversionTest) {
    auto parent_node     = std::make_shared<SafeNode<int>>(10);
    auto child_node1     = std::make_shared<SafeNode<int>>(5);
    auto child_node2     = std::make_shared<SafeNode<int>>(15);
    auto grandchild_node = std::make_shared<SafeNode<int>>(20);

    child_node2->add_child(grandchild_node);
    parent_node->add_child(child_node1);
    parent_node->add_child(child_node2);

    auto raw_node = parent_node->to_raw();

    EXPECT_EQ(raw_node->get_children().size(), 2);
    EXPECT_EQ(raw_node->get_children()[0]->get_data(), 5);
    EXPECT_EQ(raw_node->get_children()[1]->get_data(), 15);
    EXPECT_EQ(raw_node->get_children()[1]->get_children()[0]->get_data(), 20);
    EXPECT_EQ(raw_node->get_children()[0]->get_parent(), raw_node.get());
    EXPECT_EQ(raw_node->get_children()[1]->get_parent(), raw_node.get());
    EXPECT_EQ(raw_node->get_children()[1]->get_children()[0]->get_parent(), raw_node->get_children()[1].get());
}

TEST(NodeTest, StringNodeTest) {
    auto raw_node = std::make_unique<RawNode<std::string>>("hello");
    EXPECT_EQ(raw_node->get_data(), "hello");
    EXPECT_EQ(raw_node->get_parent(), nullptr);
    EXPECT_EQ(raw_node->get_children().size(), 0);

    auto safe_node = std::make_shared<SafeNode<std::string>>("world");
    EXPECT_EQ(safe_node->get_data(), "world");
    EXPECT_EQ(safe_node->get_parent().lock(), nullptr);
    EXPECT_EQ(safe_node->get_children().size(), 0);
}

TEST(NodeTest, RawNodeHasParentTest) {
    auto parent = std::make_unique<RawNode<int>>(10);
    auto child  = std::make_unique<RawNode<int>>(5);

    EXPECT_FALSE(parent->has_parent());
    EXPECT_FALSE(child->has_parent());

    parent->add_child(std::move(child));

    EXPECT_FALSE(parent->has_parent());
    EXPECT_TRUE(parent->get_children()[0]->has_parent());
}

TEST(NodeTest, SafeNodeHasParentTest) {
    auto parent = std::make_shared<SafeNode<int>>(10);
    auto child  = std::make_shared<SafeNode<int>>(5);

    EXPECT_FALSE(parent->has_parent());
    EXPECT_FALSE(child->has_parent());

    parent->add_child(child);

    EXPECT_FALSE(parent->has_parent());
    EXPECT_TRUE(child->has_parent());
}


TEST(NodeTest, SafeNodeSerializationTest) {
    auto node   = std::make_shared<SafeNode<int>>(10);
    auto child1 = std::make_shared<SafeNode<int>>(5);
    auto child2 = std::make_shared<SafeNode<int>>(15);

    node->add_child(child1);
    node->add_child(child2);

    {// Start of scope for output archive
        std::ofstream             ofs("safe_node.json");
        cereal::JSONOutputArchive oarchive(ofs);
        oarchive(node);
    }// Output archive will be flushed here

    std::shared_ptr<SafeNode<int>> deserialized_node;

    {// Start of scope for input archive
        std::ifstream            ifs("safe_node.json");
        cereal::JSONInputArchive iarchive(ifs);
        iarchive(deserialized_node);
    }// Input archive will be finished here

    EXPECT_EQ(node->get_data(), deserialized_node->get_data());
    EXPECT_EQ(node->get_children().size(), deserialized_node->get_children().size());
    for (size_t i = 0; i < node->get_children().size(); i++) {
        EXPECT_EQ(node->get_children()[i]->get_data(), deserialized_node->get_children()[i]->get_data());
    }
}
TEST(NodeTest, RawNodeSerializationTest) {
    auto node   = std::make_unique<RawNode<int>>(10);
    auto child1 = std::make_unique<RawNode<int>>(5);
    auto child2 = std::make_unique<RawNode<int>>(15);

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));

    {// Start of scope for output archive
        std::ofstream             ofs("raw_node.json");
        cereal::JSONOutputArchive oarchive(ofs);
        oarchive(node);
    }// Output archive will be flushed here

    std::unique_ptr<RawNode<int>> deserialized_node;

    {// Start of scope for input archive
        std::ifstream            ifs("raw_node.json");
        cereal::JSONInputArchive iarchive(ifs);
        iarchive(deserialized_node);
    }// Input archive will be finished here

    EXPECT_EQ(node->get_data(), deserialized_node->get_data());                      // use object, not pointer
    EXPECT_EQ(node->get_children().size(), deserialized_node->get_children().size());// use object, not pointer
    for (size_t i = 0; i < node->get_children().size(); i++) {
        EXPECT_EQ(node->get_children()[i]->get_data(), deserialized_node->get_children()[i]->get_data());// use object, not pointer
    }
}

TEST(NodeTest, SafeNodeChildrenOwnershipTest) {
    auto node   = std::make_shared<SafeNode<int>>(10);
    auto child1 = std::make_shared<SafeNode<int>>(5);
    auto child2 = std::make_shared<SafeNode<int>>(15);

    node->add_child(child1);
    node->add_child(child2);

    EXPECT_EQ(node->get_children().size(), 2);

    // Check that child nodes are owned by the parent
    EXPECT_EQ(child1.use_count(), 2);
    EXPECT_EQ(child2.use_count(), 2);
}

TEST(NodeTest, SafeNodeRemoveChildTest) {
    auto node   = std::make_shared<SafeNode<int>>(10);
    auto child1 = std::make_shared<SafeNode<int>>(5);
    auto child2 = std::make_shared<SafeNode<int>>(15);

    node->add_child(child1);
    node->add_child(child2);

    EXPECT_EQ(node->get_children().size(), 2);

    node->remove_child(child1);

    // Check that the child node is correctly removed
    EXPECT_EQ(node->get_children().size(), 1);
    EXPECT_FALSE(child1->has_parent());
    EXPECT_TRUE(child2->has_parent());
}
//
TEST(NodeTest, RawNodeRemoveChildTest) {
    auto node   = std::make_unique<RawNode<int>>(10);
    auto child1 = std::make_unique<RawNode<int>>(5);
    auto child2 = std::make_unique<RawNode<int>>(15);

    auto child1_raw = child1.get();

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));

    EXPECT_EQ(node->get_children().size(), 2);

    node->remove_child(child1_raw);

    // Check that the child node is correctly removed
    EXPECT_EQ(node->get_children().size(), 1);
}


TEST(NodeTest, SafeNodeJSONSerializationTest) {
    auto node   = std::make_shared<SafeNode<int>>(10);
    auto child1 = std::make_shared<SafeNode<int>>(5);
    auto child2 = std::make_shared<SafeNode<int>>(15);

    node->add_child(child1);
    node->add_child(child2);

    save_to_json(node, "safe_node.json");

    auto deserialized_node = load_from_json<std::shared_ptr<SafeNode<int>>>("safe_node.json");

    EXPECT_EQ(node->get_data(), deserialized_node->get_data());
    EXPECT_EQ(node->get_children().size(), deserialized_node->get_children().size());

    for (size_t i = 0; i < node->get_children().size(); i++) {
        EXPECT_EQ(node->get_children()[i]->get_data(), deserialized_node->get_children()[i]->get_data());
    }
}

TEST(NodeTest, RawNodeBinarySerializationTest) {
    auto node   = std::make_unique<RawNode<int>>(10);
    auto child1 = std::make_unique<RawNode<int>>(5);
    auto child2 = std::make_unique<RawNode<int>>(15);

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));

    save_to_binary(node, "raw_node.bin");

    auto deserialized_node = load_from_binary<std::unique_ptr<RawNode<int>>>("raw_node.bin");

    EXPECT_EQ(node->get_data(), deserialized_node->get_data());
    EXPECT_EQ(node->get_children().size(), deserialized_node->get_children().size());

    for (size_t i = 0; i < node->get_children().size(); i++) {
        EXPECT_EQ(node->get_children()[i]->get_data(), deserialized_node->get_children()[i]->get_data());
    }
}

TEST(NodeTest, SafeNodeBaseFunctionsTest) {
    auto node   = std::make_shared<SafeNode<int>>(10);
    auto child1 = std::make_shared<SafeNode<int>>(5);
    auto child2 = std::make_shared<SafeNode<int>>(15);

    node->add_child(child1);
    node->add_child(child2);

    EXPECT_EQ(node->get_node_count(), 3);
    EXPECT_EQ(node->get_max_depth(), 2);
    EXPECT_EQ(node->get_average_branch_factor(), 2);
}

TEST(NodeTest, RawNodeBaseFunctionsTest) {
    auto node   = std::make_unique<RawNode<int>>(10);
    auto child1 = std::make_unique<RawNode<int>>(5);
    auto child2 = std::make_unique<RawNode<int>>(15);

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));

    EXPECT_EQ(node->get_node_count(), 3);
    EXPECT_EQ(node->get_max_depth(), 2);
    EXPECT_EQ(node->get_average_branch_factor(), 2);
}

// Simple tree with one node
TEST(NodeTest, SingleNodeBaseFunctionsTest) {
    auto node = std::make_unique<RawNode<int>>(10);

    EXPECT_EQ(node->get_node_count(), 1);
    EXPECT_EQ(node->get_max_depth(), 1);
    EXPECT_EQ(node->get_average_branch_factor(), 0);
}

// Test with two level tree and varying number of children at each node
TEST(NodeTest, TwoLevelTreeBaseFunctionsTest) {
    auto node   = std::make_unique<RawNode<int>>(10);
    auto child1 = std::make_unique<RawNode<int>>(5);
    auto child2 = std::make_unique<RawNode<int>>(15);
    auto child3 = std::make_unique<RawNode<int>>(20);

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));
    node->add_child(std::move(child3));

    EXPECT_EQ(node->get_node_count(), 4);
    EXPECT_EQ(node->get_max_depth(), 2);
    EXPECT_DOUBLE_EQ(node->get_average_branch_factor(), 3);
}

// Test with deeper tree and varying number of children at each node
TEST(NodeTest, DeepTreeBaseFunctionsTest) {
    auto node        = std::make_unique<RawNode<int>>(10);
    auto child1      = std::make_unique<RawNode<int>>(5);
    auto child2      = std::make_unique<RawNode<int>>(15);
    auto child3      = std::make_unique<RawNode<int>>(20);
    auto grandchild1 = std::make_unique<RawNode<int>>(25);
    auto grandchild2 = std::make_unique<RawNode<int>>(30);

    child1->add_child(std::move(grandchild1));
    child2->add_child(std::move(grandchild2));
    node->add_child(std::move(child1));
    node->add_child(std::move(child2));
    node->add_child(std::move(child3));

    EXPECT_EQ(node->get_node_count(), 6);
    EXPECT_EQ(node->get_max_depth(), 3);
    EXPECT_DOUBLE_EQ(node->get_average_branch_factor(), 1.6666666666666667);
}


// Repeat the tests for SafeNode
// (It might be better to template these tests in the future to avoid repetition)
TEST(NodeTest, SafeSingleNodeBaseFunctionsTest) {
    auto node = std::make_shared<SafeNode<int>>(10);

    EXPECT_EQ(node->get_node_count(), 1);
    EXPECT_EQ(node->get_max_depth(), 1);
    EXPECT_EQ(node->get_average_branch_factor(), 0);
}

TEST(NodeTest, SafeTwoLevelTreeBaseFunctionsTest) {
    auto node   = std::make_shared<SafeNode<int>>(10);
    auto child1 = std::make_shared<SafeNode<int>>(5);
    auto child2 = std::make_shared<SafeNode<int>>(15);
    auto child3 = std::make_shared<SafeNode<int>>(20);

    node->add_child(child1);
    node->add_child(child2);
    node->add_child(child3);

    EXPECT_EQ(node->get_node_count(), 4);
    EXPECT_EQ(node->get_max_depth(), 2);
    EXPECT_DOUBLE_EQ(node->get_average_branch_factor(), 3);
}

// Repeat the tests for SafeNode
// (It might be better to template these tests in the future to avoid repetition)
TEST(NodeTest, SafeDeepTreeBaseFunctionsTest) {
    auto node        = std::make_shared<SafeNode<int>>(10);
    auto child1      = std::make_shared<SafeNode<int>>(5);
    auto child2      = std::make_shared<SafeNode<int>>(15);
    auto child3      = std::make_shared<SafeNode<int>>(20);
    auto grandchild1 = std::make_shared<SafeNode<int>>(25);
    auto grandchild2 = std::make_shared<SafeNode<int>>(30);

    child1->add_child(grandchild1);
    child2->add_child(grandchild2);
    node->add_child(child1);
    node->add_child(child2);
    node->add_child(child3);

    EXPECT_EQ(node->get_node_count(), 6);
    EXPECT_EQ(node->get_max_depth(), 3);
    EXPECT_DOUBLE_EQ(node->get_average_branch_factor(), 1.6666666666666667);
}


template<template<typename> typename NodeType, typename DataType>
struct TypeParams {
    using Node    = NodeType<DataType>;
    using NodePtr = std::shared_ptr<Node>;
    using Data    = DataType;// Add this
};

template<typename T>
struct TypedNodeTest : public ::testing::Test {
    using Node    = typename T::Node;
    using NodePtr = typename T::NodePtr;
    using Data    = typename T::Data;// And this
};

struct TestStruct {
    int a;
    int b;

    template<typename Archive>
    void serialize(Archive &ar) {
        ar(a, b);
    }
};

template<typename T>
struct DataMaker {
    static T make_data() { return T{}; }
};

template<>
struct DataMaker<int> {
    static int make_data() { return 10; }
};

template<>
struct DataMaker<float> {
    static float make_data() { return 10.0f; }
};

template<>
struct DataMaker<TestStruct> {
    static TestStruct make_data() { return {10, 20}; }
};


using MyTypes = ::testing::Types<
        TypeParams<RawNode, int>,
        TypeParams<RawNode, float>,
        TypeParams<SafeNode, int>,
        TypeParams<SafeNode, float>,
        TypeParams<RawNode, TestStruct>,
        TypeParams<SafeNode, TestStruct>>;

TYPED_TEST_SUITE(TypedNodeTest, MyTypes);

TYPED_TEST(TypedNodeTest, SingleNodeBaseFunctionsTest) {
    auto                          data = DataMaker<typename TestFixture::Data>::make_data();
    typename TestFixture::NodePtr node = std::make_shared<typename TestFixture::Node>(data);

    EXPECT_EQ(node->get_node_count(), 1);
    EXPECT_EQ(node->get_max_depth(), 1);
    EXPECT_EQ(node->get_average_branch_factor(), 0);
}


TYPED_TEST(TypedNodeTest, TwoLevelTreeBaseFunctionsTest) {
    auto data = DataMaker<typename TestFixture::Data>::make_data();

    typename TestFixture::NodePtr node   = std::make_unique<typename TestFixture::Node>(data);
    auto                          child1 = std::make_unique<typename TestFixture::Node>(data);
    auto                          child2 = std::make_unique<typename TestFixture::Node>(data);
    auto                          child3 = std::make_unique<typename TestFixture::Node>(data);

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));
    node->add_child(std::move(child3));

    EXPECT_EQ(node->get_node_count(), 4);
    EXPECT_EQ(node->get_max_depth(), 2);
    EXPECT_DOUBLE_EQ(node->get_average_branch_factor(), 3.0);
}

TYPED_TEST(TypedNodeTest, DeepTreeBaseFunctionsTest) {
    auto data        = DataMaker<typename TestFixture::Data>::make_data();
    auto node        = std::make_unique<typename TestFixture::Node>(data);
    auto child1      = std::make_unique<typename TestFixture::Node>(data);
    auto child2      = std::make_unique<typename TestFixture::Node>(data);
    auto child3      = std::make_unique<typename TestFixture::Node>(data);
    auto grandchild1 = std::make_unique<typename TestFixture::Node>(data);
    auto grandchild2 = std::make_unique<typename TestFixture::Node>(data);

    child1->add_child(std::move(grandchild1));
    child2->add_child(std::move(grandchild2));
    node->add_child(std::move(child1));
    node->add_child(std::move(child2));
    node->add_child(std::move(child3));

    EXPECT_EQ(node->get_node_count(), 6);
    EXPECT_EQ(node->get_max_depth(), 3);
    EXPECT_DOUBLE_EQ(node->get_average_branch_factor(), 1.6666666666666667);
}

// The serialization test is not included in the template test because the serialization function
// may not be available for all types of data, especially if the data type is a user-defined structure.
TEST(NodeTest, SafeNodeXMLSerializationTest) {
    auto node   = std::make_unique<SafeNode<int>>(10);
    auto child1 = std::make_unique<SafeNode<int>>(5);
    auto child2 = std::make_unique<SafeNode<int>>(15);

    node->add_child(std::move(child1));
    node->add_child(std::move(child2));

    save_to_xml(node, "safe_node.xml");

    auto deserialized_node = load_from_xml<std::unique_ptr<SafeNode<int>>>("safe_node.xml");

    EXPECT_EQ(node->get_data(), deserialized_node->get_data());
    EXPECT_EQ(node->get_children().size(), deserialized_node->get_children().size());

    for (size_t i = 0; i < node->get_children().size(); i++) {
        EXPECT_EQ(node->get_children()[i]->get_data(), deserialized_node->get_children()[i]->get_data());
    }
}