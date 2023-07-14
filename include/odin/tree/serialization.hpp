#ifndef ODIN_TREE_SERIALIZATION_HPP
#define ODIN_TREE_SERIALIZATION_HPP

#include <odin/tree/base.hpp>

#include <fstream>
#include <type_traits>

#include <cereal/access.hpp>// For LoadAndConstruct
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

using namespace odin;

namespace cereal {

// Save a RawNode by converting it to SafeNode and then serializing
template<class Archive, typename T>
void save(Archive &archive, RawNode<T> const &node) {
    auto safe_node = node.to_safe();
    archive(cereal::make_nvp("data", safe_node->get_data()),
            cereal::make_nvp("children", safe_node->get_children()),
            cereal::make_nvp("parent", safe_node->get_parent()));
}

template<class Archive, typename T>
void serialize(Archive &archive, SafeNode<T> &node) {
    archive(cereal::make_nvp("data", node.get_data()),
            cereal::make_nvp("children", node.get_children()),
            cereal::make_nvp("parent", node.get_parent()));
}

// Specialize the LoadAndConstruct template for SafeNode
template<typename T>
struct LoadAndConstruct<SafeNode<T>> {
    template<typename Archive>
    static void load_and_construct(Archive &ar, cereal::construct<SafeNode<T>> &construct) {
        T                                         data;
        std::vector<std::shared_ptr<SafeNode<T>>> children;
        std::weak_ptr<SafeNode<T>>                parent;

        ar(cereal::make_nvp("data", data),
           cereal::make_nvp("children", children),
           cereal::make_nvp("parent", parent));

        construct(data);
        for (auto &child: children) { construct->add_child(child); }
        if (auto p = parent.lock()) { construct->set_parent(p); }
    }
};

// Specialize the LoadAndConstruct template for RawNode
template<typename T>
struct LoadAndConstruct<RawNode<T>> {
    template<typename Archive>
    static void load_and_construct(Archive &ar, cereal::construct<RawNode<T>> &construct) {
        T                                         data;
        std::vector<std::shared_ptr<SafeNode<T>>> children;
        std::weak_ptr<SafeNode<T>>                parent;

        ar(cereal::make_nvp("data", data),
           cereal::make_nvp("children", children),
           cereal::make_nvp("parent", parent));

        construct(data);
        for (auto &child: children) { construct->add_child(child->to_raw()); }
        if (auto p = parent.lock()) {
            auto raw_node_ptr = p->to_raw().release();
            construct->set_parent(raw_node_ptr);
        }
    }
};

// Serialization function for the cereal library
template<class Archive, typename T>
void serialize(Archive &archive, Tree<T> &tree) {
    archive(cereal::make_nvp("root", tree.get_root()));

    // check if Derived has a member called serialized_derived, if so call it (check this at compile time)
    //        if constexpr (has_serialized_derived<Archive>::value) {
    //            static_cast<Derived &>(tree).serialized_derived(archive);
    //        }
}

}// namespace cereal

#endif// ODIN_TREE_SERIALIZATION_HPP