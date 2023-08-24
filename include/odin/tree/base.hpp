#ifndef ODIN_TREE_BASE_H
#define ODIN_TREE_BASE_H

#include <memory>
#include <vector>

#include <oneapi/tbb/task_group.h>

namespace odin {

template<typename T>
class RawNode;

template<typename T>
class SafeNode;

template<typename Derived,
         typename T,
         template<typename, typename...>
         class ParentPtrType,
         template<typename, typename...>
         class ChildPtrType>
struct NodeTrait {
    using p_parent = ParentPtrType<Derived>;
    using p_child  = ChildPtrType<Derived>;
};

template<typename T>
using RawNodeTrait = NodeTrait<RawNode<T>, T, std::add_pointer_t, std::unique_ptr>;

template<typename T>
using SafeNodeTrait = NodeTrait<SafeNode<T>, T, std::weak_ptr, std::shared_ptr>;


template<typename Derived, typename T, typename Trait, typename Float = double>
class NodeBase {
public:
    using p_parent      = typename Trait::p_parent;
    using p_child       = typename Trait::p_child;
    using children_type = std::vector<p_child>;


    explicit NodeBase(const T &data) : m_data(data) {}

    void update_levels(Derived *node, int level) {
        node->m_level = level;
        for (auto &child: node->get_children()) {
            update_levels(child, level + 1);
        }
    }

    Derived &as_derived() {
        return static_cast<Derived &>(*this);
    }

    const Derived &as_derived() const {
        return static_cast<const Derived &>(*this);
    }

    std::shared_ptr<SafeNode<T>> to_safe() {
        return as_derived().to_safe_impl();
    }

    std::shared_ptr<SafeNode<T>> to_safe() const {
        return as_derived().to_safe_impl();
    }

    std::unique_ptr<RawNode<T>> to_raw() {
        return as_derived().to_raw_impl();
    }

    void add_child(p_child child) {
        child->m_level = m_level + 1;
        as_derived().add_child_impl(std::move(child));
    }

    [[nodiscard]] bool has_parent() const {
        return as_derived().has_parent_impl();
    }

    void set_parent(p_parent parent) {
        as_derived().set_parent_impl(std::move(parent));
        m_level = as_derived().get_parent_level() + 1;
    }

    [[nodiscard]] size_t get_parent_level() const {
        return as_derived().get_parent_level_impl();
    }

    [[nodiscard]] const T &get_data() const {
        return m_data;
    }

    [[nodiscard]] p_parent get_parent() const {
        return m_parent;
    }

    [[nodiscard]] size_t get_level() const {
        return m_level;
    }

    [[nodiscard]] const children_type &get_children() const {
        return m_children;
    }

    p_child remove_child(Derived *child) {
        return as_derived().remove_child_impl(child);
    }

    p_child clone() const {
        return as_derived().clone_impl();
    }

    [[nodiscard]] size_t get_node_count() const {
        size_t count = 1;
        for (const auto &child: m_children) count += child->get_node_count();
        return count;
    }

    [[nodiscard]] size_t get_max_depth() const {
        size_t max_depth = 0;
        for (const auto &child: m_children)
            max_depth = std::max(max_depth, child->get_max_depth());
        return (m_children.empty()) ? 1 : max_depth + 1;
    }

    [[nodiscard]] Float get_average_branch_factor() const {
        size_t non_root_nodes = get_node_count() - 1;
        size_t non_leaf_nodes = get_non_leaf_node_count();
        if (non_leaf_nodes == 0) return 0;
        return static_cast<double>(non_root_nodes) / static_cast<double>(non_leaf_nodes);
    }

    [[nodiscard]] size_t get_non_leaf_node_count() const {
        if (m_children.empty()) return 0;
        size_t count = 1;
        for (const auto &child: m_children) count += child->get_non_leaf_node_count();
        return count;
    }

//    [[nodiscard]] p_parent get_siblings() {
//        return as_derived().get_siblings_impl();
//    }
//
//    [[nodiscard]] p_parent get_siblings() const {
//        return as_derived().get_siblings_impl();
//    }

    [[nodiscard]] bool is_leaf() const {
        return m_children.empty();
    }

    [[nodiscard]] bool is_root() const {
        return !has_parent();
    }

//    [[nodiscard]] bool is_valid() const {
//        return as_derived().is_valid_impl();
//    }
//
//    [[nodiscard]] bool is_valid() {
//        return as_derived().is_valid_impl();
//    }

    [[nodiscard]] T &get_data_ref() {
        return m_data;
    }

    T *get_data_ptr() {
        return &m_data;
    }


protected:
    size_t        m_level = 0;// Add a level property
    T             m_data;
    p_parent      m_parent;
    children_type m_children;
};

template<typename T>
class RawNode : public NodeBase<RawNode<T>, T, RawNodeTrait<T>> {
public:
    using base_type = NodeBase<RawNode<T>, T, RawNodeTrait<T>>;
    using base_type::base_type;

    explicit RawNode(const T &data) : base_type(data) { this->m_parent = nullptr; }

    std::shared_ptr<SafeNode<T>> to_safe_impl() const {
        auto safe_node = std::make_shared<SafeNode<T>>(this->m_data);
        for (const auto &raw_child: this->m_children) {
            auto safe_child = raw_child->to_safe();
            safe_node->add_child(safe_child);
        }
        return safe_node;
    }

    std::unique_ptr<RawNode<T>> to_raw_impl() {
        // RawNode does not need to convert to RawNode, return a copy of itself
        auto raw_copy = std::make_unique<RawNode<T>>(this->m_data);
        for (const auto &child: this->m_children) {
            auto child_copy = child->to_raw();
            raw_copy->add_child(std::move(child_copy));
        }
        return raw_copy;
    }

    void add_child_impl(typename base_type::p_child child) {
        child->m_parent = this;
        this->m_children.emplace_back(std::move(child));
    }

    [[nodiscard]] size_t get_parent_level_impl() const {
        if (this->get_parent()) {
            return this->get_parent()->m_level;
        } else {
            return 0;// or whatever value makes sense for no parent
        }
    }

    [[nodiscard]] bool has_parent_impl() const { return this->m_parent != nullptr; }

    void set_parent_impl(typename base_type::p_parent parent) { this->m_parent = parent; }

    typename base_type::p_child remove_child_impl(RawNode<T> *child) {
        auto it = std::find_if(
                this->m_children.begin(),
                this->m_children.end(),
                [child](const auto &unique_ptr) { return unique_ptr.get() == child; });

        if (it != this->m_children.end()) {
            // Unset parent from child node
            (*it)->m_parent = nullptr;
            // Save child to return to caller
            auto removed_child = std::move(*it);
            // Remove child from this node's children
            this->m_children.erase(it);
            // Return child to caller
            return removed_child;
        }
        // Return nullptr if child was not found
        return nullptr;
    }

    std::unique_ptr<RawNode<T>> clone_impl() const {
        auto cloned_node = std::make_unique<RawNode<T>>(this->m_data);
        for (const auto &child: this->m_children) cloned_node->add_child(child->clone());
        return cloned_node;
    }
};

template<typename T>
class SafeNode : public NodeBase<SafeNode<T>, T, SafeNodeTrait<T>>,
                 public std::enable_shared_from_this<SafeNode<T>> {
public:
    using base_type = NodeBase<SafeNode<T>, T, SafeNodeTrait<T>>;

    explicit SafeNode(const T &data) : base_type(data) {}

    std::shared_ptr<SafeNode<T>> to_safe_impl() {
        // SafeNode does not need to convert to SafeNode, return a copy of itself
        auto safe_copy = std::make_shared<SafeNode<T>>(this->m_data);
        for (const auto &child: this->m_children) {
            auto child_copy = child->to_safe();
            safe_copy->add_child(child_copy);
        }
        return safe_copy;
    }

    std::unique_ptr<RawNode<T>> to_raw_impl() {
        auto raw_node = std::make_unique<RawNode<T>>(this->m_data);
        for (const auto &safe_child: this->m_children) {
            auto raw_child = safe_child->to_raw();
            raw_node->add_child(std::move(raw_child));
        }
        return raw_node;
    }

    void add_child_impl(typename base_type::p_child child) {
        child->m_parent = this->weak_from_this();
        this->m_children.push_back(child);
    }

    [[nodiscard]] size_t get_parent_level_impl() const {
        if (auto parent = this->m_parent.lock()) {
            return parent->get_level();
        } else {
            return 0;// or whatever value makes sense for no parent
        }
    }

    [[nodiscard]] bool has_parent_impl() const { return !this->m_parent.expired(); }

    void set_parent_impl(typename base_type::p_parent parent) { this->m_parent = parent; }

    typename base_type::p_child remove_child_impl(SafeNode<T> *child) {
        auto it = std::find_if(
                this->m_children.begin(),
                this->m_children.end(),
                [child](const auto &shared_ptr) { return shared_ptr.get() == child; });

        if (it != this->m_children.end()) {
            // Unset parent from child node
            (*it)->m_parent.reset();
            // Save child to return to caller
            auto removed_child = std::move(*it);
            // Remove child from this node's children
            this->m_children.erase(it);
            // Return child to caller
            return removed_child;
        }
        // Return nullptr if child was not found
        return nullptr;
    }

    std::shared_ptr<SafeNode<T>> clone_impl() const {
        auto cloned_node = std::make_shared<SafeNode<T>>(this->m_data);
        for (const auto &child: this->m_children) cloned_node->add_child(child->clone());
        return cloned_node;
    }
};
//
//template<typename T>
//class NodeView {
//private:
//    RawNode<T> *m_node;
//
//    class NodeMediator {
//    };
//
//public:
//    [[nodiscard]] bool is_valid() const {
//        return m_node != nullptr;
//    }
//
//    [[nodiscard]] bool is_valid_else_throw() const {
//        if (!is_valid()) {
//            throw std::runtime_error("NodeView is not valid");
//        }
//        return true;
//    }
//
//    bool is_valid_then(NodeContainer container) const {
//        if (!is_valid()) {
//            return false;
//        }
//        return true;
//    }
//};


template<typename T>
class Tree {
public:
    // TODO: Implement a "NodeView" concept that allows for efficient traversal of the tree,
    //   without the need to copy the entire tree.
    // Aliases for easier reference
    using raw_node_type  = RawNode<T>;
    using safe_node_type = SafeNode<T>;
    using p_raw_node     = std::unique_ptr<raw_node_type>;
    using p_safe_node    = std::shared_ptr<safe_node_type>;

    // Construct tree with unique pointer to raw node
    explicit Tree(p_raw_node node) : m_root(std::move(node)) {}

    // Construct tree with shared pointer to safe node
    explicit Tree(const p_safe_node &node) : Tree(node->to_raw()) {}

    // Get root of tree as raw pointer. User should not delete this pointer.
    raw_node_type *get_root() const { return m_root.get(); }

    // Get root of tree as safe (shared) pointer
    p_safe_node get_safe_root() const { return m_root->to_safe(); }

    void add_child(raw_node_type *parent, p_raw_node child) {
        parent->add_child(std::move(child));
    }

    void add_child(p_safe_node parent, p_safe_node child) {
        parent->add_child(std::move(child));
    }

    p_raw_node remove_child(raw_node_type *parent, raw_node_type *child) {
        return parent->remove_child(child);
    }

    p_safe_node remove_child(p_safe_node parent, p_safe_node child) {
        return parent->remove_child(child);
    }

    // Add a child to a parent node. Thread-safe.
    // Transfers ownership of child from caller to Tree.
    void add_child_mt(raw_node_type *parent, p_raw_node child) {
        std::lock_guard<std::mutex> lock(m_tree_mutex);
        parent->add_child(std::move(child));
    }

    // Remove a child from a parent node. Thread-safe.
    // Transfers ownership of child from Tree to caller.
    p_raw_node remove_child_mt(raw_node_type *parent, raw_node_type *child) {
        std::lock_guard<std::mutex> lock(m_tree_mutex);
        return parent->remove_child(child);
    }

private:
    p_raw_node         m_root;      // Root of the tree
    mutable std::mutex m_tree_mutex;// Mutex to protect tree from simultaneous modifications
};


//
//template<typename T, template<typename> class Trait>
//class MCTSNode : public NodeBase<MCTSNode<T, Trait>, T, Trait<T>> {
//    // MCTS specific functions
//};
//
//using SafeMCTSNode = MCTSNode<int, SafeNodeTrait>;
//using RawMCTSNode  = MCTSNode<int, RawNodeTrait>;
//
//template<typename Derived, typename Node>
//struct TreeBase {
//
//    // Root node of the tree
//    std::unique_ptr<Node> m_root;
//
//    explicit TreeBase(const Node &node) {
//        m_root = std::make_unique<Node>(node);
//    }
//
//    TreeBase(const TreeBase &other) {
//        this->root = std::make_unique<Node>(*other.root);
//    }
//
//    TreeBase &operator=(const TreeBase &other) {
//    }
//
//    friend std::ostream &operator<<(std::ostream &os, const TreeBase &tree) {
//        os << "Tree (root): \n";
//        print_node(os, *(tree.root), "", true);
//        return os;
//    }
//
//    // Serialization function definition
//    template<class Archive>
//    void serialize(Archive &ar) {
//        ar(cereal::make_nvp("root", convert_node_to_py_node(*root)));
//
//        // Call method to serialize derived parts.
//        // This will call the appropriate method for the actual type of the object.
//        static_cast<Derived *>(this)->serialize_derived(ar);
//    }
//
//    template<class Archive>
//    static void
//    load_and_construct(Archive &ar, cereal::construct<Derived> &construct) {
//        std::shared_ptr<Node> root;
//
//        // Call deserialize on all members
//        ar(CEREAL_NVP(root));
//
//        // Call construct on the PyNode
//        construct(convert_py_node_to_node(*root));
//
//        // Call method to serialize derived parts.
//        // This will call the appropriate method for the actual type of the object.
//        //        Derived::load_and_construct_derived(ar);
//    }
//
//private:
//    static void print_node(std::ostream &os, const node_type &node, const std::string &prefix, bool is_tail) {
//        os << prefix << (is_tail ? "└─ " : "├─ ") << node << '\n';
//        const auto &children = node.children;
//        for (size_t i = 0; i < children.size(); ++i) {
//            print_node(os, *children[i], prefix + (is_tail ? "   " : "│  "), i == children.size() - 1);
//        }
//    }
//};

}// namespace odin
#endif// ODIN_TREE_BASE_H