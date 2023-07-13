#ifndef ODIN_TREE_BASE_HPP
#define ODIN_TREE_BASE_HPP

#include <memory>
#include <oneapi/tbb/task_group.h>
#include <vector>

template<typename T>
class RawNode;

template<typename T>
class SafeNode;

template<typename Derived, typename T, template<typename, typename...> class ParentPtrType, template<typename, typename...> class ChildPtrType>
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

    T                  m_data;
    p_parent           m_parent;
    children_type      m_children;
    mutable std::mutex m_mutex;

    explicit NodeBase(const T &data) : m_data(data) {}

    std::shared_ptr<SafeNode<T>> to_safe() {
        return static_cast<Derived *>(this)->to_safe_impl();
    }

    std::shared_ptr<SafeNode<T>> to_safe() const {
        return static_cast<const Derived *>(this)->to_safe_impl();
    }

    std::unique_ptr<RawNode<T>> to_raw() {
        return static_cast<Derived *>(this)->to_raw_impl();
    }

    void add_child(p_child child) {
        static_cast<Derived *>(this)->add_child_impl(std::move(child));
    }

    [[nodiscard]] bool has_parent() const {
        return static_cast<const Derived *>(this)->has_parent_impl();
    }

    void set_parent(p_parent parent) {
        static_cast<Derived *>(this)->set_parent_impl(std::move(parent));
    }

    const T &get_data() const {
        return m_data;
    }

    p_parent get_parent() const {
        return m_parent;
    }

    const children_type &get_children() const {
        return m_children;
    }

    void remove_child(p_child child) {
        static_cast<Derived *>(this)->remove_child_impl(std::move(child));
    }

    void add_child_mt(p_child child) {
        std::lock_guard<std::mutex> lock(m_mutex);
        add_child(std::move(child));
    }

    void remove_child_mt(p_child child) {
        std::lock_guard<std::mutex> lock(m_mutex);
        remove_child(std::move(child));
    }

    [[nodiscard]] size_t get_node_count() const {
        size_t count = 1;// This node itself
        for (const auto &child: m_children)
            count += child->get_node_count();
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

        size_t count = 1;// Count this node because it is non-leaf
        for (const auto &child: m_children)
            count += child->get_non_leaf_node_count();

        return count;
    }
};

template<typename T>
class RawNode : public NodeBase<RawNode<T>, T, RawNodeTrait<T>> {
public:
    using base_type = NodeBase<RawNode<T>, T, RawNodeTrait<T>>;
    using base_type::base_type;

    explicit RawNode(const T &data) : base_type(data) {
        this->m_parent = nullptr;
    }

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

    void add_child_impl(base_type::p_child child) {
        child->m_parent = this;
        this->m_children.emplace_back(std::move(child));
    }

    [[nodiscard]] bool has_parent_impl() const {
        return this->m_parent != nullptr;
    }

    void set_parent_impl(base_type::p_parent parent) {
        this->m_parent = parent;
    }

    void remove_child(RawNode<int> *child) {
        auto it = std::find_if(this->m_children.begin(), this->m_children.end(),
                               [child](const auto &unique_ptr) { return unique_ptr.get() == child; });

        if (it != this->m_children.end()) {
            (*it)->m_parent = nullptr;
            this->m_children.erase(it);
        }
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

    void add_child_impl(base_type::p_child child) {
        child->m_parent = this->weak_from_this();
        this->m_children.push_back(child);
    }

    [[nodiscard]] bool has_parent_impl() const {
        return !this->m_parent.expired();
    }

    void set_parent_impl(base_type::p_parent parent) {
        this->m_parent = parent;
    }

    void remove_child_impl(base_type::p_child child) {
        auto it = std::find(this->m_children.begin(), this->m_children.end(), child);
        if (it != this->m_children.end()) {
            // remove parent
            (*it)->m_parent.reset();
            this->m_children.erase(it);
        }
    }
};


template<typename T>
class Tree {
public:
    using raw_node_type  = RawNode<T>;
    using safe_node_type = SafeNode<T>;
    using p_raw_node     = std::unique_ptr<raw_node_type>;
    using p_safe_node    = std::shared_ptr<safe_node_type>;

    explicit Tree(const raw_node_type &node) : m_root(std::make_unique<raw_node_type>(node)) {}
    explicit Tree(const safe_node_type &node) : m_root(node.to_raw()) {}

    raw_node_type *get_root() const {
        return m_root.get();
    }

    p_safe_node get_safe_root() const {
        return m_root->to_safe();
    }

    void parallel_in_order_traversal(std::function<void(const T &)> visit) {
        tbb::task_group group;

        auto worker = [&visit](p_raw_node node) {
            if (node) {
                node->parallel_in_order_traversal(visit);
            }
        };

        for (auto &child: m_root->get_children()) {
            group.run([=]() { worker(std::move(child)); });
        }

        group.wait();
    }

    void add_child_mt(raw_node_type *parent, p_raw_node child) {
        std::lock_guard<std::mutex> lock(m_tree_mutex);
        parent->add_child(std::move(child));
    }

    void remove_child_mt(raw_node_type *parent, p_raw_node child) {
        std::lock_guard<std::mutex> lock(m_tree_mutex);
        parent->remove_child(std::move(child));
    }

private:
    p_raw_node         m_root;
    mutable std::mutex m_tree_mutex;
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

#endif//ODIN_TREE_BASE_H