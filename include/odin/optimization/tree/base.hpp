#ifndef ODIN_TREE_BASE_HPP
#define ODIN_TREE_BASE_HPP

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/variant.hpp>
#include <cereal/types/vector.hpp>
#include <optional>

#include <odin/core/policy/serialization/eigen_serialization_policy.hpp>

//#include "mcts_serialization.hpp"
#include <Eigen/Core>
#include <cereal/access.hpp>
#include <memory>
#include <vector>

// Default formatter
template<typename T, typename Enable = void>
struct formatter {
    static void print(std::ostream &os, const T &t) {
        os << t;
    }
};

// Specialized formatter for Eigen::Matrix
template<typename T>
struct formatter<T, typename std::enable_if<std::is_base_of<Eigen::MatrixBase<T>, T>::value>::type> {
    static void print(std::ostream &os, const T &t) {
        Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");
        os << t.format(CommaInitFmt);
    }
};

// Specialized formatter for std::tuple of std::vector<int> and std::vector<float>
template<>
struct formatter<std::tuple<std::vector<int>, std::vector<float>>> {
    static void print(std::ostream &os, const std::tuple<std::vector<int>, std::vector<float>> &t) {
        os << "[";
        for (const auto &val: std::get<0>(t)) {
            os << val << ", ";
        }
        os << "], [";
        for (const auto &val: std::get<1>(t)) {
            os << val << ", ";
        }
        os << "]";
    }
};

// Specialized formatter for std::vector
template<typename T>
struct formatter<std::vector<T>> {
    static void print(std::ostream &os, const std::vector<T> &t) {
        os << "[";
        for (size_t i = 0; i < t.size(); ++i) {
            formatter<T>::print(os, t[i]);
            if (i < t.size() - 1) {
                os << ", ";
            }
        }
        os << "]";
    }
};


template<typename Float = double>
class RunningStats {
public:
    explicit RunningStats(bool calc_mean = false, bool calc_variance = false, bool calc_min = false, bool calc_max = true)
        : m_n(0),
          m_calc_mean(calc_mean),
          m_calc_variance(calc_variance),
          m_calc_min(calc_min),
          m_calc_max(calc_max),
          m_min(calc_min ? std::numeric_limits<Float>::max() : 0),
          m_max(0),
          //          m_max(calc_max ? std::numeric_limits<Float>::lowest() : 0),
          m_old_mean(calc_mean || calc_variance ? 0.0 : 0),
          m_new_mean(calc_mean || calc_variance ? 0.0 : 0),
          m_old_sqr_diff(calc_variance ? 0.0 : 0), m_new_sqr_diff(calc_variance ? 0.0 : 0) {}

    friend std::ostream &operator<<(std::ostream &os, const RunningStats<Float> &stats) {
        os << "{n: " << stats.count();
        if (stats.m_calc_min) os << ", min: " << stats.min();
        if (stats.m_calc_max) os << ", max: " << stats.max();
        if (stats.m_calc_mean) os << ", x̄: " << stats.mean();
        if (stats.m_calc_variance) os << ", σ²: " << stats.variance();
        os << "}";
        return os;
    }

    void set_calc_mean(bool calc_mean) {
        m_calc_mean = calc_mean;
    }

    void set_calc_variance(bool calc_variance) {
        m_calc_variance = calc_variance;
    }

    void set_calc_min(bool calc_min) {
        m_calc_min = calc_min;
    }

    void set_calc_max(bool calc_max) {
        m_calc_max = calc_max;
    }

    [[nodiscard]] bool calc_mean() const {
        return m_calc_mean;
    }

    [[nodiscard]] bool calc_variance() const {
        return m_calc_variance;
    }

    [[nodiscard]] bool calc_min() const {
        return m_calc_min;
    }

    [[nodiscard]] bool calc_max() const {
        return m_calc_max;
    }

    void clear() {
        // TODO: Think this through more carefully. Idiomatic is the key.
        m_n = 0;
        if (m_calc_min) m_min = std::numeric_limits<Float>::max();
        if (m_calc_max) m_max = std::numeric_limits<Float>::min();
        if (m_calc_mean || m_calc_variance) {
            m_old_mean = 0;
            m_new_mean = 0;
        }
        if (m_calc_variance) {
            m_old_sqr_diff = 0;
            m_new_sqr_diff = 0;
        }
    }

    void push(Float x) {
        m_n++;

        // Update min and max
        if (m_calc_min && x < m_min) m_min = x;
        if (m_calc_max && x > m_max) m_max = x;

        if (m_calc_mean || m_calc_variance) {
            if (m_n == 1) {
                m_old_mean = m_new_mean = x;
                m_old_sqr_diff          = m_calc_variance ? 0.0 : 0;
            } else {
                m_new_mean = m_old_mean + (x - m_old_mean) / static_cast<Float>(m_n);
                if (m_calc_variance) {
                    m_new_sqr_diff = m_old_sqr_diff + (x - m_old_mean) * (x - m_new_mean);
                }

                // Prepare for next iteration
                m_old_mean     = m_new_mean;
                m_old_sqr_diff = m_new_sqr_diff;
            }
        }
    }

    [[nodiscard]] int count() const {
        return m_n;
    }

    Float mean() const {
        if (!m_calc_mean) {
            throw std::runtime_error("Mean calculation not enabled.");
        }
        return (m_n > 0) ? m_new_mean : 0.0;
    }

    Float variance() const {
        if (!m_calc_variance) {
            throw std::runtime_error("Variance calculation not enabled.");
        }
        return ((m_n > 1) ? m_new_sqr_diff / static_cast<Float>(m_n - 1) : 0.0);
    }

    Float standard_deviation() const {
        return std::sqrt(variance());
    }

    Float min() const {
        if (!m_calc_min) {
            throw std::runtime_error("Min calculation not enabled.");
        }
        return m_min;
    }

    Float max() const {
        if (!m_calc_max) {
            throw std::runtime_error("Max calculation not enabled.");
        }
        return m_max;
    }

private:
    int   m_n;
    bool  m_calc_mean;
    bool  m_calc_variance;
    bool  m_calc_min;
    bool  m_calc_max;
    Float m_min, m_max;
    Float m_old_mean, m_new_mean, m_old_sqr_diff, m_new_sqr_diff;

    // making the cereal namespace a friend
    friend class cereal::access;

    template<class Archive>
    void serialize(Archive &archive) {
        archive(cereal::make_nvp("calc_mean", m_calc_mean),
                cereal::make_nvp("calc_variance", m_calc_variance),
                cereal::make_nvp("calc_min", m_calc_min),
                cereal::make_nvp("calc_max", m_calc_max),
                cereal::make_nvp("n", m_n),
                cereal::make_nvp("min", m_min),
                cereal::make_nvp("max", m_max),
                cereal::make_nvp("old_mean", m_old_mean),
                cereal::make_nvp("new_mean", m_new_mean),
                cereal::make_nvp("old_sqr_diff", m_old_sqr_diff),
                cereal::make_nvp("new_sqr_diff", m_new_sqr_diff));
    }
};

template<typename Derived, typename State, typename Action, typename Reward>
struct NodeBase {
    State                 state;
    std::optional<Action> action;// Optional for root node.
    RunningStats<Reward>  reward_stats;
    size_t                visit_count{};
    size_t                update_count{};
    bool                  is_terminal{};
    bool                  is_leaf{};
    Derived              *parent;

    virtual ~NodeBase() = default;

    [[nodiscard]] bool is_leaf_node() const {
        return ((static_cast<const Derived &>(*this).children_size() == 0) || this->is_leaf);
    }

    NodeBase() = default;

    NodeBase(State                 state,
             std::optional<Action> action,
             Derived              *parent               = nullptr,
             RunningStats<Reward>  initial_reward_stats = RunningStats<Reward>())
        : state(state),
          action(action),
          visit_count(0),
          update_count(0),
          is_terminal(false),
          is_leaf(true),
          parent(parent),
          reward_stats(initial_reward_stats) {}

    NodeBase(const NodeBase &other)
        : state(other.state), action(other.action), reward_stats(other.reward_stats),
          visit_count(other.visit_count), update_count(other.update_count),
          is_terminal(other.is_terminal), is_leaf(other.is_leaf), parent(other.parent) {
    }

    friend std::ostream &operator<<(std::ostream &os, const NodeBase &node) {
        os << "Node(A: ";
        if (node.action) {// Check if action is present
            formatter<Action>::print(os, *node.action);
        } else {
            os << "None";
        }
        os << ", R: " << node.reward_stats;
        os << ", n: " << node.visit_count;
        //        os << ", update_count: " << node.update_count;
        os << ", is_terminal: " << std::boolalpha << node.is_terminal;
        os << ", is_leaf: " << std::boolalpha << node.is_leaf;
        os << ", child_count: " << static_cast<const Derived &>(node).children_size() << ")";
        return os;
    }

    // Serialization function definition
    template<class Archive>
    void serialize(Archive &ar) {
        ar(cereal::make_nvp("state", state),
           cereal::make_nvp("action", action),
           cereal::make_nvp("reward_stats", reward_stats),
           cereal::make_nvp("visit_count", visit_count),
           cereal::make_nvp("update_count", update_count),
           cereal::make_nvp("is_terminal", is_terminal),
           cereal::make_nvp("is_leaf", is_leaf));

        // Call method to serialize derived parts.
        // This will call the appropriate method for the actual type of the object.
        static_cast<Derived *>(this)->serialize_derived(ar);
    }

    template<class Archive>
    static void load_and_construct(Archive &ar, cereal::construct<Derived> &construct) {
        State                 state;
        std::optional<Action> action;
        RunningStats<Reward>  reward_stats;
        size_t                visit_count;
        size_t                update_count;
        bool                  is_terminal;
        bool                  is_leaf;

        // Call deserialize on all members
        ar(cereal::make_nvp("state", state),
           cereal::make_nvp("action", action),
           cereal::make_nvp("reward_stats", reward_stats),
           cereal::make_nvp("visit_count", visit_count),
           cereal::make_nvp("update_count", update_count),
           cereal::make_nvp("is_terminal", is_terminal),
           cereal::make_nvp("is_leaf", is_leaf));

        Derived::load_and_construct_derived(ar, construct, state, action, reward_stats);
    }
};

// Forward declarations
template<typename State, typename Action, typename Reward>
struct Node;

template<typename State, typename Action, typename Reward>
struct PyNode;

template<typename State, typename Action, typename Reward>
std::unique_ptr<Node<State, Action, Reward>> convert_py_node_to_node(const PyNode<State, Action, Reward> &py_node, Node<State, Action, Reward> *parent = nullptr);

template<typename State, typename Action, typename Reward>
std::shared_ptr<PyNode<State, Action, Reward>> convert_node_to_py_node(const Node<State, Action, Reward> &node, std::shared_ptr<PyNode<State, Action, Reward>> parent = nullptr);

template<typename State, typename Action, typename Reward>
struct Node : NodeBase<Node<State, Action, Reward>, State, Action, Reward> {
    using NodeBase<Node, State, Action, Reward>::NodeBase;
    std::vector<std::unique_ptr<Node<State, Action, Reward>>> children;
    Node() = default;

    Node(const Node &other)
        : NodeBase<Node, State, Action, Reward>(other) {
        for (const auto &child: other.children) {
            auto new_child    = std::make_unique<Node>(*child);
            new_child->parent = this;// update parent pointer
            this->children.push_back(std::move(new_child));
        }
    }

    [[nodiscard]] size_t children_size() const {
        return children.size();
    }

    explicit operator std::shared_ptr<PyNode<State, Action, Reward>>() const {
        return convert_node_to_py_node(*this);
    }

    //    // Serialization function definition
    //    template<class Archive>
    //    void serialize_derived(Archive &ar) {
    //        ar(cereal::make_nvp("children", children));
    //    }
    //
    //    template<class Archive>
    //    static void
    //    load_and_construct_derived(Archive                                        &ar,
    //                               cereal::construct<Node<State, Action, Reward>> &construct,
    //                               State                                          &state,
    //                               std::optional<Action>                          &action,
    //                               RunningStats<Reward>                           &reward_stats,
    //                               size_t                                          visit_count,
    //                               size_t                                          update_count,
    //                               bool                                            is_terminal,
    //                               bool                                            is_leaf) {
    //        std::vector<std::unique_ptr<Node<State, Action, Reward>>> children;
    //        std::weak_ptr<PyNode<State, Action, Reward>>              parent;
    //
    //        // Call deserialize on parent and children
    //        ar(cereal::make_nvp("children", children));
    //
    //        // Call construct on the PyNode
    //        construct(state,
    //                  action,
    //                  nullptr,
    //                  reward_stats,
    //                  visit_count,
    //                  update_count,
    //                  is_terminal,
    //                  is_leaf);
    //    }
};


template<typename Derived, typename State, typename Action, typename Reward>
struct MCTSBase {
    using node_type = Node<State, Action, Reward>;
    using p_node    = std::unique_ptr<node_type>;
    // Root node of the tree
    p_node root;

    explicit MCTSBase(State state, std::optional<Action> action = std::nullopt) {
        root = std::make_unique<node_type>(state, action, nullptr);
    }

    explicit MCTSBase(node_type *node) {
        root = std::make_unique<node_type>(*node);
    }

    MCTSBase(const MCTSBase &other) {
        this->root = std::make_unique<node_type>(*other.root);
    }

    MCTSBase &operator=(const MCTSBase &other) {
        if (this != &other) {
            this->root = std::make_unique<node_type>(*other.root);
        }
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const MCTSBase &mcts) {
        os << "MCTS (root): \n";
        print_node(os, *(mcts.root), "", true);
        return os;
    }

    // Serialization function definition
    template<class Archive>
    void serialize(Archive &ar) {
        //        ar(cereal::make_nvp("root", convert_node_to_py_node(*root)));
        ar(cereal::make_nvp("root", convert_node_to_py_node(*root)));

        // Call method to serialize derived parts.
        // This will call the appropriate method for the actual type of the object.
        static_cast<Derived *>(this)->serialize_derived(ar);
    }

    template<class Archive>
    static void
    load_and_construct(Archive &ar, cereal::construct<Derived> &construct) {
        std::shared_ptr<PyNode<State, Action, Reward>> root_py;

        //        // Call deserialize on all members
        //        ar(cereal::make_nvp("root", convert_py_node_to_node(*root)));

        // Call deserialize on all members
        ar(cereal::make_nvp("root", root_py));

        //        // Call construct on the PyNode
        //        construct(*root);

        // Call construct on the Node from PyNode
        construct(convert_py_node_to_node(root_py).get());

        // Call method to serialize derived parts.
        // This will call the appropriate method for the actual type of the object.
        //        Derived::load_and_construct_derived(ar);
    }


private:
    static void print_node(std::ostream &os, const node_type &node, const std::string &prefix, bool is_tail) {
        os << prefix << (is_tail ? "└─ " : "├─ ") << node << '\n';
        const auto &children = node.children;
        for (size_t i = 0; i < children.size(); ++i) {
            print_node(os, *children[i], prefix + (is_tail ? "   " : "│  "), i == children.size() - 1);
        }
    }
};


#endif// ODIN_TREE_BASE_HPP