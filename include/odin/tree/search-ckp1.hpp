
#include "base.hpp"
#include "include/oneapi/tbb/concurrent_priority_queue.h"
#include "include/oneapi/tbb/concurrent_vector.h"
#include "include/oneapi/tbb/parallel_sort.h"
#include "include/oneapi/tbb/task_group.h"

#include <functional>

#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_group.h>
//
//template<typename T, typename NodeEvaluator>
//class BeamSearch {
//public:
//    using tree_type      = Tree<T>;
//    using raw_node_type  = typename tree_type::raw_node_type;
//    using safe_node_type = typename tree_type::safe_node_type;
//    using p_raw_node     = typename tree_type::p_raw_node;
//    using p_safe_node    = typename tree_type::p_safe_node;
//
//    explicit BeamSearch(tree_type &tree, NodeEvaluator eval, size_t beam_width)
//        : m_tree(tree), m_eval(std::move(eval)), m_beam_width(beam_width) {}
//
//    std::vector<p_raw_node> search() {
//        std::vector<p_raw_node> frontier{m_tree.get_root()};
//        std::vector<p_raw_node> next_frontier;
//
//        while (!frontier.empty()) {
//            next_frontier.clear();
//
//            for (auto &node: frontier) {
//                auto children = node->get_children();
//                next_frontier.insert(
//                        next_frontier.end(),
//                        std::make_move_iterator(children.begin()),
//                        std::make_move_iterator(children.end()));
//            }
//
//            std::sort(next_frontier.begin(), next_frontier.end(),
//                      [this](const p_raw_node &a, const p_raw_node &b) {
//                          return m_eval(*a) > m_eval(*b);
//                      });
//
//            if (next_frontier.size() > m_beam_width) {
//                next_frontier.resize(m_beam_width);
//            }
//
//            frontier = std::move(next_frontier);
//        }
//
//        return frontier;
//    }
//
//private:
//    tree_type    &m_tree;
//    NodeEvaluator m_eval;
//    size_t        m_beam_width;
//};

#include <functional>
#include <future>

#include <stack>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_group.h>

//template<typename T>
//class BaseSearch {
//public:
//    using raw_node_type = typename Tree<T>::raw_node_type;
//
//    template<typename Derived, typename... Constraints>
//    std::vector<raw_node_type *> search(Tree<T> &tree, Constraints &...constraints) {
//        auto m_constraints = std::make_tuple(constraints...);
//
//        bool should_terminate = std::apply([](auto &...cs) {
//            return (cs.should_terminate() || ...);
//        },
//                                           m_constraints);
//
//        return static_cast<Derived *>(this)->search_impl(tree, should_terminate);
//    }
//};
//
//template<typename T, typename Comparator = std::greater<>>
//class BeamSearch {
//public:
//    using raw_node_type = typename Tree<T>::raw_node_type;
//    using NodeScorePair = std::pair<raw_node_type *, double>;
//
//    explicit BeamSearch(Tree<T> &tree, std::function<double(raw_node_type *)> eval, unsigned int beam_width)
//        : m_tree(tree), m_eval(std::move(eval)), m_beam_width(beam_width), m_comp() {}
//
//    std::vector<raw_node_type *> search() {
//        std::vector<raw_node_type *> frontier = {m_tree.get_root()};
//        while (!frontier.empty()) {
//            tbb::concurrent_vector<NodeScorePair> next_frontier;
//            tbb::task_group                       tasks;
//
//            for (auto node: frontier) {
//                tasks.run([&, node] {
//                    for (auto &child: node->get_children()) {
//                        double score = m_eval(child.get());
//                        next_frontier.push_back({child.get(), score});
//                    }
//                });
//            }
//
//            tasks.wait();
//
//            std::vector<NodeScorePair> sorted_frontier(next_frontier.begin(), next_frontier.end());
//            tbb::parallel_sort(sorted_frontier.begin(), sorted_frontier.end(),
//                               [&](const NodeScorePair &a, const NodeScorePair &b) {
//                                   return m_comp(a.second, b.second);
//                               });
//
//            frontier.clear();
//            for (size_t i = 0; i < std::min(m_beam_width, sorted_frontier.size()); ++i) {
//                frontier.push_back(sorted_frontier[i].first);
//            }
//        }
//
//        return frontier;
//    }
//
//private:
//    Tree<T>                               &m_tree;
//    std::function<double(raw_node_type *)> m_eval;
//    unsigned int                           m_beam_width;
//    Comparator                             m_comp;
//};

#include <optional>
#include <type_traits>
#include <variant>

// This is a helper struct and a detector function for the SFINAE trick
template<template<typename...> class Trait, typename Enabler, typename... Args>
struct is_detected : std::false_type {};

template<template<typename...> class Trait, typename... Args>
struct is_detected<Trait, std::void_t<Trait<Args...>>, Args...> : std::true_type {};

template<template<typename...> class Trait, typename... Args>
constexpr bool is_detected_v = is_detected<Trait, void, Args...>::value;

template<typename T>
using derived_metrics_t = typename T::DerivedMetrics;

template<typename T>
using derived_result_t = typename T::DerivedResult;

template<typename...>
using void_t = void;

template<typename, template<typename> class, typename = void_t<>>
struct detect : std::false_type {};

template<typename T, template<typename> class Op>
struct detect<T, Op, void_t<Op<T>>> : std::true_type {};

template<typename T>
using has_derived_metrics = detect<T, derived_metrics_t>;


// Empty struct used as a default argument for DerivedMetrics
struct EmptyMetrics {};


template<typename Derived, typename T, typename Float = double>
class BaseSearch {
public:
    // We make raw_node_type public so that derived classes can access it
    using raw_node_type = typename Tree<T>::raw_node_type;

    // Define a structure to hold the base search metrics
    // All search algorithms will have these metrics
    struct BaseMetrics {
        std::chrono::duration<Float> search_time;   // The total time taken for the search
        size_t                       iterations;    // The total number of iterations
        size_t                       nodes_explored;// The total number of nodes explored during the search
    };

    // Using SFINAE (Substitution Failure Is Not An Error) and std::conditional_t to check if the Derived class
    // has defined its own DerivedMetrics. If DerivedMetrics is present in the Derived class, it is used.
    // Otherwise, we use the default EmptyMetrics. This flexibility allows the derived class to define
    // its own metrics for its search algorithm, and falls back to a no-op metrics (EmptyMetrics) if not provided.
    struct DefaultDerivedMetrics {};
    using DerivedMetrics = DefaultDerivedMetrics;

    //    using DerivedMetrics = std::conditional_t<has_derived_metrics<Derived>::value, derived_metrics_t<Derived>, EmptyMetrics>;

    // SearchMetrics is a combination of BaseMetrics and DerivedMetrics.
    // BaseMetrics are metrics common to all search algorithms such as search time, number of iterations, and nodes explored.
    // DerivedMetrics are specific to the derived class and may not exist if the derived class does not define it.
    // If DerivedMetrics is not provided by the derived class, it defaults to EmptyMetrics (which adds nothing).
    //    struct SearchMetrics : BaseMetrics, DerivedMetrics {};
    //
    //    // If Derived provides a search_result_t, it's used; otherwise, it defaults to an optional vector of raw_node_type pointers
    //    using DerivedResult = std::conditional_t<is_detected_v<derived_result_t, Derived>, derived_result_t<Derived>, std::optional<std::vector<raw_node_type *>>>;
    //
    //    // The search function is the main entry point for the search algorithm
    //    using search_return_type = std::conditional_t<std::is_same_v<DerivedResult, std::nullopt_t>, SearchMetrics, std::pair<DerivedResult, SearchMetrics>>;

    struct EmptyMetrics {};
    using SearchMetrics = BaseMetrics;
    //    using DerivedMetrics = EmptyMetrics;
    using DerivedResult = std::optional<std::vector<raw_node_type *>>;

    using search_return_type = std::conditional_t<std::is_same_v<DerivedResult, std::nullopt_t>, SearchMetrics, std::pair<DerivedResult, SearchMetrics>>;


    //    // search function checks the constraints and then calls the derived class's search implementation
    //    template<typename... Constraints>
    //    search_return_type search(Tree<T> &tree, Constraints &...constraints) {
    //        auto m_constraints = std::make_tuple(constraints...);
    //
    //        bool should_terminate
    //                = std::apply([](auto &...cs) { return (cs.should_terminate() || ...); },
    //                             m_constraints);
    //
    //        return static_cast<Derived *>(this)->search_impl(tree, should_terminate);
    //    }

    template<typename Strategy, typename... Constraints>
    search_return_type search(Tree<T> &tree, Strategy strategy, Constraints... constraints) {
        std::function<bool()> should_terminate = [&]() {
            return std::apply([](auto &...cs) { return (cs.should_terminate() || ...); },
                              std::make_tuple(constraints...));
        };

        return strategy(static_cast<Derived *>(this), tree, should_terminate);
    }

    // Asynchronous search method
    // This will run the search_impl of the derived class in a separate thread
    template<typename... Constraints>
    std::future<search_return_type> search_async(Tree<T> &tree, Constraints &&...constraints) {
        return std::async(std::launch::async,
                          &Derived::search_impl,
                          static_cast<Derived *>(this),
                          std::ref(tree),
                          std::forward<Constraints>(constraints)...);
    }

    // Parallel search method using TBB
    // This will call the derived class's parallel search implementation
    template<typename... Constraints>
    search_return_type search_parallel(Tree<T> &tree, Constraints &&...constraints) {
        return static_cast<Derived *>(this)->search_impl_parallel(
                tree,
                std::forward<Constraints>(constraints)...);
    }


protected:
    // Helper function to check constraints
    // This is used by the search function to decide when to terminate the search
    template<typename... Constraints>
    bool check_constraints(std::tuple<Constraints...> &constraints) {
        return std::apply([](auto &...cs) { return (cs.should_terminate() || ...); }, constraints);
    }
};

template<typename T, typename Float = double>
class DepthFirstSearch : public BaseSearch<DepthFirstSearch<T, Float>, T, Float> {
public:
    using Base               = BaseSearch<DepthFirstSearch<T, Float>, T, Float>;
    using raw_node_type      = typename Base::raw_node_type;
    using SearchMetrics      = typename Base::SearchMetrics;
    using return_type        = std::variant<raw_node_type *, std::vector<raw_node_type *>>;
    using search_return_type = std::pair<return_type, SearchMetrics>;
    using DerivedResult      = typename Base::DerivedResult;
    // DFS does not need any additional metrics, so we use BaseMetrics
    //    using DerivedMetrics = Base::BaseMetrics;

    search_return_type search_impl(Tree<T>                             &tree,
                                   std::function<void(raw_node_type *)> process_node,
                                   const std::function<bool()>         &should_terminate) {
        std::stack<raw_node_type *>  stack;
        std::vector<raw_node_type *> visited;

        auto root = tree.get_root();
        if (root) {
            stack.push(root);
        }

        auto start_time = std::chrono::high_resolution_clock::now();
        while (!stack.empty() && !should_terminate()) {
            auto node = stack.top();
            stack.pop();

            visited.push_back(node);

            // Process the node using the lambda function
            process_node(node);

            for (auto &child: node->get_children()) {
                stack.push(child.get());
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();

        SearchMetrics metrics;
        metrics.search_time    = std::chrono::duration<Float>(end_time - start_time);
        metrics.iterations     = visited.size();
        metrics.nodes_explored = visited.size();

        // Return pair of result and metrics
        return std::make_pair(DerivedResult{visited}, metrics);
    }
    // DFS is inherently sequential and cannot be easily parallelized
    search_return_type search_impl_parallel(Tree<T> &tree, std::function<std::pair<bool, return_type>(raw_node_type *)> condition, bool should_terminate) {
        return this->search_impl(tree, condition, should_terminate);
    }
};

//
//template<typename T, typename Comparator = std::greater<>>
//class BeamSearch : public BaseSearch<T> {
//public:
//    using raw_node_type = typename Tree<T>::raw_node_type;
//    using NodeScorePair = std::pair<raw_node_type *, double>;
//
//    explicit BeamSearch(std::function<double(raw_node_type *)> eval, unsigned int beam_width)
//        : m_eval(std::move(eval)), m_beam_width(beam_width), m_comp() {}
//
//protected:
//    std::vector<raw_node_type *> search_impl(Tree<T> &tree, bool should_terminate) override {
//        std::vector<raw_node_type *> frontier = {tree.get_root()};
//        while (!frontier.empty() && !should_terminate) {
//            tbb::concurrent_vector<NodeScorePair> next_frontier;
//            tbb::task_group                       tasks;
//
//            for (auto node: frontier) {
//                tasks.run([&, node] {
//                    for (auto &child: node->get_children()) {
//                        double score = m_eval(child.get());
//                        next_frontier.push_back({child.get(), score});
//                    }
//                });
//            }
//
//            tasks.wait();
//
//            std::vector<NodeScorePair> sorted_frontier(next_frontier.begin(), next_frontier.end());
//            tbb::parallel_sort(sorted_frontier.begin(), sorted_frontier.end(),
//                               [&](const NodeScorePair &a, const NodeScorePair &b) {
//                                   return m_comp(a.second, b.second);
//                               });
//
//            frontier.clear();
//            for (size_t i = 0; i < std::min(m_beam_width, sorted_frontier.size()); ++i) {
//                frontier.push_back(sorted_frontier[i].first);
//            }
//        }
//
//        return frontier;
//    }
//
//private:
//    std::function<double(raw_node_type *)> m_eval;
//    unsigned int                           m_beam_width;
//    Comparator                             m_comp;
//};

template<typename T, typename Comparator = std::greater<>>
class BeamSearch : public BaseSearch<BeamSearch<T, Comparator>, T> {
public:
    using raw_node_type = typename Tree<T>::raw_node_type;
    using NodeScorePair = std::pair<double, raw_node_type *>;

    explicit BeamSearch(std::function<double(raw_node_type *)> eval, unsigned int beam_width)
        : m_eval(std::move(eval)), m_beam_width(beam_width), m_comp() {}

protected:
    std::vector<raw_node_type *> search_impl(Tree<T> &tree, bool should_terminate) override {
        tbb::concurrent_priority_queue<NodeScorePair, Comparator> frontier;
        frontier.push({m_eval(tree.get_root()), tree.get_root()});

        while (!frontier.empty() && !should_terminate) {
            tbb::concurrent_priority_queue<NodeScorePair, Comparator> next_frontier;
            tbb::task_group                                           tasks;

            while (!frontier.empty() && next_frontier.size() < m_beam_width) {
                auto [score, node] = frontier.top();
                frontier.pop();

                tasks.run([&, node] {
                    for (auto &child: node->get_children()) {
                        double score = m_eval(child.get());
                        next_frontier.push({score, child.get()});
                    }
                });
            }

            tasks.wait();
            frontier.swap(next_frontier);
        }

        std::vector<raw_node_type *> result;
        while (!frontier.empty()) {
            auto [_, node] = frontier.top();
            result.push_back(node);
            frontier.pop();
        }

        return result;
    }

private:
    std::function<double(raw_node_type *)> m_eval;
    unsigned int                           m_beam_width;
    Comparator                             m_comp;
};

//https://chat.openai.com/c/e135a210-61eb-4a2e-9bb1-18372bf88cb9
//enum class ParallelScheme {
//    NoParallelism,
//    EvalParallelism,
//    NestedParallelism
//};
//
//template<typename T, typename Comparator = std::greater<>>
//class BeamSearch : public BaseSearch<T> {
//    // ...
//    explicit BeamSearch(std::function<double(raw_node_type *)> eval, unsigned int beam_width, ParallelScheme scheme)
//        : m_eval(std::move(eval)), m_beam_width(beam_width), m_comp(), m_scheme(scheme) {}
//    // ...
//private:
//    // ...
//    ParallelScheme m_scheme;
//};