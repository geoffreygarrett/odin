
#include <functional>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_group.h>


#include "base.hpp"
#include "include/oneapi/tbb/concurrent_priority_queue.h"
#include "include/oneapi/tbb/concurrent_vector.h"
#include "include/oneapi/tbb/parallel_sort.h"
#include "include/oneapi/tbb/task_group.h"
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

template<typename T>
class BaseSearch {
public:
    using raw_node_type = typename Tree<T>::raw_node_type;

    template<typename... Constraints>
    std::vector<raw_node_type *> search(Tree<T> &tree, Constraints &...constraints) {
        auto m_constraints = std::make_tuple(constraints...);

        bool should_terminate = std::apply([](auto &...cs) {
            return (cs.should_terminate() || ...);
        },
                                           m_constraints);

        return search_impl(tree, should_terminate);
    }

protected:
    virtual std::vector<raw_node_type *> search_impl(Tree<T> &tree, bool should_terminate) = 0;
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
class BeamSearch : public BaseSearch<T> {
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