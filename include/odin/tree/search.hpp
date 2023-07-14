#ifndef ODIN_SEARCH_HPP
#define ODIN_SEARCH_HPP

#include <functional>
#include <future>
#include <odin/tree/base.hpp>
#include <odin/tree/resource.hpp>
#include <optional>
#include <stack>
#include <tbb/concurrent_priority_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_group.h>
#include <type_traits>
#include <variant>

namespace odin {

// Empty struct used as a default argument for DerivedMetrics
struct EmptyMetrics {};


template<typename Derived,
         typename T,
         typename Float,
         typename DerivedMetrics = EmptyMetrics,
         typename DerivedResult  = std::optional<std::vector<typename Tree<T>::raw_node_type *>>>
class BaseSearch {
public:
    struct SearchMetrics {
        std::chrono::duration<Float> search_time;
        size_t                       iterations     = 0;
        size_t                       nodes_explored = 0;

        void update_metrics(std::chrono::duration<Float> time, size_t iters, size_t nodes) {
            search_time += time;
            iterations += iters;
            nodes_explored += nodes;
        }

        friend std::ostream &operator<<(std::ostream &os, const SearchMetrics &m) {
            os << "Search Time: " << m.search_time.count() << "s, "
               << "Iterations: " << m.iterations << ", "
               << "Nodes Explored: " << m.nodes_explored;
            return os;
        }
    };

    using search_return_type = std::pair<DerivedResult, SearchMetrics>;

    template<typename Strategy, typename... Constraints>
    search_return_type search(Tree<T> &tree, Strategy strategy, Constraints &...constraints) {
        auto should_terminate = [&]() {
            return (... || constraints.should_terminate());
        };
        return static_cast<Derived *>(this)->search_impl(tree, strategy, static_cast<std::function<bool()>>(should_terminate));
    }

    SearchMetrics metrics;
};


template<typename T, typename Float = double>
class DepthFirstSearch : public BaseSearch<DepthFirstSearch<T, Float>, T, Float> {
public:
    using Base          = BaseSearch<DepthFirstSearch<T, Float>, T, Float>;
    using raw_node_type = typename Tree<T>::raw_node_type;
    using SearchMetrics = typename Base::SearchMetrics;


    using search_return_type = typename Base::search_return_type;
    // Define the return_type and Result
    using return_type = raw_node_type *;
    using Result      = return_type;

    search_return_type search_impl(Tree<T>                                                     &tree,
                                   std::function<std::pair<bool, return_type>(raw_node_type *)> process_node,
                                   const std::function<bool()>                                 &should_terminate) {

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
            auto [found, value] = process_node(node);
            if (found) {
                auto end_time = std::chrono::high_resolution_clock::now();
                this->metrics.update_metrics(
                        std::chrono::duration<Float>(end_time - start_time),
                        visited.size(),
                        visited.size());
                return std::make_pair(std::optional<std::vector<raw_node_type *>>{std::vector<raw_node_type *>{value}}, this->metrics);
            }

            for (auto &child: node->get_children()) {
                stack.push(child.get());
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();

        this->metrics.update_metrics(std::chrono::duration<Float>(end_time - start_time),
                                     visited.size(), visited.size());

        return std::make_pair(std::nullopt, this->metrics);
    }
};


template<typename T, typename Float = double, typename Comparator = std::greater<>>
class BeamSearch : public BaseSearch<BeamSearch<T, Float, Comparator>, T, Float> {
public:
    using raw_node_type = typename Tree<T>::raw_node_type;
    using NodeScorePair = std::pair<raw_node_type *, double>;
    using Base          = BaseSearch<BeamSearch<T, Float, Comparator>, T, Float>;
    using SearchMetrics = typename Base::SearchMetrics;

    using search_return_type = typename Base::search_return_type;
    // Define the return_type and Result
    using return_type = raw_node_type *;
    using Result      = std::vector<return_type>;

    BeamSearch(
            std::function<double(raw_node_type *)> eval,
            unsigned int                           beam_width,
            size_t                                 reserve_amount             = 1000,
            bool                                   mark_unselected_as_visited = false)
        : m_eval(std::move(eval)),
          m_beam_width(beam_width),
          m_comp(),
          m_reserve_amount(reserve_amount),
          m_mark_unselected_as_visited(mark_unselected_as_visited) {}


    search_return_type search_impl(
            Tree<T>                                                     &tree,
            std::function<std::pair<bool, return_type>(raw_node_type *)> process_node,
            const std::function<bool()>                                 &should_terminate) {

        tbb::concurrent_vector<raw_node_type *> frontier = {tree.get_root()};
        tbb::concurrent_vector<raw_node_type *> visited;

        size_t total_nodes_processed = 0;
        size_t iterations            = 0;
        auto   start_time            = std::chrono::high_resolution_clock::now();

        while (!frontier.empty() && !should_terminate()) {
            tbb::concurrent_vector<NodeScorePair> next_frontier;

            tbb::parallel_for_each(
                    frontier.begin(),
                    frontier.end(),
                    [this, &next_frontier](raw_node_type *node) {
                        for (auto &child: node->get_children()) {
                            double score = m_eval(child.get());
                            next_frontier.push_back({child.get(), score});
                        }
                    });

            tbb::parallel_sort(
                    next_frontier.begin(),
                    next_frontier.end(),
                    [&](const NodeScorePair &a, const NodeScorePair &b) {
                        return m_comp(a.second, b.second);
                    });

            // Get the start of the unvisited nodes in next_frontier.
            auto unvisited_start = next_frontier.begin() + std::min(static_cast<size_t>(m_beam_width), next_frontier.size());

            if (m_mark_unselected_as_visited) {
                std::for_each(unvisited_start, next_frontier.end(), [&](auto &node_score_pair) {
                    visited.push_back(node_score_pair.first);
                });
            }

            // Set frontier to the next set of nodes to visit.
            frontier.clear();
            std::transform(next_frontier.begin(),
                           unvisited_start,
                           std::back_inserter(frontier),
                           [](const auto &node_score_pair) { return node_score_pair.first; });
            total_nodes_processed += next_frontier.size();
            iterations++;
            next_frontier.clear();
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        this->metrics.update_metrics(std::chrono::duration<Float>(end_time - start_time), iterations, total_nodes_processed);

        std::vector<raw_node_type *> visited_vector(visited.begin(), visited.end());
        return std::make_pair(std::optional<std::vector<raw_node_type *>>(std::move(visited_vector)), this->metrics);
    }

private:
    std::function<double(raw_node_type *)> m_eval;
    unsigned int                           m_beam_width;
    Comparator                             m_comp;
    size_t                                 m_reserve_amount;
    bool                                   m_mark_unselected_as_visited;
};

}// namespace odin
#endif// ODIN_SEARCH_H