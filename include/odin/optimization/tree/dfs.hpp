
#include "base.hpp"
#include <stack>

template<typename State, typename Action, typename Metric>
struct DepthFirstSearch : public TreeSearchBase<DepthFirstSearch<State, Action, Metric>, State, Action, Metric> {
    using Base   = TreeSearchBase<DepthFirstSearch<State, Action, Metric>, State, Action, Metric>;
    using p_node = typename Base::p_node;
    std::stack<p_node> stack;// Stack to store the nodes for DFS

    p_node expand_node_impl(p_node node) {
        // In DFS, we would generate the successors of the current node and add them to the stack.
        for (const auto &child: node->children) {
            stack.push(child);
        }
        return node;
    }

    p_node best_child_impl(p_node node) {
        // In DFS, the best child is the last one added to the stack (LIFO).
        // The caller needs to ensure that the node has children before calling this method.
        return stack.top();
    }

    p_node search_impl() {
        // Push the root node into the stack
        stack.push(this->root);

        // Continue while there are still nodes in the stack
        while (!stack.empty()) {
            // Pop the top node from the stack
            auto node = stack.top();
            stack.pop();

            // Expand the node and add its children to the stack
            expand_node(node);

            // If this node is a goal (however that's defined), return it
            if (is_goal(node)) {
                return node;
            }
        }

        // No goal was found
        return nullptr;
    }

    // A placeholder for the actual goal check function
    // This function should be defined appropriately according to the problem at hand
    bool is_goal(p_node node) {
        // TODO: implement
        return false;
    }
};
