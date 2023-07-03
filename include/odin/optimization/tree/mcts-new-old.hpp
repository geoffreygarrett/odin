#include "base.hpp"
#include <algorithm>
#include <complex>
#include <functional>
#include <memory>
#include <random>
#include <stack>
#include <vector>


// Base Selection Policy
template<typename Derived, typename Node>
struct SelectionPolicy {

    double operator()(Node node) {
        return static_cast<Derived *>(this)->select_score(node);
    }
};

// UBC1 Selection Policy
template<typename Node, typename Float = double>
struct UBC1 : public SelectionPolicy<UBC1<Node>, Node> {

    Float cp_;

    explicit UBC1(Float cp) : cp_(cp) {}

    Node select_node(const Node &node) {
        Float parent_visits = static_cast<double>(node->parent->visits);
        return node->reward + cp_ * std::sqrt(std::log(parent_visits) / node->visits);
    }
};

template<typename Node, typename Float = double>
struct EpsilonGreedy : public SelectionPolicy<EpsilonGreedy<Node>, Node> {
    Float                                 epsilon;
    std::mt19937                          engine;
    std::uniform_real_distribution<Float> dist;

    explicit EpsilonGreedy(Float epsilon) : epsilon(epsilon), dist(0.0, 1.0), engine(std::random_device{}()) {}

    std::shared_ptr<Node> select_node(std::shared_ptr<Node> node) {
        std::shared_ptr<Node> selected_node = node;
        Float                 min_score     = std::numeric_limits<Float>::max();
        for (auto &child: node->children) {
            Float avg_reward = child->total_reward / child->visit_count;
            Float score      = avg_reward + epsilon * node->visit_count / child->visit_count;
            if (score < min_score) {
                min_score     = score;
                selected_node = child;
            }
        }
        return selected_node;
    }
};

//// Helper function to calculate UCB1 score
//template<typename Node, typename Integer>
//double fn_ubc1(const Node &node, double Cp) {
//    Integer parent_visits = static_cast<double>(node->parent->visits);
//    return node->reward + Cp * std::sqrt(std::log(parent_visits) / node->visits);
//}
//
//template<typename Node, typename Float>
//double fn_epsilon_greedy(const Node &node, double epsilon) {
//    return dist(engine) < epsilon ? std::numeric_limits<double>::max() : node.reward;
//}
//
//// Helper function to calculate UCB1-Tuned score
//template<typename Node, typename Float>
//double fn_ubc1_tuned(const Node &node) {
//    Float parent_visits  = static_cast<double>(node->parent->visits);
//    Float child_variance = node->reward_variance + std::sqrt((2.0 * std::log(parent_visits)) / node->visits);
//    return node->reward + std::sqrt((std::log(parent_visits) / node->visits) * std::min(0.25, child_variance));
//}

template<typename State, typename Action, typename Result, typename S>
class MCTS : public TreeSearchBase<MCTS<State, Action, Result>, State, Action, Result> {
public:
    using Base       = TreeSearchBase<MCTS<State, Action, Result>, State, Action, Result>;
    using p_node     = typename Base::p_node;
    using trajectory = std::vector<std::tuple<State, Action, Result>>;

    std::function<std::vector<Action>(State)> mfn_get_actions;// Function to get available actions for a state
    std::function<State(State, Action)>       mfn_transition; // Function to get new state after an action
    std::function<bool(State)>                mfn_is_terminal;// Function to check if a state is terminal
    std::function<Result(State, bool)>        mfn_reward;     // Function to get the reward of a state
    SelectionPolicy                           m_selection_policy;
    double                                    m_epsilon{}; // Epsilon for epsilon-greedy selection
    double                                    m_Cp{};      // Exploration constant for UCB1
    double                                    m_Cp_tuned{};// Exploration constant for UCB1-Tuned


    MCTS(State                                     initial_state,
         Action                                    initial_action,
         std::function<std::vector<Action>(State)> get_actions,
         std::function<State(State, Action)>       transition,
         std::function<bool(State)>                is_terminal,
         std::function<Result(State, bool)>        reward,
         SelectionPolicy                           selection_policy = EpsilonGreedy(0.1)
        : Base(initial_state, initial_action),
          mfn_get_actions(get_actions),
          mfn_transition(transition),
          mfn_is_terminal(is_terminal),
          m_selection_policy(selection_policy),
          mfn_reward(reward) {}


    p_node select_node(p_node node) {
        auto selected_child = std::max_element(
                node->children.begin(),
                node->children.end(),
                [&selection_fn](const auto &a, const auto &b) {
                    return selection_fn(a.get()) < selection_fn(b.get());
                });

        return selected_child->get();

        //        if (m_selection_policy == UCB1) {
        //            return select_node_impl(node, fn_ubc1<p_node, int>);
        //        } else if (m_selection_policy == UCB1_TUNED) {
        //            return select_node_impl(node, fn_ubc1_tuned<p_node, double>);
        //        } else if (m_selection_policy == EPSILON_GREEDY) {
        //            return select_node_impl(node, fn_epsilon_greedy<p_node, double>);
        //        } else {
        //            throw std::runtime_error("Unknown selection policy: " + std::to_string(m_selection_policy));
        //        }
    }

    p_node expand_node_impl(p_node node) {
        // Stop expanding if the node is terminal
        if (is_terminal(node->state)) {
            node->is_terminal = true;
            return node;
        }

        // Generate possible actions
        auto possible_actions = get_possible_actions(node->state);

        // Expand node by creating children for all possible actions
        for (auto const &action: possible_actions) {
            State new_state = perform_action(node->state, action);

            // Check for duplicate nodes
            if (std::none_of(node->children.begin(), node->children.end(), [&](auto &child) { return child->state == new_state; })) {
                auto child = std::make_shared<Node>(node, new_state, action, !node->move);
                node->children.push_back(child);
            }
        }

        // If the expansion resulted in no children (all actions lead to already existing states), mark the node as terminal
        if (node->children.empty()) {
            node->is_terminal = true;
        }
    }

    //    void expand_node_impl(p_node node) {
    //        auto trajectory = simulate(node->state);
    //        for (const auto &[state, action, result]: trajectory) {
    //            node->children.push_back(std::make_unique<Node<State, Action, Result>>(state, action, result, {}));
    //        }
    //    }

    p_node best_child(const p_node &node) const {
        return *std::max_element(
                node->children.begin(), node->children.end(),
                [](const p_node &a, const p_node &b) {
                    return a->reward < b->reward;
                });
    }

    std::vector<std::tuple<State, Action, Result>> best_trajectory(p_node node) const {
        std::vector<std::tuple<State, Action, Result>> trajectory;

        // Descend through the tree to the best leaf node.
        while (!node->is_leaf()) {
            node = best_child(node);

            // Add the state, action, and result associated with the node to the trajectory.
            trajectory.emplace_back(node->state, node->action, node->result);
        }

        return trajectory;
    }


    trajectory search_impl(int iterations = 1000) {
        while (iterations--) {
            // 1. Selection
            // Start from the root node and use the selection policy to descend through the tree,
            // until reaching a leaf node that hasn't been expanded yet.
            p_node node = this->root;
            while (!node->is_leaf()) {
                node = select_node(node);
            }

            // 2. Expansion
            // For the selected leaf node, compute all possible actions from its state and add a
            // new child node for each action to the tree.
            expand_node(node);

            // 3. Simulation
            // For each child node of the newly expanded node, simulate a trajectory from the child
            // node's state to a terminal state. Each state-action-result triplet in the trajectory
            // is added to the tree as a child node of the preceding state. This results in a chain
            // of nodes representing the simulated trajectory. The final state of the trajectory is
            // a terminal state, and the associated result is backpropagated up through the chain
            // of nodes.
            // The 'Simulation' step in `search_impl`
            for (auto &child_node: node->children) {
                // Perform a simulation from the child node's state and append the trajectory to the child node.
                auto trajectory = simulate(child_node->state);

                // Expand the tree along the simulated trajectory.
                auto terminal_node = expand_trajectory(child_node.get(), trajectory);

                // 4. Backpropagation
                // Backpropagate the result of the simulated trajectory up the tree through all
                // parent nodes, starting with the terminal node and up to the child node.
                for (auto current_node = terminal_node; current_node != child_node.get(); current_node = current_node->parent) {
                    backpropagate(current_node, current_node->result);
                }
            }

            // Step 5: Contraction
            // Based on an extension to the four steps of MCTS (selection, expansion, simulation, backpropagation),
            // introduced by Hennes and Izzo in "Interplanetary Trajectory Planning with Monte Carlo Tree Search"
            // (IJCAI 2015), we include a contraction step.
            //
            // This step is particularly relevant for trajectory planning problems that span a very broad tree,
            // with a high branching factor and limited depth. For such problems, the internal tree of the MCTS
            // is often able to reach final nodes. When a subtree of the internal tree fully covers the corresponding
            // problem search space (i.e., each child has been played at least once and all leaves are terminal states),
            // it can be ignored in further search iterations.
            //
            // By removing such subtrees from the internal tree, the contraction step helps to keep the tree compact
            // and manageable across multiple iterations, enhancing the overall performance of the algorithm.
            //
            // Reference:
            // Hennes, D., & Izzo, D. (2015). Interplanetary Trajectory Planning with Monte Carlo Tree Search.
            // In Proceedings of the Twenty-Fourth International Joint Conference on Artificial Intelligence (IJCAI 2015).
            // European Space Agency, Advanced Concepts Team, Noordwijk, The Netherlands.
            contract(this->root);
        }

        // After time limit is reached, return the best node from the root according to a policy (e.g., the node with the highest value).
        return best_trajectory(this->root);
    }


private:
    // Expands the tree along a simulated trajectory by creating a chain of nodes.
    // The chain starts from a given node (initially an unexpanded leaf node) and
    // extends to a terminal node representing the end of a simulated trajectory.
    p_node expand_trajectory(p_node node, const std::vector<std::tuple<State, Action, Result>> &trajectory) {
        // Iterate over each state-action-result triplet in the simulated trajectory.
        for (const auto &[state, action, result]: trajectory) {
            // Create a new node for each triplet.
            // The parent of each new node is set to the preceding node in the trajectory.
            // `std::make_unique` is used to create a new unique_ptr object that takes ownership of the new node.
            auto new_node = std::make_unique<Node<State, Action, Result>>(state, action, result, node.get());

            // Add the new node to the list of children of the preceding node.
            // `std::move` is used to transfer ownership of the new node from the unique_ptr to the child list.
            node->children.push_back(std::move(new_node));

            // Update `node` to point to the new node, making it the "current" node for the next iteration.
            node = node->children.back().get();
        }

        // After all triplets in the trajectory have been processed, `node` points to the terminal node of the trajectory.
        // Return a raw pointer to the terminal node.
        return node;
    }

    //    Result simulate(State state) {
    //        while (!is_terminal(state)) {
    //            std::vector<Action> actions = get_actions(state);
    //            Action              action  = actions[rand() % actions.size()];// Select a random action
    //            state                       = transition(state, action);       // Apply the action to get the new state
    //        }
    //        return reward(state, true);// Compute the reward for the terminal state
    //    }
    //    p_node select_node_impl(p_node node, SelectionFunction<p_node> selection_fn) {
    //        auto best_child = std::max_element(
    //                node->children.begin(),
    //                node->children.end(),
    //                [&selection_fn](const auto &a, const auto &b) {
    //                    return selection_fn(a.get()) < selection_fn(b.get());
    //                });
    //
    //        return best_child->get();
    //    }
    void contract(p_node node) {
        for (auto it = node->children.begin(); it != node->children.end();) {
            if (all_children_terminal(*it)) {
                it = node->children.erase(it);
            } else {
                ++it;
            }
        }
        for (auto &child: node->children) {
            contract(child);
        }
    }

    bool all_children_terminal(p_node node) {
        for (auto &child: node->children) {
            if (!is_terminal(child->state) || child->children.empty()) {
                return false;
            }
        }
        return true;
    }

    std::vector<std::tuple<State, Action, Result>> simulate(State state) {
        std::vector<std::tuple<State, Action, Result>> trajectory;
        while (!is_terminal(state)) {
            std::vector<Action> actions = get_actions(state);
            Action              action  = actions[rand() % actions.size()];// Select a random action
            state                       = transition(state, action);       // Apply the action to get the new state
            Result result               = reward(state, false);            // Compute the reward for the current state
            trajectory.emplace_back(state, action, result);
        }
        return trajectory;
    }

    void backpropagate(p_node node, Result result) {
        // Start from the given node and traverse up the tree
        while (node) {
            // Increment the visit count of the node
            node->visit_count++;// TODO: with the multiple leaf expansion, this will throw off the selection policies original intended behaviours
                                // .... actually it wont, we must just make sure not to have the expanded node visited for the number of parallel simulations

            // Update the value of the node based on the result
            node->value += result;

            // Move to the parent of the current node
            node = node->parent;
        }
    }
};