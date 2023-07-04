#ifndef MCTS_HPP
#define MCTS_HPP

#include "base.hpp"
#include <algorithm>
#include <complex>
#include <functional>
#include <memory>
#include <odin/logging.hpp>
#include <random>
//#include <ranges>
#include <stack>
#include <variant>
#include <vector>


template<typename Node, typename Float = double>
class SelectionPolicy {
public:
    using value_estimator_type = std::function<Float(const Node *)>;

protected:
    value_estimator_type mfn_value_estimator;

public:
    explicit SelectionPolicy(value_estimator_type value_estimator) : mfn_value_estimator(value_estimator) {}

    Node *operator()(Node *node) {
        if (node->children.empty()) {
            throw std::runtime_error("No children to select from");
        }

        auto comparator = [&](const auto &a, const auto &b) {
            if (!a || !b) return false;

            if (a->visit_count == 0) return false;
            if (b->visit_count == 0) return true;

            Float a_val = this->value(a.get(), node);
            Float b_val = this->value(b.get(), node);
            return a_val < b_val;
        };

        auto selected_child = std::max_element(node->children.begin(), node->children.end(), comparator);

        if (selected_child == node->children.end()) {
            throw std::runtime_error("Failed to select child");
        }

        return selected_child->get();
    }

    virtual Float value(Node *child, Node *parent) = 0;// Pure virtual function

    virtual ~SelectionPolicy() = default;// Virtual destructor
};

template<typename Node, typename Float = double>
class UCB1 : public SelectionPolicy<Node, Float> {
    Float m_cp;

public:
    using value_estimator_type = typename SelectionPolicy<Node, Float>::value_estimator_type;

    UCB1(Float cp, value_estimator_type value_estimator)
            : SelectionPolicy<Node, Float>(value_estimator), m_cp(cp) {}

    Float get_cp() const {
        return m_cp;
    }

    Float value(Node *child, Node *parent) override {
        return this->mfn_value_estimator(child) + m_cp * std::sqrt(std::log(parent->visit_count) / child->visit_count);
    }
};

template<typename Node, typename Float = double>
class UCB1Tuned : public SelectionPolicy<Node, Float> {
    Float m_cp;

public:
    using value_estimator_type = typename SelectionPolicy<Node, Float>::value_estimator_type;

    UCB1Tuned(Float cp, value_estimator_type value_estimator)
            : SelectionPolicy<Node, Float>(value_estimator), m_cp(cp) {}

    Float get_cp() const {
        return m_cp;
    }

    Float value(Node *child, Node *parent) override {
        Float Vi = child->reward_stats.variance() + std::sqrt(2 * std::log(parent->visit_count) / child->visit_count);
        return this->mfn_value_estimator(child) +
               m_cp * std::sqrt(std::log(parent->visit_count) / child->visit_count * std::min(0.25, Vi));
    }
};

template<typename Node, typename Float = double>
class EpsilonGreedy : public SelectionPolicy<Node, Float> {
    Float m_epsilon;

public:
    using value_estimator_type = typename SelectionPolicy<Node, Float>::value_estimator_type;

    explicit EpsilonGreedy(Float epsilon, value_estimator_type value_estimator)
            : SelectionPolicy<Node, Float>(value_estimator), m_epsilon(epsilon) {}

    Float get_epsilon() const {
        return m_epsilon;
    }

    Float value(Node *child, Node *parent) override {
        return this->mfn_value_estimator(child) + m_epsilon * parent->visit_count / child->visit_count;
    }
};


template<typename State, typename Action, typename Reward, typename Float = double>
class MCTS : public MCTSBase<MCTS<State, Action, Reward, Float>, State, Action, Reward> {
public:
    using Base = MCTSBase<MCTS<State, Action, Reward, Float>, State, Action, Reward>;
    using p_node = typename Base::p_node;
    using node_type = typename Base::node_type;
    using trajectory_type = std::vector<std::tuple<State, std::optional<Action>, Reward>>;
    using state_transition = std::function<State(State, Action)>;
    using action_generator = std::function<std::vector<Action>(State &)>;
    using selection_policy = std::shared_ptr<SelectionPolicy<node_type, Float>>;
    using is_terminal = std::function<bool(State &)>;
    using reward = std::function<Reward(State &)>;

    action_generator mfn_get_actions;// Function to get available actions for a state
    state_transition mfn_transition; // Function to get new state after an action
    is_terminal mfn_is_terminal;// Function to check if a state is terminal
    reward mfn_reward;     // Function to get the reward of a state
    selection_policy mfn_selection_policy;

    std::mt19937 engine;


    MCTS(State initial_state,
         action_generator get_actions,
         state_transition transition,
         is_terminal is_terminal,
         reward reward,
         selection_policy selection_policy = std::make_shared<UCB1<node_type, Float>>(1.0,
                                                                                      [](const node_type *node) {
                                                                                          return node->reward_stats.max();
                                                                                      }),
         size_t seed = std::random_device{}(),
         std::optional<Action> initial_action = std::nullopt)// actions are generally assigned by the environment through the process
            : Base(initial_state, std::nullopt),                 // use std::nullopt as the initial action
              mfn_get_actions(get_actions),
              mfn_transition(transition),
              mfn_is_terminal(is_terminal),
              mfn_reward(reward),
              mfn_selection_policy(std::move(selection_policy)),
              engine(seed) {}


    node_type *get_root() {
        return Base::root.get();
    }

    node_type *select_node(node_type *node) {
        return mfn_selection_policy->operator()(node);
    }


    void expand_all_actions(node_type *node) {
        try {
            // Stop expanding if the node is terminal
            if (mfn_is_terminal(node->state)) {
                node->is_terminal = true;
                return;
            }
            node->is_leaf = false;

            // Generate possible actions
            auto possible_actions = mfn_get_actions(node->state);

            // If no actions are possible, mark the node as terminal
            if (possible_actions.empty()) {
                node->is_terminal = true;
                return;
            }

            // Expand node by creating children for all possible actions
            for (auto const &action: possible_actions) {
                State new_state = mfn_transition(node->state, action);

                // Check for duplicate nodes
                if (std::none_of(node->children.begin(), node->children.end(),
                                 [&](const std::unique_ptr<node_type> &child) {
                                     return child != nullptr && child->state == new_state;
                                 })) {
                    auto reward_stats = RunningStats<Reward>{};
                    reward_stats.set_calc_max(true);
                    auto child = std::make_unique<node_type>(new_state, action, node, reward_stats);
                    node->children.push_back(std::move(child));// take ownership of the child
                }
            }

            // If the expansion resulted in no children (all actions lead to already existing states), mark the node as terminal
            if (node->children.empty()) {
                node->is_terminal = true;
            }
        } catch (const std::exception &e) {
            std::cerr << "Exception caught in expand_all_actions: " << e.what() << "\n";
        } catch (...) {
            std::cerr << "Non-standard exception caught in expand_all_actions.\n";
        }
    }


    node_type *expand_node(node_type *node, bool expand_all = false) {
        // Stop expanding if the node is terminal
        if (mfn_is_terminal(node->state)) {
            node->is_terminal = true;
            return nullptr;
        }

        // Generate possible actions
        auto possible_actions = mfn_get_actions(node->state);

        // If no actions are possible, mark the node as terminal
        if (possible_actions.empty()) {
            node->is_terminal = true;
            node->is_leaf = false;
            return nullptr;
        }


        // Filter the possible actions to exclude those that would lead to duplicate nodes
        // TODO: See if this is actually possible with LLVM 15/16. It's a bit of a hassle right now.

//        auto filtered_actions = possible_actions | std::views::filter([&](auto &action) {
//            State new_state = mfn_transition(node->state, action);
//            return std::none_of(node->children.begin(), node->children.end(),
//                                [&](auto &child) { return child->state == new_state; });
//        });
        std::vector<typename decltype(possible_actions)::value_type> filtered_actions;

        for (const auto &action: possible_actions) {
            State new_state = mfn_transition(node->state, action);
            bool is_duplicate = false;

            for (const auto &child: node->children) {
                if (child->state == new_state) {
                    is_duplicate = true;
                    break;
                }
            }

            if (!is_duplicate) {
                filtered_actions.push_back(action);
            }
        }

        // Convert the filtered range into a vector
        std::vector<Action> unique_actions(begin(filtered_actions), end(filtered_actions));

        // If all actions have been removed (i.e., all lead to existing states), mark the node as terminal
        if (unique_actions.empty()) {
            node->is_terminal = true;
            node->is_leaf = false;
            return nullptr;
        }

        // Create a uniform integer distribution to select one of the remaining possible actions at random
        std::uniform_int_distribution<int> dist(0, unique_actions.size() - 1);

        // Select one of the possible actions at random
        auto action = unique_actions[dist(engine)];

        State new_state = mfn_transition(node->state, action);

        RunningStats<Reward> reward_stats;
        auto child = std::make_unique<node_type>(new_state, action, node, reward_stats);
        //        child->is_leaf = true; // Nodes are instantiated as leaves by default.
        node->children.push_back(std::move(child));// take ownership of the child

        // Return the last child that was added, i.e., the one we just expanded
        return node->children.back().get();
    }

    node_type *best_child(node_type *node) const {
        return std::max_element(
                node->children.begin(), node->children.end(),
                [](const p_node &a, const p_node &b) {
                    return a->reward_stats.max() < b->reward_stats.max();
                })
                ->get();
    }

    trajectory_type get_best_trajectory(node_type *node) const {
        trajectory_type trajectory;

        while (node->children_size() != 0) {
            node = best_child(node);

            // Add the state, action, and result associated with the node to the trajectory.
            trajectory.emplace_back(node->state, node->action, node->reward_stats.max());
        }

        return trajectory;
    }

    // Initial root expansion
    //        if (this->root->children.empty()) {
    //            ODIN_LOG_INFO << "Root node children are empty. Expanding all actions.";
    //            expand_all_actions(this->root.get());
    //        }


    trajectory_type
    search(int iterations = 1e3, double seconds = -1.0, bool expand_all = false, bool contraction = true) {
        using namespace std::chrono;

        auto start_time = steady_clock::now();
        auto end_time = (seconds < 0.0)
                        ? time_point < steady_clock, duration < double >> (duration<double>::max())
                        : start_time + std::chrono::seconds(static_cast<int>(seconds));


        for (int i = 0; i < iterations && steady_clock::now() < end_time; ++i) {
            ODIN_VLOG(20) << "Iteration: " << iterations;

            // 1. Selection
            // Start from the root node and use the selection policy to descend through the tree,
            // until reaching a leaf node that hasn't been expanded yet.
            node_type *node = this->root.get();
            node->visit_count++;

            while (!node->is_leaf_node()) {
                node = select_node(node);
                node->visit_count++;
            }

            // 2. Expansion
            // For the selected leaf node, compute all possible actions from its state and add a
            // new child node for each action to the tree.
            auto new_node = expand_node(node, expand_all);

            // 3. Simulation
            // The child of the newly expanded node, simulate a trajectory from the child
            // node's state to a terminal state. Each state-action-result triplet in the trajectory
            // is added to the tree as a child node of the preceding state. This results in a chain
            // of nodes representing the simulated trajectory. The final state of the trajectory is
            // a terminal state, and the associated result is backpropagated up through the chain
            // of nodes.
            if (new_node) {
                // Perform a simulation from the child node's state and append the trajectory to the child node.
                const trajectory_type trajectory = simulate(new_node->state);

                // Expand the tree along the simulated trajectory.
                node_type *terminal_node = expand_trajectory(new_node, trajectory);

                terminal_node->is_terminal = true;

                // 4. Backpropagation
                // Backpropagate the result of the simulated trajectory up the tree through all
                // parent nodes, starting with the terminal node and up to the child node.
                backpropagate(terminal_node, terminal_node->reward_stats.max());
            }

            // Step 5: Contraction
            // Based on an extension to the four steps of MCTS (selection, expansion, simulation, backpropagation),
            // introduced by Hennes and Izzo in "Interplanetary Trajectory Planning with Monte Carlo Tree Search"
            if (contraction) {
                ODIN_VLOG(20) << "Contraction";
                // TODO: This contraction is removing the final leaf node of the trajectory. This is not correct. Fix.
                contract(this->root.get());
            }
        }

        // After time limit is reached, return the best node from the root according to a policy (e.g., the node with the highest value).
        return get_best_trajectory(this->root.get());
    }

    //    // Example function to save a single trajectory to a file.
    //    // This would need to be adapted based on how exactly you're representing trajectories and states.
    //    void save_trajectory_to_file(const trajectory& traj, const std::string& filename) {
    //        std::ofstream file(filename, std::ios::app); // append to the file
    //        for (const auto& [state, action, reward] : traj) {
    //            file << state << ',' << action << ',' << reward << '\n';
    //        }
    //    }
    //
    //    // Recursive function to save all admissible trajectories in a subtree to a file.
    //    void save_admissible_trajectories(node_type* node, const std::string& filename) {
    //        if (can_contract(node)) {
    //            // This node is the root of an admissible trajectory.
    //            // Retrieve the trajectory and save it.
    //            trajectory traj = get_trajectory(node);
    //            save_trajectory_to_file(traj, filename);
    //        }
    //
    //        // Recursively process all child nodes.
    //        for (auto& child : node->children) {
    //            save_admissible_trajectories(child.get(), filename);
    //        }
    //    }

private:
    // Expands the tree along a simulated trajectory by creating a chain of nodes.
    // The chain starts from a given node (initially an unexpanded leaf node) and
    // extends to a terminal node representing the end of a simulated trajectory.
    node_type *expand_trajectory(node_type *node, const trajectory_type &trajectory) {
        // Iterate over each state-action-result triplet in the simulated trajectory.
        for (const auto &[state, action, reward]: trajectory) {
            // Create a new node for each triplet.
            // The parent of each new node is set to the preceding node in the trajectory.
            // `std::make_unique` is used to create a new unique_ptr object that takes ownership of the new node.
            RunningStats<Reward> reward_stats;
            reward_stats.set_calc_max(true);
            reward_stats.push(reward);
            auto new_node = std::make_unique<node_type>(state, action, node, reward_stats);

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


    bool can_contract(node_type *node) {
        if (node->visit_count < 1) {
            // not all children have been visited
            return false;
        }

        if (!mfn_is_terminal(node->state) || node->is_leaf) {
            // this node is not a terminal node
            return false;
        }

        for (auto &child: node->children) {
            if (!can_contract(child.get())) {
                return false;
            }
        }

        // this node and all its children are terminal nodes and have been visited
        return true;
    }

    void contract(node_type *node) {
        auto it = node->children.begin();
        while (it != node->children.end()) {
            if (can_contract(it->get())) {
                it = node->children.erase(it);
            } else {
                contract(it->get());
                ++it;
            }
        }
    }

    bool all_children_terminal(node_type *node) {
        for (auto &child: node->children) {
            if (!mfn_is_terminal(child->state) || child->is_leaf) {
                return false;
            }
        }
        return true;
    }


    trajectory_type simulate(State state) {
        trajectory_type trajectory;
        while (!mfn_is_terminal(state)) {
            std::vector<Action> actions = mfn_get_actions(state);

            // If no actions are possible, mark the node as terminal
            if (actions.empty()) {
                ODIN_LOG_WARNING << "No actions possible from state during simulation.";
                break;
            }

            // Create a uniform distribution from 0 to size - 1
            std::uniform_int_distribution<size_t> distribution(0, actions.size() - 1);

            // Generate a random number using the distribution and the generator
            Action action = actions[distribution(engine)];
            state = mfn_transition(state, action);// Apply the action to get the new state
            Reward reward = mfn_reward(state);            // Compute the reward for the current state
            trajectory.emplace_back(state, action, reward);
        }
        return trajectory;
    };

    //    void backpropagate(node_type *node, Reward reward) {
    //        // Start from the given node and traverse up the tree
    //        while (node) {
    //            // Increment the visit count of the node
    //            node->visit_count++;// TODO: with the multiple leaf expansion, this will throw off the selection policies original intended behaviours
    //                                // .... actually it wont, we must just make sure not to have the expanded node visited for the number of parallel simulations
    //
    //            // Update the value of the node based on the result
    //            node->reward += reward;
    //
    //            // Move to the parent of the current node
    //            node = node->parent;
    //        }

    //    void backpropagate(node_type *node, Reward final_reward) {
    //        // Start from the given node and traverse up the tree
    //        while (node) {
    //            // Increment the visit count of the node
    //            node->visit_count++;
    //
    //            // Update the value of the node
    //            // Keep track of the max reward encountered so far
    //            node->reward = std::max(node->reward, final_reward);
    //
    //            // Move to the parent of the current node
    //            node = node->parent;
    //        }
    //    }
    void backpropagate(node_type *node, Reward final_reward) {
        // Start from the given node and traverse up the tree
        //        int         depth           = node->depth;// The depth of the current node
        //        const Float discount_factor = 0.99;       // The discount factor for future rewards
        while (node) {
            // Increment the visit count of the node
            node->visit_count++;

            // Update the value of the node
            node->reward_stats.push(final_reward);

            // Move to the parent of the current node
            node = node->parent;
        }
    }
};

#endif//MCTS_HPP