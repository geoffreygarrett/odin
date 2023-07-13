#ifndef MCTS_HPP
#define MCTS_HPP

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/common.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#include "base.hpp"
#include <algorithm>
#include <atomic>
#include <complex>
#include <functional>
#include <memory>
#include <odin/logging.hpp>
#include <oneapi/tbb/concurrent_vector.h>
#include <optional>
#include <random>
#include <stack>
#include <tbb/tbb.h>
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
        return this->mfn_value_estimator(child) + m_cp * std::sqrt(std::log(parent->visit_count) / child->visit_count * std::min(0.25, Vi));
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

struct SearchMetrics {
    double_t search_seconds;
    size_t   search_iterations;
    size_t   fevals_transitions_evaluated;
    size_t   fevals_rewards_evaluated;
    size_t   fevals_actions_generated;
    size_t   fevals_terminal_checks;
    size_t   fevals_selection_policy;

    friend std::ostream &operator<<(std::ostream &os, const SearchMetrics &metrics) {
        os << "SearchMetrics(";
        os << "search_seconds: " << metrics.search_seconds;
        os << ", search_iterations: " << metrics.search_iterations;
        os << ", fevals_transitions_evaluated: " << metrics.fevals_transitions_evaluated;
        os << ", fevals_rewards_evaluated: " << metrics.fevals_rewards_evaluated;
        os << ", fevals_actions_generated: " << metrics.fevals_actions_generated;
        os << ", fevals_terminal_checks: " << metrics.fevals_terminal_checks;
        os << ", fevals_selection_policy: " << metrics.fevals_selection_policy;
        os << ")";
        return os;
    }

    bool operator==(const SearchMetrics &other) const {
        return (search_seconds == other.search_seconds)
            && (search_iterations == other.search_iterations)
            && (fevals_transitions_evaluated == other.fevals_transitions_evaluated)
            && (fevals_rewards_evaluated == other.fevals_rewards_evaluated)
            && (fevals_actions_generated == other.fevals_actions_generated)
            && (fevals_terminal_checks == other.fevals_terminal_checks)
            && (fevals_selection_policy == other.fevals_selection_policy);
    }

    template<class Archive>
    inline void serialize(Archive &archive) {
        archive(CEREAL_NVP(search_seconds),
                CEREAL_NVP(search_iterations),
                CEREAL_NVP(fevals_transitions_evaluated),
                CEREAL_NVP(fevals_rewards_evaluated),
                CEREAL_NVP(fevals_actions_generated),
                CEREAL_NVP(fevals_terminal_checks),
                CEREAL_NVP(fevals_selection_policy));
    }
};

template<typename State, typename Action, typename Reward, typename Float = double>
class MCTS : public MCTSBase<MCTS<State, Action, Reward, Float>, State, Action, Reward> {
public:
    using Base = MCTSBase<MCTS<State, Action, Reward, Float>, State, Action, Reward>;
    //    using Base::Base;
    using p_node              = typename Base::p_node;
    using node_type           = typename Base::node_type;
    using trajectory_type     = std::vector<std::tuple<State, std::optional<Action>, Reward>>;
    using fn_state_transition = std::function<State(State, Action)>;
    using fn_action_generator = std::function<std::vector<Action>(State &)>;
    using fn_selection_policy = std::shared_ptr<SelectionPolicy<node_type, Float>>;
    using fn_is_terminal      = std::function<bool(State &)>;
    using fn_reward           = std::function<Reward(State &)>;

    std::atomic<size_t> fevals_transitions_evaluated;
    std::atomic<size_t> fevals_rewards_evaluated;
    std::atomic<size_t> fevals_actions_generated;
    std::atomic<size_t> fevals_terminal_checks;
    std::atomic<size_t> fevals_selection_policy;

    fn_action_generator mfn_get_actions;
    fn_state_transition mfn_transition;
    fn_is_terminal      mfn_is_terminal;
    fn_reward           mfn_reward;
    fn_selection_policy mfn_selection_policy;

    using search_result_type = std::tuple<trajectory_type, SearchMetrics>;
    std::mt19937 engine;


    struct SimulationRunner {
        std::mt19937                             local_engine;
        State                                    m_initial_state;
        tbb::concurrent_vector<trajectory_type> &trajectories;
        fn_action_generator                      mfn_get_actions;
        fn_state_transition                      mfn_transition;
        fn_reward                                mfn_reward;
        fn_is_terminal                           mfn_is_terminal;

        SimulationRunner(
                State                                   &initial_state,
                tbb::concurrent_vector<trajectory_type> &trajectories,
                fn_action_generator                      mfn_get_actions,
                fn_state_transition                      mfn_transition,
                fn_reward                                mfn_reward,
                fn_is_terminal                           mfn_is_terminal)
            : local_engine(std::random_device{}()),
              trajectories(trajectories),
              mfn_get_actions(mfn_get_actions),
              mfn_transition(mfn_transition),
              mfn_reward(mfn_reward),
              mfn_is_terminal(mfn_is_terminal),
              m_initial_state(initial_state) {}

        void operator()(int j) {
            SPDLOG_DEBUG("SimulationRunner: Starting simulation for thread {}", j);
            trajectory_type trajectory;
            auto            state = m_initial_state;
            while (!mfn_is_terminal(state)) {
                std::vector<Action> actions = mfn_get_actions(state);

                // If no actions are possible, mark the node as terminal
                if (actions.empty()) {
                    SPDLOG_DEBUG("No actions possible from state during simulation in thread {}", j);
                    break;
                }
                // Create a uniform distribution from 0 to size - 1
                std::uniform_int_distribution<size_t> distribution(0, actions.size() - 1);

                // Generate a random number using the distribution and the generator
                Action action = actions[distribution(local_engine)];
                state         = mfn_transition(state, action);
                Reward reward = mfn_reward(state);
                trajectory.emplace_back(state, action, reward);
            }
            trajectories[j] = trajectory;
        }
    };

    MCTS(State                 initial_state,
         fn_action_generator   get_actions,
         fn_state_transition   transition,
         fn_is_terminal        is_terminal,
         fn_reward             reward,
         fn_selection_policy   selection_policy = std::make_shared<UCB1<node_type, Float>>(1.0,
                                                                                         [](const node_type *node) {
                                                                                             return node->reward_stats.max();
                                                                                         }),
         size_t                seed             = std::random_device{}(),
         std::optional<Action> initial_action   = std::nullopt)// actions are generally assigned by the environment through the process
        : Base(initial_state, std::nullopt),
          mfn_get_actions([this, get_actions](State s) {
              fevals_actions_generated++;
              return get_actions(s);
          }),
          mfn_transition([this, transition](State s, Action a) {
              fevals_transitions_evaluated++;
              return transition(s, a);
          }),
          mfn_is_terminal([this, is_terminal](State s) {
              fevals_terminal_checks++;
              return is_terminal(s);
          }),
          mfn_reward([this, reward](State s) {
              fevals_rewards_evaluated++;
              return reward(s);
          }),
          mfn_selection_policy(selection_policy),
          engine(seed),
          fevals_transitions_evaluated(0),
          fevals_rewards_evaluated(0),
          fevals_actions_generated(0),
          fevals_terminal_checks(0),
          fevals_selection_policy(0) {}

    //    // default constructor
    //    MCTS() : Base(State(), std::nullopt),
    //             engine(std::random_device{}()),
    //             fevals_transitions_evaluated(0),
    //             fevals_rewards_evaluated(0),
    //             fevals_actions_generated(0),
    //             fevals_terminal_checks(0),
    //             fevals_selection_policy(0) {}

    //    // constructor without the functions, for serialization
    //    explicit MCTS(State initial_state, size_t seed = std::random_device{}())
    //        : Base(initial_state, std::nullopt),
    //          engine(seed),
    //          fevals_transitions_evaluated(0),
    //          fevals_rewards_evaluated(0),
    //          fevals_actions_generated(0),
    //          fevals_terminal_checks(0),
    //          fevals_selection_policy(0) {}

    bool is_ready() {
        return (mfn_get_actions != nullptr)
            && (mfn_transition != nullptr)
            && (mfn_is_terminal != nullptr)
            && (mfn_reward != nullptr)
            && (mfn_selection_policy != nullptr);
    }

//    explicit MCTS(std::shared_ptr<node_type> root) : Base(root),
//                                                     engine(std::random_device{}()),
//                                                     fevals_transitions_evaluated(0),
//                                                     fevals_rewards_evaluated(0),
//                                                     fevals_actions_generated(0),
//                                                     fevals_terminal_checks(0),
//                                                     fevals_selection_policy(0) {}

    explicit MCTS(node_type *node) : Base(node),
                                     engine(std::random_device{}()),
                                     fevals_transitions_evaluated(0),
                                     fevals_rewards_evaluated(0),
                                     fevals_actions_generated(0),
                                     fevals_terminal_checks(0),
                                     fevals_selection_policy(0) {}

    MCTS(MCTS &&other) noexcept
        : Base(other.root.get()),
          fevals_transitions_evaluated(other.fevals_transitions_evaluated.load()),
          fevals_rewards_evaluated(other.fevals_rewards_evaluated.load()),
          fevals_actions_generated(other.fevals_actions_generated.load()),
          fevals_terminal_checks(other.fevals_terminal_checks.load()),
          fevals_selection_policy(other.fevals_selection_policy.load()),
          mfn_get_actions(std::move(other.mfn_get_actions)),
          mfn_transition(std::move(other.mfn_transition)),
          mfn_is_terminal(std::move(other.mfn_is_terminal)),
          mfn_reward(std::move(other.mfn_reward)),
          mfn_selection_policy(std::move(other.mfn_selection_policy)),
          engine(std::move(other.engine)) {}
    //
    //    MCTS &operator=(MCTS &&other) noexcept {
    //        Base::operator=(std::move(other));
    //        fevals_transitions_evaluated = other.fevals_transitions_evaluated.load();
    //        fevals_rewards_evaluated     = other.fevals_rewards_evaluated.load();
    //        fevals_actions_generated     = other.fevals_actions_generated.load();
    //        fevals_terminal_checks       = other.fevals_terminal_checks.load();
    //        fevals_selection_policy      = other.fevals_selection_policy.load();
    //        mfn_get_actions              = std::move(other.mfn_get_actions);
    //        mfn_transition               = std::move(other.mfn_transition);
    //        mfn_is_terminal              = std::move(other.mfn_is_terminal);
    //        mfn_reward                   = std::move(other.mfn_reward);
    //        mfn_selection_policy         = std::move(other.mfn_selection_policy);
    //        engine                       = std::move(other.engine);
    //        return *this;
    //    }

    void set_fn_is_terminal(fn_is_terminal is_terminal) {
        mfn_is_terminal = is_terminal;
    }

    void set_fn_reward(fn_reward reward) {
        mfn_reward = reward;
    }

    void set_fn_transition(fn_state_transition transition) {
        mfn_transition = transition;
    }

    void set_fn_get_actions(fn_action_generator get_actions) {
        mfn_get_actions = get_actions;
    }

    void set_fn_selection_policy(fn_selection_policy selection_policy) {
        mfn_selection_policy = selection_policy;
    }

    node_type *get_root() {
        return Base::root.get();
    }

    void set_root(node_type *root) {
        Base::root.reset(root);
    }

    node_type *select_node(node_type *node) {
        fevals_selection_policy++;
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
        spdlog::debug("Expanding node");
        //        std::cout << node->state << std::endl;
        // Stop expanding if the node is terminal
        if (mfn_is_terminal(node->state)) {
            spdlog::debug("Expanding node");

            node->is_terminal = true;
            return nullptr;
        }
        spdlog::debug("Node is not terminal.");
        // Generate possible actions
        auto possible_actions = mfn_get_actions(node->state);

        // If no actions are possible, mark the node as terminal
        if (possible_actions.empty()) {
            node->is_terminal = true;
            node->is_leaf     = false;
            return nullptr;
        }
        spdlog::debug("Possible actions: {}", possible_actions.size());

        // Filter the possible actions to exclude those that would lead to duplicate nodes
        // TODO: See if this is actually possible with LLVM 15/16. It's a bit of a hassle right now.

        //        auto filtered_actions = possible_actions | std::views::filter([&](auto &action) {
        //            State new_state = mfn_transition(node->state, action);
        //            return std::none_of(node->children.begin(), node->children.end(),
        //                                [&](auto &child) { return child->state == new_state; });
        //        });
        std::vector<typename decltype(possible_actions)::value_type> filtered_actions;

        for (const auto &action: possible_actions) {
            State new_state    = mfn_transition(node->state, action);
            bool  is_duplicate = false;

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
            node->is_leaf     = false;
            return nullptr;
        }

        // Create a uniform integer distribution to select one of the remaining possible actions at random
        std::uniform_int_distribution<int> dist(0, unique_actions.size() - 1);

        // Select one of the possible actions at random
        auto action = unique_actions[dist(engine)];

        State new_state = mfn_transition(node->state, action);

        RunningStats<Reward> reward_stats;
        auto                 child = std::make_unique<node_type>(new_state, action, node, reward_stats);
        //        child->is_leaf = true; // Nodes are instantiated as leaves by default.
        node->children.push_back(std::move(child));// take ownership of the child

        // Return the last child that was added, i.e., the one we just expanded
        return node->children.back().get();
    }

    node_type *best_child(node_type *node) const {
        return std::max_element(node->children.begin(),
                                node->children.end(),
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

    search_result_type search(
            int    iterations       = 1e3,
            double seconds          = -1.0,
            bool   expand_all       = false,
            bool   contraction      = true,
            int    leaf_parallelism = 4) {
        using namespace std::chrono;

        if (!is_ready()) {
            throw std::runtime_error("The search tree is not ready. You need the state-action-reward (SAR) functions to be defined.");
        }

        auto start_time = steady_clock::now();
        auto end_time   = (seconds < 0.0)
                                ? time_point<steady_clock, duration<double>>(duration<double>::max())
                                : start_time + std::chrono::seconds(static_cast<int>(seconds));


        for (int i = 0; i < iterations && steady_clock::now() < end_time; ++i) {
            spdlog::debug("Starting iteration: {}", i);

            // Phase 1: Selection
            // Descend from the root node to a leaf node that has not been expanded yet.
            node_type *node = this->root.get();
            node->visit_count++;

            while (!node->is_leaf_node()) {
                node = select_node(node);
                node->visit_count++;
            }
            spdlog::debug("Selected node for expansion");

            // PHASE 2: Expansion
            // Add a new child node for each possible action from the state of the selected leaf node.
            auto new_node = expand_node(node, expand_all);
            if (new_node) {
                spdlog::debug("Expanded node");

                // PHASE 3: Simulation
                tbb::concurrent_vector<trajectory_type> trajectories(leaf_parallelism);
                std::vector<SimulationRunner>           runners;

                for (int j = 0; j < leaf_parallelism; ++j) {
                    runners.emplace_back(new_node->state,
                                         trajectories,
                                         mfn_get_actions,
                                         mfn_transition,
                                         mfn_reward,
                                         mfn_is_terminal);
                }

                tbb::parallel_for(0, leaf_parallelism, [&](int j) {
                    spdlog::debug("(leaf_parallelism) START    leaf {}/{}", j + 1, leaf_parallelism);
                    runners[j](j);
                    spdlog::debug("(leaf_parallelism) COMPLETE leaf {}/{}", j + 1, leaf_parallelism);
                });

                // For each simulated trajectory, add a chain of nodes to the tree.
                std::vector<node_type *> terminal_nodes(leaf_parallelism);
                for (int j = 0; j < leaf_parallelism; ++j) {
                    terminal_nodes[j]              = expand_trajectory(new_node, trajectories[j]);
                    terminal_nodes[j]->is_terminal = true;
                    spdlog::debug("Expanded trajectory for thread {}", j);
                }

                // PHASE 4: Backpropagation
                for (int j = 0; j < leaf_parallelism; ++j) {
                    backpropagate(terminal_nodes[j], terminal_nodes[j]->reward_stats.max());
                    spdlog::debug("Backpropagated result for thread {}", j);
                }

            } else {
                spdlog::debug("Skipped simulation and backpropagation because no new node was expanded");
            }
            // Step 5: Contraction
            // Based on an extension to the four steps of MCTS (selection, expansion, simulation, backpropagation),
            // introduced by Hennes and Izzo in "Interplanetary Trajectory Planning with Monte Carlo Tree Search"
            if (contraction) {
                spdlog::debug("Contracting tree");
                // TODO: This contraction is removing the final leaf node of the trajectory. This is not correct. Fix.
                contract(this->root.get());
            }
        }

        // After time limit is reached, return the best node from the root according to a policy (e.g., the node with the highest value).
        SearchMetrics metrics{};
        metrics.search_iterations            = iterations;
        metrics.search_seconds               = duration_cast<duration<double>>(steady_clock::now() - start_time).count();
        metrics.fevals_actions_generated     = this->fevals_actions_generated;
        metrics.fevals_rewards_evaluated     = this->fevals_rewards_evaluated;
        metrics.fevals_selection_policy      = this->fevals_selection_policy;
        metrics.fevals_terminal_checks       = this->fevals_terminal_checks;
        metrics.fevals_transitions_evaluated = this->fevals_transitions_evaluated;
        return std::make_tuple(get_best_trajectory(this->root.get()), metrics);
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

    // Serialization function definition
    template<class Archive>
    void serialize_derived(Archive &ar) {
        //                ar(CEREAL_NVP(fevals_transitions_evaluated),
        //                   CEREAL_NVP(fevals_rewards_evaluated),
        //                   CEREAL_NVP(fevals_actions_generated),
        //                   CEREAL_NVP(fevals_terminal_checks),
        //                   CEREAL_NVP(fevals_selection_policy));
    }

    //    template<class Archive>
    //    static void
    //    load_and_construct_derived(Archive &ar, cereal::construct<MCTS<State, Action, Reward, Float>> &construct) {
    //        std::atomic<size_t> fevals_transitions_evaluated;
    //        std::atomic<size_t> fevals_rewards_evaluated;
    //        std::atomic<size_t> fevals_actions_generated;
    //        std::atomic<size_t> fevals_terminal_checks;
    //        std::atomic<size_t> fevals_selection_policy;
    //
    //        ar(CEREAL_NVP(fevals_transitions_evaluated),
    //           CEREAL_NVP(fevals_rewards_evaluated),
    //           CEREAL_NVP(fevals_actions_generated),
    //           CEREAL_NVP(fevals_terminal_checks),
    //           CEREAL_NVP(fevals_selection_policy));
    //
    //        construct(fevals_transitions_evaluated,
    //                  fevals_rewards_evaluated,
    //                  fevals_actions_generated,
    //                  fevals_terminal_checks,
    //                  fevals_selection_policy);
    //    }

private:
    // Expands the tree along a simulated trajectory by creating a chain of nodes.
    // The chain starts from a given node (initially an unexpanded leaf node) and
    // extends to a terminal node representing the end of a simulated trajectory.
    node_type *
    expand_trajectory(node_type *node, const trajectory_type &trajectory) {
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

    template<typename SimulationRunner>
    void simulate(SimulationRunner runner, int leaf_parallelism) {
        tbb::parallel_for(0, leaf_parallelism, runner);
    }


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