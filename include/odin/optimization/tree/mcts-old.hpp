#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <unordered_set>
#include <vector>

template<typename T>
std::string to_string(const T &value);

template<>
std::string to_string<int>(const int &value);

template<>
std::string to_string<double>(const double &value);

template<>
std::string to_string<std::vector<int>>(const std::vector<int> &value);

template<>
std::string to_string<std::vector<double>>(const std::vector<double> &value);

template<typename T>
std::string vector_to_string(const std::vector<T> &value) {
    std::ostringstream oss;
    for (const auto &v: value) {
        oss << to_string(v) << ' ';
    }
    return oss.str();
}

template<>
std::string to_string<int>(const int &value) {
    return std::to_string(value);
}

template<>
std::string to_string<double>(const double &value) {
    return std::to_string(value);
}

template<>
std::string to_string<std::vector<int>>(const std::vector<int> &value) {
    return vector_to_string(value);
}

template<>
std::string to_string<std::vector<double>>(const std::vector<double> &value) {
    return vector_to_string(value);
}

template<>
std::string to_string<Eigen::VectorX<int>>(const Eigen::VectorX<int> &value) {
    std::ostringstream oss;
    oss << value.transpose();
    return oss.str();
}

template<>
std::string to_string<Eigen::VectorXd>(const Eigen::VectorXd &value) {
    std::ostringstream oss;
    oss << value.transpose();
    return oss.str();
}

template<typename T>
std::string to_string(const T &value) {
    return value.to_string();
}

template<typename T>
std::string to_string(const std::vector<T> &value) {
    return vector_to_string(value);
}


template<typename State, typename Action, typename Float = double>
struct Node {
    std::vector<std::shared_ptr<Node<State, Action, Float>>> children;     // Child nodes
    std::shared_ptr<Node<State, Action, Float>>              parent;       // Parent node
    Float                                                    total_score;  // Win score for the node
    int                                                      visit_count;  // Number of times node has been visited
    bool                                                     move;         // "Player" who made the move
    State                                                    state;        // State represented by this node
    Action                                                   action;       // Action that resulted in this state
    bool                                                     is_terminal{};// Whether this node is terminal

    Node() : parent(nullptr), total_score(0), visit_count(0), move(false) {}
    Node(std::shared_ptr<Node<State, Action, Float>> parent, State state, Action action, bool move)
        : parent(parent),
          total_score(0),
          visit_count(0),
          move(move),
          state(state),
          action(action) {}
    void print(const std::string &prefix = "", bool is_tail = true, int max_depth = -1, bool best_branch_only = false, bool is_root = true) const {
        if (max_depth == 0) {
            return;
        }


        std::cout << prefix << (is_tail ? "└── " : "├── ")
                  << "State: " << to_string(state)
                  << ", Action: " << to_string(action)
                  << ", Score: " << total_score
                  << ", Visits: " << visit_count
                  << ", Children: " << children.size()
                  << "\n";

        if (!children.empty()) {
            std::string new_prefix = prefix + (is_tail ? "    " : "│   ");

            auto best_child_it = std::max_element(children.begin(), children.end(),
                                                  [](const auto &a, const auto &b) {
                                                      return a->total_score < b->total_score;
                                                  });

            for (size_t i = 0; i < children.size(); ++i) {
                const auto &child         = children[i];
                bool        is_child_tail = (i == children.size() - 1);

                if (is_root) {
                    if (child == *best_child_it) {
                        // this is the best child, we print it fully
                        child->print(new_prefix, is_child_tail, max_depth - 1, false, false);
                    } else {
                        // we only print direct children of root, not their children
                        child->print(new_prefix, is_child_tail, 1, false, false);
                    }
                } else if (child == *best_child_it) {
                    // this is the best branch, we print it fully
                    child->print(new_prefix, is_child_tail, max_depth - 1, false, false);
                }
            }
        }
    }
};

template<typename State, typename Action, typename Float = double>
class MCTS {
private:
    std::shared_ptr<Node<State, Action, Float>> root;// Use shared_ptr here.

    std::function<std::vector<Action>(State)> get_actions;// Function to get available actions for a state
    std::function<State(State, Action)>       transition; // Function to get new state after an action
    std::function<bool(State)>                is_terminal;// Function to check if a state is terminal
    std::function<Float(State, bool)>         reward;     // Function to get the reward of a state

    std::shared_ptr<Node<State, Action, Float>> select_node(std::shared_ptr<Node<State, Action, Float>> node) {
        std::shared_ptr<Node<State, Action, Float>> selected_node = node;
        Float                                       max_ucb       = std::numeric_limits<Float>::lowest();

        for (auto child: node->children) {
            Float ucb = 0;
            if (child->visit_count == 0) {
                ucb = std::numeric_limits<Float>::max();
            } else {
                ucb = (child->total_score / child->visit_count) + std::sqrt(2.0 * std::log(node->visit_count) / child->visit_count);
            }
            if (ucb > max_ucb) {
                max_ucb       = ucb;
                selected_node = child;
            }
        }

        return selected_node;
    }

    void expand_node(std::shared_ptr<Node<State, Action, Float>> node) {
        if (node->is_terminal) {
            return;
        }
        std::vector<Action> possible_actions      = get_actions(node->state);
        bool                all_children_terminal = true;// Assume initially all children are terminal
        for (auto action: possible_actions) {
            State new_state = transition(node->state, action);

            // Check if a child with the same state and action already exists
            bool duplicate = false;
            for (auto &child: node->children) {
                if (child->state == new_state && child->action == action) {
                    duplicate = true;
                    break;
                }
            }

            // If no such child exists, create a new one
            if (!duplicate) {
                auto child = std::make_shared<Node<State, Action, Float>>(node, new_state, action, !node->move);
                // Check if the state of the child is terminal
                if (is_terminal(new_state)) {
                    child->is_terminal = true;
                } else {
                    all_children_terminal = false;
                }
                node->children.push_back(child);
            }
        }
        // If all children are terminal, set parent node as terminal too
        if (all_children_terminal) {
            node->is_terminal = true;
        }
    }

    Float simulate(std::shared_ptr<Node<State, Action, Float>> node) {
        static std::mt19937 gen(std::random_device{}());// Create the generator once

        State temp_state   = node->state;
        bool  current_move = node->move;

        while (!is_terminal(temp_state)) {
            std::vector<Action>             possible_actions = get_actions(temp_state);
            std::uniform_int_distribution<> dis(0, possible_actions.size() - 1);
            Action                          action = possible_actions[dis(gen)];
            temp_state                             = transition(temp_state, action);
            current_move                           = !current_move;
        }

        return reward(temp_state, !current_move);// Reward for the move that led to the terminal state
    }

    void backpropagate(std::shared_ptr<Node<State, Action, Float>> node, Float reward) {
        while (node) {
            ++(node->visit_count);
            node->total_score += reward;
            node = node->parent;
        }
    }


public:
    MCTS(State initial_state, Action initial_action,
         std::function<std::vector<Action>(State)> get_actions,
         std::function<State(State, Action)>       transition,
         std::function<bool(State)>                is_terminal,
         std::function<Float(State, bool)>         reward)
        : root(std::make_shared<Node<State, Action, Float>>(nullptr, initial_state, initial_action, false)),
          get_actions(get_actions),
          transition(transition),
          is_terminal(is_terminal),
          reward(reward) {}

    Action best_action() {
        std::shared_ptr<Node<State, Action, Float>> best_node = nullptr;
        Float                                       max_score = std::numeric_limits<Float>::lowest();

        for (auto child: root->children) {
            if (child->total_score > max_score) {
                max_score = child->total_score;
                best_node = child;
            }
        }

        return best_node->action;
    }

    std::shared_ptr<Node<State, Action, Float>> select_newly_expanded_node(std::shared_ptr<Node<State, Action, Float>> node) {
        std::random_device              rd;
        std::mt19937                    gen(rd());
        std::uniform_int_distribution<> dis(0, node->children.size() - 1);
        return node->children[dis(gen)];
    }

    void print_tree(int max_depth = -1, bool best_branch_only = false) const {
        root->print("", true, max_depth, best_branch_only);
    }

    void execute(int iterations) {
        for (int i = 0; i < iterations; ++i) {
            std::shared_ptr<Node<State, Action, Float>> node = root;
            // Keep selecting and expanding nodes until a terminal node is reached
            while (!is_terminal(node->state)) {
                node = select_node(node);
                expand_node(node);
            }
            node->is_terminal  = true;
            Float reward_value = simulate(node);
            backpropagate(node, reward_value);
        }
    }
};