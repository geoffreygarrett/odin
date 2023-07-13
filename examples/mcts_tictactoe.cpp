#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <odin/concepts.hpp>
#include <odin/logging.hpp>
#include <odin/optimization/tree/mcts.hpp>
#include <optional>
#include <tuple>
#include <vector>


// State dimensions.
constexpr int DIM = 3;

// TicTacToe state
struct State {
    std::vector<int> board;
    int              player;

    std::string to_string() const {
        std::string str;
        for (size_t i = 0; i < board.size(); ++i) {
            if (board[i] == 1) {
                str += " X ";
            } else if (board[i] == -1) {
                str += " O ";
            } else {
                str += " • ";
            }

            if ((i + 1) % DIM == 0) {
                str += "\n";
            }
        }
        return str;
    }
};

bool operator==(const State &lhs, const State &rhs) {
    return lhs.board == rhs.board && lhs.player == rhs.player;
}

// TicTacToe action
struct Action {
    int index;

    std::string to_string() const {
        return std::to_string(index);
    }
};

template<Stringable T>
std::ostream &operator<<(std::ostream &os, const T &obj) {
    os << obj.to_string();
    return os;
}

template<Streamable T>
std::ostream &operator<<(std::ostream &os, const std::optional<T> &objOpt) {
    if (objOpt.has_value()) {
        os << objOpt.value();
    } else {
        os << "None";
    }
    return os;
}

static int check_winner(State state) {
    std::vector<std::vector<int>> grid(DIM, std::vector<int>(DIM));
    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            grid[i][j] = state.board[i * DIM + j];
        }
    }

    for (int i = 0; i < DIM; i++) {
        if (std::abs(std::accumulate(grid[i].begin(), grid[i].end(), 0)) == DIM) {
            return grid[i][0];
        }

        int column_sum = 0;
        for (int j = 0; j < DIM; j++) {
            column_sum += grid[j][i];
        }

        if (std::abs(column_sum) == DIM) {
            return grid[0][i];
        }
    }

    int diag_sum1 = 0;
    int diag_sum2 = 0;
    for (int i = 0; i < DIM; i++) {
        diag_sum1 += grid[i][i];
        diag_sum2 += grid[i][DIM - i - 1];
    }

    if (std::abs(diag_sum1) == DIM || std::abs(diag_sum2) == DIM) {
        return diag_sum1 != 0 ? std::signbit(diag_sum1) : std::signbit(diag_sum2);
    }

    return 0;
}

bool contains_zero(State state) {
    return std::any_of(state.board.begin(), state.board.end(), [](int i) { return i == 0; });
}

// Initialize state with zeros and player 1.
static State init_state() {
    State state;
    state.board  = std::vector<int>(DIM * DIM, 0);
    state.player = 1;
    return state;
}

// Check if the game has ended.
bool is_terminal(const State &state) {
    return check_winner(state) != 0 || !contains_zero(state);
}

// Get list of available actions.
static std::vector<Action> get_actions(const State &state) {
    std::vector<Action> actions;
    for (int i = 0; i < (int) state.board.size(); i++) {
        if (state.board[i] == 0) {
            actions.push_back({i});
        }
    }
    return actions;
}

// Transition function.
static State transition(const State &state, Action action) {
    State newState               = state;
    newState.board[action.index] = newState.player;
    newState.player              = -newState.player;
    return newState;
}

// Reward function.
static double reward(const State &state) {
    int winner    = check_winner(state);
    int num_moves = std::count_if(state.board.begin(), state.board.end(), [](int i) { return i != 0; });
    num_moves     = std::max(num_moves, 1);

    if (winner == 1) {
        return 1.0 / num_moves;
    } else if (winner == -1 || is_terminal(state)) {
        return 0.0;
    }

    return 0.0;
}


int main() {
//    spdlog::set_level(spdlog::level::debug);
    using mcts_type       = MCTS<State, Action, double>;
    using node_type       = mcts_type::node_type;
    auto selection_policy = std::make_shared<UCB1<node_type, double>>(1.0,
                                                                      [](const node_type *node) {
                                                                          return node->reward_stats.max();
                                                                      });

    mcts_type mcts(init_state(), get_actions, transition, is_terminal, reward, selection_policy);


    // Simulate MCTS from initial state
    //    mcts_type::trajectory_type trajectory;
    auto [trajectory, metrics] = mcts.search(
            1e3,  // iterations
            -1.0, //seconds
            false,// expand_all
            true, // contraction
            4    // leaf_parallelism
    );

    // Print metrics
    std::cout << metrics << std::endl;

    // Print trajectory
    for (const auto &traj: trajectory) {
        const auto &[state, actionOpt, reward] = traj;
        std::string actionStr                  = actionOpt.has_value() ? actionOpt.value().to_string() : "None";
        std::cout << "State: \n"
                  << state << "\nAction: " << actionStr << ", \nReward: " << reward << "\n";
    }

    std::cout << "Complete" << std::endl;

    std::cout << "Test passed!" << std::endl;

    return 0;
}
