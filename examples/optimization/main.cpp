#include <iostream>
#include <odin/core/optimization/mcts.hpp>
#include <vector>

int main() {
    auto get_actions = [](int state) {
        std::vector<int> actions;
        for (int i = 1; i <= 3; ++i) {
            if (state + i <= 5) {
                actions.push_back(i);
            }
        }
        return actions;
    };

    auto transition = [](int state, int action) {
        return state + action;
    };

    auto is_terminal = [](int state) {
        return state >= 5;
    };

    // Reward function modified to encourage reaching the goal faster
    auto reward = [](int state, bool terminal) {
        if (terminal && state == 5)
            return 100.0; // Large reward for reaching goal
        else if (state > 5)
            return -20.0; // Penalty for going over goal
        else
            return -10.0; // Small penalty for not reaching goal
    };

    MCTS<int, int> mcts(0, 0, get_actions, transition, is_terminal, reward);
    mcts.execute(10000);
    std::cout << "Best Action: " << mcts.best_action() << std::endl;
    mcts.print_tree();

    return 0;
}
