#ifndef ODIN_MCTS_SERIALIZATION_H
#define ODIN_MCTS_SERIALIZATION_H

#include "base.hpp"// Include the header file of your struct
#include "mcts.hpp"// Include the header file of your struct
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/common.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>

template<class Archive, typename Float>
void serialize(Archive &archive, RunningStats<Float> &rs) {
    archive(cereal::make_nvp("calc_mean", rs.m_calc_mean),
            cereal::make_nvp("calc_variance", rs.m_calc_variance),
            cereal::make_nvp("calc_min", rs.m_calc_min),
            cereal::make_nvp("calc_max", rs.m_calc_max),
            cereal::make_nvp("n", rs.m_n),
            cereal::make_nvp("min", rs.m_min),
            cereal::make_nvp("max", rs.m_max),
            cereal::make_nvp("old_mean", rs.m_old_mean),
            cereal::make_nvp("new_mean", rs.m_new_mean),
            cereal::make_nvp("old_sqr_diff", rs.m_old_sqr_diff),
            cereal::make_nvp("new_sqr_diff", rs.m_new_sqr_diff));
}

template<class Archive>
inline void serialize(Archive &archive, SearchMetrics &sm) {
    archive(cereal::make_nvp("search_seconds", sm.search_seconds),
            cereal::make_nvp("search_iterations", sm.search_iterations),
            cereal::make_nvp("fevals_transitions_evaluated", sm.fevals_transitions_evaluated),
            cereal::make_nvp("fevals_rewards_evaluated", sm.fevals_rewards_evaluated),
            cereal::make_nvp("fevals_actions_generated", sm.fevals_actions_generated),
            cereal::make_nvp("fevals_terminal_checks", sm.fevals_terminal_checks),
            cereal::make_nvp("fevals_selection_policy", sm.fevals_selection_policy));
}

#endif// ODIN_MCTS_SERIALIZATION_H