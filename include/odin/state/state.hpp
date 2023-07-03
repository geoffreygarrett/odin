/**
 * @file state.hpp
 * @brief This file contains various state classes.
 * @note Currently, state classes are designed for stack allocation, thus they do not use the EIGEN_MAKE_ALIGNED_OPERATOR_NEW macro.
 *       This macro ensures proper memory alignment for dynamic allocations of Eigen types, enabling efficient SIMD optimizations.
 *       However, it's only necessary if these classes are dynamically allocated with the new operator.
 *
 *       The use of the EIGEN_MAKE_ALIGNED_OPERATOR_NEW macro doesn't affect current usage, as instances are stack-allocated.
 *       Nevertheless, if dynamic allocation is introduced in the future, each class with fixed-size Eigen types will need
 *       to include the EIGEN_MAKE_ALIGNED_OPERATOR_NEW macro to ensure proper alignment and optimal performance.
 *
 * @todo Should dynamic allocation be required in the future, include EIGEN_MAKE_ALIGNED_OPERATOR_NEW in relevant classes.
 */
#ifndef STATE_HPP
#define STATE_HPP

#include <Eigen/Core>
#include <utility>

#define DEFAULT_SCALAR double
#define DEFAULT_DIMENSION 3

#include "state_types.hpp"
#include "state_handler.hpp"

#endif // STATE_HPP