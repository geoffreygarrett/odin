#ifndef ODIN_GRID_HPP
#define ODIN_GRID_HPP

#include "parallel_concepts.hpp"
#include <Eigen/Dense>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/blocked_range3d.h>
#include <tbb/parallel_for.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace odin {

    template<int InputDims, int OutputDims = 1>
    struct GridCalculator;

    // Default for any dimension
    template<int InputDims, int OutputDims>
    struct GridCalculator {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                               &model,
                const Eigen::Tensor<IType, InputDims> &input,
                Func                                 &&calculate_fn,
                Eigen::Tensor<OType, OutputDims>      &output) {

            const auto &dims = input.dimensions();

            tbb::parallel_for(
                    tbb::blocked_range<int>(0, input.size()),
                    [&](const tbb::blocked_range<int> &r) {
                        thread_local Derived tl_model = get_thread_local_instance(model);
                        for (int i = r.begin(); i != r.end(); ++i) {
                            Eigen::array<int, InputDims> indices;
                            int                          index = i;
                            for (int d = InputDims - 1; d >= 0; --d) {
                                indices[d] = index % dims[d];
                                index /= dims[d];
                            }

                            auto var_set = input(indices);
                            auto result  = calculate_fn(tl_model, var_set);

                            // Assuming output has same size and dimensions
                            if constexpr (OutputDims == 1)
                                output(i) = result;
                            else
                                output(indices) = result;
                        }
                    });
        }
    };

    // Specialization for 2D grid
    template<int OutputDims>
    struct GridCalculator<2, OutputDims> {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                          &model,
                const Eigen::Tensor<IType, 2>    &input,
                Func                            &&calculate_fn,
                Eigen::Tensor<OType, OutputDims> &output) {

            const auto &dims = input.dimensions();
            tbb::parallel_for(
                    tbb::blocked_range2d<int>(0, dims[0], 0, dims[1]),
                    [&](const tbb::blocked_range2d<int> &r) {
                        thread_local Derived tl_model = get_thread_local_instance(model);
                        for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
                            for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
                                Eigen::array<int, 2> indices{i, j};
                                auto                 var_set = input(indices);
                                auto                 result  = calculate_fn(tl_model, var_set);
                                if constexpr (OutputDims == 1)
                                    output(i, j) = result;
                                else
                                    output(indices) = result;
                            }
                        }
                    });
        }
    };

    // Specialization for 3D grid
    template<int OutputDims>
    struct GridCalculator<3, OutputDims> {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                          &model,
                const Eigen::Tensor<IType, 3>    &input,
                Func                            &&calculate_fn,
                Eigen::Tensor<OType, OutputDims> &output) {

            const auto &dims = input.dimensions();
            tbb::parallel_for(
                    tbb::blocked_range3d<int>(0, dims[0], 0, dims[1], 0, dims[2]),
                    [&](const tbb::blocked_range3d<int> &r) {
                        thread_local Derived tl_model = get_thread_local_instance(model);
                        for (int i = r.pages().begin(); i != r.pages().end(); ++i) {
                            for (int j = r.rows().begin(); j != r.rows().end(); ++j) {
                                for (int k = r.cols().begin(); k != r.cols().end(); ++k) {
                                    Eigen::array<int, 3> indices{i, j, k};
                                    auto                 var_set = input(indices);
                                    auto                 result  = calculate_fn(tl_model, var_set);
                                    if constexpr (OutputDims == 1)
                                        output(i, j, k) = result;
                                    else
                                        output(indices) = result;
                                }
                            }
                        }
                    });
        }
    };

    template<ThreadSafe Derived, typename Func, typename IType, typename OType, int InputDims, int OutputDims>
    void calculate_on_grid(
            Derived                               &model,
            const Eigen::Tensor<IType, InputDims> &input,
            Func                                 &&calculate_fn,
            Eigen::Tensor<OType, OutputDims>      &output) {
        GridCalculator<InputDims, OutputDims>::calculate(model, input, std::forward<Func>(calculate_fn), output);
    }

}// namespace odin

#endif// ODIN_GRID_HPP
