#ifndef ODIN_SERIES_HPP
#define ODIN_SERIES_HPP

#include "parallel_concepts.hpp"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace odin {

    // Pre-declaration of SeriesCalculator struct template
    template<int InputDims, int OutputDims = 1>
    struct SeriesCalculator;

    template<int InputDims, int OutputDims>
    struct SeriesCalculator {
        template<ThreadSafe Derived, typename Func, typename IMatrix, typename OMatrix>
        static void calculate(
                Derived       &model,
                const IMatrix &input,
                Func         &&calculate_fn,
                OMatrix       &output) {
            static_assert(IMatrix::ColsAtCompileTime == InputDims || InputDims == Eigen::Dynamic, "Incorrect number of input dimensions.");

            if constexpr (OutputDims != Eigen::Dynamic) {
                static_assert(OMatrix::ColsAtCompileTime == OutputDims, "Incorrect number of output dimensions.");
            }

            using vector_type = Eigen::Matrix<typename IMatrix::Scalar, 1, IMatrix::ColsAtCompileTime>;

            const auto &dims = input.rows();

            tbb::parallel_for(
                    tbb::blocked_range<int>(0, dims),
                    [&](const tbb::blocked_range<int> &r) {
                        Derived tl_model = model.thread_local_copy();
                        for (int i = r.begin(); i != r.end(); ++i) {
                            vector_type var_set = input.row(i);
                            for (int j = 0; j < output.cols(); ++j) {
                                auto result = calculate_fn(tl_model, var_set);
                                if constexpr (OutputDims == 1) {
                                    output(i, 0) = result;
                                } else {
                                    output(i, j) = result;
                                }
                            }
                        }
                    });
        }
    };


    template<ThreadSafe Derived, typename Func, typename IMatrix, typename OMatrix>
    void calculate_in_series(
            Derived       &model,       // A model object that is thread-safe
            const IMatrix &input,       // The input data
            Func         &&calculate_fn,// The function to calculate output
            OMatrix       &output) {          // The output data
        SeriesCalculator<IMatrix::ColsAtCompileTime, OMatrix::ColsAtCompileTime>::calculate(
                model,
                input,
                std::forward<Func>(calculate_fn),
                output);
    }

}// namespace odin

#endif// ODIN_SERIES_HPP
