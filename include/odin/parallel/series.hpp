#ifndef ODIN_SERIES_HPP
#define ODIN_SERIES_HPP

#include "parallel_concepts.hpp"
#include <spdlog/spdlog.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>


template<typename Derived>
std::string to_string(const Eigen::MatrixBase<Derived> &mat) {
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

namespace odin {

    template<typename Scalar, int Dim>
    using vector_type = Eigen::Vector<Scalar, Dim>;

    template<typename Scalar>
    using scalar_type = Scalar;

    template<typename Scalar, int Dim>
    using series_vector_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Dim>;

    template<typename Scalar>
    using series_scalar_type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    template<typename ISingle, typename OSingle, typename ISeries, typename OSeries>
    struct SeriesCalculator;

    // Declaration of the SeriesCalculator
    template<typename ISingle, typename OSingle, typename ISeries, typename OSeries>
    struct SeriesCalculator;

    // Base version
    template<typename ISingle, typename OSingle, typename ISeries, typename OSeries>
    struct SeriesCalculator {
        template<ThreadSafe Derived, typename Func>
        static void calculate(
                Derived       &model,
                const ISeries &input,
                Func         &&calculate_fn,
                OSeries       &output) {

            SPDLOG_DEBUG("Calculating in series");

            const auto &dims = input.rows();

            tbb::parallel_for(
                    tbb::blocked_range<int>(0, dims),
                    [&](const tbb::blocked_range<int> &r) {
                        thread_local Derived tl_model = model.thread_local_copy();
                        for (int i = r.begin(); i != r.end(); ++i) {
                            try {
                                SPDLOG_DEBUG("{} of {}", i, dims);
                                ISingle var_set = input.row(i);
                                SPDLOG_DEBUG("input: {}", to_string(var_set));
                                OSingle result = calculate_fn(tl_model, var_set);
                                output.row(i)  = result;
                            } catch (const std::exception &e) {
                                spdlog::error("Exception caught in parallel loop: {}", e.what());
                                throw;// re-throw the exception after logging the error
                            }
                        }
                    });
        }
    };

    template<typename Scalar, int Dim>
    struct SeriesCalculator<vector_type<Scalar, Dim>, scalar_type<Scalar>, series_vector_type<Scalar, Dim>, series_scalar_type<Scalar>> {
        template<ThreadSafe Derived, typename Func>
        static void calculate(
                Derived                               &model,
                const series_vector_type<Scalar, Dim> &input,
                Func                                 &&calculate_fn,
                series_scalar_type<Scalar>            &output) {

            const auto &dims = input.rows();

            tbb::parallel_for(
                    tbb::blocked_range<int>(0, dims),
                    [&](const tbb::blocked_range<int> &r) {
                        thread_local Derived tl_model = model.thread_local_copy();
                        for (int i = r.begin(); i != r.end(); ++i) {
                            try {
                                SPDLOG_DEBUG("Processing row: {}", i);
                                vector_type<Scalar, Dim> var_set = input.row(i);
                                SPDLOG_DEBUG("Input vector: {}", to_string(var_set));
                                scalar_type<Scalar> result = calculate_fn(tl_model, var_set);
                                SPDLOG_DEBUG("Calculated result: {}", result);
                                output(i) = result;
                            } catch (const std::exception &e) {
                                spdlog::error("Exception caught in parallel loop at iteration {}: {}", i, e.what());
                                throw;// re-throw the exception after logging the error
                            }
                        }
                    });
        }
    };

//    template<typename Func, typename OutputType>
//    OutputType calculate_series(const series_vector_type &positions, Func &&func)
//        requires ThreadSafe<Derived>
//    {
//        OutputType results(positions.rows(), positions.cols());
//
//        SeriesCalculator<
//                vector_type,
//                decltype(func(std::declval<Derived &>(), std::declval<const vector_type &>())),
//                series_vector_type,
//                OutputType>::calculate(as_derived(), positions, std::forward<Func>(func), results);
//
//        return results;
//    }

}// namespace odin

#endif// ODIN_SERIES_HPP
