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
    // String conversion for Eigen::array
    template<typename T, std::size_t N>
    std::string to_string(const Eigen::array<T, N> &arr) {
        std::stringstream ss;
        for (const auto &v: arr) ss << v << ' ';
        return ss.str();
    }

    // String conversion for Eigen::Tensor
    template<typename T, int D>
    std::string to_string(const Eigen::Tensor<T, D> &tensor) {
        std::stringstream ss;
        ss << tensor;
        return ss.str();
    }

    // String conversion for Eigen::Matrix
    template<typename T, int R, int C>
    std::string to_string(const Eigen::Matrix<T, R, C> &matrix) {
        std::stringstream ss;
        ss << matrix;
        return ss.str();
    }

    template<typename Scalar, int Dim>
    using vector_type = Eigen::Vector<Scalar, Dim>;

    template<typename Scalar>
    using scalar_type = Scalar;

    template<typename Scalar, int Dim>
    using tensor_type = Eigen::Tensor<Scalar, Dim>;

    template<typename Scalar>
    using scalar_tensor_type = Eigen::Tensor<Scalar, 1>;

    template<typename Scalar, int Dim>
    using vector_tensor_type = Eigen::Tensor<Scalar, Dim + 1>;


    template<int InputDims, int OutputDims = 1>
    struct GridCalculator {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                                                      &model,
                const std::array<Eigen::Tensor<IType, InputDims>, InputDims> &inputTensors,
                Func                                                        &&calculate_fn,
                Eigen::Tensor<OType, OutputDims>                             &output) {

            // Verify that all input tensors have the same shape
            const auto &dims = inputTensors[0].dimensions();
            for (const auto &tensor: inputTensors) {
                assert(tensor.dimensions() == dims);
            }

            Eigen::array<int, InputDims> indices;

            tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, inputTensors[0].size()),
                    [&](const tbb::blocked_range<size_t> &range) {
                        thread_local Derived tl_model = model.thread_local_copy();
                        for (size_t i = range.begin(); i != range.end(); ++i) {
                            try {
                                // Convert the flat index to multidimensional indices
                                size_t idx = i;
                                for (int d = InputDims - 1; d >= 0; --d) {
                                    indices[d] = idx % dims[d];
                                    idx /= dims[d];
                                }

                                SPDLOG_DEBUG("Processing tensor indices: {}", to_string(indices));

                                // Gather input values from each tensor at the current indices
                                Eigen::array<IType, InputDims> var_set;
                                for (int d = 0; d < InputDims; ++d) {
                                    var_set[d] = inputTensors[d](indices);
                                }

                                auto result     = calculate_fn(tl_model, var_set);
                                output(indices) = result;
                            } catch (const std::exception &e) {
                                spdlog::error("Exception caught in parallel loop at indices {}: {}", to_string(indices), e.what());
                                throw;// re-throw the exception after logging the error
                            }
                        }
                    });
        }
    };

    // Specialization for 3D input of three tensors
    template<>
    struct GridCalculator<3, 1> {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                                      &model,
                const std::array<Eigen::Tensor<IType, 3>, 3> &input,
                Func                                        &&calculate_fn,
                Eigen::Tensor<OType, 3>                      &output) {

            auto &input1 = input[0];
            auto &input2 = input[1];
            auto &input3 = input[2];

            if (input1.dimensions() != input2.dimensions() || input1.dimensions() != input3.dimensions()) {
                throw std::invalid_argument("Input tensors must have the same dimensions");
            }

            const auto &dims = input1.dimensions();

            tbb::parallel_for(
                    tbb::blocked_range3d<size_t>(0, dims[0], 0, dims[1], 0, dims[2]),
                    [&](const tbb::blocked_range3d<size_t> &r) {
                        thread_local Derived tl_model = model.thread_local_copy();
                        for (size_t i = r.pages().begin(); i != r.pages().end(); ++i) {
                            for (size_t j = r.rows().begin(); j != r.rows().end(); ++j) {
                                for (size_t k = r.cols().begin(); k != r.cols().end(); ++k) {
                                    Eigen::array<size_t, 3> indices{i, j, k};

                                    try {
                                        SPDLOG_DEBUG("Processing tensor indices: {}", to_string(indices));
                                        Eigen::Matrix<IType, 3, 1> var_set;
                                        var_set << input1(indices), input2(indices), input3(indices);
                                        OType result    = calculate_fn(tl_model, var_set);
                                        output(indices) = result;
                                    } catch (const std::exception &e) {
                                        spdlog::error("Exception caught in parallel loop at indices {}: {}", to_string(indices), e.what());
                                        throw;// re-throw the exception after logging the error
                                    }
                                }
                            }
                        }
                    });
        }
    };

    // Specialization for 2D input of three tensors
    // Specialization for 2D input of three matrices
    template<>
    struct GridCalculator<2, 1> {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                                                                   &model,
                const std::array<Eigen::Matrix<IType, Eigen::Dynamic, Eigen::Dynamic>, 3> &input,
                Func                                                                     &&calculate_fn,
                Eigen::Matrix<OType, Eigen::Dynamic, Eigen::Dynamic>                      &output) {

            auto &input1 = input[0];
            auto &input2 = input[1];
            auto &input3 = input[2];

            if (input1.rows() != input2.rows() || input1.cols() != input2.cols() || input1.rows() != input3.rows() || input1.cols() != input3.cols()) {
                throw std::invalid_argument("Input matrices must have the same dimensions");
            }

            const size_t rows = input1.rows();
            const size_t cols = input1.cols();

            tbb::parallel_for(
                    tbb::blocked_range2d<size_t>(0, rows, 0, cols),
                    [&](const tbb::blocked_range2d<size_t> &r) {
                        thread_local Derived tl_model = model.thread_local_copy();
                        for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
                            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                                Eigen::array<size_t, 2> indices{i, j};

                                try {
                                    SPDLOG_DEBUG("Processing tensor indices: {}", to_string(indices));
                                    Eigen::Matrix<IType, 3, 1> var_set;
                                    var_set << input1(i, j), input2(i, j), input3(i, j);
                                    OType result = calculate_fn(tl_model, var_set);
                                    output(i, j) = result;
                                } catch (const std::exception &e) {
                                    spdlog::error("Exception caught in parallel loop at indices {}: {}", to_string(indices), e.what());
                                    throw;// re-throw the exception after logging the error
                                }
                            }
                        }
                    });
        }
    };

    // Specialization for 3D input of three tensors
    template<>
    struct GridCalculator<3, 3> {
        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
        static void calculate(
                Derived                                      &model,
                const std::array<Eigen::Tensor<IType, 3>, 3> &input,
                Func                                        &&calculate_fn,
                std::array<Eigen::Tensor<OType, 3>, 3>       &output) {

            auto &input1 = input[0];
            auto &input2 = input[1];
            auto &input3 = input[2];

            if (input1.dimensions() != input2.dimensions() || input1.dimensions() != input3.dimensions()) {
                throw std::invalid_argument("Input tensors must have the same dimensions");
            }

            const auto &dims = input1.dimensions();

            tbb::parallel_for(
                    tbb::blocked_range3d<size_t>(0, dims[0], 0, dims[1], 0, dims[2]),
                    [&](const tbb::blocked_range3d<size_t> &r) {
                        thread_local Derived tl_model = model.thread_local_copy();
                        for (size_t i = r.pages().begin(); i != r.pages().end(); ++i) {
                            for (size_t j = r.rows().begin(); j != r.rows().end(); ++j) {
                                for (size_t k = r.cols().begin(); k != r.cols().end(); ++k) {
                                    Eigen::array<size_t, 3> indices{i, j, k};

                                    try {
                                        SPDLOG_DEBUG("Processing tensor indices: {}", to_string(indices));
                                        Eigen::Matrix<IType, 3, 1> var_set;
                                        var_set << input1(indices), input2(indices), input3(indices);
                                        Eigen::Matrix<OType, 3, 1> result = calculate_fn(tl_model, var_set);
                                        for (int l = 0; l < 3; ++l) {
                                            output[l](indices) = result(l);
                                        }
                                    } catch (const std::exception &e) {
                                        spdlog::error("Exception caught in parallel loop at indices {}: {}", to_string(indices), e.what());
                                        throw;// re-throw the exception after logging the error
                                    }
                                }
                            }
                        }
                    });
        }
    };
    //
    //    // Specialization for 3D input of two tensors
    //    template<>
    //    struct GridCalculator<3, 1> {
    //        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
    //        static void calculate(
    //                Derived                                      &model,
    //                const std::array<Eigen::Tensor<IType, 3>, 2> &input,
    //                Func                                        &&calculate_fn,
    //                Eigen::Tensor<OType, 3>                      &output) {
    //
    //            auto &input1 = input[0];
    //            auto &input2 = input[1];
    //
    //            if (input1.dimensions() != input2.dimensions()) {
    //                throw std::invalid_argument("Input tensors must have the same dimensions");
    //            }
    //
    //            const auto &dims = input1.dimensions();
    //
    //            tbb::parallel_for(
    //                    tbb::blocked_range3d<size_t>(0, dims[0], 0, dims[1], 0, dims[2]),
    //                    [&](const tbb::blocked_range3d<size_t> &r) {
    //                        thread_local Derived tl_model = model.thread_local_copy();
    //                        for (size_t i = r.pages().begin(); i != r.pages().end(); ++i) {
    //                            for (size_t j = r.rows().begin(); j != r.rows().end(); ++j) {
    //                                for (size_t k = r.cols().begin(); k != r.cols().end(); ++k) {
    //                                    Eigen::array<size_t, 3> indices{i, j, k};
    //
    //                                    try {
    //                                        SPDLOG_DEBUG("Processing tensor indices: {}", to_string(indices));
    //                                        Eigen::array<IType, 3> var_set{input1(indices), input2(indices)};
    //                                        OType                  result = calculate_fn(tl_model, var_set);
    //                                        output(indices)               = result;
    //                                    } catch (const std::exception &e) {
    //                                        spdlog::error("Exception caught in parallel loop at indices {}: {}", to_string(indices), e.what());
    //                                        throw;// re-throw the exception after logging the error
    //                                    }
    //                                }
    //                            }
    //                        }
    //                    });
    //        }
    //    };
    //
    //
    //    // Specialization for 2D grid
    //    template<int OutputDims>
    //    struct GridCalculator<2, OutputDims> {
    //        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
    //        static void calculate(
    //                Derived                          &model,
    //                const Eigen::Tensor<IType, 2>    &input,
    //                Func                            &&calculate_fn,
    //                Eigen::Tensor<OType, OutputDims> &output) {
    //
    //            const auto &dims = input.dimensions();
    //            tbb::parallel_for(
    //                    tbb::blocked_range2d<int>(0, dims[0], 0, dims[1]),
    //                    [&](const tbb::blocked_range2d<int> &r) {
    //                        thread_local Derived tl_model = get_thread_local_instance(model);
    //                        for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
    //                            for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
    //                                Eigen::array<int, 2> indices{i, j};
    //                                auto                 var_set = input(indices);
    //                                auto                 result  = calculate_fn(tl_model, var_set);
    //                                if constexpr (OutputDims == 1)
    //                                    output(i, j) = result;
    //                                else
    //                                    output(indices) = result;
    //                            }
    //                        }
    //                    });
    //        }
    //    };
    //
    //    // Specialization for 3D grid
    //    template<int OutputDims>
    //    struct GridCalculator<3, OutputDims> {
    //        template<ThreadSafe Derived, typename Func, typename IType, typename OType>
    //        static void calculate(
    //                Derived                          &model,
    //                const Eigen::Tensor<IType, 3>    &input,
    //                Func                            &&calculate_fn,
    //                Eigen::Tensor<OType, OutputDims> &output) {
    //
    //            const auto &dims = input.dimensions();
    //            tbb::parallel_for(
    //                    tbb::blocked_range3d<int>(0, dims[0], 0, dims[1], 0, dims[2]),
    //                    [&](const tbb::blocked_range3d<int> &r) {
    //                        thread_local Derived tl_model = get_thread_local_instance(model);
    //                        for (int i = r.pages().begin(); i != r.pages().end(); ++i) {
    //                            for (int j = r.rows().begin(); j != r.rows().end(); ++j) {
    //                                for (int k = r.cols().begin(); k != r.cols().end(); ++k) {
    //                                    Eigen::array<int, 3> indices{i, j, k};
    //                                    auto                 var_set = input(indices);
    //                                    auto                 result  = calculate_fn(tl_model, var_set);
    //                                    if constexpr (OutputDims == 1)
    //                                        output(i, j, k) = result;
    //                                    else
    //                                        output(indices) = result;
    //                                }
    //                            }
    //                        }
    //                    });
    //        }
    //    };
    //
    //    template<ThreadSafe Derived, typename Func, typename IType, typename OType, int InputDims, int OutputDims>
    //    void calculate_on_grid(
    //            Derived                               &model,
    //            const Eigen::Tensor<IType, InputDims> &input,
    //            Func                                 &&calculate_fn,
    //            Eigen::Tensor<OType, OutputDims>      &output) {
    //        GridCalculator<InputDims, OutputDims>::calculate(model, input, std::forward<Func>(calculate_fn), output);
    //    }

}// namespace odin

#endif// ODIN_GRID_HPP
