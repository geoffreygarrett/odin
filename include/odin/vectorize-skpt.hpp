#ifndef ODIN_VECTORIZE_HPP
#define ODIN_VECTORIZE_HPP

#include "include/oneapi/tbb/blocked_range2d.h"
#include "include/oneapi/tbb/blocked_range3d.h"
#include "include/tbb/parallel_for.h"
#include <Eigen/Dense>
#include <functional>
#include <tbb/parallel_for.h>
#include <type_traits>
#include <unsupported/Eigen/CXX11/Tensor>

namespace odin {

#include <tuple>


    namespace stdlib {

        template<typename Float = double>
        std::vector<Float> linspace(Float start, Float end, std::size_t num, bool endpoint = true) {
            std::vector<Float> linspaced(num);

            if (num > 0) {
                Float delta = endpoint ? (end - start) / (num - 1) : (end - start) / num;

                for (std::size_t i = 0; i < num; ++i) {
                    linspaced[i] = start + delta * i;
                }

                if (!endpoint) {
                    linspaced[num - 1] = start + delta * (num - 1);
                }
            }
            return linspaced;
        }


        template<typename Float = double>
        std::vector<Float> linspace_mt(Float start, Float end, std::size_t num, bool endpoint = true) {
            std::vector<Float> linspaced(num);

            if (num > 0) {
                Float delta = endpoint ? (end - start) / (num - 1) : (end - start) / num;

                tbb::parallel_for(std::size_t(0), num, [&](std::size_t i) {
                    linspaced[i] = start + delta * i;
                });

                if (!endpoint) {
                    linspaced[num - 1] = start + delta * (num - 1);
                }
            }
            return linspaced;
        }

    }// namespace stdlib

    namespace eigen {

        template<typename Float = double>
        Eigen::VectorXd linspace(Float start, Float end, std::size_t num, bool endpoint = true) {
            Eigen::VectorXd linspaced(num);

            if (num > 0) {
                Float delta = endpoint ? (end - start) / (num - 1) : (end - start) / num;

                for (std::size_t i = 0; i < num; ++i) {
                    linspaced[i] = start + delta * i;
                }

                if (!endpoint) {
                    linspaced[num - 1] = start + delta * (num - 1);
                }
            }
            return linspaced;
        }

        template<typename Float = double>
        Eigen::VectorXd linspace_mt(Float start, Float end, std::size_t num, bool endpoint = true) {
            Eigen::VectorXd linspaced(num);

            if (num > 0) {
                Float delta = endpoint ? (end - start) / (num - 1) : (end - start) / num;

                tbb::parallel_for(std::size_t(0), num, [&](std::size_t i) {
                    linspaced[i] = start + delta * i;
                });

                if (!endpoint) {
                    linspaced[num - 1] = start + delta * (num - 1);
                }
            }
            return linspaced;
        }

    }// namespace eigen

    enum class Indexing {
        xy,
        ij
    };

    template<typename T, int D>
    Eigen::Tensor<T, D, 0, int> make_tensor(const std::array<int, D> &dims) {
        Eigen::Tensor<T, D, 0, int> tensor;
        tensor.resize(dims);
        return tensor;
    }


    template<typename T, int N>
    void fill_tensor(Eigen::Tensor<T, N, 0, int>               &tensor,
                     const Eigen::Matrix<T, Eigen::Dynamic, 1> &x,
                     const std::array<int, N>                  &dims,
                     int                                        tensor_index,
                     Indexing                                   indexing) {
        std::array<Eigen::DenseIndex, N> index;
        auto                             total_size = tensor.size();
        for (Eigen::DenseIndex i = 0; i < total_size; ++i) {
            Eigen::DenseIndex temp = i;
            for (int d = 0; d < N; ++d) {
                index[d] = temp % dims[d];
                temp /= dims[d];
            }
            tensor(index) = (indexing == Indexing::xy) ? x(index[(N - 1) - tensor_index]) : x(index[tensor_index]);
        }
    }

    template<typename... Ts, std::size_t... Is>
    auto meshgrid_impl(Indexing indexing,
                       std::integer_sequence<std::size_t, Is...>,
                       const Eigen::Matrix<Ts, Eigen::Dynamic, 1> &...xs) {
        constexpr int N  = sizeof...(Ts);
        using TensorType = Eigen::Tensor<double, N, 0, int>; // Use int index type

        // Create a tuple from xs
        auto xs_tuple = std::make_tuple(xs...);

        std::array<int, N>        dims;
        std::array<TensorType, N> tensors;

        // Use fold expression to initialize dims and tensors
        (..., (dims[Is] = std::get<Is>(xs_tuple).size(), tensors[Is] = make_tensor<double, N>(std::array<int, N>{dims[Is]})));

        // Use fold expression to fill the tensors
        (..., fill_tensor(tensors[Is], std::get<Is>(xs_tuple), dims, static_cast<int>(Is), indexing)); // cast std::size_t to int

        return tensors;
    }


    template<typename... Ts>
    auto meshgrid(const Eigen::Matrix<Ts, Eigen::Dynamic, 1> &...xs) {
        return meshgrid_impl(Indexing::xy, std::make_integer_sequence<int, sizeof...(Ts)>(), xs...);
    }

    template<typename... Ts, std::size_t... Is>
    auto meshgrid(Indexing indexing, std::integer_sequence<std::size_t, Is...>, const Eigen::Matrix<Ts, Eigen::Dynamic, 1> &...xs) {
        return meshgrid_impl(indexing, std::make_integer_sequence<int, sizeof...(Ts)>(), xs...);
    }

    template<typename... Ts>
    auto meshgrid(const Eigen::Matrix<Ts, Eigen::Dynamic, 1> &...xs, Indexing indexing) {
        return meshgrid_impl(indexing, std::make_integer_sequence<int, sizeof...(Ts)>(), xs...);
    }

    //    template<typename T, int D, std::size_t... Is>
    //    void fill_tensor(Eigen::Tensor<T, D> &tensor, const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, std::index_sequence<Is...>) {
    //        auto reshaped = tensor.reshape(std::array<Eigen::DenseIndex, D>{x.rows(), 1});
    //        (void)std::initializer_list<int>{(call_and_discard([&] {
    //                                               Eigen::Tensor<T, 1> chip = reshaped.template chip<Is>(0);
    //                                               auto filled_chip = fill_with_value_parallel(chip, x(Is), Is);
    //                                               reshaped.template chip<Is>(0) = filled_chip;
    //                                           }), 0)...};
    //        tensor = reshaped.reshape(tensor.dimensions());
    //    }

    //    template<typename T, int D, std::size_t... Is>
    //    void fill_tensor(Eigen::Tensor<T, D> &tensor, const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, std::index_sequence<Is...>) {
    //        auto reshaped = tensor.reshape(std::array<Eigen::DenseIndex, D>{x.rows(), 1});
    //        (void(std::initializer_list<int>{(call_and_discard([&]{
    //                                              Eigen::Tensor<T, 1> chip = reshaped.template chip<Is>(0);
    //                                              return fill_with_value_parallel(chip, x(Is), Is);
    //                                          }()), 0)...}));
    //    }

    template<typename I, typename O, typename M, typename... Args>
    using method_type = std::function<O(M &, const I &, Args...)>;
    // For Eigen::Matrix
    // New tensorize function overload for Eigen::Matrix as input
    template<typename I, typename O, int Rank>
    void tensorize(
            Eigen::Tensor<I, Rank>                                       &input,
            Eigen::Tensor<O, Rank>                                       &output,
            std::function<O(const Eigen::Matrix<I, Eigen::Dynamic, 1> &)> func) {
        tensor_for_each(input, [&func, &output](I &value, const size_t &index) {
            output(index) = func(Eigen::Matrix<I, Eigen::Dynamic, 1>::Constant(1, 1, value));
        });
    }

    // Corrected tensorize function overload for Eigen::Tensor as input and Eigen::Matrix as std::function input
    template<typename I, typename O, int Rank>
    void tensorize(
            Eigen::Tensor<I, Rank>                             &input,
            Eigen::Tensor<O, Rank>                             &output,
            std::function<O(const Eigen::Matrix<I, Rank, 1> &)> func) {
        tensor_for_each(input, [&func, &output](I &value, const size_t &index) {
            Eigen::Matrix<I, Rank, 1> matrix;
            matrix.setConstant(value);
            output(index) = func(matrix);
        });
    }

    // Existing functions remain the same
    template<typename I, typename O>
    void tensorize(
            Eigen::Matrix<I, Eigen::Dynamic, Eigen::Dynamic> &input,
            Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic> &output,
            std::function<O(const I &)>                       func) {
        output = input.unaryExpr(func);
    }

    template<typename I, typename O, int Rank>
    void tensorize(
            Eigen::Tensor<I, Rank>     &input,
            Eigen::Tensor<O, Rank>     &output,
            std::function<O(const I &)> func) {
        tensor_for_each(input, [&func, &output](I &value, const size_t &index) {
            output(index) = func(value);
        });
    }

    template<typename I, typename O, int Rank, typename M, typename... Args>
    void tensorize(
            Eigen::Tensor<I, Rank>       &input,
            Eigen::Tensor<O, Rank>       &output,
            method_type<I, O, M, Args...> func,
            M                            &model,
            Args... args) {
        tensor_for_each(input, [&func, &model, &output, args...](I &value, const size_t &index) {
            thread_local M tl_model = model;
            output(index)           = func(tl_model, value, args...);
        });
    }

    template<typename Tensor, typename Function>
    void tensor_for_each(
            Tensor  &tensor,
            Function func) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, tensor.size()), [&](const tbb::blocked_range<size_t> &range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                func(tensor.data()[i], i);
            }
        });
    }

    using Vec3d = Eigen::Matrix<double, 3, 1>;

    template<typename I, typename O, typename M, typename... Args>
    void vectorize(const Eigen::Tensor<I, 3> &input, Eigen::Tensor<O, 4> &output, method_type<Vec3d, O, M, Args...> func, M &model, bool use_thread_local_copy = false, Args &&...args) {
        if (input.dimensions() != output.dimensions()) {
            throw std::invalid_argument("Input and output dimensions do not match");
        }

        tbb::parallel_for(tbb::blocked_range3d<int>(0, input.dimension(0), 0, input.dimension(1), 0, input.dimension(2)),
                          [&](const tbb::blocked_range3d<int> &r) {
                              if (use_thread_local_copy) {
                                  thread_local M tl_model = model;
                                  for (auto x = r.pages().begin(); x != r.pages().end(); ++x)
                                      for (auto y = r.rows().begin(); y != r.rows().end(); ++y)
                                          for (auto z = r.cols().begin(); z != r.cols().end(); ++z)
                                              output(x, y, z) = func(tl_model, Vec3d(x, y, z), std::forward<Args>(args)...);
                              } else {
                                  for (auto x = r.pages().begin(); x != r.pages().end(); ++x)
                                      for (auto y = r.rows().begin(); y != r.rows().end(); ++y)
                                          for (auto z = r.cols().begin(); z != r.cols().end(); ++z)
                                              output(x, y, z) = func(model, Vec3d(x, y, z), std::forward<Args>(args)...);
                              }
                          });
    }
}// namespace odin

#endif//ODIN_VECTORIZE_HPP