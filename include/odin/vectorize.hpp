#ifndef ODIN_VECTORIZE_HPP
#define ODIN_VECTORIZE_HPP

#include "include/oneapi/tbb/blocked_range2d.h"
#include "include/oneapi/tbb/blocked_range3d.h"
#include "include/tbb/parallel_for.h"
#include "spdlog/sinks/stdout_color_sinks.h"//support for stdout logging
#include "spdlog/spdlog.h"
#include <Eigen/Dense>
#include <functional>
#include <tbb/parallel_for.h>
#include <type_traits>
#include <unsupported/Eigen/CXX11/Tensor>

#include <omp.h>
#include <span>


namespace odin {

    enum class ParallelBackend {
        NONE,// For no parallelization
        OMP,
        TBB
    };


    namespace stdlib {

        template<int N = 0, typename Scalar = double>
        auto create_vector() {
            if constexpr (N == 0) {
                return std::vector<Scalar>();
            } else {
                return std::array<Scalar, N>();
            }
        }

    }// namespace stdlib

    namespace eigen {

        template<int N = Eigen::Dynamic, typename Scalar = double>
        auto create_vector() {
            return Eigen::Vector<Scalar, N>();
        }

    }// namespace eigen

    namespace detail {

        // A space type enumerator to differentiate between linspace and logspace
        enum class SpaceType {
            LINSPACE,
            LOGSPACE,
            GEOMSPACE
        };

        template<typename Scalar>
        void fill_space(Scalar *space, std::size_t num, Scalar start, Scalar end, bool endpoint, bool parallel, ParallelBackend backend, SpaceType type) {
            if (num == 0) {
                SPDLOG_WARN("The number of points requested is 0.");
                return;
            }

            if (num == 1) {
                space[0] = start;
                return;
            }

            Scalar inv_num = 1.0 / (endpoint ? static_cast<Scalar>(num - 1) : static_cast<Scalar>(num));
            Scalar ratio;
            switch (type) {
                case SpaceType::LINSPACE:
                    ratio = (end - start) * inv_num;
                    break;
                case SpaceType::LOGSPACE:
                    start = std::pow(10, start);// adjust start
                    end   = std::pow(10, end);  // adjust end
                    ratio = std::pow(end / start, inv_num);
                    break;
                case SpaceType::GEOMSPACE:
                    ratio = std::pow(end / start, inv_num);
                    break;
                default:
                    throw std::runtime_error(fmt::format("Invalid space type: {}", static_cast<int>(type)));
            }

            auto linear_fill = [&](Scalar *s, std::size_t size, Scalar start, Scalar ratio) {
                for (std::size_t i = 0; i < size; ++i) {
                    s[i] = start + i * ratio;
                }
            };

            auto exponential_fill = [&](Scalar *s, std::size_t size, Scalar start, Scalar ratio) {
                Scalar val = start;
                for (std::size_t i = 0; i < size; ++i) {
                    s[i] = val;
                    val *= ratio;
                }
            };

            if (!parallel) {
                switch (type) {
                    case SpaceType::LINSPACE:
                        linear_fill(space, num, start, ratio);
                        break;
                    case SpaceType::LOGSPACE:
                    case SpaceType::GEOMSPACE:
                        exponential_fill(space, num, start, ratio);
                        break;
                }
            } else {
                switch (backend) {
                    case ParallelBackend::TBB:
                        SPDLOG_DEBUG("Executing in parallel mode with TBB.");
                        switch (type) {
                            case SpaceType::LINSPACE:
                                tbb::parallel_for(std::size_t(0), num, [&](std::size_t i) {
                                    space[i] = start + i * ratio;
                                });
                                break;
                            case SpaceType::LOGSPACE:
                            case SpaceType::GEOMSPACE:
                                tbb::parallel_for(std::size_t(0), num, [&](std::size_t i) {
                                    space[i] = start * std::pow(ratio, i);
                                });
                                break;
                        }
                        break;
                    case ParallelBackend::OMP:
                        SPDLOG_DEBUG("Executing in parallel mode with OpenMP.");
                        switch (type) {
                            case SpaceType::LINSPACE:
#pragma omp parallel for
                                for (std::size_t i = 0; i < num; i++) {
                                    space[i] = start + i * ratio;
                                }
                                break;
                            case SpaceType::LOGSPACE:
                            case SpaceType::GEOMSPACE:
#pragma omp parallel for
                                for (std::size_t i = 0; i < num; i++) {
                                    space[i] = start * std::pow(ratio, i);
                                }
                                break;
                        }
                        break;
                    case ParallelBackend::NONE:
                        break;
                }
            }

            if (!endpoint) {
                switch (type) {
                    case SpaceType::LINSPACE:
                        space[num - 1] = start + ratio * (num - 1);
                        break;
                    case SpaceType::LOGSPACE:
                    case SpaceType::GEOMSPACE:
                        space[num - 1] = start * std::pow(ratio, num - 1);
                        break;
                }
            }
        }

    }// namespace detail

    namespace stdlib {
        template<detail::SpaceType spaceType, int N = 0, typename Scalar = double>
        auto make_space(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            auto space = create_vector<N, Scalar>();
            if constexpr (N == 0) {
                space.resize(num);
            }
            ParallelBackend backend = parallel ? ParallelBackend::OMP : ParallelBackend::NONE;
            detail::fill_space(space.data(), space.size(), start, end, endpoint, parallel, backend, spaceType);
            return space;
        }
    }// namespace stdlib

    namespace stdlib {
        template<int N = 0, typename Scalar = double>
        auto linspace(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            return make_space<detail::SpaceType::LINSPACE, N, Scalar>(start, end, num, endpoint, parallel);
        }

        template<int N = 0, typename Scalar = double>
        auto logspace(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            return make_space<detail::SpaceType::LOGSPACE, N, Scalar>(start, end, num, endpoint, parallel);
        }

        template<int N = 0, typename Scalar = double>
        auto geomspace(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            return make_space<detail::SpaceType::GEOMSPACE, N, Scalar>(start, end, num, endpoint, parallel);
        }
    }// namespace stdlib

    namespace eigen {
        template<detail::SpaceType spaceType, int N = Eigen::Dynamic, typename Scalar = double>
        auto make_space(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            auto space = create_vector<N, Scalar>();
            if constexpr (N == Eigen::Dynamic) {
                space.resize(num);
            }
            ParallelBackend backend = parallel ? ParallelBackend::OMP : ParallelBackend::NONE;
            detail::fill_space(space.data(), space.size(), start, end, endpoint, parallel, backend, spaceType);
            return space;
        }
    }// namespace eigen

    namespace eigen {
        template<int N = Eigen::Dynamic, typename Scalar = double>
        auto linspace(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            return make_space<detail::SpaceType::LINSPACE, N, Scalar>(start, end, num, endpoint, parallel);
        }

        template<int N = Eigen::Dynamic, typename Scalar = double>
        auto logspace(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            return make_space<detail::SpaceType::LOGSPACE, N, Scalar>(start, end, num, endpoint, parallel);
        }

        template<int N = Eigen::Dynamic, typename Scalar = double>
        auto geomspace(Scalar start, Scalar end, std::size_t num, bool endpoint = true, bool parallel = false) {
            return make_space<detail::SpaceType::GEOMSPACE, N, Scalar>(start, end, num, endpoint, parallel);
        }
    }// namespace eigen

    // Enumeration for indexing options.
    enum class Indexing {
        xy,// x-first and y-second (typical for Cartesian coordinates)
        ij // i-first and j-second (typical for matrix indices)
    };

    // Creates a tensor with dimensions specified in 'dims'.
    template<typename T, std::size_t D, int Options, typename IndexType>
    Eigen::Tensor<T, D, Options, IndexType> make_tensor(const std::array<IndexType, D> &dims) {
        Eigen::Tensor<T, D, Options, IndexType> tensor;// Declare the tensor.
        tensor.resize(dims);                           // Resize it to the provided dimensions.
        return tensor;                                 // Return the tensor.
    }

    // Fills a tensor with values from a given vector, according to the specified indexing option.
    template<typename T, std::size_t N, typename IndexType>
    void fill_tensor(Eigen::Tensor<T, N, 0, IndexType>         &tensor,
                     const Eigen::Matrix<T, Eigen::Dynamic, 1> &x,
                     const std::array<IndexType, N>            &dims,
                     int                                        tensor_index,
                     Indexing                                   indexing,
                     bool                                       parallel) {
        // Lambda to calculate the index in the tensor and assign the value from 'x'.
        auto calc_index_and_assign_value = [&](IndexType i) {
            std::array<IndexType, N> index;
            IndexType                temp = i;
            for (IndexType d = 0; d < N; ++d) {
                index[d] = temp % tensor.dimension(d);
                temp /= tensor.dimension(d);
            }
            // Assign x value according to indexing and tensor_index.
            auto idx      = (indexing == Indexing::xy) ? index[(N - 1) - tensor_index] : index[tensor_index];
            tensor(index) = x(idx);
        };
        // Total size of the tensor.
        auto total_size = tensor.size();

        // Parallelize the operation if specified, otherwise do it serially.
        if (parallel) {
            tbb::parallel_for(IndexType(0), total_size, calc_index_and_assign_value);
        } else {
            for (IndexType i = 0; i < total_size; ++i) {
                calc_index_and_assign_value(i);
            }
        }
    }

    // Fill collection of tensors using provided input vectors.
    template<std::size_t... Is, typename... Ts, typename IndexType>
    void fill_tensors(const std::tuple<Ts...>                                                       &tuple,
                      std::array<Eigen::Tensor<double, sizeof...(Is), 0, IndexType>, sizeof...(Is)> &tensors,
                      std::array<IndexType, sizeof...(Is)> &dims, Indexing indexing, bool parallel) {
        // For each tensor, call fill_tensor with the corresponding vector from the tuple.
        (fill_tensor<double, sizeof...(Is), IndexType>(tensors[Is], std::get<Is>(tuple), dims, static_cast<int>(Is), indexing, parallel), ...);
    }

    // Meshgrid implementation.
    template<typename Scalar, typename IndexType, int Options, typename... Args, std::size_t... Is>
    auto meshgrid_impl(Indexing indexing, bool parallel,
                       std::index_sequence<Is...>,
                       const std::tuple<Args...> &xs) {
        // Obtain the number of dimensions from the input tuple.
        constexpr IndexType N = sizeof...(Is);

        // Define the type of tensors that we are going to work with.
        using tensor_type = Eigen::Tensor<Scalar, N, Options, IndexType>;

        // Gather the dimensions of the input vectors.
        std::array<IndexType, N> dims = {std::get<Is>(xs).size()...};

        // Initialize the array to store the output tensors.
        std::array<tensor_type, N> tensors;

        // Generate and store initial tensors using the dimensions gathered above.
        for (IndexType i = 0; i < N; ++i) {
            tensors[i] = make_tensor<Scalar, N, Options, IndexType>(dims);
        }

        // For each tensor, resize it according to the indexing option (xy or ij).
        for (IndexType i = 0; i < N; ++i) {
            std::array<IndexType, N> new_dims = dims;

            // If 'xy' indexing is specified and there are more than one dimensions,
            // swap the first two dimensions.
            if (indexing == Indexing::xy && N > 1) {
                std::swap(new_dims[0], new_dims[1]);
            }
            tensors[i].resize(new_dims);
        }

        // Fill the tensors with the appropriate values from the input vectors.
        fill_tensors<Is...>(xs, tensors, dims, indexing, parallel);

        // Return the array of filled tensors.
        return tensors;
    }

    // Meshgrid function taking an Indexing enum, a boolean for parallel execution, and a variadic list of arguments.
    template<typename... Args>
    auto meshgrid(Indexing indexing = Indexing::xy, bool parallel = true, const Args &...xs) {
        // The function passes the inputs to meshgrid_impl for actual processing.
        return meshgrid_impl<double, Eigen::Index, Eigen::ColMajor>(indexing, parallel, std::index_sequence_for<Args...>(), std::make_tuple(xs...));
    }

    // Overload of the meshgrid function only taking the variadic list of arguments.
    // This function assumes default values for the indexing and parallel execution.
    template<typename... Args>
    auto meshgrid(const Args &...xs) {
        // Here, the indexing defaults to Indexing::xy, and parallel execution is set to true.
        return meshgrid(Indexing::xy, true, xs...);
    }

    // Overload of the meshgrid function taking a boolean for parallel execution and a variadic list of arguments.
    // This function assumes a default value for the indexing.
    template<typename... Args>
    auto meshgrid(bool parallel, const Args &...xs) {
        // Here, the indexing defaults to Indexing::xy.
        return meshgrid(Indexing::xy, parallel, xs...);
    }

    // Overload of the meshgrid function taking an Indexing enum and a variadic list of arguments.
    // This function assumes a default value for the parallel execution.
    template<typename... Args>
    auto meshgrid(Indexing indexing, const Args &...xs) {
        // Here, parallel execution defaults to true.
        return meshgrid(indexing, true, xs...);
    }


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