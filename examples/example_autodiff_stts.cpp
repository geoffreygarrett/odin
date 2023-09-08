
// include for size_t
#include <Eigen/Core>
#include <tuple>
//#include <unsupported/Eigen/CXX11/Tensor>

//#ifdef ODIN_AUTODIFF
//#define ODIN_CONST
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
using namespace autodiff;
//#else
////#define ODIN_CONST const
//#endif

#include <Eigen/LU>
#include <Eigen/QR>
#include <iomanip>
#include <unsupported/Eigen/CXX11/Tensor>

//// FOR PRINTING TENSORS
template<typename Scalar, int Dimensions>
std::ostream &operator<<(std::ostream &os, const Eigen::Tensor<Scalar, Dimensions> &tensor) {
    if (Dimensions < 1) return os;// Invalid for tensors with less than one dimension

    Eigen::array<Eigen::Index, Dimensions> dims;
    for (int i = 0; i < Dimensions; ++i) { dims[i] = tensor.dimension(i); }

    std::vector<Eigen::Index> index(Dimensions, 0);// Initialize index vector to zeros

    if (Dimensions == 1) {
        // Specific case for 1D tensor (vector)
        for (index[0] = 0; index[0] < dims[0]; ++index[0]) {
            Eigen::array<Eigen::Index, Dimensions> index_array = {index[0]};
            os << std::setw(10) << tensor.coeff(index_array);
        }
        os << "\n";
        return os;
    }

    auto print_recursive
            = [&](auto &self, int dim, std::vector<Eigen::Index> index_so_far) -> void {
        if (dim == Dimensions - 1) {
            for (index[dim] = 0; index[dim] < dims[dim]; ++index[dim]) {
                Eigen::array<Eigen::Index, Dimensions> index_array;
                std::copy(index.begin(), index.end(), index_array.begin());
                os << std::setw(10) << tensor.coeff(index_array);
            }
            os << "\n";
        } else {
            for (index[dim] = 0; index[dim] < dims[dim]; ++index[dim]) {
                index_so_far.push_back(index[dim]);
                if (Dimensions > 2
                        && dim == Dimensions - 3) {// Only print "Slice" for 3D and higher
                    std::string einstein_notation = "Slice (";
                    for (size_t i = 0; i < index_so_far.size(); ++i) {
                        einstein_notation += std::to_string(index_so_far[i]);
                        if (i < index_so_far.size() - 1) { einstein_notation += ", "; }
                    }
                    einstein_notation += ", :, ..., :)";
                    os << einstein_notation << ":\n";
                }
                self(self, dim + 1, index_so_far);
                index_so_far.pop_back();// Remove the last element to revert the state
            }
        }
    };

    print_recursive(print_recursive, 0, {});

    return os;
}

// Function to fill 1st-order tensor
template<typename Fun, typename... Vars, typename... Args, typename U, typename T2>
void sst1(const Fun &eom, const Wrt<Vars...> &wrt, const At<Args...> &at, U &u, T2 &t2) {
    long   n = wrt_total_length(wrt);
    size_t m = 0;
    ForEachWrtVar(wrt, [&](auto &&j, auto &&xj) constexpr {
        u = eval(eom, at, detail::wrt(xj));// Seed xi, xj, xk and evaluate u
        if (m == 0) {
            m = u.size();
            t2.resize(m, n);
        };
        auto du = derivative<1>(u);
        for (size_t i = 0; i < m; ++i) t2(i, j) = du[i];
    });
}

// Function to fill 2nd-order state transition tensor
template<typename Fun, typename... Vars, typename... Args, typename U, typename T3, typename T2>
void sst2(const Fun &eom, const Wrt<Vars...> &wrt, const At<Args...> &at, U &u, T3 &t3, T2 &t2) {
    using namespace autodiff::detail;
    long   n = wrt_total_length(wrt);
    size_t m = 0;
    ForEachWrtVar(wrt, [&](auto &&j, auto &&xj) constexpr {
        ForEachWrtVar(wrt, [&](auto &&k, auto &&xk) constexpr {
            static_assert(!isConst<decltype(xj)> && !isConst<decltype(xk)>,
                    "Expecting non-const autodiff numbers in wrt list because these need to be "
                    "seeded, and thus altered!");
            if (j < j) return;// this take advantage of the fact the Hessian matrix is symmetric
            // evaluate u with xi and xj seeded to produce u0, du/dxi, d2u/dxidxj
            u = eval(eom, at, detail::wrt(xj, xk));
            if (m == 0) {
                m = u.size();
                t3.resize(m, n, n);
                t2.resize(m, n);
            };
            if (j == k) {
                auto du = derivative<1>(u);
                for (size_t i = 0; i < m; ++i) t2(i, j) = du[i];
            }
            auto d2u = derivative<2>(u);
            for (size_t i = 0; i < m; ++i) t3(i, j, k) = t3(i, k, j) = d2u[i];
        });
    });
}

// Function to fill 3rd-order state transition tensor
template<typename Fun, typename... Vars, typename... Args, typename U, typename T4>
void sst3(const Fun &eom, const Wrt<Vars...> &wrt, const At<Args...> &at, U &u, T4 &t4) {
    using namespace autodiff::detail;
    long   n = wrt_total_length(wrt);
    size_t m = 0;
    ForEachWrtVar(wrt, [&](auto &&j, auto &&xj) constexpr {
        ForEachWrtVar(wrt, [&](auto &&k, auto &&xk) constexpr {
            ForEachWrtVar(wrt, [&](auto &&l, auto &&xl) constexpr {
                // TODO: Improve efficiency with symmetry
                static_assert(
                        !isConst<decltype(xj)> && !isConst<decltype(xk)> && !isConst<decltype(xl)>,
                        "Expecting non-const autodiff numbers in wrt list because these need to be "
                        "seeded, and thus altered!");
                if (j < k || k < l)
                    return;// this take advantage of the fact the Hessian matrix is symmetric
                u = eval(eom, at, detail::wrt(xj, xk, xl));
                if (m == 0) {
                    m = u.size();
                    t4.resize(m, n, n, n);
                };
                auto d3u = derivative<3>(u);
                // Fill symmetric entries
                for (int i = 0; i < m; ++i)
                    t4(i, j, k, l) = t4(i, k, j, l) = t4(i, j, l, k) = t4(i, l, j, k)
                            = t4(i, k, l, j) = t4(i, l, k, j) = d3u[i];
            });
        });
    });
}

// Function to fill 4th-order state transition tensor
template<typename Fun, typename... Vars, typename... Args, typename U, typename T5>
void sst4(const Fun &eom, const Wrt<Vars...> &wrt, const At<Args...> &at, U &u, T5 &t5) {
    long   dim_n = wrt_total_length(wrt);
    size_t dim_m = 0;
    ForEachWrtVar(wrt, [&](auto &&j, auto &&xj) constexpr {
        ForEachWrtVar(wrt, [&](auto &&k, auto &&xk) constexpr {
            ForEachWrtVar(wrt, [&](auto &&l, auto &&xl) constexpr {
                ForEachWrtVar(wrt, [&](auto &&m, auto &&xm) constexpr {
                    // TODO: Improve efficiency with symmetry
                    if (j < k || k < l || l < m)
                        return;// this take advantage of the fact the Hessian matrix is symmetric
                    u = eval(eom, at, detail::wrt(xj, xk, xl, xm));
                    if (dim_m == 0) {
                        dim_m = u.size();
                        t5.resize(dim_m, dim_n, dim_n, dim_n, dim_n);
                    };
                    auto d4u = derivative<4>(u);
                    for (int i = 0; i < dim_m; ++i)
                        t5(i, j, k, l, m) = t5(i, j, k, m, l) = t5(i, j, l, k, m)
                                = t5(i, j, l, m, k) = t5(i, j, m, k, l) = t5(i, j, m, l, k)
                                = t5(i, k, j, l, m) = t5(i, k, j, m, l) = t5(i, k, l, j, m)
                                = t5(i, k, l, m, j) = t5(i, k, m, j, l) = t5(i, k, m, l, j)
                                = t5(i, l, j, k, m) = t5(i, l, j, m, k) = t5(i, l, k, j, m)
                                = t5(i, l, k, m, j) = t5(i, l, m, j, k) = t5(i, l, m, k, j)
                                = t5(i, m, j, k, l) = t5(i, m, j, l, k) = t5(i, m, k, j, l)
                                = t5(i, m, k, l, j) = t5(i, m, l, j, k) = t5(i, m, l, k, j)
                                = d4u[i];
                });
            });
        });
    });
}


// Function to retrieve the Nth-order state transition tensor series
template<int N, typename Fun, typename... Vars, typename... Args, typename... Tensors, typename U>
auto stts(const Fun             &eom,
        const Wrt<Vars...>      &wrt,
        const At<Args...>       &at,
        U                       &u,
        std::tuple<Tensors &...> tensors_tuple) {

    using namespace autodiff::detail;
    static_assert(sizeof...(Vars) >= 1, "Expecting at least one wrt variable");
    static_assert(sizeof...(Args) >= 1, "Expecting at least one at variable");
    static_assert(sizeof...(Tensors) >= 1, "Expecting at least one tensor");
    static_assert(N >= 1 && N <= 4, "Expecting N to be between 1 and 3");

    size_t n = wrt_total_length(wrt);
    u.resize(n);

    // For STM (1st order)
    //    sst1(eom, wrt, at, u, std::get<0>(tensors_tuple));

    // For 2nd order
    if constexpr (N >= 2) {
        sst2(eom, wrt, at, u, std::get<1>(tensors_tuple), std::get<0>(tensors_tuple));
    }

    // For 3rd order
    if constexpr (N >= 3) { sst3(eom, wrt, at, u, std::get<2>(tensors_tuple)); }

    // For 4th order
    if constexpr (N >= 4) { sst4(eom, wrt, at, u, std::get<3>(tensors_tuple)); }
}


template<typename Scalar>
struct ScalarDivision {
    Scalar divisor;
    explicit ScalarDivision(const Scalar &divisor)
        : divisor(divisor) {}

    const Scalar operator()(const Scalar &x) const {
        return x / divisor;
    }
};

int main() {

    // Define types
    using scalar_type  = HigherOrderDual<4, double>;// dual1st
    using vector_type  = Eigen::Matrix<scalar_type, 3, 1>;
    using state_type   = Eigen::Matrix<scalar_type, 6, 1>;
    using tensor1_type = Eigen::Tensor<scalar_type, 1>;
    using tensor2_type = Eigen::Tensor<scalar_type, 2>;
    using tensor3_type = Eigen::Tensor<scalar_type, 3>;
    using tensor4_type = Eigen::Tensor<scalar_type, 4>;
    using tensor5_type = Eigen::Tensor<scalar_type, 5>;

    // Initial conditions & parameters
    scalar_type gravity_parameter = 1.0;
    vector_type central_position(0.0, 0.0, 0.0);
    //    scalar_type dt = 1.0;
    state_type state;
    state << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    Eigen::Vector<scalar_type, 4> params;
    params << central_position(0), central_position(1), central_position(2), gravity_parameter;

    // Counter for eom evaluations
    int fevals = 0;

    // Arbitrary point mass acceleration
    auto point_mass_acceleration = [](auto &state, auto &params) -> vector_type {
        auto r = state.template segment<3>(0) - params.template segment<3>(0);
        return -params(3) * r / r.norm() / r.norm() / r.norm();
    };

    // Arbitrary drag acceleration
    scalar_type drag_coefficient  = 1.0;
    scalar_type air_density       = 1.0;
    scalar_type sc_area           = 1.0;
    scalar_type sc_mass           = 1.0;
    auto        drag_acceleration = [&drag_coefficient, &sc_mass, &sc_area, &air_density](
                                     auto &state, auto &params) -> vector_type {
        auto v = state.template segment<3>(3);
        auto partial_dynamic_pressure = 0.5 * air_density * v.norm();
        return -drag_coefficient * partial_dynamic_pressure * sc_area / sc_mass * v;
    };

    // EOM (point mass + drag)
    auto eom = [&](auto &state, auto &params) -> state_type {
        state_type result;
        result.template segment<3>(0) = state.template segment<3>(3);
        result.template segment<3>(3)
                = (point_mass_acceleration(state, params) + drag_acceleration(state, params));
        fevals++;
        return result;
    };

    // State transition tensors
    tensor2_type stt1;
    tensor3_type stt2;
    tensor4_type stt3;
    tensor5_type stt4;

    // State derivative
    state_type state_derivative;

    // Perform differentiation for 1st, 2nd, 3rd, and 4th order
    stts<4>(eom, wrt(state), at(state, params), state_derivative, std::tie(stt1, stt2, stt3, stt4));

    // State derivative
    std::cout << "State derivative: \n" << state_derivative << std::endl;

    // Perturb state
    tensor1_type perturbation;
    perturbation.resize(6);
    perturbation.setZero();
    perturbation(0) = 0.1;
    perturbation(1) = 0.1;
    perturbation(2) = 0.1;
    perturbation(3) = 0.1;
    perturbation(4) = 0.1;
    perturbation(5) = 0.1;

    // Define contraction dimensions for each tensor
    std::cout << "\nSS1: dx/dx0 " << stt1.dimensions() << "\n" << stt1 << std::endl;
    std::cout << "\nSS2: d2x/dx0dx0 " << stt2.dimensions() << "\n" << stt2 << std::endl;
    std::cout << "\nSS3: d3x/dx0dx0dx0 " << stt3.dimensions() << "\n" << stt3 << std::endl;

    // Perform the contractions for terms 2, 3, and 4 (Contraction is general multiplication between tensors)
    // I'm not convinced the indexing is 100% correct yet.
    constexpr auto contraction_1_0
            = Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(1, 0)};
    constexpr auto contraction_2_0
            = Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(2, 0)};
    constexpr auto contraction_3_0
            = Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(3, 0)};
    constexpr auto contraction_4_0
            = Eigen::array<Eigen::IndexPair<int>, 1>{Eigen::IndexPair<int>(4, 0)};

    // [6, 6] -> [m, n] -> [i, j]
    tensor1_type term1 = stt1.contract(perturbation, contraction_1_0);

    // [6, 6, 6] -> [m, n, n] -> [i, j, k]
    tensor1_type term2 = stt2.contract(perturbation, contraction_2_0)
                                 .contract(perturbation, contraction_1_0)
                                 .unaryExpr(ScalarDivision<scalar_type>(2.0));

    // [6, 6, 6, 6] -> [m, n, n, n] -> [i, j, j, k]
    tensor1_type term3 = stt3.contract(perturbation, contraction_3_0)
                                 .contract(perturbation, contraction_2_0)
                                 .contract(perturbation, contraction_1_0)
                                 .unaryExpr(ScalarDivision<scalar_type>(6.0));

    // [6, 6, 6, 6, 6] -> [m, n, n, n, n] -> [i, j, k, l, m]
    tensor1_type term4 = stt4.contract(perturbation, contraction_4_0)
                                 .contract(perturbation, contraction_3_0)
                                 .contract(perturbation, contraction_2_0)
                                 .contract(perturbation, contraction_1_0)
                                 .unaryExpr(ScalarDivision<scalar_type>(24.0));

    // Print the terms
    std::cout << "term1: \n" << term1 << std::endl;
    std::cout << "term2: \n" << term2 << std::endl;
    std::cout << "term3: \n" << term3 << std::endl;
    std::cout << "term4: \n" << term4 << std::endl;

    // Dx = dx/dx0 * Dx0 + d2x/dx0dx0 * Dx0 * Dx0 / 2! + d3x/dx0dx0dx0 * Dx0 * Dx0 * Dx0 / 3! + d4x/dx0dx0dx0dx0 * Dx0 * Dx0 * Dx0 * Dx0 / 4!
    tensor1_type delta = term1 + term2 + term3 + term4;
    std::cout << "delta (sum of terms): \n" << delta << std::endl;

    // TODO: can be reduced with symmetry checks in the ss1,ss2 derivative functions
    std::cout << "\nfevals: " << fevals << std::endl;
    return 0;
}
