
// include for size_t
#include <Eigen/Core>
#include <cstddef>
#include <tuple>

#ifdef ODIN_AUTODIFF
#define ODIN_CONST
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
using namespace autodiff;
#else
#define ODIN_CONST const
#endif

#include <Eigen/LU>
#include <Eigen/QR>


// macro to check for C++20
#if __cplusplus > 201703L
// Two argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_2(TYPE1, TYPE2) requires std::is_same_v<TYPE1, TYPE2>
// Three argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_3(TYPE1, TYPE2, TYPE3) \
    requires(std::is_same_v<TYPE1, TYPE2> || std::is_same_v<TYPE1, TYPE3>)
#else
// Two argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_2(TYPE1, TYPE2) \
    std::enable_if_t<std::is_same_v<TYPE1, TYPE2>, TYPE2>
// Three argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_3(TYPE1, TYPE2, TYPE3)                                       \
    std::enable_if_t<std::disjunction_v<std::is_same<TYPE1, TYPE2>, std::is_same<TYPE1, TYPE3>>, \
            TYPE1>
#endif

// Using helper macros to differentiate between the number of arguments
#define GET_MACRO(_1, _2, _3, NAME, ...) NAME
#define ODIN_ENABLE_IF_TYPE_MATCHES(...)                                                 \
    GET_MACRO(__VA_ARGS__, ODIN_ENABLE_IF_TYPE_MATCHES_3, ODIN_ENABLE_IF_TYPE_MATCHES_2) \
    (__VA_ARGS__)


template<typename Derived>
class crtp_base {
public:
    Derived &as_derived() {
        return static_cast<Derived &>(*this);
    }

    const Derived &as_derived() const {
        return static_cast<const Derived &>(*this);
    }
};


template<typename Derived, typename Scalar, size_t Dim = 3, typename... Params>
class model_base {
public:
    using params_type                     = std::tuple<Params...>;
    using scalar_type                     = Scalar;
    using vector_type                     = Eigen::Matrix<scalar_type, Dim, 1>;
    constexpr static std::size_t n_params = std::tuple_size_v<params_type>;

    explicit model_base(Params... params)
        : m_params(params...) {}

    // getter for params
    template<std::size_t I>
    auto &get_param() {
        return std::get<I>(this->m_params);
    }

    // setter for params
    template<std::size_t I>
    void set_param(auto &&value) {
        std::get<I>(this->m_params) = std::forward<decltype(value)>(value);
    }

    // getter for params
    auto &get_params() {
        return this->m_params;
    }

    // setter for params
    void set_params(auto &params) {
        this->m_params = params;
    }

    template<typename O, typename Func, typename WrtArgs, typename AtArgs, typename U, typename G>
    auto _calc_partial(Func lambda, WrtArgs wrt_args, AtArgs at_args, U &eval, G &grad) {
        if constexpr (std::is_same_v<O, scalar_type>) {
            return autodiff::gradient(lambda, wrt_args, at_args, eval, grad);
        } else if constexpr (std::is_same_v<O, vector_type>) {
            return autodiff::jacobian(lambda, wrt_args, at_args, eval, grad);
        } else {
            static_assert(sizeof(O) == -1, "Invalid return type for partial deducing template");
        }
    }

    template<typename O, typename Func, typename WrtArgs, typename AtArgs, typename U>
    auto _calc_partial(Func lambda, WrtArgs wrt_args, AtArgs at_args, U &eval) {
        if constexpr (std::is_same_v<O, scalar_type>) {
            return autodiff::gradient(lambda, wrt_args, at_args, eval);
        } else if constexpr (std::is_same_v<O, vector_type>) {
            return autodiff::jacobian(lambda, wrt_args, at_args, eval);
        } else {
            static_assert(sizeof(O) == -1, "Invalid return type for partial deducing template");
        }
    }

protected:
    params_type m_params;
};


// Using named constants for better code readability
constexpr std::size_t CENTRAL_POSITION  = 1;
constexpr std::size_t GRAVITY_PARAMETER = 0;

template<typename Derived, typename Scalar, std::size_t Dim = 3>
class gravity_base : public model_base<Derived, Scalar, Dim, Scalar, Eigen::Matrix<Scalar, Dim, 1>>,
                     public crtp_base<Derived> {
public:
    using vector_type  = Eigen::Matrix<Scalar, Dim, 1>;
    using scalar_type  = Scalar;
    using derived_type = Derived;

    template<typename... Params>
    explicit gravity_base(Params... params)
        : model_base<Derived, Scalar, Dim, Params...>(params...) {}

    // Thread-local copy of the derived class
    derived_type thread_local_copy() {
        return this->as_derived().thread_local_copy_impl();
    }

    // Wrapper function for potential
    scalar_type potential(ODIN_CONST vector_type &position) {
        return std::apply(
                [&](auto &&...params) {
                    return Derived::potential_static(
                            position, std::forward<decltype(params)>(params)...);
                },
                this->m_params);
    }

    // Wrapper function for acceleration
    vector_type acceleration(ODIN_CONST vector_type &position) {
        return std::apply(
                [&](auto &&...params) {
                    return Derived::acceleration_static(
                            position, std::forward<decltype(params)>(params)...);
                },
                this->m_params);
    }

    // getter for params
    template<std::size_t I>
    auto &get() {
        return std::get<I>(this->m_params);
    }

    // setter for params
    template<std::size_t I>
    void set(auto &&value) {
        std::get<I>(this->m_params) = std::forward<decltype(value)>(value);
    }

    //    // Generalized calc_partial, pass differential operator as a parameter
    //    template<typename DiffOp, typename Func, typename WrtArgs, typename AtArgs, typename U>
    //    auto calc_partial(DiffOp diff_operator, Func lambda, WrtArgs wrt_args, AtArgs at_args, U &eval) {
    //        // Compile-time check
    //        static_assert(std::is_invocable_v<DiffOp, Func, WrtArgs, AtArgs, U>,
    //                "Invalid differential operator signature");
    //
    //        return diff_operator(lambda, wrt_args, at_args, eval);
    //    }

    template<typename O, typename Func, typename WrtArgs, typename U, typename G>
    void calc_partial(Func lambda, vector_type &position, WrtArgs wrt_args, U &eval, G &grad) {

        auto ref_params_tuple = std::apply(
                [](auto &...args) { return std::make_tuple(std::ref(args)...); }, this->m_params);

        auto params_tuple = std::tuple_cat(std::make_tuple(std::ref(position)), ref_params_tuple);

        auto at_args = std::apply(
                [](auto &&...args) { return at(std::forward<decltype(args)>(args)...); },
                params_tuple);

        this->template _calc_partial<O>(lambda, wrt_args, at_args, eval, grad);
    }

    template<typename O, typename Func, typename WrtArgs, typename U>
    auto calc_partial(Func lambda, vector_type &position, WrtArgs wrt_args, U &eval) {

        auto ref_params_tuple = std::apply(
                [](auto &...args) { return std::make_tuple(std::ref(args)...); }, this->m_params);

        auto params_tuple = std::tuple_cat(std::make_tuple(std::ref(position)), ref_params_tuple);

        auto at_args = std::apply(
                [](auto &&...args) { return at(std::forward<decltype(args)>(args)...); },
                params_tuple);

        return this->template _calc_partial<O>(lambda, wrt_args, at_args, eval);
    }
};

template<typename Scalar, std::size_t Dim = 3>
class point_mass_gravity : public gravity_base<point_mass_gravity<Scalar, Dim>, Scalar, Dim> {
public:
    using scalar_type = Scalar;
    using vector_type = Eigen::Matrix<Scalar, Dim, 1>;

    explicit point_mass_gravity(
            scalar_type &gravity_parameter, vector_type &central_position = vector_type::Zero())
        : gravity_base<point_mass_gravity<Scalar, Dim>, Scalar, Dim>(
                gravity_parameter, central_position) {}

    static scalar_type potential_static(
            ODIN_CONST vector_type &position, scalar_type GM, vector_type &central_position) {
        vector_type r;
        r.noalias() = position - central_position;
        return -GM / r.norm();
    }

    static vector_type acceleration_static(
            ODIN_CONST vector_type &position, scalar_type GM, vector_type &central_position) {
        vector_type r;
        r.noalias()        = position - central_position;
        scalar_type r_norm = r.norm();
        return -GM * r / (r_norm * r_norm * r_norm);
    }

    template<typename T>
    auto thread_local_copy_impl() {
        return std::apply(
                [&](auto &&...params) {
                    return point_mass_gravity<T, Dim>(std::forward<decltype(params)>(params)...);
                },
                this->m_params);
    }
};

using scalar_type = real;
using vector_type = Eigen::Matrix<scalar_type, 3, 1>;

static scalar_type test_potential_static(
        ODIN_CONST vector_type &position, scalar_type GM, vector_type &central_position) {
    auto r = position - central_position;
    return -GM / r.norm();
}

// Forward declaration
template<typename T>
int count_elements(const T &arg);

// Special case for Eigen::Matrix
template<typename T, int Rows, int Cols>
int count_elements(const Eigen::Matrix<T, Rows, Cols> &matrix) {
    return matrix.size();
}

// General function to count elements in a tuple
template<typename... Types>
int count_tuple_elements(const std::tuple<Types...> &tuple) {
    return (... + count_elements(std::get<Types>(tuple)));
}


template<typename, typename = std::void_t<>>
struct has_size_method : std::false_type {};

template<typename T>
struct has_size_method<T, std::void_t<decltype(std::declval<T>().size())>> : std::true_type {};

template<class... Ts>
struct overload : Ts... {
    using Ts::operator()...;
};

template<class... Ts>
overload(Ts...) -> overload<Ts...>;

#include <iostream>
#include <type_traits>
#include <utility>


namespace odin::measurement {

// CRTP Base class for measurement types
template<typename Derived, typename Scalar, std::size_t Dim, typename Data = void>
struct measurement {
    using size        = std::integral_constant<std::size_t, Dim>;
    using scalar_type = Scalar;
    using data_type   = Data;

    // Common utilities for all measurements
    void print_size() const {
        std::cout << "Dimension: " << size::value << std::endl;
    }

    // Function to be specialized
    template<typename... Args>
    void compute(Args &&...args) {
        static_cast<Derived *>(this)->compute_impl(std::forward<Args>(args)...);
    }
};

// Type traits to detect whether a type is a measurement
template<typename T>
struct is_measurement : std::false_type {};

template<typename Derived, typename Scalar, std::size_t Dim, typename Data>
struct is_measurement<measurement<Derived, Scalar, Dim, Data>> : std::true_type {};

// Derived measurement: pseudorange
template<typename Scalar>
struct pseudorange : measurement<pseudorange<Scalar>, Scalar, 1, double> {
    void compute_impl(double input_data) {
        std::cout << "Computing pseudorange with data: " << input_data << std::endl;
    }
};

// Generic N-way radar range measurements
template<typename Scalar>
struct range : measurement<range<Scalar>, Scalar, 1, std::vector<double>> {
    using scalar_type = Scalar;
    using data_type   = std::vector<double>;

    explicit range(size_t n_way, scalar_type c = scalar_type(299792458.0))
        : n_way(n_way),
          c(c) {}

    void compute_impl(const std::vector<scalar_type> &timestamps) {
        if (timestamps.size() != n_way) {
            std::cerr << "Error: Number of timestamps must equal n_way" << std::endl;
            return;
        }

        std::cout << "Computing N-way radar range..." << std::endl;

        // Your computation logic here...
    }

    size_t      n_way;
    scalar_type c;
};

// Derived measurement: range_rate
template<typename Scalar>
struct range_rate : measurement<range_rate<Scalar>, Scalar, 1, double> {
    void compute_impl(double input_data) {
        std::cout << "Computing range rate with data: " << input_data << std::endl;
    }
};

// Derived measurement: doppler
template<typename Scalar>
struct doppler : measurement<doppler<Scalar>, Scalar, 1, double> {
    void compute_impl(double input_data) {
        std::cout << "Computing doppler with data: " << input_data << std::endl;
    }
};

// Function to handle any measurement
template<typename T, typename... Args>
void handle_measurement(T &&m, Args &&...args) {
    static_assert(is_measurement<std::decay_t<T>>::value, "Not a measurement type");
    m.compute(std::forward<Args>(args)...);
}

}// namespace odin::measurement


namespace odin::acceleration {

// CRTP Base class for acceleration types
template<typename Derived, typename Scalar, std::size_t Dim, typename Data = void>
struct acceleration {
    using size        = std::integral_constant<std::size_t, Dim>;
    using scalar_type = Scalar;
    using data_type   = Data;

    // Common utilities for all accelerations
    void print_size() const {
        std::cout << "Dimension: " << size::value << std::endl;
    }

    // Function to be specialized
    template<typename... Args>
    void compute(Args &&...args) {
        static_cast<Derived *>(this)->compute_impl(std::forward<Args>(args)...);
    }
};


}// namespace odin::acceleration

//int main() {
//    pseudorange<double> pr;
//    range_rate<double>  rr;
//    doppler<double>     dp;
//
//    pr.print_size();
//    rr.print_size();
//    dp.print_size();
//
//    handle_measurement(pr, 5.0);
//    handle_measurement(rr, 10.0);
//    handle_measurement(dp, 15.0);
//
//    return 0;
//}


int main() {
    using scalar_type             = dual;
    using vector_type             = Eigen::Matrix<scalar_type, 3, 1>;
    scalar_type gravity_parameter = 200.0;
    vector_type central_position(0.0, 0.0, 0.0);

    // Instantiate a point_mass_gravity object
    point_mass_gravity<scalar_type> gravity(gravity_parameter, central_position);

    //    // Test the potential function
    //    vector_type test_position(20.0, 0.0, 0.0);
    //    scalar_type potential_energy = gravity.potential(test_position);
    //    std::cout << "Potential energy at (1, 0, 0): " << potential_energy << std::endl;
    //
    //    // Test the acceleration function
    //    vector_type acceleration = gravity.acceleration(test_position);
    //    std::cout << "Acceleration at (1, 0, 0): " << acceleration.transpose() << std::endl;
    //    scalar_type potential_result;
    //    //    vector_type potential_partial_result(0.0, 0.0, 0.0);
    //    auto position = vector_type(1.0, 0.0, 0.0);
    //    //    auto gravity_parameter = gravity_parameter;
    //
    //
    //    //    auto at_args = at(position, gravity_parameter, central_position);
    //    //    auto at_args = std::make_tuple(position, gravity_parameter, central_position);
    //    //    auto at_args = std::apply(at, std::make_tuple(position, gravity_parameter, central_position));
    //
    //    // Define your tuple of parameters
    //    //    auto params_tuple1 = std::make_tuple(position, gravity_parameter, central_position);
    //
    //    auto params_tuple1 = std::make_tuple(
    //            std::ref(position), std::ref(gravity_parameter), std::ref(central_position));
    //
    //    //    std::tuple<scalar_type, vector_type> params_tuple = gravity.params();
    //
    //    //    auto params_tuple = std::tuple_cat(std::make_tuple(position), gravity.params());
    //    auto ref_position_tuple = std::make_tuple(std::ref(position));
    //
    //    // Suppose gravity.params() returns a std::tuple<scalar_type, vector_type>
    //    auto original_params  = gravity.get_params();
    //    auto ref_params_tuple = std::make_tuple(
    //            std::ref(std::get<0>(original_params)), std::ref(std::get<1>(original_params)));
    //
    //    auto params_tuple = std::tuple_cat(ref_position_tuple, ref_params_tuple);
    //
    //
    //    std::cout << typeid(params_tuple).name() << std::endl;
    //    std::cout << typeid(params_tuple1).name() << std::endl;
    //
    //    // Unpack tuple into `At` using `at` function
    //    const auto at_args = std::apply(
    //            [](auto &&...args) { return at(std::forward<decltype(args)>(args)...); }, params_tuple);
    //
    //
    //    // Using lambda for partial potential
    //    auto potential_partial_result = gradient(                  // return type
    //            point_mass_gravity<scalar_type>::potential_static, // lambda
    //            wrt(position, gravity_parameter, central_position),// wrt_args
    //            //            at(position, gravity_parameter, central_position),                                  // result
    //            at_args,        // result
    //            potential_result//,       // lambda
    //                            //            potential_partial_result// lambda
    //    );
    //
    //    const auto wrt_args
    //            = std::apply([](auto &&...args) { return wrt(std::forward<decltype(args)>(args)...); },
    //                    params_tuple);
    //
    //
    //    gravity.calc_partial<scalar_type>(                        // return type
    //            point_mass_gravity<scalar_type>::potential_static,// lambda
    //            position,
    //            //            wrt(position,
    //            //                    gravity.get<GRAVITY_PARAMETER>(),
    //            //                    gravity.get<CENTRAL_POSITION>()),// wrt_args
    //            wrt_args,
    //            potential_result,       // result
    //            potential_partial_result// result
    //    );
    //
    //    std::cout << "Test: " << gravity.get<GRAVITY_PARAMETER>() << std::endl;
    //    std::cout << "Partial derivatives of potential: " << potential_partial_result.transpose()
    //              << std::endl;
    //    std::cout << "Potential energy at (1, 0, 0): " << potential_result << std::endl;
    //    vector_type acceleration_result;
    //    //
    //    //    //    // Using lambda for partial acceleration
    //
    //    //    Eigen::Matrix<scalar_type, 3, 3> jacobian_result;


    auto position1 = vector_type(0.0, 10.0, 0.0);
    auto velocity1 = vector_type(0.0, 10.0, 0.0);
    auto dt        = scalar_type(0.2);

    // Initialize state and matrices
    Eigen::Matrix<scalar_type, 6, 1> eom_state_1, eom_state_2, eom_state_derivative;

    eom_state_1 << position1, velocity1;

    std::cout << "State: " << eom_state_1.transpose() << std::endl;

    using jacobian_type = Eigen::Matrix<scalar_type, 3, 3>;
    using vector_type   = Eigen::Matrix<scalar_type, 3, 1>;

    jacobian_type jacobian_eval;
    vector_type   acceleration_eval;
    //
    //    gravity.calc_partial<vector_type>(                           //<-- return type
    //            point_mass_gravity<scalar_type>::acceleration_static,// function
    //            position1,
    //            autodiff::wrt(position1),// with respect to
    //                                     //                    gravity.get<GRAVITY_PARAMETER>(),
    //                                     //                    gravity.get<CENTRAL_POSITION>()),
    //            acceleration_eval,       // eval result for F
    //            jacobian_eval);          // eval result for dF/dx

    std::cout << "Acceleration: " << acceleration_eval.transpose() << std::endl;
    std::cout << "Jacobian: " << std::endl << jacobian_eval << std::endl;


    //    eom_state_1.segment(0, 3) = position1;
    //    eom_state_1.segment(3, 3) = velocity1;

    std::cout << "Initial state: " << eom_state_1.transpose() << std::endl;
    eom_state_2.setZero();
    eom_state_derivative.setZero();

    Eigen::Matrix<scalar_type, 6, 6> eom_jacobian_eval, Phi;
    Phi.setIdentity();
    eom_jacobian_eval.setZero();

    // Equation of motion function
    const auto equation_of_motion = [&](const Eigen::Matrix<scalar_type, 6, 1> &state) {
        Eigen::Matrix<scalar_type, 6, 1> dstate;
        Eigen::Matrix<scalar_type, 3, 1> position = state.segment(0, 3);
        //        std::cout << "Position: " << position.transpose() << std::endl;
        //        std::cout << "Velocity: " << state.segment(3, 3).transpose() << std::endl;
        dstate << state(3), state(4), state(5),
                point_mass_gravity<scalar_type>::acceleration_static(position,
                        gravity.get<GRAVITY_PARAMETER>(),
                        gravity.get<CENTRAL_POSITION>());
        return dstate;
    };

    Eigen::Matrix<scalar_type, 6, 1> state_derivative;
    Eigen::Matrix<scalar_type, 6, 6> state_jacobian;

    // Runge-Kutta with autodiff step function
    auto derivative_and_jacobian
            = [&](auto &current_state, auto &state_derivative, auto &state_jacobian) {
                  jacobian(equation_of_motion,
                          autodiff::wrt(current_state),
                          autodiff::at(current_state),
                          state_derivative,
                          state_jacobian);
                  std::cout << "State derivative: " << state_derivative.transpose() << std::endl;
                  std::cout << "State jacobian: " << std::endl << state_jacobian << std::endl;
              };

    Eigen::Matrix<scalar_type, 6, 1> state_derivative_1;
    Eigen::Matrix<scalar_type, 6, 6> state_jacobian_1;

    derivative_and_jacobian(eom_state_1, state_derivative_1, state_jacobian_1);
    auto state_int = eom_state_1 + state_derivative_1 * dt;
    Phi += (state_jacobian_1 * dt) * Phi;
    auto state_stm = Phi * eom_state_1;

    std::cout << "State 1: " << std::endl << eom_state_1.transpose() << std::endl;
    std::cout << "State 2 (INT): " << std::endl << state_int.transpose() << std::endl;
    std::cout << "State 2 (STM): " << std::endl << state_stm.transpose() << std::endl;
    std::cout << "STM: " << std::endl << Phi << std::endl;


    //    gravity.calc_partial<vector_type>(                           //<-- return type
    //            point_mass_gravity<scalar_type>::acceleration_static,// function
    //            position,
    //            autodiff::wrt(position),
    //            acceleration_eval,// eval result for F
    //            jacobian_eval);   // eval result for dF/dx

    //    gravity.calc_partial<vector_type>(                           //<-- return type
    //            point_mass_gravity<scalar_type>::acceleration_static,// function
    //            position,
    //            autodiff::wrt(position),
    //            acceleration_eval,// eval result for F
    //            jacobian_eval);   // eval result for dF/dx


    //    Eigen::Matrix<scalar_type, 6, 6> stm;
    //    stm.setIdentity();
    //
    //    Eigen::Matrix<scalar_type, 6, 6> full_jacobian;
    //    full_jacobian.setZero();
    //
    //    // Set the upper-left 3x3 block to identity (representing dx/dx usually)
    //    full_jacobian.block<3, 3>(0, 0) = Eigen::Matrix<scalar_type, 3, 3>::Identity().eval();
    //
    //    // Set the bottom-right 3x3 block to identity (representing dv/dv usually)
    //    full_jacobian.block<3, 3>(3, 3) = Eigen::Matrix<scalar_type, 3, 3>::Identity().eval();
    //
    //    // Set the bottom-left 3x3 block to the Jacobian of the acceleration function (representing dv/dx usually)
    //    full_jacobian.block<3, 3>(3, 0) = jacobian_eval.block<3, 3>(0, 0);
    //
    //    // Set the upper-right 3x3 block to the Jacobian of the acceleration function (representing dx/dv usually)
    //    full_jacobian.block<3, 3>(0, 3) = velocity1.asDiagonal();

    //
    //    std::cout << "Full Jacobian: \n" << full_jacobian << std::endl;
    //    stm = full_jacobian * stm;
    //    std::cout << "STM: \n" << stm << std::endl;
    //
    //    std::cout << "Partial derivatives of acceleration: \n" << jacobian_eval << std::endl;
    //    std::cout << "Acceleration at (1, 0, 0): \n" << acceleration_eval.transpose() << std::endl;


    // test estimation process
    using param_type      = std::variant<scalar_type, vector_type>;
    auto current_position = vector_type(1.0, 0.0, 0.0);
    auto current_gm       = scalar_type(1.0);
    auto true_params = std::make_tuple(scalar_type(1.0), vector_type(-0.1987, 1.2345675999, 1.2));
    auto curr_params = std::make_tuple(current_gm, current_position);

    auto true_gravity
            = point_mass_gravity<scalar_type>(std::get<0>(true_params), std::get<1>(true_params));

    auto curr_gravity
            = point_mass_gravity<scalar_type>(std::get<0>(curr_params), std::get<1>(curr_params));

    // sample positions in a circle
    auto num_samples        = 300;
    auto radius             = 3.0;
    auto true_positions     = std::vector<vector_type>(num_samples);
    auto states             = std::vector<vector_type>(num_samples);
    auto observations       = std::vector<vector_type>(num_samples);
    auto curr_accelerations = std::vector<vector_type>(num_samples);
    auto true_errors        = std::vector<vector_type>(num_samples);


    // estimated parameters
    auto &current_params = curr_gravity.get_params();
    auto  est_wrt_params
            = std::apply([](auto &&...args) { return wrt(std::forward<decltype(args)>(args)...); },
                    current_params);

    const auto n_params = wrt_total_length(est_wrt_params);
    const auto dim_obs  = 3;

    auto delta_params = Eigen::Vector<scalar_type, 4>(1.0, 1.0, 1.0, 1.0);

    // sample observations
    for (auto i = 0; i < num_samples; ++i) {
        auto angle         = 2.0 * M_PI * i / num_samples;
        auto test_position = vector_type(radius * std::cos(angle), radius * std::sin(angle), 0.0);
        states[i]          = test_position;
        observations[i]    = true_gravity.acceleration(test_position);
    }

    using h_matrix_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
    using y_vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    // Initialize H matrix and y vector
    h_matrix_type H(dim_obs * num_samples, n_params);// [(dim_obs * num_samples) x n_params]
    y_vector_type y(dim_obs * num_samples);          // [(dim_obs * num_samples) x 1]

    // Initialize R matrix and its inverse
    scalar_type noise_variance = 1e-3;
    auto R = h_matrix_type::Identity(dim_obs * num_samples, dim_obs * num_samples) * noise_variance;
    auto invR = R.inverse();

    // For the current sample, compute Jacobian and acceleration using calc_partial
    vector_type eval;// Replace this type if necessary

    while (delta_params.norm() > 1e-8) {

        // Filling H matrix and y vector
        for (int i = 0; i < num_samples; ++i) {

            // Create lvalues to hold the temporary segments and blocks
            auto H_block = H.block(dim_obs * i, 0, dim_obs, n_params);
            curr_gravity.calc_partial<vector_type>(
                    point_mass_gravity<scalar_type>::acceleration_static,
                    states[i],
                    est_wrt_params,
                    eval,
                    H_block);

            // Storing the observables error as a 3x1 block in y
            y.segment<dim_obs>(dim_obs * i) = (observations[i] - eval);
        }

        // Normalize entire H matrix along the rows
        //        H.rowwise().normalize();

        // Estimate parameter adjustments using Least Squares Estimation
        delta_params = (H.transpose() * H).inverse() * H.transpose() * y;
        // Estimate parameter adjustments using Weighted Least Squares Estimation
        //        delta_params = (H.transpose() * invR * H).inverse() * H.transpose() * invR * y;

        //        auto A = H.transpose() * invR * H;
        //        auto b = H.transpose() * invR * y;
        //        delta_params = A.colPivHouseholderQr().solve(b);

        // Offset for indexing into delta_params
        int offset = 0;

        // Function to update scalar values
        auto update_scalar
                = [&delta_params, &offset](auto &scalar) { scalar += delta_params[offset++]; };

        // Function to update vector values
        auto update_vector = [&delta_params, &offset](auto &vec) {
            int size = vec.size();
            vec += delta_params.segment(offset, size);
            offset += size;
        };

        // Function to update values within the variant param_type
        auto update_param = [&delta_params, &offset, &update_scalar, &update_vector](auto &param) {
            using T = std::decay_t<decltype(param)>;
            if constexpr (std::is_same_v<T, scalar_type>) {
                update_scalar(param);
            } else if constexpr (std::is_same_v<T, vector_type>) {
                update_vector(param);
            }
        };

        // Apply the update function to the current parameters
        std::apply([&update_param](auto &&...args) { (..., update_param(args)); }, current_params);
        std::cout << "delta_params: " << delta_params.transpose() << std::endl;
    }

    return 0;
}
