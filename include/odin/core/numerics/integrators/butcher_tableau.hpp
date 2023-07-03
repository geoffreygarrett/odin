#ifndef BUTCHER_TABLEAU_HPP
#define BUTCHER_TABLEAU_HPP

#include <array>
#include <tuple>
#include <vector>

/**
 * @class ButcherTableau
 * @brief Contains the definition of the Butcher Tableau for the representation of a Runge-Kutta method.
 *
 * This struct is based on the work of John C. Butcher, a New Zealand mathematician who specializes
 * in numerical methods for the solution of ordinary differential equations. Butcher works on multistage
 * methods for initial value problems, such as Runge-Kutta and general linear methods. The Butcher tableau
 * is named after him.
 *
 * Reference:
 * A stability property of implicit Runge-Kutta methods - J. C. Butcher, BIT Numerical Mathematics volume 15,
 * pages 358–361 (1975) https://doi.org/10.1007/BF01931672
 *
 * John C. Butcher - https://en.wikipedia.org/wiki/John_C._Butcher
 *
 * @tparam T The type of elements in the arrays.
 * @tparam N The size of arrays a, b, and c.
 * @tparam M The size of array b_star in the extended version. Defaults to 0.
 * @tparam is_extended A boolean to check if M > 0 and therefore the extended version is used. Defaults to false.
 */
//template<typename T, std::size_t N, std::size_t M = 0, bool is_extended = (M > 0)>
//struct ButcherTableau;

// Non-extended version
/*
 * @tparam S The number of stages.
 */

template<typename Scalar, size_t Stages, typename ExtensionType = void>
struct ButcherTableau {
    static constexpr size_t stages      = Stages;
    static constexpr bool   is_extended = false;

    std::array<std::array<Scalar, Stages>, Stages> a;// Runge-Kutta matrix
    std::array<Scalar, Stages>                     b;// weights
    std::array<Scalar, Stages>                     c;// nodes
};

template<typename Scalar, size_t Stages>
struct ButcherTableau<Scalar, Stages, std::array<Scalar, Stages>> {
    static constexpr size_t stages      = Stages;
    static constexpr bool   is_extended = true;

    std::array<std::array<Scalar, Stages - 1>, Stages> a;     // Runge-Kutta matrix
    std::array<Scalar, Stages>                         b;     // weights
    std::array<Scalar, Stages>                         c;     // nodes
    std::array<Scalar, Stages>                         b_star;// lower-order weights
};

namespace butcher_tableau {
    // A non-adaptive tableau
    // Forward Euler method
    template<typename T = double>
    struct Euler : public ButcherTableau<T, 1> {
        constexpr static std::array<std::array<T, 1>, 1> a = {{{0.0}}};
        constexpr static std::array<T, 1>                c = {0.0};
        constexpr static std::array<T, 1>                b = {1.0};
    };

    template<typename T = double>
    struct Heun : public ButcherTableau<T, 2> {
        constexpr static std::array<std::array<T, 2>, 2> a = {
                {{0.0, 0.0}, {1.0, 0.0}}
        };
        constexpr static std::array<T, 2> c = {0.0, 1.0};
        constexpr static std::array<T, 2> b = {0.5, 0.5};
    };

    // Heun's method with Euler method as its lower order method
    template<typename T = double>
    struct HeunEuler : public ButcherTableau<T, 2, std::array<T, 2>> {
        constexpr static auto             a      = Heun<T>::a;
        constexpr static auto             c      = Heun<T>::c;
        constexpr static auto             b      = Heun<T>::b;
        constexpr static std::array<T, 2> b_star = {1.0, 0.0};
    };

    /*
     *
     *
     * a.k.a. Runge–Kutta–Fehlberg, RKF45, Fehlberg 4(5), RK4(5), RK45
     *
     */
    //    template<typename T, std::size_t N, std::size_t M = 0, bool is_extended = (M > 0)>
    //    struct ButcherTableau;
    //
    //    // Non-extended version
    //    template<typename T, std::size_t N, std::size_t M>
    //    struct ButcherTableau<T, N, M, false> {
    //        static constexpr std::size_t N_value = N;
    //        static constexpr std::size_t M_value = M;
    //        static constexpr bool is_extended_value = false;
    //
    //        std::array<std::array<T, N_value>, N_value> a;
    //        std::array<T, N_value> b;
    //        std::array<T, N_value> c;
    //    };


    template<typename T = double>
    struct Fehlberg : public ButcherTableau<T, 6> {
        constexpr static std::array<std::array<T, 6>, 6> a = {
                {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {1.0 / 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0, 0.0},
                 {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0, 0.0},
                 {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0, 0.0},
                 {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0, 0.0}}
        };
        constexpr static std::array<T, 6> c = {
                0.0,
                1.0 / 4.0,
                3.0 / 8.0,
                12.0 / 13.0,
                1.0,
                1.0 / 2.0};
        constexpr static std::array<T, 6> b = {
                16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};
        constexpr static std::array<T, 6> b_star = {
                25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0};
    };

    /**
     * @class DormandPrince
     * @brief Implements the Dormand-Prince method for solving ordinary differential equations.
     *
     * a.k.a. RKDP, DOPRI, default method in MATLAB's ode45
     *
     * This class encapsulates the Dormand-Prince method, an explicit Runge-Kutta method that
     * was originally described in a paper by J.R. Dormand and P.J. Prince published in Journal of
     * Computational and Applied Mathematics in 1980. The paper was titled "A Family of Embedded
     * Runge-Kutta Formulae" and was centered on the derivation of a family of Runge-Kutta formulae
     * RK5(4). These formulae offer 'small' principal truncation terms in the fifth order and
     * extended regions of absolute stability.
     *
     * The Dormand-Prince method is recognized for its efficient handling of non-stiff initial value
     * problems using the Runge-Kutta embedding technique. This technique uses two RK formulas of
     * orders p and q (q > p), which share the same function evaluations. The method has been
     * extensively used for numerical simulations in areas such as dynamical astronomy.
     *
     * @b Reference: Dormand, J.R., Prince, P.J. (1980). A Family of Embedded Runge-Kutta Formulae.
     *               Journal of Computational and Applied Mathematics, 6(1).
     * @b URL: https://www.sciencedirect.com/science/article/pii/0771050X80900133
     */
    template<typename T = double>
    struct DormandPrince : public ButcherTableau<T, 7> {
        constexpr static std::array<std::array<T, 7>, 7> a = {
                {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
                 {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
                 {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
                 {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}}
        };
        constexpr static std::array<T, 7> c = {
                0.0,
                0.2,
                0.3,
                0.8,
                8.0 / 9.0,
                1.0,
                1.0};
        constexpr static std::array<T, 7> b = {
                35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
        constexpr static std::array<T, 7> b_star = {
                5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
    };

    /**
     * @class CashKarp
     * @brief Implements the Cash-Karp method for solving ordinary differential equations.
     *
     * This class encapsulates the Cash-Karp method, which is an explicit Runge-Kutta method that
     * was originally described in a paper by J.R. Cash and A.H. Karp published in ACM Transactions
     * on Mathematical Software in 1990. This paper was titled "A Variable Order Runge-Kutta Method
     * for Initial Value Problems with Rapidly Varying Right-Hand Sides" and was focused on the
     * development of a family of explicit Runge-Kutta formulas for the efficient numerical
     * integration of non-stiff, initial value problems.
     *
     * Each member of the Cash-Karp family consists of a fifth-order formula that contains embedded
     * formulas of all orders 1 through 4. By computing solutions at several different orders, it is
     * possible to detect sharp fronts or discontinuities before all the function evaluations defining
     * the full Runge-Kutta step have been computed.
     *
     * @b Reference: Cash, J.R., Karp, A.H. (1990). A Variable Order Runge-Kutta Method for Initial
     * Value Problems with Rapidly Varying Right-Hand Sides. ACM Transactions on Mathematical
     * Software, 16(3), 201-222.
     */

    template<typename T = double>
    struct CashKarp : public ButcherTableau<T, 6> {
        constexpr static std::array<std::array<T, 6>, 6> a = {
                {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.2, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
                 {0.3, -0.9, 1.2, 0.0, 0.0, 0.0},
                 {-11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0, 0.0, 0.0},
                 {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0, 0.0}}
        };
        constexpr static std::array<T, 6> c = {
                0.0,
                0.2,
                0.3,
                0.6,
                1.0,
                7.0 / 8.0};
        constexpr static std::array<T, 6> b = {
                37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0};
        constexpr static std::array<T, 6> b_star = {
                2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 0.25};
    };


    // Fourth-order Runge-Kutta method
    template<typename T = double>
    struct RK4 : public ButcherTableau<T, 4> {
        constexpr static std::array<std::array<T, 4>, 4> a = {
                {{0, 0, 0, 0},
                 {0.5, 0, 0, 0},
                 {0, 0.5, 0, 0},
                 {0, 0, 1, 0}}
        };
        constexpr static std::array<T, 4> c = {1 / 6., 1 / 3., 1 / 3., 1 / 6.};
        constexpr static std::array<T, 4> b = {0, 0.5, 0.5, 1};
    };
}// namespace butcher_tableau


template<typename T, size_t N>
struct DefaultButcherTableau;

template<typename T>
struct DefaultButcherTableau<T, 1> {
    using type = butcher_tableau::Euler<T>;
};

template<typename T>
struct DefaultButcherTableau<T, 2> {
    using type = butcher_tableau::Euler<T>;
};

template<typename T>
struct DefaultButcherTableau<T, 4> {
    using type = butcher_tableau::RK4<T>;
};


#endif//BUTCHER_TABLEAU_HPP