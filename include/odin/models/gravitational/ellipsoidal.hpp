#ifndef TRI_AXIAL_ELLIPSOIDAL_HPP
#define TRI_AXIAL_ELLIPSOIDAL_HPP

#include <Eigen/Core>
#include <algorithm>
#include <complex>
#include <iomanip>// std::setprecision, std::setw
#include <odin/logging.hpp>
#include <vector>

#include "gravitational_base.hpp"
#include "gsl.hpp"


namespace odin {

template<typename S,                                // Friday before: typename Scalar,
        typename F        = gsl_root_fsolver_type,  // bracketing solver type
        typename FDF      = gsl_root_fdfsolver_type,// polishing/derivative solver type
        int GSL_PREC_MODE = PrecisionSelector<S>::value>
class TriAxialEllipsoid
    : public GravitationalModelBase<TriAxialEllipsoid<S, F, FDF, GSL_PREC_MODE>, S, 3> {
    //        static_assert(
    //                GravitationalConcept<TriAxialEllipsoid<S, F, FDF, GSL_PREC_MODE>, 3>,
    //                "TriAxialEllipsoid must satisfy the GravitationalModel concept");

private:
    // Define solver traits
    SolverTraits<F>   bracketing_solver_;
    SolverTraits<FDF> derivative_solver_;

public:
    using Scalar = S;

    static_assert(((GSL_PREC_MODE == GSL_PREC_DOUBLE) || (GSL_PREC_MODE == GSL_PREC_SINGLE)
                          || (GSL_PREC_MODE == GSL_PREC_APPROX)),
            "Invalid GSL precision mode");

    // destructor must free the solvers
    ~TriAxialEllipsoid() {
        bracketing_solver_.free();
        derivative_solver_.free();
    }

    //    using scalar = Scalar;

    [[nodiscard]] TriAxialEllipsoid thread_local_copy_impl() const {
        TriAxialEllipsoid copy
                = *this;// The SolverTraits copy constructor will allocate new solvers
        return copy;
    }

    TriAxialEllipsoid(Scalar a, Scalar b, Scalar c, Scalar mu)
        : bracketing_solver_(gsl_root_fsolver_brent),
          derivative_solver_(gsl_root_fdfsolver_steffenson),
          mu_(mu) {

        //        error_occurred = false;
        // Initialization flag
        bool swaps_occurred = false;
        spdlog::debug("Initializing TriAxialEllipsoid with a = {}, b = {}, c = {}", a, b, c);

        // Check and enforce a >= b
        if (a < b) {
            std::swap(a, b);
            swaps_occurred = true;
            spdlog::debug("Swapped a and b to enforce a >= b. New a = {}, b = {}", a, b);
        }

        // Check and enforce b >= c
        if (b < c) {
            std::swap(b, c);
            swaps_occurred = true;
            spdlog::debug("Swapped b and c to enforce b >= c. New b = {}, c = {}", b, c);
        }

        // Check and enforce a >= b again, in case the first swap affected this condition
        if (a < b) {
            std::swap(a, b);
            swaps_occurred = true;
            spdlog::debug("Swapped a and b to enforce a >= b. New a = {}, b = {}", a, b);
        }

        // a >= b >= c is now guaranteed to be true
        // Assign the values to member variables
        a_ = a;
        b_ = b;
        c_ = c;

        // Log a warning if swaps occurred
        if (swaps_occurred) {
            spdlog::warn(
                    "Input axes do not follow the order a >= b >= c. "
                    "The inputs have been swapped to satisfy this condition. "
                    "This might be a result of an estimation process where two axes have very "
                    "similar values. "
                    "Ensure that the provided axes satisfy the condition a >= b >= c.");
        }

        // Use spdlog's built-in formatting for easy-to-read logging
        spdlog::debug("Ellipsoid initialized with a = {}, b = {}, c = {}, mu = {}", a, b, c, mu);

        if constexpr (!std::is_same_v<Scalar, float> && !std::is_same_v<Scalar, double>) {
            spdlog::warn(
                    "Scalar type is neither float nor double. "
                    "Note that the GSL special functions only support up to double precision. "
                    "See https://www.gnu.org/software/gsl/doc/html/specfunc.html#modes for more "
                    "details.");
        }
    }

    [[nodiscard]] Scalar potential_impl(const Eigen::Vector3<Scalar> &position) const {
        Scalar kappa_0 = compute_confocal_ellipsoid_root(position.x(), position.y(), position.z());
        Scalar potential = (// clang-format off
                -((3.0 / 2.0) * mu_ * R_F(a_ * a_ + kappa_0, b_ * b_ + kappa_0, c_ * c_ + kappa_0)
                  - 0.5 * mu_ *
                    (position.x() * position.x() * R_D(b_ * b_ + kappa_0, c_ * c_ + kappa_0, a_ * a_ + kappa_0)
                     + position.y() * position.y() * R_D(a_ * a_ + kappa_0, c_ * c_ + kappa_0, b_ * b_ + kappa_0)
                     + position.z() * position.z() * R_D(a_ * a_ + kappa_0, b_ * b_ + kappa_0, c_ * c_ + kappa_0)
                    ))// clang-format on
        );
        return potential;
    }

    [[nodiscard]] Eigen::Vector3<Scalar> acceleration_impl(
            const Eigen::Vector3<Scalar> &position) const {
        Scalar kappa_0 = compute_confocal_ellipsoid_root(position.x(), position.y(), position.z());
        Eigen::Vector3<Scalar> acceleration = {
                -mu_ * position.x() * R_D(b_ * b_ + kappa_0, c_ * c_ + kappa_0, a_ * a_ + kappa_0),
                -mu_ * position.y() * R_D(a_ * a_ + kappa_0, c_ * c_ + kappa_0, b_ * b_ + kappa_0),
                -mu_ * position.z() * R_D(a_ * a_ + kappa_0, b_ * b_ + kappa_0, c_ * c_ + kappa_0)};
        return acceleration;
    }

    //    TriAxialEllipsoid(Scalar a, Scalar b, Scalar c, Scalar mu) : TriAxialEllipsoid(
    //            a, b, c, mu,
    //            gsl_root_fsolver_brent,
    //            gsl_root_fdfsolver_steffenson) {}

    //    template<typename Scalar, size_t Dim>
    //    struct Params {
    //        std::array<Scalar, Dim> position;
    //        std::array<Scalar, Dim> axes;
    //    };
    //
    //    template<typename Scalar, size_t Dim, size_t... Is>
    //    Scalar C_static_impl(Scalar kappa, const Params<Scalar, Dim>& params, std::index_sequence<Is...>) {
    //        return ((std::pow(params.position[Is], 2) / (std::pow(params.axes[Is], 2) + kappa)) + ...);
    //    }
    //
    //    template<typename Scalar, size_t Dim>
    //    Scalar C_static(Scalar kappa, const Params<Scalar, Dim>& params) {
    //        return C_static_impl(kappa, params, std::make_index_sequence<Dim>{}) - 1;
    //    }
    //
    //    template<typename Scalar, size_t Dim, size_t... Is>
    //    Scalar dC_dkappa_static_impl(Scalar kappa, const Params<Scalar, Dim>& params, std::index_sequence<Is...>) {
    //        return (-std::pow(params.position[Is], 2) / std::pow(std::pow(params.axes[Is], 2) + kappa, 2) + ...);
    //    }
    //
    //    template<typename Scalar, size_t Dim>
    //    Scalar dC_dkappa_static(Scalar kappa, const Params<Scalar, Dim>& params) {
    //        return dC_dkappa_static_impl(kappa, params, std::make_index_sequence<Dim>{});
    //    }
    //
    //    template<typename Scalar, size_t Dim, size_t... Is>
    //    void C_and_dC_dkappa_static_impl(Scalar kappa, const Params<Scalar, Dim>& params, Scalar* C, Scalar* dC, std::index_sequence<Is...>) {
    //        *C  = C_static_impl(kappa, params, std::make_index_sequence<Dim>{});
    //        *dC = dC_dkappa_static_impl(kappa, params, std::make_index_sequence<Dim>{});
    //    }
    //
    //    template<typename Scalar, size_t Dim>
    //    void C_and_dC_dkappa_static(Scalar kappa, const Params<Scalar, Dim>& params, Scalar* C, Scalar* dC) {
    //        C_and_dC_dkappa_static_impl(kappa, params, C, dC, std::make_index_sequence<Dim>{});
    //    }
    struct C_params {
        double x, y, z, a, b, c;
    };

private:
    // GSL expects double precision functions
    static double C_static(double kappa, void *p) {
        auto  *params = (struct C_params *) p;
        double x      = params->x;
        double y      = params->y;
        double z      = params->z;
        double a      = params->a;
        double b      = params->b;
        double c      = params->c;

        double a_sqr_plus_kappa = a * a + kappa;
        double b_sqr_plus_kappa = b * b + kappa;
        double c_sqr_plus_kappa = c * c + kappa;

        return (x * x) / a_sqr_plus_kappa + (y * y) / b_sqr_plus_kappa + (z * z) / c_sqr_plus_kappa
             - 1;
    }

    // GSL expects double precision functions
    static double dC_dkappa_static(double kappa, void *p) {
        auto  *params = (struct C_params *) p;
        double x      = params->x;
        double y      = params->y;
        double z      = params->z;
        double a      = params->a;
        double b      = params->b;
        double c      = params->c;

        double a_sqr_plus_kappa = a * a + kappa;
        double b_sqr_plus_kappa = b * b + kappa;
        double c_sqr_plus_kappa = c * c + kappa;

        return -(x * x) / (a_sqr_plus_kappa * a_sqr_plus_kappa)
             - (y * y) / (b_sqr_plus_kappa * b_sqr_plus_kappa)
             - (z * z) / (c_sqr_plus_kappa * c_sqr_plus_kappa);
    }

    // GSL expects double precision functions
    static void C_and_dC_dkappa_static(double kappa, void *p, double *C, double *dC) {
        auto  *params = (struct C_params *) p;
        double x      = params->x;
        double y      = params->y;
        double z      = params->z;
        double a      = params->a;
        double b      = params->b;
        double c      = params->c;

        double a_sqr_plus_kappa = a * a + kappa;
        double b_sqr_plus_kappa = b * b + kappa;
        double c_sqr_plus_kappa = c * c + kappa;

        *C = (x * x) / a_sqr_plus_kappa + (y * y) / b_sqr_plus_kappa + (z * z) / c_sqr_plus_kappa
           - 1;
        *dC = -(x * x) / (a_sqr_plus_kappa * a_sqr_plus_kappa)
            - (y * y) / (b_sqr_plus_kappa * b_sqr_plus_kappa)
            - (z * z) / (c_sqr_plus_kappa * c_sqr_plus_kappa);
    }


    Scalar a_, b_, c_, mu_;

    static void handler(const char *reason, const char *file, int line, int gsl_errno) {
        SPDLOG_ERROR("GSL error number: {}, File: {}, Line: {}, Reason: {}",
                gsl_errno,
                file,
                line,
                reason);
        throw std::runtime_error("GSL error");
    }

    Scalar compute_confocal_ellipsoid_root(Scalar x,
            Scalar                                y,
            Scalar                                z,
            int                                   max_iter         = 100,
            double                                interval_abs_eps = 0.0,
            double                                interval_rel_eps = 1e-3,
            double                                delta_abs_eps    = 0.0,
            double                                delta_rel_eps    = 1e-3,
            double                                residual_abs_eps = 0.0,
            double                                residual_rel_eps = 0.0) const {
        //        There is only one root that applies to the confocal ellipsoid, the other two
        //        roots belong to the families of the 1 and 2 confocal hyperboloids. The root is
        //        the largest one, so we need to find the largest, therefore I say "The" root.
        //        ODIN_VLOG(50) << "Computing largest root for x = " << x << ", y = " << y << ", z = " << z;
        //        The above is incorrect if we are inside the ellipsoid. I'm not 100% sure how to
        //        handle this case. I think the root is simply 0.0, but I'm not sure.


        // Setting up the function and its derivative for the solver
        C_params         params;
        gsl_function_fdf function_fdf;
        function_fdf.f      = &C_static;
        function_fdf.df     = &dC_dkappa_static;
        function_fdf.fdf    = &C_and_dC_dkappa_static;
        params.x            = x;
        params.y            = y;
        params.z            = z;
        params.a            = a_;
        params.b            = b_;
        params.c            = c_;
        function_fdf.params = &params;
        gsl_function function_f;
        function_f.function = &C_static;
        function_f.params   = &params;

        //        spdlog::error("a: {0}, b: {1}, c: {2}, x: {3}, y: {4}, z: {5}", a_, b_, c_, x, y, z);
        Scalar lower_bound, upper_bound;

        // Check if we are inside the ellipsoid
        static const Scalar perturbation = 1e-16;
        Scalar              C            = C_static(0.0, &params);
        if (C < 0) {

            return 0.0;
        } else {
            // Lower bound is a little > than -c^2. We add a small but meaningful value to -c^2
            // to make sure we don't start from exactly -c^2 (which is undefined).
            lower_bound = -std::pow(c_, 2) + perturbation * std::pow(c_, 2);

            // Upper bound is the square of the maximum value among |a*x|, |b*y|, and |c*z|.
            // We take the absolute values to ensure we deal with positive numbers as we are squaring it.
            // This upper limit is based on an approximation of the potential field outside the ellipsoid.
            upper_bound = std::pow(
                    std::max({std::abs(a_ * x), std::abs(b_ * y), std::abs(c_ * z)}), 2.5);
        }

        //        ODIN_LOG_INFO << "upper_bound: " << upper_bound << " lower_bound: " << lower_bound << std::endl;
        //        Scalar upper_bound = std::pow(std::max({std::abs(a_ * x), std::abs(b_ * y), std::abs(c_ * z)}), 2.5);
        //        Scalar upper_bound = std::max({std::pow(x, 2), std::pow(z, 2), std::pow(z, 2)});
        //        SPDLOG_INFO("upper_bound: {0}, lower_bound: {1}", upper_bound, lower_bound);

        //        SetContext(x, y, z, a_, b_, c_);
        gsl_set_error_handler(&TriAxialEllipsoid::handler);
        //        ODIN_LOG_INFO << "a: " << a_ << " b: " << b_ << " c: " << c_ << " x: " << x << " y: " << y << " z: " << z << " lower_bound: " << lower_bound << " upper_bound: " << upper_bound;
        //        ODIN_LOG_INFO << "C(lower_bound): " << C_static(lower_bound, &params) << " C(upper_bound): " << C_static(upper_bound, &params);
        // Setting up the bracketing solver
        bracketing_solver_.set(&function_f, lower_bound, upper_bound);


        //        if (error_occurred) {
        //            // print out the a,b,c,x,y,z values
        //            ODIN_LOG_ERROR << "Error occurred in GSL solver. Printing out the values of a,b,c,x,y,z";
        //            ODIN_LOG_ERROR << "a: " << a_ << " b: " << b_ << " c: " << c_ << " x: " << x << " y: " << y << " z: " << z;
        //            std::exit(1);
        //        }

        // Setting up the hybrid solver
        int                  status;
        int                  iterations        = 0;
        bool                 use_bracketing    = true;
        static constexpr int STATUS_BRACKETING = 10;
        static constexpr int STATUS_CONVERGED  = 0;

        Scalar root_0;
        Scalar root_1;
        int    total_iters      = 0;
        int    bracketing_iters = 0;
        int    polishing_iters  = 0;

        // Hybrid solver
        do {// This should only run once, but just in case the polishing solver fails, we could bracket for the
            // final desired tolerance.
            //            ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES([&]() { return Blue(bracketing_solver_.log_method()); }(), 60);
            //            ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES([&]() { return bracketing_solver_.log_header(); }(), 60);
            // Bracketing phase
            do {
                root_0 = root_1;
                bracketing_iters++, total_iters++;
                status      = bracketing_solver_.iterate();
                upper_bound = bracketing_solver_.x_upper();
                lower_bound = bracketing_solver_.x_lower();

                root_1 = bracketing_solver_.root();

                status = gsl_root_test_interval(
                        lower_bound, upper_bound, interval_abs_eps, interval_rel_eps);

                if (status == GSL_CONTINUE) {
                    //                    ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES(
                    //                            [&]() { return bracketing_solver_.log_status(bracketing_iters); }(), 60);
                }

            } while (status == GSL_CONTINUE);
            //            ODIN_VLOG(60)
            //                    << ODIN_IF_VERBOSITY_MATCHES([&]() { return Green(bracketing_solver_.log_convergence()); }(), 60);
            //            ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES(
            //                    [&]() { return Green(bracketing_solver_.log_status(bracketing_iters)); }(), 60);


            // Gradient (polish) phase
            Scalar midpoint = (lower_bound + upper_bound) / 2;
            derivative_solver_.set(&function_fdf, midpoint);
            //            ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES([&]() { return Blue(derivative_solver_.log_method()); }(), 60);
            //            ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES([&]() { return derivative_solver_.log_header(); }(), 60);
            //            do {
            //                root_0 = root_1;
            //                polishing_iters++, total_iters++;
            //                status = derivative_solver_.iterate();
            //                root_1 = derivative_solver_.root();
            //                status = gsl_root_test_delta(root_1,
            //                                             root_0,
            //                                             delta_abs_eps,
            //                                             delta_rel_eps);
            //
            //                if (status == GSL_CONTINUE) {
            ////                    ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES(
            ////                            [&]() { return derivative_solver_.log_status(polishing_iters); }(), 60);
            //                }
            //            } while (status == GSL_CONTINUE);
            //            ODIN_VLOG(60)
            //                    << ODIN_IF_VERBOSITY_MATCHES([&]() { return Green(derivative_solver_.log_convergence()); }(), 60);
            //            ODIN_VLOG(60) << ODIN_IF_VERBOSITY_MATCHES(
            //                    [&]() { return Green(derivative_solver_.log_status(polishing_iters)); }(), 60);
        } while (status == GSL_CONTINUE && iterations < max_iter);
        return root_1;
    }

    [[nodiscard]] Scalar R_D(Scalar x, Scalar y, Scalar z) const {
        return gsl_sf_ellint_RD(x, y, z, GSL_PREC_MODE);
    }

    [[nodiscard]] Scalar R_F(Scalar x, Scalar y, Scalar z) const {
        return gsl_sf_ellint_RF(x, y, z, GSL_PREC_MODE);
    }
};

// Verify that it adheres to the GravitationalModel concept
static_assert(is_gravitational_model_v<TriAxialEllipsoid<double>>,
        "TriAxialEllipsoid does not comply with GravitationalConcept.");

}// namespace odin

#endif//TRI_AXIAL_ELLIPSOIDAL_HPP
