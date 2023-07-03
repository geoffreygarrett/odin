#include <gsl/gsl_errno.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_ellint.h>

template<typename Scalar>
struct PrecisionSelector;

template<>
struct PrecisionSelector<float> {
    static constexpr int value = GSL_PREC_SINGLE;
};

template<>
struct PrecisionSelector<double> {
    static constexpr int value = GSL_PREC_DOUBLE;
};

template<>
struct PrecisionSelector<long double> {
    static constexpr int value = GSL_PREC_DOUBLE;
};


[[nodiscard]] static auto format_value(const double val, const int width,
                                       const int significant_digits   = 10,
                                       const int fixed_threshold      = 7,
                                       const int decimal_point_length = 1) {
    int int_length  = (val != 0) ? static_cast<int>(std::log10(std::abs(val))) + 1 : 1;
    int sign_length = (val < 0) ? 1 : 0;// Account for negative sign

    std::stringstream tmp_stream;

    if (int_length + sign_length <= fixed_threshold) {
        int fixed_precision = std::max(0, 2 + significant_digits - int_length - sign_length - decimal_point_length);
        tmp_stream << std::fixed << std::setprecision(fixed_precision) << val;
    } else {
        tmp_stream << std::scientific << std::setprecision(significant_digits) << val;
    }

    return tmp_stream.str();
};


template<typename SolverType>
struct SolverTraits;

template<>
struct SolverTraits<gsl_root_fsolver_type> {
    gsl_root_fsolver            *s;
    const gsl_root_fsolver_type *T;

    explicit SolverTraits(const gsl_root_fsolver_type *T) : s(gsl_root_fsolver_alloc(T)), T(T) {}

    SolverTraits(const SolverTraits &other) : s(gsl_root_fsolver_alloc(other.T)), T(other.T) {}

    int set(gsl_function *F, double x_lo, double x_hi) const {
        return gsl_root_fsolver_set(s, F, x_lo, x_hi);
    }

    [[nodiscard]] const char *name() const {
        return gsl_root_fsolver_name(s);
    }

    [[nodiscard]] int iterate() const {
        return gsl_root_fsolver_iterate(s);
    }

    [[nodiscard]] double root() const {
        return gsl_root_fsolver_root(s);
    }

    [[nodiscard]] double x_lower() const {
        return gsl_root_fsolver_x_lower(s);
    }

    [[nodiscard]] double x_upper() const {
        return gsl_root_fsolver_x_upper(s);
    }

    void free() const {
        gsl_root_fsolver_free(s);
    }

    constexpr static const std::array<std::size_t, 5> column_widths   = {5, 18, 18, 18, 18};// Define the column widths
    constexpr static const char                      *separator       = "   ";
    constexpr static const char                      *separator_left  = "  [";
    constexpr static const char                      *separator_right = " ] ";
    constexpr static const char                      *separator_comma = ",  ";

    [[nodiscard]] std::string log_method() const {
        return "Using " + std::string(name()) + " method";
    }

    [[nodiscard]] static std::string log_header() {
        std::stringstream header_stream;
        header_stream << std::setw(column_widths[0]) << "iter" << separator_left
                      << std::setw(column_widths[1]) << "lower_bound" << separator_comma
                      << std::setw(column_widths[2]) << "upper_bound" << separator_right
                      << std::setw(column_widths[3]) << "root_estimate" << separator
                      << std::setw(column_widths[4]) << "error_estimate"
                      << "\n";
        return header_stream.str();
    }


    [[nodiscard]] std::string log_status(int iter) const {
        double            x_lo      = x_lower();
        double            x_hi      = x_upper();
        double            est_error = x_hi - x_lo;
        std::stringstream log_stream;
        log_stream << std::setw(column_widths[0]) << iter << separator_left
                   << std::setw(column_widths[1]) << format_value(x_lo, column_widths[1]) << separator_comma
                   << std::setw(column_widths[2]) << format_value(x_hi, column_widths[2]) << separator_right
                   << std::setw(column_widths[3]) << format_value(root(), column_widths[3]) << separator
                   << std::setw(column_widths[4]) << format_value(est_error, column_widths[4]);
        return log_stream.str();
    }

    [[nodiscard]] static std::string log_convergence() {
        return "Converged:";
    }
};

template<>
struct SolverTraits<gsl_root_fdfsolver_type> {
    gsl_root_fdfsolver            *s;
    const gsl_root_fdfsolver_type *T;

    explicit SolverTraits(const gsl_root_fdfsolver_type *T) : s(gsl_root_fdfsolver_alloc(T)), T(T) {}

    SolverTraits(const SolverTraits &other) : s(gsl_root_fdfsolver_alloc(other.T)), T(other.T) {}

    int set(gsl_function_fdf *FDF, double x0) const {
        return gsl_root_fdfsolver_set(s, FDF, x0);
    }

    [[nodiscard]] const char *name() const {
        return gsl_root_fdfsolver_name(s);
    }

    [[nodiscard]] int iterate() const {
        return gsl_root_fdfsolver_iterate(s);
    }

    [[nodiscard]] double root() const {
        return gsl_root_fdfsolver_root(s);
    }

    void free() const {
        gsl_root_fdfsolver_free(s);
    }

    [[nodiscard]] std::string log_method() const {
        return "Using " + std::string(name()) + " method";
    }

    constexpr static const std::array<std::size_t, 5> column_widths   = {5, 18, 18, 18, 18};// Define the column widths
    constexpr static const char                      *separator_iter  = " | ";
    constexpr static const char                      *separator       = "   ";
    constexpr static const char                      *separator_left  = "  [";
    constexpr static const char                      *separator_right = " ] ";
    constexpr static const char                      *separator_comma = ",  ";

    [[nodiscard]] static std::string log_header() {
        std::stringstream header_stream;
        header_stream << std::setw(column_widths[0]) << "iter" << separator
                      << std::setw(column_widths[3]) << "root_estimate" << separator
                      << std::setw(column_widths[4]) << "error_estimate"
                      << "\n";
        return header_stream.str();
    }

    [[nodiscard]] std::string log_status(int iter) const {
        std::stringstream log_stream;
        log_stream << std::setw(column_widths[0]) << iter << separator_left
                   << std::setw(column_widths[3]) << format_value(root(), column_widths[3]) << separator;
        return log_stream.str();
    }

    [[nodiscard]] static std::string log_convergence() {
        return "Converged:";
    }
};
