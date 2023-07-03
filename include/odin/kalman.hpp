#ifndef KALMAN_HPP
#define KALMAN_HPP

#include <Eigen/Dense>
#include <glog/logging.h>
#include <iostream>
#include <odin/utilities/sigma_point_calculator.hpp>
#include <optional>
#include <utility>
#include <vector>

typedef Eigen::VectorXd Parameter;

class ParameterEstimator {
public:
    virtual ~ParameterEstimator() = default;

    //    virtual void estimate_parameters(const Eigen::VectorXd &x, const
    //    Eigen::VectorXd &y) = 0;

    [[nodiscard]] Parameter get_parameters() const { return parameters; }

protected:
    Parameter parameters;
};

template<typename StateType, typename MeasurementType,
         typename StateCovarianceType, typename MeasurementCovarianceType>
class KalmanFilter : public ParameterEstimator {
public:
    static_assert(
            static_cast<int>(StateType::RowsAtCompileTime) ==
                            static_cast<int>(StateCovarianceType::RowsAtCompileTime) &&
                    static_cast<int>(StateType::RowsAtCompileTime) ==
                            static_cast<int>(StateCovarianceType::ColsAtCompileTime),
            "StateType dimensions must match with StateCovarianceType");
    static_assert(
            static_cast<int>(MeasurementType::RowsAtCompileTime) ==
                            static_cast<int>(MeasurementCovarianceType::RowsAtCompileTime) &&
                    static_cast<int>(MeasurementType::RowsAtCompileTime) ==
                            static_cast<int>(MeasurementCovarianceType::ColsAtCompileTime),
            "MeasurementType dimensions must match with MeasurementCovarianceType");

    // Various function pointer types used in KalmanFilter
    using scalar_t = StateType::Scalar;
    using state_t = StateType;
    using measurement_t = MeasurementType;
    using StateFunc = std::function<StateType(const state_t &, scalar_t)>;
    using MeasurementFunc = std::function<MeasurementType(const state_t &)>;
    using ProcessCovarianceFunc =
            std::function<StateCovarianceType(const state_t &, scalar_t)>;
    using MeasurementCovarianceFunc =
            std::function<MeasurementCovarianceType(const measurement_t &)>;

    KalmanFilter(StateFunc state_transition, MeasurementFunc measurement_model,
                 ProcessCovarianceFunc process_covariance,
                 MeasurementCovarianceFunc measurement_covariance)
        : state_transition_(std::move(state_transition)),
          measurement_model_(std::move(measurement_model)),
          process_covariance_(std::move(process_covariance)),
          measurement_covariance_(std::move(measurement_covariance)) {}

    ~KalmanFilter() override = default;

    virtual void predict(const scalar_t dt) {
        //        x_ = state_transition_(x_, dt);
        //        P_ = F_ * P_ * F_.transpose() + process_covariance_(x_);
    }

    virtual void update(const MeasurementType &measurement,
                        const MeasurementCovarianceType &measurement_noise) {
        //        MeasurementType y = measurement - measurement_model_(x_);
        //        Eigen::MatrixXd S = H_ * P_ * H_.transpose() +
        //        measurement_covariance_(measurement); Eigen::MatrixXd K = P_ *
        //        H_.transpose() * S.inverse();
        //
        //        x_ = x_ + K * y;
        //        P_ = (Eigen::MatrixXd::Identity(P_.rows(), P_.cols()) - K * H_) *
        //        P_;
    }

protected:
    StateType x_;          // state
    StateCovarianceType P_;// covariance
    Eigen::MatrixXd F_;    // state transition matrix
    Eigen::MatrixXd H_;    // observation model

    StateFunc state_transition_;
    MeasurementFunc measurement_model_;
    ProcessCovarianceFunc process_covariance_;
    MeasurementCovarianceFunc measurement_covariance_;

    // private:
};

template<
        typename StateType, typename MeasurementType, typename StateCovarianceType,
        typename MeasurementCovarianceType,
        typename StateTransitionJacobianFunc = Eigen::MatrixXd (*)(
                const StateType &, const typename KalmanFilter<
                                           StateType, MeasurementType, StateCovarianceType,
                                           MeasurementCovarianceType>::StateFunc &),
        typename MeasurementModelJacobianFunc = Eigen::MatrixXd (*)(
                const StateType &, const typename KalmanFilter<
                                           StateType, MeasurementType, StateCovarianceType,
                                           MeasurementCovarianceType>::MeasurementFunc &)>
class ExtendedKalmanFilter
    : public KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                          MeasurementCovarianceType> {

    static_assert(std::is_arithmetic<StateType>::value,
                  "StateType must be a numeric type");
    static_assert(std::is_arithmetic<MeasurementType>::value,
                  "MeasurementType must be a numeric type");

public:
    using typename KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                                MeasurementCovarianceType>::StateFunc;
    using typename KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                                MeasurementCovarianceType>::MeasurementFunc;
    using typename KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                                MeasurementCovarianceType>::ProcessCovarianceFunc;
    using typename KalmanFilter<
            StateType, MeasurementType, StateCovarianceType,
            MeasurementCovarianceType>::MeasurementCovarianceFunc;

    using JacobianFunc = std::function<Eigen::MatrixXd(const StateType &)>;

    using scalar_t = StateType::Scalar;

    ExtendedKalmanFilter(StateFunc state_transition,
                         MeasurementFunc measurement_model,
                         ProcessCovarianceFunc process_covariance,
                         MeasurementCovarianceFunc measurement_covariance,
                         StateTransitionJacobianFunc state_transition_jacobian =
                                 numerical_estimate_jacobian,
                         MeasurementModelJacobianFunc measurement_model_jacobian =
                                 numerical_estimate_jacobian)
        : KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                       MeasurementCovarianceType>(
                  std::move(state_transition), std::move(measurement_model),
                  std::move(process_covariance), std::move(measurement_covariance)),
          state_transition_jacobian_(std::move(state_transition_jacobian)),
          measurement_model_jacobian_(std::move(measurement_model_jacobian)) {}

    void predict(const scalar_t dt) override {
        this->x_ = this->state_transition_(this->x_, dt);
        this->F_ = state_transition_jacobian_(this->x_, this->state_transition_);
        this->P_ = this->F_ * this->P_ * this->F_.transpose() +
                   this->process_covariance_(this->x_);
    }

    //    template<int MeasurementSize>
    //    void update(const Eigen::Matrix<double, MeasurementSize, 1>
    //    &measurement,
    //                const Eigen::Matrix<double, MeasurementSize,
    //                MeasurementSize> &measurement_noise) {

    void update(const MeasurementType &z) override {
        MeasurementType y = z - this->measurement_model_(this->x_);
        this->H_ = measurement_model_jacobian_(this->x_, this->measurement_model_);
        Eigen::MatrixXd S = this->H_ * this->P_ * this->H_.transpose() +
                            this->measurement_covariance_(z);
        Eigen::MatrixXd K = this->P_ * this->H_.transpose() * S.inverse();

        this->x_ = this->x_ + K * y;
        this->P_ = (Eigen::MatrixXd::Identity(this->P_.rows(), this->P_.cols()) -
                    K * this->H_) *
                   this->P_;
    }

    static Eigen::MatrixXd
    numerical_estimate_jacobian(const StateType &x, const StateFunc &state_func) {
        using state_t = typename StateType::Scalar;
        Eigen::MatrixXd jacobian(x.size(), x.size());
        StateType x_plus_dx, x_minus_dx;
        constexpr state_t dx = 1e-6;

        for (int i = 0; i < x.size(); ++i) {
            x_plus_dx = x_minus_dx = x;
            x_plus_dx(i) += dx;
            x_minus_dx(i) -= dx;

            Eigen::VectorXd derivative =
                    (state_func(x_plus_dx) - state_func(x_minus_dx)) / (2.0 * dx);
            jacobian.col(i) = derivative;
        }

        return jacobian;
    }

private:
    // define your Jacobian calculations here
    Eigen::MatrixXd calculate_jacobian_f(const StateType &x);

    Eigen::MatrixXd calculate_jacobian_h(const StateType &x);

    StateTransitionJacobianFunc state_transition_jacobian_;
    MeasurementModelJacobianFunc measurement_model_jacobian_;
};

struct InversionPolicy {
    template<typename MatrixType>
    static MatrixType invert(const MatrixType &matrix) {
        static_assert(always_false<MatrixType>::value,
                      "InversionPolicy::invert() not implemented");
        return matrix;// this will never be executed
    }

private:
    template<typename T>
    struct always_false : std::false_type {};
};

struct DirectInversePolicy : public InversionPolicy {
    template<typename MatrixType>
    static MatrixType invert(const MatrixType &matrix) {
        return matrix.inverse();
    }
};

//
// struct RegularizedInversePolicy : public InversionPolicy {
//    template<typename MatrixType>
//    static MatrixType invert(const MatrixType &matrix) {
//        double lambda = 1e-3; // Regularization factor
//        MatrixType identity = MatrixType::Identity(matrix.rows(),
//        matrix.cols()); return (matrix.transpose() * matrix + lambda *
//        identity).ldlt().solve(matrix.transpose());
//    }
//};
//
// struct LDLTInversePolicy : public InversionPolicy {
//    template<typename MatrixType>
//    static MatrixType invert(const MatrixType &matrix) {
//        return matrix.ldlt().solve(MatrixType::Identity(matrix.rows(),
//        matrix.cols()));
//    }
//};

/**
 * @brief Checks for any numerical issues within a given matrix.
 *
 * This function checks if a matrix has NaN or infinity values which could cause
 * computational errors.
 *
 * @tparam Derived The Eigen::Matrix derived type.
 * @param mat A const Eigen::MatrixBase<Derived> reference to the matrix to be
 * checked.
 * @return Boolean value indicating the presence of any numerical issues.
 */
template<typename Derived>
bool check_numerical_issues(const Eigen::MatrixBase<Derived> &mat) {
    if (mat.hasNaN() ||
        (mat.array() == std::numeric_limits<typename Derived::Scalar>::infinity())
                .any() ||
        (mat.array() ==
         -std::numeric_limits<typename Derived::Scalar>::infinity())
                .any()) {
        return true;
    }
    return false;
}

/**
 * @brief Checks if the given matrix is well-conditioned for computations.
 *
 * This function performs a FullPivLU decomposition on the input matrix and
 * checks if it is full-rank. If the matrix is not full-rank, it may lead to
 * issues in computations like inverse calculations.
 *
 * @tparam Derived The Eigen::Matrix derived type.
 * @param mat A const Eigen::MatrixBase<Derived> reference to the matrix to be
 * checked.
 * @return Boolean value indicating whether the matrix is well-conditioned.
 */

template<typename MatrixType>
struct ConditioningResult {
    using scalar_t = typename MatrixType::Scalar;
    using complex_scalar_t = std::complex<scalar_t>;
    static constexpr int Rows = MatrixType::RowsAtCompileTime;
    static constexpr int Cols = MatrixType::ColsAtCompileTime;
    static constexpr int Dimension =
            (Rows != Eigen::Dynamic && Cols != Eigen::Dynamic) ? Rows
                                                               : Eigen::Dynamic;
    bool is_well_conditioned;
    scalar_t condition_number;
    std::optional<scalar_t> determinant;
    std::optional<Eigen::Matrix<complex_scalar_t, Dimension, 1>> eigenvalues;
};

template<typename Derived>
bool is_symmetric(const Eigen::MatrixBase<Derived> &matrix) {
    return matrix.isApprox(matrix.transpose(), 1e-10);
}

// #include <glog/logging.h>

template<typename MatrixType, bool ThrowOnError = false, bool Verbose = true,
         bool ReturnFullResult = true>
ConditioningResult<MatrixType>
check_conditioning(const MatrixType &matrix,
                   typename MatrixType::Scalar condition_threshold = 1.0e14,
                   const std::string &matrix_name = "") {
    Eigen::JacobiSVD<MatrixType> svd(matrix);
    using scalar_t = typename MatrixType::Scalar;
    scalar_t smallest_singular_value =
            svd.singularValues()(svd.singularValues().size() - 1);
    scalar_t largest_singular_value = svd.singularValues()(0);
    scalar_t cond = largest_singular_value / smallest_singular_value;

    ConditioningResult<MatrixType> result;
    result.condition_number = cond;
    result.is_well_conditioned = cond <= condition_threshold;

    if constexpr (MatrixType::RowsAtCompileTime != Eigen::Dynamic &&
                  MatrixType::ColsAtCompileTime != Eigen::Dynamic &&
                  MatrixType::RowsAtCompileTime ==
                          MatrixType::ColsAtCompileTime) {
        if (matrix.isApprox(matrix.transpose(),
                            1e-10)) {// Check if the matrix is symmetric
            result.determinant = matrix.determinant();
            result.eigenvalues = matrix.eigenvalues();
        }
    } else if (matrix.rows() == matrix.cols()) {
        if (matrix.isApprox(matrix.transpose(),
                            1e-10)) {// Check if the matrix is symmetric
            result.determinant = matrix.determinant();
            result.eigenvalues = matrix.eigenvalues();
        }
    }

    if (!result.is_well_conditioned) {
        std::stringstream error_message;
        error_message << "Matrix is ill-conditioned. Condition number is: " << cond
                      << std::endl;
        error_message << "Smallest and largest singular values are: "
                      << smallest_singular_value << ", " << largest_singular_value
                      << ". ";
        if (result.determinant.has_value()) {
            error_message << "Determinant is: " << *result.determinant << ". ";
            if (result.eigenvalues.has_value()) {
                std::stringstream eigenvalues_ss;
                for (const auto &eigenvalue: *result.eigenvalues) {
                    eigenvalues_ss << eigenvalue << ", ";
                }
                std::string eigenvalues_str = eigenvalues_ss.str();
                eigenvalues_str = eigenvalues_str.substr(
                        0, eigenvalues_str.length() - 2);// get rid of the last comma
                error_message << "Eigenvalues are: " << eigenvalues_str << ". ";
            } else {
                error_message << "Eigenvalues are: N/A. ";
            }
        }
        if (!matrix_name.empty()) {
            std::stringstream final_error_message;
            final_error_message << "Error occurred with matrix '" << matrix_name
                                << "'. " << error_message.str();
            error_message.swap(final_error_message);
        }

        if (ThrowOnError) {
            throw std::runtime_error(error_message.str());
        } else if (Verbose) {
            LOG(ERROR) << error_message.str();
        }
    }

    if constexpr (!ReturnFullResult) {
        ConditioningResult<MatrixType> boolean_result;
        boolean_result.is_well_conditioned = result.is_well_conditioned;
        return boolean_result;
    } else {
        return result;
    }
}

// template<typename MatrixType, bool ThrowOnError = false, bool Verbose = true,
// bool ReturnFullResult = true> ConditioningResult<MatrixType>
// check_conditioning(const MatrixType &matrix, double condition_threshold
// = 1.0e14, const std::string &matrix_name = "") {
//     Eigen::JacobiSVD<MatrixType> svd(matrix);
//     double smallest_singular_value =
//     svd.singularValues()(svd.singularValues().size() - 1); double
//     largest_singular_value = svd.singularValues()(0); double cond =
//     largest_singular_value / smallest_singular_value;
//
//     ConditioningResult<MatrixType> result;
//     result.condition_number = cond;
//     result.is_well_conditioned = cond <= condition_threshold;
//
//     if constexpr (MatrixType::RowsAtCompileTime != Eigen::Dynamic &&
//                   MatrixType::ColsAtCompileTime != Eigen::Dynamic &&
//                   MatrixType::RowsAtCompileTime ==
//                   MatrixType::ColsAtCompileTime) {
//         if (matrix.isApprox(matrix.transpose(), 1e-10)) { // Check if the
//         matrix is symmetric
//             result.determinant = matrix.determinant();
//             result.eigenvalues = matrix.eigenvalues();
//         }
//     } else if (matrix.rows() == matrix.cols()) {
//         if (matrix.isApprox(matrix.transpose(), 1e-10)) { // Check if the
//         matrix is symmetric
//             result.determinant = matrix.determinant();
//             result.eigenvalues = matrix.eigenvalues();
//         }
//     }
//
//     if (!result.is_well_conditioned) {
//         std::stringstream error_message;
//         error_message << "Matrix is ill-conditioned. Condition number is: "
//         << cond << std::endl; error_message << "Smallest and largest singular
//         values are: " << smallest_singular_value << ", "
//                       << largest_singular_value << ". ";
//         if (result.determinant.has_value()) {
//             error_message << "Determinant is: " << *result.determinant << ".
//             "; if (result.eigenvalues.has_value()) {
//                 std::stringstream eigenvalues_ss;
//                 for (const auto &eigenvalue: *result.eigenvalues) {
//                     eigenvalues_ss << eigenvalue << ", ";
//                 }
//                 std::string eigenvalues_str = eigenvalues_ss.str();
//                 eigenvalues_str = eigenvalues_str.substr(0,
//                 eigenvalues_str.length() - 2);  // get rid of the last comma
//                 error_message << "Eigenvalues are: " << eigenvalues_str << ".
//                 ";
//             } else {
//                 error_message << "Eigenvalues are: N/A. ";
//             }
//         }
//         if (!matrix_name.empty()) {
//             std::stringstream final_error_message;
//             final_error_message << "Error occurred with matrix '" <<
//             matrix_name << "'. " << error_message.str();
//             error_message.swap(final_error_message);
//         }
//
//         if (ThrowOnError) {
//             throw std::runtime_error(error_message.str());
//         } else if (Verbose) {
//             LOG(ERROR) << error_message.str();
//         }
//     }
//
//     if constexpr (!ReturnFullResult) {
//         ConditioningResult<MatrixType> boolean_result;
//         boolean_result.is_well_conditioned = result.is_well_conditioned;
//         return boolean_result;
//     } else {
//         return result;
//     }
// }

template<typename MatrixType, typename Scalar = typename MatrixType::Scalar>
bool is_well_conditioned(const MatrixType &matrix,
                         Scalar condition_threshold = 1.0e12,
                         const std::string &matrix_name = "") {
    auto result = check_conditioning<MatrixType, false, false, false>(
            matrix, condition_threshold, matrix_name);
    return result.is_well_conditioned;
}

/**
 * @brief Computes the condition number of eigenvalues of a given matrix.
 *
 * This function calculates the condition number of eigenvalues using the
 * formula: κ(λ)=∥x∥2∥y∥2/|x∗y|. It returns the maximum condition number across
 * all eigenvalues as the most sensitive eigenvalue would determine the overall
 * sensitivity.
 *
 * @param mat A const Eigen::MatrixXd reference to the matrix whose eigenvalues'
 * condition numbers are to be computed.
 * @return The maximum condition number across all eigenvalues of the input
 * matrix.
 */
// double condition_number_eigenvalues(const Eigen::MatrixXd &mat) {
//     Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(mat);
//     Eigen::VectorXd eigenvalues = eigensolver.eigenvalues().real();
//     Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors().real();
//
//     double max_condition_number = 0.0;
//     for (int i = 0; i < eigenvalues.size(); ++i) {
//         double norm_x = eigenvectors.col(i).norm();
//         double norm_y = eigenvectors.row(i).norm();
//         double xy = std::abs(eigenvectors.col(i).dot(eigenvectors.row(i)));
//         double condition_number = (norm_x * norm_x * norm_y * norm_y) / xy;
//         if (condition_number > max_condition_number) {
//             max_condition_number = condition_number;
//         }
//     }
//     return max_condition_number;
// }

/**
 * @class UnscentedKalmanFilter
 *
 * @tparam StateType Represents the state vector type, should be a column vector
 * of dimension L.
 * @tparam MeasurementType Represents the measurement vector type, should be a
 * column vector of dimension M.
 * @tparam StateCovarianceType Represents the state covariance matrix type,
 * should be a square matrix of dimension LxL.
 * @tparam MeasurementCovarianceType Represents the measurement covariance
 * matrix type, should be a square matrix of dimension MxM.
 *
 * @brief Implements the Unscented Kalman Filter algorithm.
 *
 * @details The Unscented Kalman Filter (UKF) is a state estimation algorithm
 * that uses a deterministic sampling technique to capture the mean and
 * covariance estimates of the state distribution. It propagates the state's
 * mean and covariance through nonlinear transformations without requiring
 * linearization.
 *
 *          The parameters alpha, beta, and kappa control the generation of
 * sigma points, where alpha determines the spread of the sigma points around
 * the mean state, beta incorporates any prior knowledge about the distribution
 * of the state, and kappa is a secondary scaling parameter.
 *
 *          References:
 *          [1] Julier, S.J. and Uhlmann, J.K., 1997. A new extension of the
 * Kalman filter to nonlinear systems. In AeroSense'97 (pp. 182-193).
 * International Society for Optics and Photonics. [2] Wan, E.A. and Van Der
 * Merwe, R., 2000. The unscented Kalman filter for nonlinear estimation. In
 * Adaptive Systems for Signal Processing, Communications, and Control Symposium
 * 2000. AS-SPCC. The IEEE 2000 (pp. 153-158). IEEE.
 */
#include <odin/eigen/inversion_policy.hpp>

template<typename StateType, typename MeasurementType,
         typename StateCovarianceType, typename MeasurementCovarianceType,
         typename InversionPolicy = policy::pinv>
class UnscentedKalmanFilter
    : public KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                          MeasurementCovarianceType> {
public:
    using scalar_t = typename StateType::Scalar;
    using state_t = StateType;
    const static int state_dim = state_t::RowsAtCompileTime;
    using measurement_t = MeasurementType;
    using state_covariance_t = StateCovarianceType;
    using measurement_covariance_t = MeasurementCovarianceType;

    static constexpr int SigmaPointSizeValue =
            SigmaPointSize<state_t::RowsAtCompileTime>::value;

    // Define a method selection variable
    //    static constexpr int inversionMethod = 0; // 0 for normal inverse, 1 for
    //    pseudo-inverse, 2 for regularized inverse

    // Perform static dimension checks at compile time for static matrices.
    // Compile-time checks to ensure matrix/vector dimensions are compatible.
    static_assert(
            static_cast<int>(StateType::RowsAtCompileTime) ==
                            static_cast<int>(StateCovarianceType::RowsAtCompileTime) &&
                    static_cast<int>(StateType::RowsAtCompileTime) ==
                            static_cast<int>(StateCovarianceType::ColsAtCompileTime),
            "StateType dimensions must match with StateCovarianceType");
    static_assert(
            static_cast<int>(MeasurementType::RowsAtCompileTime) ==
                            static_cast<int>(MeasurementCovarianceType::RowsAtCompileTime) &&
                    static_cast<int>(MeasurementType::RowsAtCompileTime) ==
                            static_cast<int>(MeasurementCovarianceType::ColsAtCompileTime),
            "MeasurementType dimensions must match with MeasurementCovarianceType");

    // Various function pointer types used in UKF
    using typename KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                                MeasurementCovarianceType>::StateFunc;
    using typename KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                                MeasurementCovarianceType>::MeasurementFunc;
    using typename KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                                MeasurementCovarianceType>::ProcessCovarianceFunc;
    using typename KalmanFilter<
            StateType, MeasurementType, StateCovarianceType,
            MeasurementCovarianceType>::MeasurementCovarianceFunc;

    /**
   * @brief Constructs an instance of UnscentedKalmanFilter.
   *
   * @param state_transition A function representing the state transition model.
   * @param measurement_model A function representing the measurement model.
   * @param process_covariance A function representing the process covariance
   * model.
   * @param measurement_covariance A function representing the measurement
   * covariance model.
   * @param alpha A parameter for the UKF, controlling the spread of the sigma
   * points around the mean state. Default is 1.0.
   * @param beta A parameter for the UKF, used to incorporate prior knowledge of
   * the distribution of the state. Default is 2.0.
   * @param kappa A parameter for the UKF, used to control the secondary scaling
   * of the sigma points. Default is 0.0.
   *
   */
    UnscentedKalmanFilter(StateFunc state_transition,
                          MeasurementFunc measurement_model,
                          ProcessCovarianceFunc process_covariance,
                          MeasurementCovarianceFunc measurement_covariance,
                          scalar_t alpha = 1.0, scalar_t beta = 2.0,
                          scalar_t kappa = 0.0)
        : KalmanFilter<StateType, MeasurementType, StateCovarianceType,
                       MeasurementCovarianceType>(
                  std::move(state_transition), std::move(measurement_model),
                  std::move(process_covariance), std::move(measurement_covariance)),
          sigma_point_calculator_(alpha, beta, kappa) {}

    /**
   * @brief Destructor for UnscentedKalmanFilter
   */
    ~UnscentedKalmanFilter() override = default;

    //    void estimate_parameters(const Eigen::VectorXd &x, const Eigen::VectorXd
    //    &y) override {
    //        // Implement method here
    //    }

    void set_initial_state(const state_t &initial_state) {
        state_ = initial_state;
    }

    state_t &get_state() { return state_; }

    state_covariance_t &get_state_covariance() { return state_covariance_; }

    void set_initial_covariance(const StateCovarianceType &initial_covariance) {
        state_covariance_ = initial_covariance;
    }

    /**
   * @brief Performs the prediction step in the Unscented Kalman Filter (UKF).
   *
   * @param dt The time step to be used for the prediction.
   *
   * @details This function generates sigma points (matrix X, dimension
   * Lx(2L+1)), propagates them through the state transition function, and
   * recombines them to provide the a priori estimate. This process follows the
   * equations: X_k = chol((L+kappa)*P_k), where chol() is the Cholesky
   * decomposition, and L is the dimension of the state vector, x_{k|k-1} = Σ
   * W_j * x_j, where x_j are the propagated sigma points and W_j are the
   * weights, P_{k|k-1} = Σ W_j * (x_j - x_{k|k-1}) * (x_j - x_{k|k-1})^T + Q_k,
   *          where Q_k is the process noise covariance (dimension LxL).
   */
    void predict(const scalar_t dt) override {
        // Generate sigma points (matrix X, dimension Lx(2L+1))
        sigma_point_calculator_.calculateSigmaPoints(state_, state_covariance_);
        auto sigma_points = sigma_point_calculator_.get_sigma_points();

        // Propagate the sigma points through the state transition function,
        // outputting propagated sigma points (dimension Lx(2L+1)).
        for (int i = 0; i < sigma_points.cols(); ++i) {
            // Note: In template-derived classes, we need to use 'this->' to access
            //       protected members of the base class when the base class is a
            //       dependent base class (i.e., it depends on the template parameters
            //       of the derived class). This is because the C++ compiler does not
            //       look in dependent base classes when resolving non-dependent
            //       names. By using 'this->', we explicitly instruct the compiler
            //       that these members belong to the current object, so it properly
            //       checks the dependent base classes and finds these members. I
            //       noted this as it took me a while to remember this and debug a
            //       "bad function call".
            sigma_points.col(i) = this->state_transition_(sigma_points.col(i), dt);
        }

        // Recombine the propagated sigma points to generate the a priori state
        // estimate (dimension Lx1).
        state_.setZero();
        for (int i = 0; i < sigma_points.cols(); ++i) {
            state_ +=
                    sigma_point_calculator_.get_weights_mean()(i) * sigma_points.col(i);
        }

        // Compute the a priori estimate of the state covariance (dimension LxL).
        state_covariance_.setZero();
        for (int i = 0; i < sigma_points.cols(); ++i) {
            StateType diff = sigma_points.col(i) - state_;
            state_covariance_ += sigma_point_calculator_.get_weights_covariance()(i) *
                                 diff * diff.transpose();
        }

        // Calculate process covariance once (dimension LxL)
        StateCovarianceType process_covariance =
                this->process_covariance_(state_, dt);// add dt here

        // Validate the dimensions of the state_covariance_ and process_covariance
        // for compatibility
        check_dimensions(
                state_covariance_, process_covariance, MatrixOperation::ADD,
                "Inconsistent dimensions of state covariance and process covariance.");

        // Add the process noise covariance to the a priori state covariance
        // estimate.
        state_covariance_ += process_covariance;
    }

    /**
   * @brief Template function to perform the update step in UKF
   * @details Transforms the sigma points through the measurement function,
   *          calculates predicted measurement and its covariance,
   *          calculates cross covariance, Kalman gain and updates the state and
   * covariance This version allows for measurements of dynamic sizes. For
   * fixed-size measurements, the size can be provided as a template argument
   * for more efficient computations. The measurement noise covariance is also
   * provided as a parameter.
   * @tparam MeasurementSize - The size of the measurement. Can be a fixed size
   * or Eigen::Dynamic for dynamic size (i.e. some measurement sources aren't
   * available).
   * @param measurement - The measurement (dimension Mx1)
   * @param measurement_noise - The measurement noise covariance (dimension MxM)
   */
    template<int MeasurementSize>
    void update(const Eigen::Vector<scalar_t, MeasurementSize> &measurement,
                const Eigen::Matrix<scalar_t, MeasurementSize, MeasurementSize>
                        &measurement_noise) {

        // Ensure that the number of sigma points is as expected
        auto sigma_points = sigma_point_calculator_.get_sigma_points();
        assert(sigma_points.cols() == 2 * state_.size() + 1 &&
               "Inconsistent number of sigma points");

        /*
     *  Transform the sigma points through the measurement function h (dimension
     * Mx(2L+1)). Here, compile-time conditions (constexpr) are employed to
     * optimize the initialization or resizing of the matrix Z, depending on
     *  whether MeasurementSize and SigmaPointSizeValue are known at compile
     * time or are dynamic.
     */
        Eigen::Matrix<scalar_t, MeasurementSize, SigmaPointSizeValue> Z;
        if constexpr (MeasurementSize == Eigen::Dynamic &&
                      SigmaPointSizeValue == Eigen::Dynamic) {
            /*
       *  Both MeasurementSize and SigmaPointSizeValue are dynamic, implying
       * that their values are not known at compile time. Therefore, the matrix
       * Z's dimensions are determined based on the runtime dimensions of the
       *  input measurement and sigma points.
       */
            Z.resize(measurement.rows(), sigma_points.cols());
        } else if constexpr (MeasurementSize == Eigen::Dynamic) {
            /*
       *  Only MeasurementSize is dynamic. Hence, Z's row count is set based on
       * the measurement's runtime row count. The column count is set to
       * SigmaPointSizeValue, which is known at compile time.
       */
            Z.resize(measurement.rows(), SigmaPointSizeValue);
        } else if constexpr (SigmaPointSizeValue == Eigen::Dynamic) {
            /*
       *  Only SigmaPointSizeValue is dynamic. Therefore, the number of Z's
       * columns is set based on the sigma points' runtime column count. The row
       * count is set to MeasurementSize, known at compile time.
       */
            Z.resize(MeasurementSize, sigma_points.cols());
        } else {
            /*
       *  Both MeasurementSize and SigmaPointSizeValue are known at compile
       * time. Therefore, we do not need to resize Z at runtime, as its size is
       * fixed.
       */
        }
        /*
     *  TECHNICAL NOTE:
     *
     *  Sigma points are transformed via the measurement model in this loop,
     * with the results stored in the 'Z' matrix. The Eigen operation employed
     * depends on whether 'MeasurementSize' is known at compile-time or runtime:
     *
     *  - Z.col(i): For a runtime-determined 'MeasurementSize', accesses the ith
     * column of 'Z'. May be less efficient if 'MeasurementSize' is a
     * compile-time constant.
     *
     *  - Z.block<MeasurementSize, 1>(0, i): Optimized for a compile-time
     * constant 'MeasurementSize'. Accesses a block of 'Z' starting from the top
     * of the ith column, including 'MeasurementSize' rows and 1 column.
     *
     *  An 'if constexpr' statement selects the appropriate operation based on
     * the properties of 'MeasurementSize'. If 'measurement_model_' can be
     * vectorized, this loop could potentially be eliminated, leveraging
     * parallel computing resources for significant performance gains.
     */
        for (int i = 0; i < sigma_points.cols(); ++i) {
            if constexpr (MeasurementSize == Eigen::Dynamic) {
                Z.col(i) = this->measurement_model_(sigma_points.col(i));
            } else {
                Z.template block<MeasurementSize, 1>(0, i) =
                        this->measurement_model_(sigma_points.col(i));
                //                LOG(INFO) << "Z:
                //                --------------------------------------"; LOG(INFO) <<
                //                "\n" << Z; LOG(INFO) << "Sig:
                //                ------------------------------------"; LOG(INFO) <<
                //                "\n" << sigma_points; LOG(INFO) << "State:
                //                ----------------------------------"; LOG(INFO) << "\n"
                //                << state_; LOG(INFO) << "State Cov:
                //                ------------------------------"; LOG(INFO) << "\n" <<
                //                state_covariance_;
            }
        }
        //        LOG(INFO) <<
        //        "\n\033[0;34m████████████████████████████████████████████████████████\033[0m";

        // (1) Predicted measurement: z = Sum(W_j * Z_j) for all j
        Eigen::Matrix<scalar_t, MeasurementSize, 1> z_pred =
                Z * sigma_point_calculator_.get_weights_mean().asDiagonal() *
                Eigen::Matrix<scalar_t, SigmaPointSizeValue, 1>::Ones();

        // (2) Predicted measurement covariance: S = Sum(W_j * (Z_j - z) * (Z_j -
        // z)^T) + R for all j
        Eigen::Array<scalar_t, MeasurementSize, SigmaPointSizeValue> diff_Z =
                Z.array().colwise() - z_pred.array();
        Eigen::Array<scalar_t, SigmaPointSizeValue, 1> weights =
                sigma_point_calculator_.get_weights_covariance().array();
        Eigen::Array<scalar_t, MeasurementSize, SigmaPointSizeValue>
                diff_Z_weighted =
                        (weights.transpose().replicate(MeasurementSize, 1) * diff_Z)
                                .matrix();
        Eigen::Matrix<scalar_t, MeasurementSize, MeasurementSize> S =
                (diff_Z_weighted.matrix() * diff_Z_weighted.matrix().transpose()) +
                measurement_noise;

        // Check if S is well-conditioned and has no numerical issues
        //        if (!is_well_conditioned(S) || check_numerical_issues(S)) {
        //            LOG(ERROR)
        //                    << "Measurement covariance is singular,
        //                    ill-conditioned, or contains NaN/Inf values. Skipping
        //                    update.";
        //            return;
        //        }

        // Cross covariance matrix: C_xz = Sum(W_j * (X_j - x) * (Z_j - z)^T) for
        // all j

        Eigen::Array<scalar_t, state_dim, SigmaPointSizeValue> diff_X =
                sigma_points.array().colwise() - state_.array();
        Eigen::Matrix<scalar_t, Eigen::Dynamic, MeasurementSize> cross_covariance =
                diff_X.matrix() *
                sigma_point_calculator_.get_weights_covariance().asDiagonal() *
                diff_Z.matrix().transpose();

        //        if (!is_well_conditioned(cross_covariance)) {
        //            LOG(ERROR) << "Cross covariance is ill-conditioned";
        //            LOG(INFO) <<
        //            "sigma_point_calculator_.getWeightsCovariance():\n" <<
        //            sigma_point_calculator_.getWeightsCovariance(); LOG(INFO) <<
        //            "sigma_points: \n" << sigma_points; LOG(INFO) << "state_: \n"
        //            << state_; LOG(INFO) << "diff_Z: \n" << diff_Z; LOG(INFO) <<
        //            "diff_X: \n" << diff_X; LOG(INFO) << "cross_covariance: \n" <<
        //            cross_covariance; return;
        //        }

        // Kalman gain: K = C_xz * S^-1
        Eigen::Matrix<scalar_t, Eigen::Dynamic, MeasurementSize> K =
                cross_covariance * InversionPolicy::invert(S);

        // Check Kalman gain for numerical issues
        //        if (check_numerical_issues(K)) {
        //            LOG(ERROR) << "Kalman gain has NaN or Inf values. Skipping
        //            update."; return;
        //        }

        // Update state
        // x = x + K * (measurement - z_pred) (dimension Lx1)
        //        std::cout << "\nstate_i:\n" << state_ << std::endl;
        state_ += K * (measurement - z_pred);
        //        std::cout << "\nstate_i+1:\n" << state_ << std::endl;

        // Update state covariance
        //         P = P - K * S * K^T (dimension LxL)
        state_covariance_ -= K * S * K.transpose();

        // Calculate state covariance update
        //        Eigen::MatrixXd updated_state_covariance = state_covariance_ - K *
        //        S * K.transpose();

        // Check if conditioning of updated_state_covariance is acceptable, and
        // check for numerical issues
        //        if (!is_well_conditioned(updated_state_covariance) ||
        //        check_numerical_issues(updated_state_covariance)) {
        //            LOG(ERROR)
        //                    << "Updated state covariance is singular,
        //                    ill-conditioned, or contains NaN/Inf values. Skipping
        //                    covariance update.";
        //        } else {
        //            state_covariance_ = updated_state_covariance;
        //        }

        /* Uncomment this block to use full diagnostic compilation

    // Get the full conditioning result
    ConditioningResult cov_check_result =
    check_conditioning(updated_state_covariance, "updated_state_covariance");

    // Check if conditioning of updated_state_covariance is acceptable, and
    check for numerical issues if (!cov_check_result.is_well_conditioned ||
    check_numerical_issues(updated_state_covariance)) {

        // Throw an error if necessary, otherwise log it
        if (cov_check_result.should_throw) {
            throw std::runtime_error(cov_check_result.error_message.str());
        } else {
            LOG(ERROR) << cov_check_result.error_message.str();
        }

    } else {
        state_covariance_ = updated_state_covariance;
    }

    */
    }

private:
    state_t state_;                      // state
    state_covariance_t state_covariance_;// covariance of the state
    UKFSigmaPointCalculator<StateType::RowsAtCompileTime, scalar_t>
            sigma_point_calculator_;// sigma point calculator
};

#endif// KALMAN_HPP
