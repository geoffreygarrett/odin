/**
 * @file sigma_point_calculator.hpp
 * @brief Provides Sigma Point calculation for Kalman Filters.
 * @details This file contains the interfaces and implementation for calculating Sigma Points, which is a crucial
 * part of the Unscented Kalman Filter (UKF). The Sigma Points are representative points of the state probability
 * distribution and are used for propagating the mean and covariance through the non-linear transformation.
 * @author Geoffrey
 * @date 2023-05-31
 */

#ifndef SIGMA_POINT_CALCULATOR_HPP
#define SIGMA_POINT_CALCULATOR_HPP

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <odin/eigen.hpp>

/**
 * @brief Specialization structure to define the size of Sigma Points
 *
 * @tparam StateSize State Dimension
 */
template<int StateSize>
struct SigmaPointSize {
    static constexpr int value = 2 * StateSize + 1;
};

template<>
struct SigmaPointSize<Eigen::Dynamic> {
    static constexpr int value = Eigen::Dynamic;
};

/**
 * @class SigmaPointCalculator
 * @brief Abstract base class for Sigma Point calculators.
 *
 * This class provides a generic interface for Sigma Point calculators, using the Curiously Recurring Template Pattern (CRTP)
 * for efficient polymorphism.
 *
 * CRTP is used here to avoid the performance penalty of traditional dynamic polymorphism while still providing
 * a similar interface structure. Derived classes override the calculateSigmaPoints method to provide
 * their own implementation, and the get* methods call these implementations.
 *
 * @tparam Derived Derived class type for CRTP
 * @tparam StateSize Size of the state (default: Dynamic)
 */
template<typename Derived, int StateSize, typename Scalar = double>
class SigmaPointCalculator {
public:
    static constexpr int SigmaPointSizeValue = SigmaPointSize<StateSize>::value;

    using scalar_t = Scalar;
    using state_t = Eigen::Vector<scalar_t, StateSize>;
    using covariance_t = Eigen::Matrix<scalar_t, StateSize, StateSize>;
    using sigma_points_t = Eigen::Matrix<scalar_t, StateSize, SigmaPointSizeValue>;
    using weights_t = Eigen::Vector<scalar_t, SigmaPointSizeValue>;

    virtual ~SigmaPointCalculator() = default;

    /**
     * @brief Pure virtual method for calculating the Sigma Points based on the state vector and covariance matrix.
     *
     * @param x State Vector
     * @param P Covariance Matrix
     */
    virtual void calculateSigmaPoints(const state_t &x, const covariance_t &P) = 0;

    /**
     * @brief Retrieve the calculated sigma points.
     *
     * This method utilizes the CRTP to call the derived class' implementation.
     */
    [[nodiscard]] sigma_points_t get_sigma_points() const {
        return static_cast<const Derived *>(this)->getSigmaPointsImpl();
    }

    /**
     * @brief Retrieve the weights for mean calculation.
     *
     * This method utilizes the CRTP to call the derived class' implementation.
     */
    [[nodiscard]] weights_t get_weights_mean() const {
        return static_cast<const Derived *>(this)->getWeightsMeanImpl();
    }

    /**
     * @brief Retrieve the weights for covariance calculation.
     *
     * This method utilizes the CRTP to call the derived class' implementation.
     */
    [[nodiscard]] weights_t get_weights_covariance() const {
        return static_cast<const Derived *>(this)->getWeightsCovarianceImpl();
    }
};


/**
 * @class UKFSigmaPointCalculator
 * @brief A class for calculating Sigma points for the Unscented Kalman Filter (UKF).
 *
 * The class is templated on the state size, with a default of dynamic sizing.
 * Three parameters, alpha, beta, and kappa, control the dispersion of the sigma points and the state distribution.
 * Typical values are alpha = 1e-3, kappa = 1, beta = 2, but these may need to be adjusted for non-linear distributions.
 * For Gaussian distributions, beta = 2 is optimal.
 *
 * @tparam StateSize The dimensionality of the state.
 */
template<int StateSize = Eigen::Dynamic, typename Scalar = double>
class UKFSigmaPointCalculator
        : public SigmaPointCalculator<UKFSigmaPointCalculator<StateSize, Scalar>, StateSize, Scalar> {
public:
    using Base = SigmaPointCalculator<UKFSigmaPointCalculator<StateSize, Scalar>, StateSize, Scalar>;
    using scalar_t = typename Base::scalar_t;
    using typename Base::state_t;
    using typename Base::covariance_t;
    using typename Base::sigma_points_t;
    using typename Base::weights_t;

    explicit UKFSigmaPointCalculator(scalar_t alpha = 1.0, scalar_t beta = 2.0, scalar_t kappa = 0.0)
            : alpha_(alpha), beta_(beta), kappa_(kappa), lambda_(alpha * alpha * (StateSize + kappa) - StateSize) {
        initializeWeights();
    }

    void calculateSigmaPoints(const state_t &x, const covariance_t &P) override {
        calculateSigmaPointsImpl(x, P);
    }

    sigma_points_t getSigmaPointsImpl() const {
        return sigma_points_;
    }

    weights_t getWeightsMeanImpl() const {
        return W_m;
    }

    weights_t getWeightsCovarianceImpl() const {
        return W_c;
    }

private:
    void initializeWeights() {
        scalar_t weight_0_m = lambda_ / (StateSize + lambda_);
        scalar_t weight_0_c = lambda_ / (StateSize + lambda_) + (1 - alpha_ * alpha_ + beta_);
        scalar_t weight_i = scalar_t(0.5) / (StateSize + lambda_);

        W_m = weights_t::Constant(weight_i);
        W_c = weights_t::Constant(weight_i);

        W_m(0) = weight_0_m;
        W_c(0) = weight_0_c;
    }

    void calculateSigmaPointsImpl(const state_t &x, const covariance_t &P) {
        sigma_points_ = sigma_points_t::Zero();
        covariance_t A = P.llt().matrixL();
        sigma_points_.col(0) = x;

        scalar_t sqrt_val = std::sqrt(StateSize + lambda_);

        for (int i = 0; i < StateSize; ++i) {
            sigma_points_.col(i + 1) = x + sqrt_val * A.col(i);
            sigma_points_.col(i + 1 + StateSize) = x - sqrt_val * A.col(i);
        }
    }

    scalar_t alpha_, beta_, lambda_;
    [[maybe_unused]] scalar_t kappa_;
    sigma_points_t sigma_points_;
    weights_t W_m;
    weights_t W_c;
};


/**
 * @class UKFSigmaPointCalculator<Eigen::Dynamic>
 * @brief Specialization of UKFSigmaPointCalculator for the case of dynamic state size.
 *
 * @details This specialization handles the case where the size of the state vector isn't known at compile time.
 * It performs the unscented transformation for the UKF, as explained in Julier, Uhlmann, and Durrant-Whyte's paper.
 * The dimensionality of the state is passed to the constructor and used in the calculation of Sigma Points and weights.
 */
template<>
class UKFSigmaPointCalculator<Eigen::Dynamic>
        : public SigmaPointCalculator<UKFSigmaPointCalculator<Eigen::Dynamic>, Eigen::Dynamic> {
public:
    /**
     * @brief Construct a new UKFSigmaPointCalculator object for dynamic state size.
     *
     * @param state_dimension Dimensionality of the state
     * @param alpha Dispersion parameter, determines the spread of the sigma points around the mean state. Usually set to a small positive value (e.g., 1 <= alpha <= 1e-4).
     * @param beta Distribution parameter, used to incorporate prior knowledge of the distribution of x (for Gaussian distributions, beta=2 is optimal).
     * @param kappa Secondary scaling parameter, usually set to 0 or 3 - L.
     */
    explicit UKFSigmaPointCalculator(int state_dimension, double alpha = 1.0, double beta = 2.0, double kappa = 0.0)
            : state_dimension_(state_dimension), alpha_(alpha), beta_(beta), kappa_(kappa),
              lambda_(alpha * alpha * (state_dimension + kappa) - state_dimension) {
        initializeWeights();
    }

    /**
     * @brief Calculate the Sigma Points based on the given state vector and covariance matrix.
     *
     * @param x State Vector (mean state)
     * @param P Covariance Matrix
     */
    void calculateSigmaPoints(const Eigen::VectorXd &x, const Eigen::MatrixXd &P) override {
        calculateSigmaPointsImpl(x, P);
    }

    /**
     * @brief Retrieve the calculated Sigma Points.
     * @return Calculated Sigma Points
     */
    [[nodiscard]] Eigen::MatrixXd getSigmaPointsImpl() const {
        return sigma_points_;
    }

    /**
     * @brief Retrieve the weights for mean calculation.
     * @return Weights for mean calculation
     */
    [[nodiscard]] Eigen::VectorXd getWeightsMeanImpl() const {
        return W_m;
    }

    /**
     * @brief Retrieve the weights for covariance calculation.
     * @return Weights for covariance calculation
     */
    [[nodiscard]] Eigen::VectorXd getWeightsCovarianceImpl() const {
        return W_c;
    }


private:
    /**
     * @brief Initialize the weights for Sigma Points calculation.
     *
     * @details According to the unscented transformation algorithm, the weights are calculated as follows:
     * W_0^(m) = lambda / (L + lambda)
     * W_0^(c) = lambda / (L + lambda) + (1 - alpha^2 + beta)
     * W_i^(m) = W_i^(c) = 1 / {2(L + lambda)} for i = 1, ..., 2L
     */
    void initializeWeights() {
        double weight_0_m = lambda_ / (state_dimension_ + lambda_);
        double weight_0_c = lambda_ / (state_dimension_ + lambda_) + (1 - alpha_ * alpha_ + beta_);
        double weight_i = 0.5 / (state_dimension_ + lambda_);

        W_m = Eigen::VectorXd::Constant(2 * state_dimension_ + 1, weight_i);
        W_c = Eigen::VectorXd::Constant(2 * state_dimension_ + 1, weight_i);

        W_m(0) = weight_0_m;
        W_c(0) = weight_0_c;
    }

    /**
     * @brief Calculate the Sigma Points based on the given state vector and covariance matrix.
     *
     * @details According to the unscented transformation algorithm, the sigma points are calculated as follows:
     * X_0 = x_bar
     * X_i = x_bar + sqrt((L + lambda) * P_x)_i for i = 1, ..., L
     * X_i = x_bar - sqrt((L + lambda) * P_x)_{i - L} for i = L + 1, ..., 2L
     *
     * @param x State Vector (mean state)
     * @param P Covariance Matrix
     */
    void calculateSigmaPointsImpl(const Eigen::VectorXd &x, const Eigen::MatrixXd &P) {
        sigma_points_ = Eigen::MatrixXd::Zero(state_dimension_, 2 * state_dimension_ + 1);
        Eigen::MatrixXd A = P.llt().matrixL();
        sigma_points_.col(0) = x;

        for (int i = 0; i < state_dimension_; ++i) {
            sigma_points_.col(i + 1) = x + std::sqrt(state_dimension_ + lambda_) * A.col(i);
            sigma_points_.col(i + 1 + state_dimension_) = x - std::sqrt(state_dimension_ + lambda_) * A.col(i);
        }
    }

    int state_dimension_;
    double alpha_, beta_, lambda_;
    [[maybe_unused]] double kappa_;
    Eigen::MatrixXd sigma_points_;
    Eigen::VectorXd W_m;
    Eigen::VectorXd W_c;
};

#endif // SIGMA_POINT_CALCULATOR_HPP