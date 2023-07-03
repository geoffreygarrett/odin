
#include "include/odin/estimation/filters/base_filter.hpp"
#include <eigen/inversion_policy.hpp>
#include <utilities/sigma_point_calculator.hpp>

template<typename X, typename M, typename P, typename R, typename InversionPolicy = policy::pinv<typename X::Scalar>>
class UnscentedKalmanFilter : public Filter<UnscentedKalmanFilter<X, M, P, R, InversionPolicy>, X, M, P> {
public:
    static constexpr int SigmaPointSizeValue = SigmaPointSize<X::RowsAtCompileTime>::value;

    using scalar_t = X::Scalar;
    using StateFunc = std::function<X(const X &, scalar_t)>;
    using MeasurementFunc = std::function<M(const X &)>;
    using ProcessCovarianceFunc = std::function<P(const X &, scalar_t)>;
    using MeasurementCovarianceFunc = std::function<R(const X &)>;

    UnscentedKalmanFilter(
            StateFunc state_transition,
            MeasurementFunc measurement_model,
            ProcessCovarianceFunc process_covariance,
            MeasurementCovarianceFunc measurement_covariance,
            scalar_t alpha = 1.0,
            scalar_t beta = 2.0,
            scalar_t kappa = 0.0) : Filter<UnscentedKalmanFilter<X, M, P, R, InversionPolicy>, X, M, P>(
            std::move(state_transition),
            std::move(measurement_model),
            std::move(process_covariance),
            std::move(measurement_covariance)),
                                  sigma_point_calculator_(alpha, beta, kappa) {}


    /**
     * @brief Destructor for UnscentedKalmanFilter
     */
    ~UnscentedKalmanFilter() override = default;

    void estimate_parameters(const Eigen::VectorXd &x, const Eigen::VectorXd &y) override {
        // Implement method here
    }

    void set_initial_state(const StateType &initial_state) { state_ = initial_state; }

    StateType &get_state() { return state_; }

    StateCovarianceType &get_state_covariance() { return state_covariance_; }


    void
    set_initial_covariance(const StateCovarianceType &initial_covariance) { state_covariance_ = initial_covariance; }


    void predict(const double dt) override {
        // Generate sigma points (matrix X, dimension Lx(2L+1))
        sigma_point_calculator_.calculateSigmaPoints(state_, state_covariance_);
        auto sigma_points = sigma_point_calculator_.getSigmaPoints();

        for (int i = 0; i < sigma_points.cols(); ++i) {

            sigma_points.col(i) = this->state_transition_(sigma_points.col(i), dt);
        }

        // Recombine the propagated sigma points to generate the a priori state estimate (dimension Lx1).
        state_.setZero();
        for (int i = 0; i < sigma_points.cols(); ++i) {
            state_ += sigma_point_calculator_.getWeightsMean()(i) * sigma_points.col(i);
        }

        // Compute the a priori estimate of the state covariance (dimension LxL).
        state_covariance_.setZero();
        for (int i = 0; i < sigma_points.cols(); ++i) {
            StateType diff = sigma_points.col(i) - state_;
            state_covariance_ += sigma_point_calculator_.getWeightsCovariance()(i) * diff * diff.transpose();
        }

        // Calculate process covariance once (dimension LxL)
        StateCovarianceType process_covariance = this->process_covariance_(state_, dt);  // add dt here

        // Validate the dimensions of the state_covariance_ and process_covariance for compatibility
        check_dimensions(state_covariance_, process_covariance, MatrixOperation::ADD,
                         "Inconsistent dimensions of state covariance and process covariance.");

        // Add the process noise covariance to the a priori state covariance estimate.
        state_covariance_ += process_covariance;
    }


    template<int MeasurementSize>
    void update(const Eigen::Matrix<double, MeasurementSize, 1> &measurement,
                const Eigen::Matrix<double, MeasurementSize, MeasurementSize> &measurement_noise) {

        // Ensure that the number of sigma points is as expected
        auto sigma_points = sigma_point_calculator_.getSigmaPoints();
        assert(sigma_points.cols() == 2 * state_.size() + 1 && "Inconsistent number of sigma points");

        /*
         *  Transform the sigma points through the measurement function h (dimension Mx(2L+1)). Here, compile-time
         *  conditions (constexpr) are employed to optimize the initialization or resizing of the matrix Z, depending on
         *  whether MeasurementSize and SigmaPointSizeValue are known at compile time or are dynamic.
         */
        Eigen::Matrix<double, MeasurementSize, SigmaPointSizeValue> Z;
        if constexpr (MeasurementSize == Eigen::Dynamic && SigmaPointSizeValue == Eigen::Dynamic) {
            /*
             *  Both MeasurementSize and SigmaPointSizeValue are dynamic, implying that their values are not known at
             *  compile time. Therefore, the matrix Z's dimensions are determined based on the runtime dimensions of the
             *  input measurement and sigma points.
             */
            Z.resize(measurement.rows(), sigma_points.cols());
        } else if constexpr (MeasurementSize == Eigen::Dynamic) {
            /*
             *  Only MeasurementSize is dynamic. Hence, Z's row count is set based on the measurement's runtime row
             *  count. The column count is set to SigmaPointSizeValue, which is known at compile time.
             */
            Z.resize(measurement.rows(), SigmaPointSizeValue);
        } else if constexpr (SigmaPointSizeValue == Eigen::Dynamic) {
            /*
             *  Only SigmaPointSizeValue is dynamic. Therefore, the number of Z's columns is set based on the sigma
             *  points' runtime column count. The row count is set to MeasurementSize, known at compile time.
             */
            Z.resize(MeasurementSize, sigma_points.cols());
        } else {
            /*
             *  Both MeasurementSize and SigmaPointSizeValue are known at compile time.
             *  Therefore, we do not need to resize Z at runtime, as its size is fixed.
             */
        }
        /*
         *  TECHNICAL NOTE:
         *
         *  Sigma points are transformed via the measurement model in this loop, with the results stored in the 'Z'
         *  matrix. The Eigen operation employed depends on whether 'MeasurementSize' is known at compile-time or
         *  runtime:
         *
         *  - Z.col(i): For a runtime-determined 'MeasurementSize', accesses the ith column of 'Z'.
         *               May be less efficient if 'MeasurementSize' is a compile-time constant.
         *
         *  - Z.block<MeasurementSize, 1>(0, i): Optimized for a compile-time constant 'MeasurementSize'.
         *                                        Accesses a block of 'Z' starting from the top of the ith column,
         *                                        including 'MeasurementSize' rows and 1 column.
         *
         *  An 'if constexpr' statement selects the appropriate operation based on the properties of 'MeasurementSize'.
         *  If 'measurement_model_' can be vectorized, this loop could potentially be eliminated, leveraging parallel
         *  computing resources for significant performance gains.
         */
        for (int i = 0; i < sigma_points.cols(); ++i) {
            if constexpr (MeasurementSize == Eigen::Dynamic) {
                Z.col(i) = this->measurement_model_(sigma_points.col(i));
            } else {
                Z.template block<MeasurementSize, 1>(0, i) = this->measurement_model_(sigma_points.col(i));
            }
        }

        // (1) Predicted measurement: z = Sum(W_j * Z_j) for all j
        Eigen::Matrix<double, MeasurementSize, 1> z_pred =
                Z * sigma_point_calculator_.getWeightsMean().asDiagonal() *
                Eigen::Matrix<double, SigmaPointSizeValue, 1>::Ones();

        // (2) Predicted measurement covariance: S = Sum(W_j * (Z_j - z) * (Z_j - z)^T) + R for all j
        Eigen::Array<double, MeasurementSize, SigmaPointSizeValue> diff_Z = Z.array().colwise() - z_pred.array();
        Eigen::Array<double, SigmaPointSizeValue, 1> weights = sigma_point_calculator_.getWeightsCovariance().array();
        Eigen::Array<double, MeasurementSize, SigmaPointSizeValue> diff_Z_weighted = (
                weights.transpose().replicate(MeasurementSize, 1) * diff_Z).matrix();
        Eigen::Matrix<double, MeasurementSize, MeasurementSize> S =
                (diff_Z_weighted.matrix() * diff_Z_weighted.matrix().transpose()) + measurement_noise;

        // Cross covariance matrix: C_xz = Sum(W_j * (X_j - x) * (Z_j - z)^T) for all j
        Eigen::ArrayXXd diff_X = sigma_points.array().colwise() - state_.array();
        Eigen::Matrix<double, Eigen::Dynamic, MeasurementSize> cross_covariance =
                diff_X.matrix() * sigma_point_calculator_.getWeightsCovariance().asDiagonal() *
                diff_Z.matrix().transpose();

        // Kalman gain: K = C_xz * S^-1
        Eigen::Matrix<double, Eigen::Dynamic, MeasurementSize> K = cross_covariance * InversionPolicy::invert(S);

        // Update state
        // x = x + K * (measurement - z_pred) (dimension Lx1)
        state_ += K * (measurement - z_pred);

        // Update state covariance
        // P = P - K * S * K^T (dimension LxL)
        state_covariance_ -= K * S * K.transpose();
    }


private:
    StateType state_; // state
    StateCovarianceType state_covariance_; // covariance of the state
    UKFSigmaPointCalculator<StateType::RowsAtCompileTime> sigma_point_calculator_; // sigma point calculator
};
