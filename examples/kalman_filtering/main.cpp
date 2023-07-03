#include <fstream>
#include <odin/core/numerics/integrators/runge_kutta.hpp>
#include <odin/kalman.hpp>
#include <odin/library.hpp>
#include <random>

#define SCALAR double

// ------------------- Constants -------------------
static constexpr SCALAR R = 6371;       // Earth's radius in kilometers
static constexpr SCALAR G = 6.67430e-20;// Gravitational constant (in km^3 kg^-1 s^-2)
static constexpr SCALAR m = 5.972e24;   // Earth's mass (in kg)

// ------------------- Initial State -------------------
Eigen::Vector<SCALAR, 4> true_state{10000, 0, 0,
                                    7.5};// Initial state of the object (position in km and velocity in km/s)

// ------------------- State Dynamics Function -------------------
Eigen::Vector<SCALAR, 4> state_dynamics(const Eigen::Vector<SCALAR, 4> &x, SCALAR) {
    SCALAR r = std::sqrt(std::pow(x[0], 2) + std::pow(x[1], 2));// distance from origin
    SCALAR a_x = -G * m * x[0] / std::pow(r, 3);                // Calculate accelerations due to gravity
    SCALAR a_y = -G * m * x[1] / std::pow(r, 3);
    return {x[2], x[3], a_x, a_y};// Return state derivative
}

auto state_dynamics_bound = std::bind(state_dynamics, std::placeholders::_1, 0);

// ------------------- Dynamics Function Structure -------------------
struct DynamicsFunction {
    Eigen::Vector<SCALAR, 4> operator()(const Eigen::Vector<SCALAR, 4> &state, SCALAR t) const {
        return state_dynamics(state, t);
    }
};


auto rk4 = RungeKuttaIntegrator<
        Eigen::Vector<SCALAR, 4>,
        DynamicsFunction,
        butcher_tableau::RK4<SCALAR>>();
//        butcher_tableau::DormandPrince<SCALAR>>();

// ------------------- State Transition Function -------------------
Eigen::Vector<SCALAR, 4> state_transition(const Eigen::Vector<SCALAR, 4> &x, SCALAR dt) {
    return rk4.step(x, dt);
}

// helper to calculate range rate between two
/**
 * The measurement model. Computes the simulated measurement for the given state.
 *
 * @param x The state vector.
 * @return The simulated measurement vector.
 */
// Helper function to calculate range
namespace observables {
    struct range {
    };
    struct range_rate {
    };
}// namespace observables

// Measurement model struct
template<typename T, typename Scalar, int Dim>
struct MeasurementModel {
};

// Specialization for Range
template<typename Scalar, int Dim>
struct MeasurementModel<observables::range, Scalar, Dim> {
public:
    static SCALAR
    calculate(const Eigen::Matrix<Scalar, Dim, 1> &position1, const Eigen::Matrix<Scalar, Dim, 1> &position2) {
        return (position1 - position2).norm();
    }
};

// Specialization for RangeRate
template<typename Scalar, int Dim>
struct MeasurementModel<observables::range_rate, Scalar, Dim> {
public:
    static SCALAR calculate(const Eigen::Matrix<Scalar, Dim, 1> &state1,
                            const Eigen::Matrix<Scalar, Dim, 1> &state2) {
        // Extract the position and velocity vectors
        Eigen::Matrix<Scalar, Dim / 2, 1> position1 = state1.head(Dim / 2);
        Eigen::Matrix<Scalar, Dim / 2, 1> position2 = state2.head(Dim / 2);
        Eigen::Matrix<Scalar, Dim / 2, 1> velocity1 = state1.tail(Dim / 2);
        Eigen::Matrix<Scalar, Dim / 2, 1> velocity2 = state2.tail(Dim / 2);

        // Compute the relative position and velocity
        Eigen::Matrix<Scalar, Dim / 2, 1> rel_position = position1 - position2;
        Eigen::Matrix<Scalar, Dim / 2, 1> rel_velocity = velocity1 - velocity2;

        // Compute the unit vector from observer to target
        Eigen::Matrix<Scalar, Dim / 2, 1> unit_vector = rel_position.normalized();

        // Compute the range rate
        SCALAR range_rate = rel_velocity.dot(unit_vector);

        return range_rate;
    }
};


// Positions of four tracking stations
static Eigen::Vector<SCALAR, 2> ts1{-R, 0};
static Eigen::Vector<SCALAR, 2> ts2{R, 0};
static Eigen::Vector<SCALAR, 2> ts3{0, R};
static Eigen::Vector<SCALAR, 2> ts4{0, -R};


Eigen::Vector<SCALAR, 8> measurement_model(const Eigen::Vector<SCALAR, 4> &x) {
    Eigen::Vector<SCALAR, 8> measurement;

    Eigen::Vector<SCALAR, 2> spacecraft_position(x[0], x[1]);
    Eigen::Vector<SCALAR, 2> spacecraft_velocity(x[2], x[3]);
    Eigen::Vector<SCALAR, 2> station_velocity(0, 0);// Assuming stationary ground stations

    measurement << MeasurementModel<observables::range, SCALAR, 2>::calculate(spacecraft_position, ts1),
            MeasurementModel<observables::range_rate, SCALAR, 2>::calculate(spacecraft_velocity, station_velocity),
            MeasurementModel<observables::range, SCALAR, 2>::calculate(spacecraft_position, ts2),
            MeasurementModel<observables::range_rate, SCALAR, 2>::calculate(spacecraft_velocity, station_velocity),
            MeasurementModel<observables::range, SCALAR, 2>::calculate(spacecraft_position, ts3),
            MeasurementModel<observables::range_rate, SCALAR, 2>::calculate(spacecraft_velocity, station_velocity),
            MeasurementModel<observables::range, SCALAR, 2>::calculate(spacecraft_position, ts4),
            MeasurementModel<observables::range_rate, SCALAR, 2>::calculate(spacecraft_velocity, station_velocity);

    return measurement;
}

/**
 * Computes the measurement noise covariance matrix R.
 *
 * @return The measurement noise covariance matrix R.
 */
template<typename T, int Dim>
struct MeasurementNoise {
};

// Specialization for Range
template<int Dim>
struct MeasurementNoise<observables::range, Dim> {
public:
    static Eigen::Matrix<SCALAR, Dim, Dim> get_noise() {
        return Eigen::Matrix<SCALAR, Dim, Dim>::Identity() * 0.02;
    }
};

// Specialization for RangeRate
template<int Dim>
struct MeasurementNoise<observables::range_rate, Dim> {
public:
    static Eigen::Matrix<SCALAR, Dim, Dim> get_noise() {
        return Eigen::Matrix<SCALAR, Dim, Dim>::Identity() * 0.01;
    }
};

Eigen::Matrix<SCALAR, 8, 8> measurement_covariance(const Eigen::Vector<SCALAR, 8> &) {
    Eigen::Matrix<SCALAR, 8, 8> noise_matrix = Eigen::Matrix<SCALAR, 8, 8>::Zero();

    noise_matrix.block<2, 2>(0, 0) = MeasurementNoise<observables::range, 2>::get_noise();
    noise_matrix.block<2, 2>(2, 2) = MeasurementNoise<observables::range_rate, 2>::get_noise();
    noise_matrix.block<2, 2>(4, 4) = MeasurementNoise<observables::range, 2>::get_noise();
    noise_matrix.block<2, 2>(6, 6) = MeasurementNoise<observables::range_rate, 2>::get_noise();

    return noise_matrix;
}


/**
 * Computes the process noise covariance matrix Q.
 *
 * @return The process noise covariance matrix Q.
 */
Eigen::Matrix<SCALAR, 4, 4> process_covariance(const Eigen::Vector<SCALAR, 4> &, SCALAR) {
    SCALAR pos_noise_stddev = 0.1;   // meters
                                     //    SCALAR pos_noise_stddev = 50.0;  // meters
    SCALAR vel_noise_stddev = 0.0001;// meters/second
                                     //    SCALAR vel_noise_stddev = 1;  // meters/second

    // Convert position noise from m^2 to km^2
    SCALAR pos_noise_km = std::pow(pos_noise_stddev / 1000.0, 2);// kilometers^2

    // Convert velocity noise from (m/s)^2 to (km/s)^2
    SCALAR vel_noise_km_per_s = std::pow(vel_noise_stddev / 1000.0, 2);// (kilometers/second)^2

    // Define the process noise covariance matrix
    Eigen::Matrix<SCALAR, 4, 4> Q = Eigen::Matrix<SCALAR, 4, 4>::Zero();// Initialize process noise covariance matrix

    // Diagonal entries
    Q(0, 0) = pos_noise_km;      // Variance for x position in km^2
    Q(1, 1) = pos_noise_km;      // Variance for y position in km^2
    Q(2, 2) = vel_noise_km_per_s;// Variance for x velocity in (km/s)^2
    Q(3, 3) = vel_noise_km_per_s;// Variance for y velocity in (km/s)^2
    return Q;
}

// Column widths - adjust as necessary for your data
const int iter_width = 15;
const int state_width = 20;
const int cov_width = 20;


/**
 * @class SundmanTransformPolicy
 * @brief A policy class to adjust the time step size according to the Sundman transformation.
 *
 * The Sundman transformation regularizes the singular three-body problem by introducing a new time variable.
 * The transformation guarantees an asymptotically slower flow near singularities and is useful in trajectory
 * propagation if equally spaced segments are considered in the new time variable's domain.
 * This policy class can be used in a generic simulation library to apply the Sundman transformation.
 *
 * @tparam StateType The type of the state variable in the simulation. Defaults to Eigen::VectorXd.
 */
//template<typename StateType = Eigen::VectorXd>
//struct SundmanTransformPolicy {
//    using state_type = StateType;
//    using step_size_type = typename StateType::Scalar;
//
//    /**
//     * @brief Adjusts the step size according to the Sundman transformation.
//     *
//     * The transformation is defined by ds = dt / r, where r is the distance from the singularity,
//     * ds is the transformed time step, and dt is the original time step. This method calculates
//     * the transformed time step.
//     *
//     * @param x The current state of the simulation.
//     * @param ds The original step size.
//     * @return The adjusted step size according to the Sundman transformation.
//     */
//    step_size_type adjust_step_size(const state_type &x, step_size_type ds) const {
//        SCALAR r = sqrt(x(0) * x(0) + x(1) * x(1) + x(2) * x(2));
//        return ds * r;
//    }
//
//    /**
//     * @brief Calculates the initial step size ds for a given r and dt.
//     *
//     * This function helps to establish a compatible flow rate between the original and transformed time variable
//     * at the start of the simulation.
//     *
//     * @param r The distance from the singularity.
//     * @param dt The original time step.
//     * @return The initial step size ds.
//     */
//    step_size_type initial_step_size(step_size_type r, step_size_type dt) const {
//        return dt / r;
//    }
//};

// Define a template struct for parameter sets
template<typename T>
struct ParameterSet {
    T alpha;
    T beta;
    T kappa;

    ParameterSet(T alpha, T beta, T kappa)
        : alpha(alpha), beta(beta), kappa(kappa) {}
};

// Define different parameter sets
ParameterSet<SCALAR> set1(0.003, 2.0, 0.0);
ParameterSet<SCALAR> set2(1.0, 2.0, 3.0);
ParameterSet<SCALAR> set3(0.001, 2.0, 1.0);

// Function to return the required parameter set
template<typename T>
ParameterSet<T> get_parameter_set(int set_num) {
    switch (set_num) {
        case 1:
            return set1;
        case 2:
            return set2;
        case 3:
            return set3;
        default:
            throw std::invalid_argument("Invalid set number");
    }
}

#include <memory>
#include <tuple>


template<typename T>
std::pair<std::shared_ptr<std::ofstream>, ParameterSet<T>> initiate_file_and_set(int set_num) {
    std::string type;

    if (typeid(T) == typeid(float)) {
        type = "float";
    } else if (typeid(T) == typeid(double)) {
        type = "double";
    } else if (typeid(T) == typeid(long double)) {
        type = "long_double";
    }

    auto output_file = std::make_shared<std::ofstream>(
            "trajectory_" + type + "_set" + std::to_string(set_num) + ".csv");
    return std::make_pair(output_file, get_parameter_set<T>(set_num));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

class CSVWriter {
    std::ofstream &output_file;
    std::vector<int> column_widths;
    int precision;

public:
    CSVWriter(std::ofstream &output_file, int precision)
        : output_file(output_file), precision(precision) {}

    void write_header(const std::vector<std::string> &headers) {
        column_widths.resize(headers.size());

        for (size_t i = 0; i < headers.size(); ++i) {
            column_widths[i] = headers[i].size();
        }

        write_row(headers);
    }

    void write_row(const std::vector<std::string> &row) {
        for (size_t i = 0; i < row.size(); ++i) {
            column_widths[i] = std::max(column_widths[i], static_cast<int>(row[i].size()));
        }

        write_row_to_file(row);
    }

    template<typename T>
    void write_row(const std::vector<T> &row) {
        std::vector<std::string> row_strings(row.size());

        for (size_t i = 0; i < row.size(); ++i) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(precision) << row[i];
            row_strings[i] = ss.str();
        }

        write_row(row_strings);
    }

private:
    void write_row_to_file(const std::vector<std::string> &row) {
        for (size_t i = 0; i < row.size(); ++i) {
            output_file << std::setw(column_widths[i]) << row[i];
            if (i < row.size() - 1) {
                output_file << ',';
            }
        }
        output_file << '\n';
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// typename StateType, typename MeasurementType, typename StateCovarianceType, typename MeasurementCovarianceType
int main(int argc, char *argv[]) {
#ifdef USE_GLOG
    google::InitGoogleLogging(argv[0]);
#endif

    // Setup initial conditions
    auto [output_file, set] = initiate_file_and_set<SCALAR>(1);
    const SCALAR delta_time = 100.0;
    SCALAR current_time = 0.0;// Initialize current time

    // Setup UKF parameters
    const SCALAR ALPHA = set.alpha;
    const SCALAR BETA = set.beta;
    const SCALAR KAPPA = set.kappa;

    // Initialize the UKF
    UnscentedKalmanFilter<Eigen::Vector<SCALAR, 4>, Eigen::Vector<SCALAR, 8>, Eigen::Matrix<SCALAR, 4, 4>, Eigen::Matrix<SCALAR, 8, 8>> ukf(
            state_transition,
            measurement_model,
            process_covariance,
            measurement_covariance,
            ALPHA,// alpha
            BETA, // beta
            KAPPA // kappa
    );

    ukf.set_initial_state(true_state);
    Eigen::Matrix<SCALAR, 4, 4> initial_covariance = process_covariance(true_state, 0.0);
    ukf.set_initial_covariance(initial_covariance);

    //    // Write header to the output file
    //    (*output_file)
    //            << "iteration,actual_x,actual_y,actual_vx,actual_vy,estimate_x,estimate_y,estimate_vx,estimate_vy,cov_xx,cov_yy,cov_xy,position_error_norm,velocity_error_norm,condition_number\n";

    CSVWriter csv_writer((*output_file), 5);// Set precision to 5
    std::vector<std::string> headers = {"iteration", "actual_x", "actual_y", "actual_vx", "actual_vy",
                                        "estimate_x", "estimate_y", "estimate_vx", "estimate_vy",
                                        "cov_xx", "cov_yy", "cov_xy", "position_error_norm",
                                        "velocity_error_norm", "condition_number"};

    csv_writer.write_header(headers);

    // Simulation Loop
    for (int i = 0; i < 40000; ++i) {
        // Update current time
        current_time += delta_time;

        // Propagate state
        true_state = state_transition(true_state, delta_time);

        // Simulate measurement
        //        auto measurement = measurement_model(true_state);

        // Perform UKF prediction and update every 400 steps
        ukf.predict(delta_time);

        //        if (i % 400 == 0) {
        //            ukf.update(measurement, measurement_covariance(measurement));
        //        }


        // Get the state and covariance
        Eigen::Vector<SCALAR, 4> pred_state = ukf.get_state();
        Eigen::Matrix<SCALAR, 4, 4> pred_state_covariance = ukf.get_state_covariance();

        // Compute the condition number of the covariance
        Eigen::JacobiSVD<Eigen::Matrix<SCALAR, 4, 4>> svd(pred_state_covariance);
        SCALAR condition_number = svd.singularValues()(0) / svd.singularValues()(2);

        // Compute position error norm
        SCALAR position_error_norm = (true_state.head(2) - pred_state.head(2)).norm();

        // Compute velocity error norm
        SCALAR velocity_error_norm = (true_state.tail(2) - pred_state.tail(2)).norm();

        // Write simulation data to the output file
        std::vector<SCALAR> row = {
                current_time,
                true_state[0],
                true_state[1],
                true_state[2],
                true_state[3],
                pred_state[0],
                pred_state[1],
                pred_state[2],
                pred_state[3],
                pred_state_covariance(0, 0),
                pred_state_covariance(1, 1),
                pred_state_covariance(0, 1),
                position_error_norm,
                velocity_error_norm,
                condition_number};

        csv_writer.write_row(row);
    }

    // Close output file
    (*output_file).close();
    return 0;
}
