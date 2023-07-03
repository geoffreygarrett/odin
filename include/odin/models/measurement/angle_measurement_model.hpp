//#include "measurement_model.hpp"
//
//
//struct AngleMeasurementModel
//        : MeasurementModel<AngleMeasurementModel, Eigen::Vector3d, Eigen::VectorXd, Eigen::Vector3d> {
//    Eigen::Vector3d ground_station_position;
//
//    AngleMeasurementModel(Eigen::Vector3d ground_station_position)
//            : ground_station_position(std::move(ground_station_position)) {}
//
//    Eigen::VectorXd
//    compute_measurement_impl(const Eigen::Vector3d &spacecraft_position, const Eigen::Vector3d &params) const {
//        // Assume that the parameters are [difference signal, sum signal, speed_of_light]
//        double delta_signal = params[0];
//        double sum_signal = params[1];
//        double speed_of_light = params[2]; // speed of light
//
//        Eigen::Vector3d B = spacecraft_position - ground_station_position; // baseline vector
//
//        // tau_g for VLBI
//        double tau_g = B.dot(spacecraft_position.normalized()) / speed_of_light;
//
//        // Compute angle offset for monopulse technique
//        double angle_offset = atan2(delta_signal, sum_signal);
//
//        // Pack the two measurements into a vector
//        Eigen::VectorXd measurement(2);
//        measurement << tau_g, angle_offset;
//
//        return measurement;
//    }
//
//    Eigen::Matrix<double, 2, 3>
//    compute_analytical_jacobian(const Eigen::Vector3d &spacecraft_position, const Eigen::Vector3d &params) const {
//        // Compute the partial derivatives of tau_g and angle_offset with respect to the spacecraft position
//        Eigen::Matrix<double, 2, 3> jacobian;
//        double speed_of_light = params[2]; // speed of light
//
//        // Derivative of tau_g with respect to spacecraft_position
//        Eigen::Vector3d dtau_g = ground_station_position.normalized() / speed_of_light;
//
//        // Derivative of angle_offset with respect to spacecraft_position
//        // This depends on how angle_offset changes with respect to spacecraft_position, which is not straightforward
//        // Here we make an assumption that the angle offset doesn't change with spacecraft position, which is a simplification
//        Eigen::Vector3d dangle_offset = Eigen::Vector3d::Zero();
//
//        jacobian.row(0) = dtau_g;
//        jacobian.row(1) = dangle_offset;
//
//        return jacobian;
//    }
//};