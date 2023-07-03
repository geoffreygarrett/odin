//#include "measurement_model.hpp"
//#include <Eigen/Dense>
//#include <utility>
//
//struct TwoWayRangeModel : MeasurementModel<TwoWayRangeModel, Eigen::Vector3d, Eigen::MatrixXd> {
//    Eigen::Vector3d ground_station_position;
//    double t_send, t_receive;
//    static constexpr double speed_of_light = 299792.458; // speed of light in km/s
//
//    TwoWayRangeModel(Eigen::Vector3d ground_station_position, double t_send, double t_receive)
//            : ground_station_position(std::move(ground_station_position)), t_send(t_send), t_receive(t_receive) {}
//
////    [[nodiscard]] double compute_measurement_impl(const Eigen::Vector3d &spacecraft_position) const {
//        double travel_time = t_receive - t_send; // two-way travel time
//        double range_measured = travel_time * speed_of_light;
//        return range_measured;
//    }
//
//    [[nodiscard]] Eigen::Matrix<double, 1, 3>
//    compute_analytical_jacobian(const Eigen::Vector3d &spacecraft_position) const {
//        Eigen::Matrix<double, 1, 3> jacobian;
//        Eigen::Vector3d diff = spacecraft_position - ground_station_position;
//        double distance = diff.norm();
//
//        // The derivative of the range with respect to the position is the unit vector pointing from the ground station to the spacecraft.
//        jacobian << diff(0) / distance, diff(1) / distance, diff(2) / distance;
//
//        return jacobian;
//    }
//
//};
