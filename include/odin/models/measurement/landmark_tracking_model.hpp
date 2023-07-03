//#ifndef LANDMARK_TRACKING_MODEL_HPP
//#define LANDMARK_TRACKING_MODEL_HPP
//
//#include <Eigen/Dense>
//#include <utility>
//#include "measurement_model.hpp"
//
///**
// * @brief Measurement model function for landmark tracking with a pin-hole camera.
// *
// * This model maps 3D positions to 2D image points for simulation and navigational purposes. The 3D position in the camera-fixed frame (x^c, y^c, z^c) is mapped to pixel points on the image plane (u,v).
// * The full procedure for landmark tracking can be found in the following reference:
// * Shuang Liu, Guowei Cai, Bingheng Li, "Landmark Tracking Algorithm Based on EKF and UKF", 2008, https://doi.org/10.1109/ROBIO.2009.4913042.
// *
// * @param position_3d Eigen::Vector3d input representing the position in the camera-fixed frame.
// * @param f double input representing the focal length of the pin-hole camera.
// *
// * @return Eigen::Vector2d output representing the projected position on the image plane.
// */
//
//struct PinholeCameraModel : public MeasurementModel<PinholeCameraModel, Eigen::Vector3d, Eigen::Vector2d> {
//    double focal_length;
//    Eigen::Vector2d optical_center;
//
//    PinholeCameraModel(double focal_length, Eigen::Vector2d optical_center)
//            : focal_length(focal_length), optical_center(std::move(optical_center)) {}
//
//    [[nodiscard]] Eigen::Vector2d compute_measurement_impl(const Eigen::Vector3d &state) const {
//        const double u = focal_length * ((state[0] - optical_center[0]) / state[2]);
//        const double v = focal_length * ((state[1] - optical_center[1]) / state[2]);
//        return Eigen::Vector2d{u, v};
//    }
//};
//
//
//#endif // LANDMARK_TRACKING_MODEL_HPP