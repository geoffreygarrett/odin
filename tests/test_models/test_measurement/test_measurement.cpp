//////
////// Created by geoffrey on 6/4/23.
//////
////
////#include <measurement/measurement.hpp>
////#include <state/state_handler.hpp>
////#include <state/state.hpp>
////#include <gtest/gtest.h>
////
////
////// The bodies are identified by their types.
////struct EARTH {
////};
////struct DELFI {
////};
////
////using GlobalStates = StateHandler<Position<EARTH>, Velocity<EARTH>, Position<DELFI>, Velocity<DELFI>>;
////
////
//template<typename E1, typename E2, typename GlobalStateType, typename MeasurementType, typename JacobianType>
//struct SimpleDistanceMeasurement
//        : MeasurementModel<SimpleDistanceMeasurement<E1, E2, GlobalStateType, MeasurementType, JacobianType>, GlobalStateType, MeasurementType, JacobianType> {
//
//    [[nodiscard]] constexpr MeasurementType compute_measurement_impl(const GlobalStateType &global_state) const {
//        const auto &state1 = GlobalStates::get<Position<E1>>(global_state).value;
//        const auto &state2 = GlobalStates::get<Position<E2>>(global_state).value;
//
//        auto difference = state1 - state2;
//        MeasurementType measurement;
//        measurement(0, 0) = difference.norm();
//        return measurement;
//    }
//};
////
////// Define the types as in the previous example
////using StateType = Eigen::Vector3d;
////using MeasurementType = Eigen::Matrix<double, 1, 1>;
////using JacobianType = Eigen::Matrix<double, 1, StateType::RowsAtCompileTime>;
////
////// Alias for brevity
////using SimpleDist = SimpleDistanceMeasurement<EARTH, DELFI, GlobalStates::state_tuple, MeasurementType, JacobianType>;
////
////class SimpleDistanceMeasurementTest : public ::testing::Test {
////protected:
////    SimpleDist model;
////};
////
////TEST_F(SimpleDistanceMeasurementTest, ZeroDistance) {
////    Position<EARTH> earth_pos{0.0, 0.0, 0.0};
////    Position<DELFI> delfi_pos{0.0, 0.0, 0.0};
////    Velocity<EARTH> earth_vel{0.0, 0.0, 0.0};
////    Velocity<DELFI> delfi_vel{0.0, 0.0, 0.0};
////
////
////    GlobalStates::state_tuple global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);
////
////    MeasurementType expected(0.0);
////    EXPECT_EQ(model.compute_measurement(global_states), expected);
////}
////
////TEST_F(SimpleDistanceMeasurementTest, UnitDistance) {
////    Position<EARTH> earth_pos{1.0, 0.0, 0.0};
////    Position<DELFI> delfi_pos{0.0, 0.0, 0.0};
////    Velocity<EARTH> earth_vel{0.0, 0.0, 0.0};
////    Velocity<DELFI> delfi_vel{0.0, 0.0, 0.0};
////    GlobalStates::state_tuple global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);
////
////    MeasurementType expected(1.0);
////    EXPECT_EQ(model.compute_measurement(global_states), expected);
////}
////
////TEST_F(SimpleDistanceMeasurementTest, GeneralCase) {
////    Position<EARTH> earth_pos{1.0, 2.0, 3.0};
////    Position<DELFI> delfi_pos{4.0, 5.0, 6.0};
////    Velocity<EARTH> earth_vel{0.0, 0.0, 0.0};
////    Velocity<DELFI> delfi_vel{0.0, 0.0, 0.0};
////    GlobalStates::state_tuple global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);
////
////    MeasurementType expected(std::sqrt(27.0));
////    EXPECT_EQ(model.compute_measurement(global_states), expected);
////}
////
////TEST_F(SimpleDistanceMeasurementTest, ZeroDistance) {
////    // Setting up the state
////    Position<EARTH> earth_pos{0.0, 0.0, 0.0};
////    Position<DELFI> delfi_pos{0.0, 0.0, 0.0};
////    Velocity<EARTH> earth_vel{0.0, 0.0, 0.0};
////    Velocity<DELFI> delfi_vel{0.0, 0.0, 0.0};
////    GlobalStates::state_tuple global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);
////
////    // Compute and print the Jacobian
////    auto jacobian = model.compute_jacobian(global_states);
////    std::cout << "Jacobian at zero distance:\n" << jacobian << "\n\n";
////}
////
////TEST_F(SimpleDistanceMeasurementTest, UnitDistance) {
////    // Setting up the state
////    Position<EARTH> earth_pos{1.0, 0.0, 0.0};
////    Position<DELFI> delfi_pos{0.0, 0.0, 0.0};
////    Velocity<EARTH> earth_vel{0.0, 0.0, 0.0};
////    Velocity<DELFI> delfi_vel{0.0, 0.0, 0.0};
////    GlobalStates::state_tuple global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);
////
////    // Compute and print the Jacobian
////    auto jacobian = model.compute_jacobian(global_states);
////    std::cout << "Jacobian at unit distance:\n" << jacobian << "\n\n";
////}
////
////TEST_F(SimpleDistanceMeasurementTest, GeneralCase) {
////    // Setting up the state
////    Position<EARTH> earth_pos{1.0, 2.0, 3.0};
////    Position<DELFI> delfi_pos{4.0, 5.0, 6.0};
////    Velocity<EARTH> earth_vel{0.0, 0.0, 0.0};
////    Velocity<DELFI> delfi_vel{0.0, 0.0, 0.0};
////    GlobalStates::state_tuple global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);
////
////    // Compute and print the Jacobian
////    auto jacobian = model.compute_jacobian(global_states);
////    std::cout << "Jacobian in the general case:\n" << jacobian << "\n\n";
////}
