//#include <state.hpp>
//#include <gtest/gtest.h>
//
//TEST(StateTest, TranslationalStateSplitCombine) {
//    TranslationalState<3> state;
//    state.pos = Eigen::Vector3d(1, 2, 3);
//    state.vel = Eigen::Vector3d(4, 5, 6);
//
//    auto split_state = TranslationalState<3>::split(state);
//    ASSERT_EQ(std::get<0>(split_state), state.pos);
//    ASSERT_EQ(std::get<1>(split_state), state.vel);
//
//    TranslationalState<3> state2;
//    TranslationalState<3>::combine(state2, state.pos, state.vel);
//    ASSERT_EQ(state2.pos, state.pos);
//    ASSERT_EQ(state2.vel, state.vel);
//}
//
//TEST(StateTest, RotationalStateSplitCombine) {
//    RotationalState<Eigen::Quaterniond> state;
//    state.rot = Eigen::Quaterniond(1, 0, 0, 0);
//    state.ang_vel = Eigen::Vector3d(1, 2, 3);
//
//    auto split_state = RotationalState<Eigen::Quaterniond>::split(state);
//    ASSERT_EQ(std::get<0>(split_state), state.rot);
//    ASSERT_EQ(std::get<1>(split_state), state.ang_vel);
//
//    RotationalState<Eigen::Quaterniond> state2;
//    RotationalState<Eigen::Quaterniond>::combine(state2, state.rot, state.ang_vel);
//    ASSERT_EQ(state2.rot.coeffs(), state.rot.coeffs());
//    ASSERT_EQ(state2.ang_vel, state.ang_vel);
//}
//
//TEST(StateTest, RigidBodyStateSplitCombine) {
//    RigidBodyState state;
//    state.center_of_gravity = Eigen::Vector3d(1, 2, 3);
//    state.mass = 4;
//    state.inertia = Eigen::Matrix3d::Identity();
//    state.mass_distribution = Eigen::Matrix3d::Identity();
//
//    auto split_state = RigidBodyState::split(state);
//    ASSERT_EQ(std::get<0>(split_state), state.center_of_gravity);
//    ASSERT_EQ(std::get<1>(split_state), state.mass);
//    ASSERT_EQ(std::get<2>(split_state), state.inertia);
//    ASSERT_EQ(std::get<3>(split_state), state.mass_distribution);
//
//    RigidBodyState state2;
//    RigidBodyState::combine(state2, state.center_of_gravity, state.mass, state.inertia, state.mass_distribution);
//    ASSERT_EQ(state2.center_of_gravity, state.center_of_gravity);
//    ASSERT_EQ(state2.mass, state.mass);
//    ASSERT_EQ(state2.inertia, state.inertia);
//    ASSERT_EQ(state2.mass_distribution, state.mass_distribution);
//}
//
//TEST(StateTest, StateHandlerSplitCombine) {
//    EntityState<TranslationalState<3>, RotationalState<Eigen::Quaterniond>, RigidBodyState> state;
//    state.pos = Eigen::Vector3d(1, 2, 3);
//    state.vel = Eigen::Vector3d(4, 5, 6);
//    state.rot = Eigen::Quaterniond(1, 0, 0, 0);
//    state.ang_vel = Eigen::Vector3d(1, 2, 3);
//    state.center_of_gravity = Eigen::Vector3d(1, 2, 3);
//    state.mass = 4;
//    state.inertia = Eigen::Matrix3d::Identity();
//    state.mass_distribution = Eigen::Matrix3d::Identity();
//
//    auto split_state = StateHandler<TranslationalState<3>, RotationalState<Eigen::Quaterniond>, RigidBodyState>::split(state);
//    ASSERT_EQ(std::get<0>(split_state), state.pos);
//    ASSERT_EQ(std::get<1>(split_state), state.vel);
//    ASSERT_EQ(std::get<2>(split_state), state.rot);
//    ASSERT_EQ(std::get<3>(split_state), state.ang_vel);
//    ASSERT_EQ(std::get<4>(split_state), state.center_of_gravity);
//    ASSERT_EQ(std::get<5>(split_state), state.mass);
//    ASSERT_EQ(std::get<6>(split_state), state.inertia);
//    ASSERT_EQ(std::get<7>(split_state), state.mass_distribution);
//
//    EntityState<TranslationalState<3>, RotationalState<Eigen::Quaterniond>, RigidBodyState> state2;
//    StateHandler<TranslationalState<3>, RotationalState<Eigen::Quaterniond>, RigidBodyState>::combine(state2, std::make_tuple(state.pos, state.vel), std::make_tuple(state.rot, state.ang_vel), std::make_tuple(state.center_of_gravity, state.mass, state.inertia, state.mass_distribution));    ASSERT_EQ(state2.vel, state.vel);
//    ASSERT_EQ(state2.rot.coeffs(), state.rot.coeffs());
//    ASSERT_EQ(state2.ang_vel, state.ang_vel);
//    ASSERT_EQ(state2.center_of_gravity, state.center_of_gravity);
//    ASSERT_EQ(state2.mass, state.mass);
//    ASSERT_EQ(state2.inertia, state.inertia);
//    ASSERT_EQ(state2.mass_distribution, state.mass_distribution);
//}
//
//int main(int argc, char **argv) {
//    testing::InitGoogleTest(&argc, argv);
//    return RUN_ALL_TESTS();
//}
