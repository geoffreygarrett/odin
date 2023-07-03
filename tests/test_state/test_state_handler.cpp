#include <include/gtest/gtest.h>

#include <state/state.hpp>
#include <state/entity.hpp>

// Static entities
DEFINE_ENTITY(EARTH);
DEFINE_ENTITY(DELFI);

// Define StateHandler
using GlobalStates2D = StateHandler<Position<EARTH, 2>, Velocity<EARTH, 2>, Position<DELFI, 2>, Velocity<DELFI, 2>>;
using GlobalStates = StateHandler<Position<EARTH>, Velocity<EARTH>, Position<DELFI>, Velocity<DELFI>>;

TEST(StateHandlerTest, GetAndModifyStateTest) {
    // Create entity instances with unique IDs
    Position<EARTH> earth_pos{10.0, 20.0, 30.0};
    Position<DELFI> delfi_pos{5.0, -2.0, 7.0};
    Velocity<EARTH> earth_vel{1.0, 2.0, 3.0};
    Velocity<DELFI> delfi_vel{-1.0, 0.5, -0.3};
    GlobalStates::state_tuple_t global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);

    // Access and modify state using StateHandler
    GlobalStates::get<Position<EARTH >>(global_states).value.x() += 10.0;
    ASSERT_EQ(GlobalStates::get<Position<EARTH >>(global_states).value.x(), 20.0);
}

TEST(StateHandlerTest, SplitAndModifyStateTest) {
    // Create entity instances with unique IDs
    Position<EARTH> earth_pos{10.0, 20.0, 30.0};
    Position<DELFI> delfi_pos{5.0, -2.0, 7.0};
    Velocity<EARTH> earth_vel{1.0, 2.0, 3.0};
    Velocity<DELFI> delfi_vel{-1.0, 0.5, -0.3};
    GlobalStates::state_tuple_t global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);

    // Access via get
    auto earth_pos_ptr = std::get<Position<EARTH> *>(GlobalStates::split(global_states));

    // Modify through pointer
    (*earth_pos_ptr).value[0] += 10.0;

    ASSERT_EQ((*earth_pos_ptr).value[0], 20.0);
    ASSERT_EQ(GlobalStates::get<Position<EARTH >>(global_states).value[0], 20.0);
}

TEST(StateHandlerTest, RecombineStateTest) {
    // Create entity instances with unique IDs
    Position<EARTH, 2> earth_pos{10.0, 20.0};
    Position<DELFI, 2> delfi_pos{5.0, -2.0};
    Velocity<EARTH, 2> earth_vel{1.0, 2.0};
    Velocity<DELFI, 2> delfi_vel{-1.0, 0.5};
    GlobalStates2D::state_tuple_t global_states = std::make_tuple(earth_pos, earth_vel, delfi_pos, delfi_vel);

    // Define entirely new position states and recombine them
    Position<EARTH, 2> new_earth_pos{30.0, 40.0};
    Position<DELFI, 2> new_delfi_pos{15.0, -7.0};
    GlobalStates2D::recombine(global_states, std::make_tuple(new_earth_pos, new_delfi_pos));
    using EARTH_POSITION = Position<EARTH, 2>;
    using DELFI_POSITION = Position<DELFI, 2>;
    ASSERT_EQ(GlobalStates2D::get<EARTH_POSITION>(global_states).value.x(), 30.0);
    ASSERT_EQ(GlobalStates2D::get<EARTH_POSITION>(global_states).value.y(), 40.0);
    ASSERT_EQ(GlobalStates2D::get<DELFI_POSITION>(global_states).value.x(), 15.0);
    ASSERT_EQ(GlobalStates2D::get<DELFI_POSITION>(global_states).value.y(), -7.0);
}
