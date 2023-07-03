#include <include/gtest/gtest.h>
#include <limits>
#include <odin/core/tracker/tracker.hpp>

// Tests for the `Tracker` class
TEST(TrackerTests, TracksMinimumValueCorrectlyWithInitialization) {
    Tracker<int, criteria::min<>> tracker(100);
    tracker.update(10);
    tracker.update(20);
    tracker.update(5);
    EXPECT_EQ(tracker.get_tracked_data().value(), 5);
    tracker.update(3);
    EXPECT_EQ(tracker.get_tracked_data().value(), 3);
    tracker.update(std::numeric_limits<int>::max());
    EXPECT_EQ(tracker.get_tracked_data().value(), 3);
}

TEST(TrackerTests, TracksMinimumValueCorrectlyWithoutInitialization) {
    Tracker<int, criteria::min<>> tracker;
    tracker.update(10);
    tracker.update(20);
    tracker.update(5);
    EXPECT_EQ(tracker.get_tracked_data().value(), 5);
    tracker.update(3);
    EXPECT_EQ(tracker.get_tracked_data().value(), 3);
    tracker.update(std::numeric_limits<int>::max());
    EXPECT_EQ(tracker.get_tracked_data().value(), 3);
}

TEST(TrackerTests, TracksMaximumValueCorrectlyWithInitialization) {
    Tracker<int, criteria::max<>> tracker(std::numeric_limits<int>::min());
    tracker.update(10);
    tracker.update(20);
    tracker.update(5);
    EXPECT_EQ(tracker.get_tracked_data().value(), 20);
    tracker.update(30);
    EXPECT_EQ(tracker.get_tracked_data().value(), 30);
    tracker.update(std::numeric_limits<int>::min());
    EXPECT_EQ(tracker.get_tracked_data().value(), 30);
}

TEST(TrackerTests, TracksMaximumValueCorrectlyWithoutInitialization) {
    Tracker<int, criteria::max<>> tracker;
    tracker.update(10);
    tracker.update(20);
    tracker.update(5);
    EXPECT_EQ(tracker.get_tracked_data().value(), 20);
    tracker.update(30);
    EXPECT_EQ(tracker.get_tracked_data().value(), 30);
    tracker.update(std::numeric_limits<int>::min());
    EXPECT_EQ(tracker.get_tracked_data().value(), 30);
}

TEST(TrackerTests, GetTrackedDataBeforeInitialization) {
    Tracker<int, criteria::min<>> tracker;
    EXPECT_FALSE(tracker.get_tracked_data().has_value());
    tracker.update(10);
    EXPECT_EQ(tracker.get_tracked_data().value(), 10);
}

TEST(TrackerTests, ResetsCorrectly) {
    Tracker<int, criteria::min<>> tracker(100);
    tracker.update(10);
    tracker.update(20);
    tracker.update(5);
    tracker.clear();
    EXPECT_FALSE(tracker.get_tracked_data().has_value());
    tracker.update(10);
    EXPECT_TRUE(tracker.get_tracked_data().has_value());
    EXPECT_EQ(tracker.get_tracked_data().value(), 10);
}


TEST(CriteriaTests, IdentityReturnsInputUnchanged) {
    criteria::identity id;
    EXPECT_EQ(id(42), 42);
    EXPECT_EQ(id("test"), "test");
}

TEST(CriteriaTests, FirstReturnsFirstOfPair) {
    criteria::first first;
    EXPECT_EQ(first(std::make_pair(42, "test")), 42);
}

TEST(CriteriaTests, SecondReturnsSecondOfPair) {
    criteria::second second;
    EXPECT_EQ(second(std::make_pair(42, "test")), "test");
}

TEST(CriteriaTests, MinReturnsMinimumValue) {
    criteria::min<> min{};
    EXPECT_EQ(min(10, 20), 10);
    EXPECT_EQ(min(-10, -20), -20);
    EXPECT_EQ(min(0, 0), 0);
    EXPECT_EQ(min(std::numeric_limits<int>::max(), 10), 10);
}

TEST(CriteriaTests, MaxReturnsMaximumValue) {
    criteria::max<> max{};
    EXPECT_EQ(max(10, 20), 20);
    EXPECT_EQ(max(-10, -20), -10);
    EXPECT_EQ(max(0, 0), 0);
    EXPECT_EQ(max(std::numeric_limits<int>::min(), 10), 10);
}
