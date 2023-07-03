#include <iostream>
#include <deque>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <functional>
#include <variant>
#include "tracker/container_policy.hpp"
#include "tracker/comparators.hpp"

#pragma once


template<typename T>
struct Helper {
    template<typename... Args>
    static T Construct(Args &&... args) {
        if constexpr (sizeof...(args) == 1)
            return T(std::forward<Args>(args)...);
        else
            return T{std::forward<Args>(args)...};
    }
};

template<typename T, typename Compare = std::less<T>>
struct TrackerIdentifier {
    Compare compare;

    explicit TrackerIdentifier(Compare comp = Compare{})
            : compare(std::move(comp)) {}
};

template<typename T, typename Container = void, typename Policy = ContainerPolicy<T, Container>, typename Compare = std::less<T>, typename Allocator = std::allocator<T>>
        class StaticTracker {
public:
    using ContainerT = typename std::conditional<!std::is_void<Container>::value, Container, std::monostate>::type;
    using AllocatorType = typename std::allocator_traits<Allocator>::template rebind_alloc<T>;
    AllocatorType allocator_;

    explicit StaticTracker(Compare comp = Compare{}) : compare_(std::move(comp)) {}

    template<typename PolicyArg>
    StaticTracker(Compare comp, PolicyArg &&arg)
            : compare_(std::move(comp)), policy_(std::forward<PolicyArg>(arg)) {}

    virtual ~StaticTracker() = default;

    template<typename... Args>
    void update(Args &&...args) {
        T data = Helper<T>::Construct(std::forward<Args>(args)...);

        if (!bestData_.has_value() || compare_(data, bestData_.value())) {
            bestData_ = data;
        }

        if constexpr (!std::is_void_v<Container>) {
            policy_.add(container_, std::forward<Args>(args)...);
        } else {
            policy_.add(std::forward<Args>(args)...);
        }
    }

    void clear() {
        if constexpr (!std::is_void_v<Container>) {
            policy_.clear(container_);
        }

        bestData_.reset();
    }

    [[nodiscard]] const std::optional<T> &get_best() const {
        return bestData_;
    }

    template<typename C = Container>
    [[nodiscard]] std::enable_if_t<!std::is_same_v<C, void>, C>
    get_container() const {
        if constexpr (std::is_array_v<C>) {
            return policy_.get_ordered_container(container_);
        } else {
            return container_;
        }
    }

    template<typename C = Container>
    std::enable_if_t<std::is_same_v<C, void>, std::monostate>
    get_container() const {
        return std::monostate{};
    }


    [[nodiscard]] size_t size_in_bytes() const {
        size_t size = sizeof(T); // Always account for the size of the bestData_ object

        if constexpr (!std::is_void_v<Container>) {
            size += policy_.size_in_bytes(container_);  // Add the size of the container
        }

        return size;
    }

    template<typename C = Container>
    [[nodiscard]] std::enable_if_t<!std::is_same_v<C, void>, C>
    move_container() {
        if constexpr (std::is_array_v<C>) {
            return policy_.get_ordered_container(std::move(container_));
        } else {
            return std::move(container_);
        }
    }

protected:
    Compare compare_;        ///< The comparator function for determining the best value.
    Policy policy_;             ///< The policy for manipulating the container.
    ContainerT container_;      ///< The container for keeping a history of all values.
    std::optional<T> bestData_; ///< The current best value.
};




template<typename... Trackers>
class TrackerManager {
public:
    template<typename T, typename... Args>
    void update(TrackerIdentifier<T> id, Args &&... args) {
        T &tracker = std::get<T>(trackers_);
        tracker.update(std::forward<Args>(args)...);
    }

    template<typename T>
    T &get_tracker(TrackerIdentifier<T> id) {
        return std::get<T>(trackers_);
    }

    void clear_all() {
        std::apply([](auto &... tracker) { (tracker.clear(), ...); }, trackers_);
    }

private:
    std::tuple<Trackers...> trackers_;
};


namespace tracker {
    template<typename T>
    using no_history = StaticTracker<T>;

    template<typename T, typename Compare>
    using no_history_custom = StaticTracker<T, void, ContainerPolicy<T, void>, Compare>;

    template<typename T>
    using vector_history = StaticTracker<T, std::vector<T>>;

    template<typename T, typename Compare>
    using vector_history_custom = StaticTracker<T, std::vector<T>, ContainerPolicy<T, std::vector<T>>, Compare>;

    template<typename T>
    using deque_history = StaticTracker<T, std::deque<T>>;

    template<typename T, typename Compare>
    using deque_history_custom = StaticTracker<T, std::deque<T>, ContainerPolicy<T, std::deque<T>>, Compare>;

    template<typename T, typename U>
    using map_history = StaticTracker<std::pair<T, U>, std::map<T, U>>;

    template<typename T, typename U, typename Compare>
    using map_history_custom = StaticTracker<std::pair<T, U>, std::map<T, U>, ContainerPolicy<std::pair<T, U>, std::map<T, U>>, Compare>;

    template<typename T, std::size_t N>
    using array_history = StaticTracker<T, std::array<T, N>, ContainerPolicy<T, std::array<T, N>>>;

    template<typename T, std::size_t N, typename Compare>
    using array_history_custom = StaticTracker<T, std::array<T, N>, ContainerPolicy<T, std::array<T, N>>, Compare>;


    template<typename T>
    TrackerIdentifier<T> make_identifier() {
        return TrackerIdentifier<T>{};
    }
}
