#ifndef TRACKER_HPP
#define TRACKER_HPP

#include <optional>
#include "criteria.hpp"

/**
 * @brief The Tracker class tracks data using a given criteria.
 * 
 * This class accepts data and applies a given criteria to determine
 * the 'tracked' value. This can be used to track values such as 
 * maximum, minimum, mean, etc.
 * 
 * @tparam T The type of data to be tracked.
 * @tparam Criteria The criteria to be used to track the data. 
 *                  The criteria class should have an operator() 
 *                  which takes two std::optional<T> and returns T, 
 *                  and a reset method to reset its state.
 */
template <typename T, typename Criteria>
class Tracker {
public:
    /**
     * @brief Default constructor.
     */
    Tracker() = default;

    /**
     * @brief Construct a new Tracker object with initial data and a criterion.
     *
     * @param initial_data The initial data.
     * @param criterion The criterion to be used to track the data.
     */
    Tracker(const T& initial_data, const Criteria& criterion)
            : tracked_data_(initial_data), criteria_(criterion) {}

    /**
     * @brief Construct a new Tracker object with initial data and a criterion.
     *
     * @param initial_data The initial data.
     * @param criterion The criterion to be used to track the data.
     */
    explicit Tracker(const T& initial_data)
            : tracked_data_(initial_data) {}

    /**
    * @brief Update the tracker with new data.
    *
    * If this is the first update (i.e., if the tracked data is uninitialized),
    * this method simply stores the new data as the current tracked value. In
    * subsequent calls, it updates the tracked value based on the new data and
    * the established criteria.
    *
    * @param data The new data.
    */
    void update(const T& data) {
        if (tracked_data_.has_value()) {
            // Use the criteria to determine the new tracked value.
            tracked_data_ = criteria_(tracked_data_.value(), data);
        } else {
            // If tracked_data_ is uninitialized, initialize it with the new data.
            tracked_data_ = data;
        }
    }

    /**
     * @brief Get the currently tracked data.
     *
     * This method returns an std::optional containing the currently tracked data.
     * If update has not yet been called, the returned std::optional will not contain a value.
     *
     * @return An std::optional containing the currently tracked data.
     */
    const std::optional<T>& get_tracked_data() const {
        return tracked_data_;
    }

    /**
     * @brief Clear the tracked data.
     *
     * This method clears the currently tracked data and resets
     * the criteria.
     */
    void clear() {
        tracked_data_.reset();
    }

private:
    std::optional<T> tracked_data_; ///< The currently tracked data.
    Criteria criteria_; ///< The criteria to be used to track the data.
};

/**
 * @brief Template deduction guide for the Tracker class.
 *
 * This deduction guide is used for Class Template Argument Deduction (CTAD)
 * when constructing a Tracker object directly. The C++ compiler will use this
 * to infer the correct template parameters based on the arguments provided to
 * the Tracker constructor.
 *
 * For instance, when we write `Tracker(0.0, criteria::mean<double>{})`,
 * the compiler uses this guide to infer that we want a
 * `Tracker<double, criteria::mean<double>>`.
 *
 * However, this deduction guide is not necessary if we use the `make_tracker`
 * factory function, which automatically deduces the correct template parameters
 * based on the type of the initial value.
 *
 * @tparam T The type of data to track.
 */
//template <typename T>
//Tracker(T, criteria::mean<T>) -> Tracker<T, criteria::mean<T>>;

/**
 * @brief Factory function to create a Tracker object with an extractor.
 *
 * This function creates a Tracker that uses the specified criteria and extractor to track data.
 * It uses automatic template argument deduction to infer the correct types
 * for the Tracker.
 *
 * This function provides a way to create a Tracker object without having
 * to explicitly specify the template parameters.
 *
 * For instance, we can write `auto tracker = make_tracker<double, criteria::max>(0.0);`
 * instead of `Tracker<double, criteria::max<double, criteria::identity>> tracker(0.0, criteria::max<double, criteria::identity>{});`.
 *
 * @tparam T The type of data to be tracked. This is deduced automatically
 * from the type of the initial_value.
 * @tparam Criteria The type of the criteria to be used.
 * @tparam Extractor The type of the extractor to be used. Default is criteria::identity.
 *
 * @param initial_value The initial value of the data to be tracked. This will also
 * determine the type of data the Tracker object will be capable of tracking.
 *
 * @return A Tracker object configured to track data according to the specified criteria and extractor.
 */
template<typename T, template <typename, typename> typename Criteria, typename Extractor = criteria::identity>
auto make_tracker(const T& initial_value) {
    return Tracker<T, Criteria<T, Extractor>>(initial_value, Criteria<T, Extractor>{});
}

template<typename T, template <typename> typename Criteria>
auto make_tracker(const T& initial_value) {
    return Tracker<T, Criteria<T>>(initial_value, Criteria<T>{});
}

/**
 * @brief Factory function to create a Tracker object with window size and an optional extractor.
 *
 * This function creates a Tracker that uses the specified criteria and optionally extractor to track data over a certain window size.
 *
 * For instance, we can write `auto tracker = make_tracker<double, criteria::window_mean, 3>(0.0);`
 * instead of `Tracker<double, criteria::window_mean<double, 3>> tracker(0.0, criteria::window_mean<double, 3>{});`.
 *
 * @tparam T The type of data to be tracked. This is deduced automatically
 * from the type of the initial_value.
 * @tparam Criteria The type of the criteria to be used.
 * @tparam WindowSize The size of the window for window-based criteria.
 * @tparam Extractor The type of the extractor to be used. Default is criteria::identity.
 *
 * @param initial_value The initial value of the data to be tracked. This will also
 * determine the type of data the Tracker object will be capable of tracking.
 *
 * @return A Tracker object configured to track data according to the specified criteria, window size and extractor.
 */
template<typename T, template <typename, size_t, typename> typename Criteria, size_t WindowSize, typename Extractor = criteria::identity>
auto make_tracker(const T& initial_value) {
    return Tracker<T, Criteria<T, WindowSize, Extractor>>(initial_value, Criteria<T, WindowSize, Extractor>{WindowSize});
}

#endif // TRACKER_HPP

