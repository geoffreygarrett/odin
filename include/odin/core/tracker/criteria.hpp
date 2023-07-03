#ifndef CRITERIA_HPP
#define CRITERIA_HPP
#pragma once

#include <type_traits>
#include <deque>
#include <utility>

/**
 * @file criteria.hpp
 * @brief Contains a set of classes and structures which can be used to extract, compare, and compute data.
 */

namespace criteria {

    /**
     * @brief Extracts a value without any modifications.
     */
    struct identity {
        /**
         * @brief Returns the same value.
         * @tparam T Any type.
         * @param value A value of type T.
         * @return The same value.
         */
        template<typename T>
        const T& operator()(const T& value) const {
            return value;
        }
    };

    /**
     * @brief Extracts the first value of a pair.
     */
    struct first {
        /**
         * @brief Returns the first value of a pair.
         * @tparam T Any type.
         * @tparam U Any type.
         * @param pair A pair of two values.
         * @return The first value of the pair.
         */
        template<typename T, typename U>
        const T& operator()(const std::pair<T, U>& pair) const {
            return pair.first;
        }
    };

    /**
     * @brief Extracts the second value of a pair.
     */
    struct second {
        /**
         * @brief Returns the second value of a pair.
         * @tparam T Any type.
         * @tparam U Any type.
         * @param pair A pair of two values.
         * @return The second value of the pair.
         */
        template<typename T, typename U>
        const U& operator()(const std::pair<T, U>& pair) const {
            return pair.second;
        }
    };

    /**
     * @brief Extracts a value using an extractor and returns the minimum of two values.
     * @tparam Extractor A functor that extracts a value from an object.
     */
    template<typename Extractor = identity>
    struct min {
        Extractor extractor;

        /**
         * @brief Returns the minimum of two values.
         * @tparam T Any type.
         * @param lhs A value of type T.
         * @param rhs A value of type T.
         * @return The minimum value.
         */
        template<typename T>
        T operator()(const T& lhs, const T& rhs) const {
            return extractor(lhs) < extractor(rhs) ? lhs : rhs;
        }
    };

    /**
     * @brief Extracts a value using an extractor and returns the maximum of two values.
     * @tparam Extractor A functor that extracts a value from an object.
     */
    template<typename Extractor = identity>
    struct max {
        Extractor extractor;

        /**
         * @brief Returns the maximum of two values.
         * @tparam T Any type.
         * @param lhs A value of type T.
         * @param rhs A value of type T.
         * @return The maximum value.
         */
        template<typename T>
        T operator()(const T& lhs, const T& rhs) const {
            return extractor(lhs) > extractor(rhs) ? lhs : rhs;
        }
    };

//#if __cplusplus >= 202002L // C++20
//
//    // Arithmetic concept for C++20 and above
//    template<typename T>
//    concept Arithmetic = requires(T a, T b) {
//        { T{} } -> std::same_as<T>;    // Default constructible
//        { a += b } -> std::same_as<T&>; // Supports addition assignment
//        { a / b } -> std::same_as<T>;   // Supports division
//    };
//
//    #define MEAN_TYPE Arithmetic
//
//#else // C++17 or older
//
//#define MEAN_TYPE typename = std::enable_if_t<std::is_arithmetic_v<T>>
//
//#endif
//
///**
// * @brief Computes the running mean of a sequence of values.
// * @tparam T Any arithmetic type.
// */
//    template<typename T, MEAN_TYPE>
//    struct mean {
//        T sum = T{};  ///< The sum of all values added so far.
//        size_t count = 0; ///< The count of all values added so far.
//
//        /**
//         * @brief Adds a value and returns the current mean.
//         * @param lhs The current mean.
//         * @param rhs A value of type T.
//         * @return The updated mean.
//         */
//        T operator()(const std::optional<T>& /*lhs*/, const T& rhs) {
//            sum += rhs;
//            count++;
//            return sum / static_cast<T>(count);
//        }
//
//        /**
//         * @brief Resets the state.
//         */
//        void reset() {
//            sum = T{};
//            count = 0;
//        }
//    };
/////
////    template <typename T>
////    mean(T) -> mean<T>;
//
//#undef MEAN_TYPE
//
//    /**
//     * @brief Computes the running mean of a window of values.
//     * @tparam T Any arithmetic type.
//     * @tparam N The size of the window.
//     */
//    template<typename T, size_t N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
//    struct window_mean {
//        T sum = T{};  ///< The sum of the values in the current window.
//        std::deque<T> window; ///< The current window of values.
//
//        /**
//         * @brief Adds a value and returns the current mean of the window.
//         * @param lhs The current mean of the window.
//         * @param rhs A value of type T.
//         * @return The updated mean of the window.
//         */
//        T operator()(const T& /*lhs*/, const T& rhs) {
//            sum += rhs;
//            window.push_back(rhs);
//
//            if (window.size() > N) {
//                sum -= window.front();
//                window.pop_front();
//            }
//
//            return sum / static_cast<T>(window.size());
//        }
//
//        /**
//         * @brief Resets the state.
//         */
//        void reset() {
//            sum = T{};
//            window.clear();
//        }
//    };
}

#endif // CRITERIA_HPP
