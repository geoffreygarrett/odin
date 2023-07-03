#ifndef CRITERIA_HPP
#define CRITERIA_HPP
#pragma once

#include <type_traits>
#include <deque>
#include <utility>

namespace criteria {

    struct identity {
        template<typename T>
        const T &operator()(const T &value) const {
            return value;
        }
    };

    struct first {
        template<typename T, typename U>
        const T &operator()(const std::pair<T, U> &pair) const {
            return pair.first;
        }
    };

    struct second {
        template<typename T, typename U>
        const U &operator()(const std::pair<T, U> &pair) const {
            return pair.second;
        }
    };

    template<typename Extractor = identity>
    struct min {
        Extractor extractor;

        template<typename T>
        T operator()(const T &lhs, const T &rhs) const {
            return extractor(lhs) < extractor(rhs) ? lhs : rhs;
        }
    };

    template<typename Extractor = identity>
    struct max {
        Extractor extractor;

        template<typename T>
        T operator()(const T &lhs, const T &rhs) const {
            return extractor(lhs) > extractor(rhs) ? lhs : rhs;
        }
    };

#if __cplusplus >= 202002L // C++20

    template<typename T>
    concept Arithmetic = requires(T a, T b) {
        { T{} } -> std::same_as<T>;    // Default constructible
        { a += b } -> std::same_as<T&>; // Supports addition assignment
        { a / b } -> std::same_as<T>;   // Supports division
    };

    template<Arithmetic T>
    struct mean {
        T sum = T{};
        size_t count = 0;

        T operator()(const T &rhs) {
            sum += rhs;
            count++;
            return sum / static_cast<T>(count);
        }

        void reset() {
            sum = T{};
            count = 0;
        }
    };

#else // C++17 or older

    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    struct mean {
        T sum = T{};
        size_t count = 0;

        T operator()(const T &rhs) {
            sum += rhs;
            count++;
            return sum / static_cast<T>(count);
        }

        void reset() {
            sum = T{};
            count = 0;
        }
    };

#endif

    template<typename T, size_t N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    struct window_mean {
        T sum = T{};
        std::deque<T> window;

        T operator()(const T &rhs) {
            sum += rhs;
            window.push_back(rhs);

            if (window.size() > N) {
                sum -= window.front();
                window.pop_front();
            }

            return sum / static_cast<T>(window.size());
        }

        void reset() {
            sum = T{};
            window.clear();
        }
    };
}

#endif // CRITERIA_HPP
