
#ifndef ODIN_CONCEPTS_HPP
#define ODIN_CONCEPTS_HPP

#include <concepts>
#include <optional>
#include <ostream>
#include <string>
#include <vector>


template<typename T>
concept Streamable = requires(T a, std::ostream &os) {
    { os << a } -> std::same_as<std::ostream &>;
};

template<typename T>
concept Comparable = requires(T a, T b) {
    { a == b } -> std::same_as<bool>;
};

template<typename T>
concept Stringable = requires(T a) {
    { a.to_string() } -> std::same_as<std::string>;
};

template<typename T, typename Archive>
concept Serializable = requires(T a, Archive &ar) {
    { serialize(ar, a) } -> std::same_as<void>;
} || requires(T a, Archive &ar) {
    { save(ar, a) } -> std::same_as<void>;
    { load(ar, a) } -> std::same_as<void>;
} || requires(T a, Archive &ar) {
    { ar(a) } -> std::same_as<void>;
};

#endif//ODIN_CONCEPTS_HPP