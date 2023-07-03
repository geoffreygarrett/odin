#ifndef BASIC_SERIALIZATION_POLICY_H
#define BASIC_SERIALIZATION_POLICY_H

#include <sstream>
#include <concepts>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/string.hpp>
#include <type_traits>
#include "serialization_policy.hpp"

template<typename T>
concept SequentialContainer = requires(T a) {
    { a.begin() } -> std::same_as<typename T::iterator>;
    { a.end() } -> std::same_as<typename T::iterator>;
    { a.size() } -> std::same_as<typename T::size_type>;
    { a.empty() } -> std::same_as<bool>;
};

template<typename T>
concept AssociativeContainer = requires(T a) {
    { a.begin() } -> std::same_as<typename T::iterator>;
    { a.end() } -> std::same_as<typename T::iterator>;
    { a.empty() } -> std::same_as<bool>;
    { a.find(typename T::key_type{}) } -> std::same_as<typename T::iterator>;
};

#endif // BASIC_SERIALIZATION_POLICY_H
