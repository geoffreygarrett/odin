#ifndef CACHE_POLICY_HPP
#define CACHE_POLICY_HPP

#include <tuple>
#include <functional>


template<typename Derived, typename... Args>
class CachePolicy {
public:
    void store(const Args &... args) {
        static_cast<Derived *>(this)->store_impl(args...);
    }

    [[nodiscard]] bool compare(const Args &... args) const {
        return static_cast<const Derived *>(this)->compare_impl(args...);
    }

    [[nodiscard]] bool compare_and_store_if_invalid(const Args &... args) const {
        return static_cast< Derived *>(this)->compare_and_store_if_invalid_impl(args...);
    }

};

template<typename... Args>
class NoCachePolicy : public CachePolicy<NoCachePolicy<Args...>> {
public:
    void store_impl(Args... args) {}

    bool compare_impl(Args... args) const {
        return false;
    }

    bool compare_and_store_if_invalid_impl(Args... args) const {
        return false;
    }
};

template<>
class NoCachePolicy<> : public CachePolicy<NoCachePolicy<>> {

    template<typename... Args>
    void store_impl(Args... args) {}

    template<typename... Args>
    bool compare_impl(Args... args) const {
        return false;
    }

};

template<typename... Args>
class DirectCachePolicy : public CachePolicy<DirectCachePolicy<Args...>> {
public:
    void store_impl(Args... args) {
        store_tuple = std::make_tuple(args...);
    }

    bool compare_impl(Args... args) const {
        return store_tuple == std::make_tuple(args...);
    }

    bool compare_and_store_if_invalid_impl(Args... args) {
        if (compare_impl(args...)) {
            return true;
        } else {
            store_impl(args...);
            return false;
        }
    }

private:
    std::tuple<Args...> store_tuple;
};


template<typename... Args>
class HashCachePolicy : public CachePolicy<HashCachePolicy<Args...>> {
public:
    auto hash_values(Args... args) const {
        return std::make_tuple(hash_args(args)...);
    }

    void store_impl(Args... args) {
        hash_values_ = hash_values(args...);
    }

    bool compare_impl(Args... args) const {
        return hash_values_ == hash_values(args...);
    }

    bool compare_and_store_if_invalid_impl(Args... args) {
        auto hash_values = hash_args(args...);
        if (hash_values_ == hash_values) {
            return true;
        } else {
            hash_values_ = hash_values;
            return false;
        }
    }

private:
    std::array<std::size_t, sizeof...(Args)> hash_values_;

    // Hash function for an individual argument
    template<typename T>
    std::size_t hash_arg(const T &arg) const {
        std::size_t seed = std::hash<T>{}(arg);
        return seed ^ (0x9e3779b9 + (seed << 6) + (seed >> 2));
    }

    // Hash function for a collection (like std::vector)
    // TODO: Look into parallelising this and benchmarking the result.
    template<typename Collection>
    std::size_t hash_collection(const Collection &coll) const {
        std::size_t seed = coll.size();
        for (const auto &elem: coll) {
            seed ^= hash_arg(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    // Hash function for a collection (like std::vector)
    template<typename Collection>
    std::size_t hash_args(const Collection &coll) const {
        // Check if T is a scalar
        if constexpr (!std::is_scalar<typename Collection::value_type>::value) {
            return hash_collection(coll);
        } else {
            return hash_arg(coll);
        }
    }
};

#endif //CACHE_POLICY_HPP