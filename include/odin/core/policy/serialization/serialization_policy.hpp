/**
 * @file
 * @brief Defines the SerializationPolicy class for serializing and deserializing data.
 */

#ifndef SERIALIZATION_POLICY_H
#define SERIALIZATION_POLICY_H


#include <sstream>
#include <cereal/archives/binary.hpp>
#include <type_traits>
#include <list>
#include <concepts>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <Eigen/Core>
#include "eigen_serialization_policy.hpp"
#include "basic_serialization_policy.hpp"

// Here we define two traits to identify STL sequential and associative containers.
// This code is currently commented out but can be used as a substitute for C++20 Concepts in C++17 environment.

/*
template<typename T>
struct is_stl_sequential_container : std::false_type {};

template<typename T>
struct is_stl_associative_container : std::false_type {};

// Specialize the traits for the containers we want to support
template<typename... Args>
struct is_stl_sequential_container<std::vector<Args...>> : std::true_type {};

template<typename... Args>
struct is_stl_sequential_container<std::list<Args...>> : std::true_type {};

template<typename... Args>
struct is_stl_associative_container<std::map<Args...>> : std::true_type {};

template<typename... Args>
struct is_stl_associative_container<std::unordered_map<Args...>> : std::true_type {};
*/

// The following Concepts are a part of C++20 and replace the type traits defined above.
// These Concepts check whether a type is a sequential or an associative container by checking the presence of
// the member functions begin(), end(), and others.

namespace archive {

    /**
     * @brief Archive traits for JSON
     *
     * Pros:
     *   - Human-readable format.
     *   - Wide support across different languages and platforms.
     *
     * Cons:
     *   - Larger file size due to text-based format.
     *   - Slower serialization/deserialization compared to binary formats.
     */
    struct json {
        using output = cereal::JSONOutputArchive;
        using input = cereal::JSONInputArchive;
    };

    /**
     * @brief Archive traits for XML
     *
     * Pros:
     *   - Human-readable format.
     *   - Wide support across different languages and platforms.
     *
     * Cons:
     *   - Larger file size due to text-based format.
     *   - Slower serialization/deserialization compared to binary formats.
     *   - More verbose compared to JSON.
     */
    struct xml {
        using output = cereal::XMLOutputArchive;
        using input = cereal::XMLInputArchive;
    };

    /**
     * @brief Archive traits for Binary
     *
     * Pros:
     *   - Fast serialization/deserialization.
     *   - Compact file size.
     *
     * Cons:
     *   - Not human-readable.
     *   - Not portable across systems with different endianness or word size.
     */
    struct binary {
        using output = cereal::BinaryOutputArchive;
        using input = cereal::BinaryInputArchive;
    };

    /**
     * @brief Archive traits for Portable Binary
     *
     * Pros:
     *   - Fast serialization/deserialization.
     *   - Compact file size.
     *   - Portable across systems with different endianness or word size.
     *
     * Cons:
     *   - Not human-readable.
     */
    struct portable_binary {
        using output = cereal::PortableBinaryOutputArchive;
        using input = cereal::PortableBinaryInputArchive;
    };

} // namespace archive



/**
 * @brief Concept defining the requirements for a type to be serializable with Cereal.
 *
 * A type is considered to be Cereal supported if it is one of the following:
 * - `std::complex<float>`
 * - `std::string`
 * - A sequential container
 * - An associative container
 * - An Eigen matrix
 * - A type with `serialize` and `deserialize` member functions
 *
 * @tparam T The type to check.
 */
template<typename T>
concept CerealSupported = requires(T a) {
    requires std::is_same_v<T, std::complex<float>> ||
             std::is_same_v<T, std::string> ||
             SequentialContainer<T> ||
             AssociativeContainer<T> ||
             EigenMatrix<T> ||
             (requires {
                 {
                 a.serialize(std::declval<typename std::remove_reference_t<T>::output &>())
                 } -> std::same_as<void>;
             } && requires {
                 {
                 a.deserialize(std::declval<typename std::remove_reference_t<T>::input &>())
                 } -> std::same_as<void>;
             });
};

/**
 * @brief SerializationPolicy class for serializing and deserializing data.
 *
 * This class uses Cereal to serialize and deserialize data. The actual serialization
 * format is determined by the `Archive` template parameter.
 *
 * @tparam T The type of data to serialize/deserialize.
 * @tparam Archive The Cereal archive to use for serialization/deserialization.
 *                 Default is `archive::portable_binary`.
 */
template<typename T, typename Archive = archive::portable_binary> requires CerealSupported<T>
class SerializationPolicy {
public:

    /**
     * @brief Serializes an object into a string.
     *
     * @param object The object to serialize.
     * @return A string representing the serialized data.
     */
    std::string serialize(const T &object) {
        std::ostringstream os;
        {
            typename Archive::output oarchive(os);
            oarchive(object);
        }
        return os.str();
    }

    /**
     * @brief Deserializes an object from a string.
     *
     * @param data The serialized data.
     * @return The deserialized object.
     */
    T deserialize(const std::string &data) {
        T object;
        std::istringstream is(data);
        typename Archive::input iarchive(is);
        iarchive(object);
        return object;
    }
};

#endif // SERIALIZATION_POLICY_H