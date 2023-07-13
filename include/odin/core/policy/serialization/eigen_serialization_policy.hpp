/**
 * @file
 * @brief Defines serialization and deserialization for Eigen matrices with Cereal.
 */
#ifndef EIGEN_SERIALIZATION_POLICY_H
#define EIGEN_SERIALIZATION_POLICY_H

#include <Eigen/Dense>
#include <cereal/archives/binary.hpp>
#include <sstream>
#include <type_traits>
#include <concepts>
//#include "serialization_policy.hpp"

/**
 * @brief Concept defining the requirements for an Eigen matrix.
 *
 * A type is considered an Eigen matrix if it satisfies the following requirements:
 * - It has a `data` member function that returns a pointer to the scalar values.
 * - It has `rows` and `cols` member functions that return the number of rows and columns respectively.
 * - It is derived from `Eigen::MatrixBase`.
 *
 * @tparam T The type to check.
 */
template<typename T>
concept EigenMatrix = requires(T a) {
    { a.data() } -> std::same_as<typename T::Scalar *>;
    { a.rows() } -> std::same_as<typename T::Index>;
    { a.cols() } -> std::same_as<typename T::Index>;
    requires std::is_base_of_v<Eigen::MatrixBase<T>, T>;
};

template<typename MatrixType, std::enable_if_t<
        MatrixType::RowsAtCompileTime == Eigen::Dynamic && MatrixType::ColsAtCompileTime == Eigen::Dynamic> * = nullptr>
void resize_matrix(MatrixType &matrix, int rows, int cols) {
    matrix.resize(rows, cols);
}

template<typename MatrixType, std::enable_if_t<
        MatrixType::RowsAtCompileTime == Eigen::Dynamic && MatrixType::ColsAtCompileTime != Eigen::Dynamic> * = nullptr>
void resize_matrix(MatrixType &matrix, int rows, int cols) {
    matrix.resize(rows);
}

template<typename MatrixType, std::enable_if_t<
        MatrixType::RowsAtCompileTime != Eigen::Dynamic && MatrixType::ColsAtCompileTime == Eigen::Dynamic> * = nullptr>
void resize_matrix(MatrixType &matrix, int rows, int cols) {
    matrix.resize(1, cols); // Here matrix is like RowVector
}

template<typename MatrixType, std::enable_if_t<
        MatrixType::RowsAtCompileTime != Eigen::Dynamic && MatrixType::ColsAtCompileTime != Eigen::Dynamic> * = nullptr>
void resize_matrix(MatrixType &matrix, int rows, int cols) {
    // Fixed size, no need to resize
}

namespace cereal {

    /// Tag struct for distinguishing binary archives.
    struct BinaryTag {
    };

    /// Tag struct for distinguishing non-binary archives.
    struct NonBinaryTag {
    };

    /// Function returning the correct tag for binary archives.
    template<class Archive>
    typename std::enable_if<std::is_same<Archive, cereal::BinaryInputArchive>::value ||
                            std::is_same<Archive, cereal::BinaryOutputArchive>::value, BinaryTag>::type
    archive_tag() {
        return BinaryTag{};
    }

    /// Function returning the correct tag for non-binary archives.
    template<class Archive>
    typename std::enable_if<!std::is_same<Archive, cereal::BinaryInputArchive>::value &&
                            !std::is_same<Archive, cereal::BinaryOutputArchive>::value, NonBinaryTag>::type
    archive_tag() {
        return NonBinaryTag{};
    }

    /**
     * @brief Serialization function for Eigen matrices.
     *
     * @param ar The archive to serialize into.
     * @param m The matrix to serialize.
     */
    template<class Archive, class Derived>
    inline void save(Archive &ar, Eigen::PlainObjectBase<Derived> const &m) requires EigenMatrix<Derived> {
        save(ar, m, archive_tag<Archive>());
    }

    /**
     * @brief Serialization function for Eigen matrices for binary archives.
     *
     * @param ar The archive to serialize into.
     * @param m The matrix to serialize.
     */
    template<class Archive, class Derived>
    inline void save(Archive &ar, Eigen::PlainObjectBase<Derived> const &m, BinaryTag) requires EigenMatrix<Derived> {
        typedef Eigen::PlainObjectBase<Derived> ArrT;
        if (ArrT::RowsAtCompileTime == Eigen::Dynamic) ar(m.rows());
        if (ArrT::ColsAtCompileTime == Eigen::Dynamic) ar(m.cols());
        ar(binary_data(m.data(), m.size() * sizeof(typename Derived::Scalar)));
    }

    /**
     * @brief Serialization function for Eigen matrices for non-binary archives.
     *
     * @param ar The archive to serialize into.
     * @param m The matrix to serialize.
     */
    template<class Archive, class Derived>
    inline void
    save(Archive &ar, Eigen::PlainObjectBase<Derived> const &m, NonBinaryTag) requires EigenMatrix<Derived> {
        int rows = m.rows(), cols = m.cols();
        ar(cereal::make_nvp("rows", static_cast<size_type>(rows)));
        ar(cereal::make_nvp("cols", static_cast<size_type>(cols)));
        std::vector<std::vector<typename Derived::Scalar>> data(rows, std::vector<typename Derived::Scalar>(cols));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                data[i][j] = m(i, j);
            }
        }
        ar(cereal::make_nvp("data", data));
    }

    /**
     * @brief Deserialization function for Eigen matrices.
     *
     * @param ar The archive to deserialize from.
     * @param m The matrix to populate with deserialized data.
     */
    template<class Archive, class Derived>
    inline void load(Archive &ar, Eigen::PlainObjectBase<Derived> &m) requires EigenMatrix<Derived> {
        load(ar, m, archive_tag<Archive>());
    }

    /**
     * @brief Deserialization function for Eigen matrices for binary archives.
     *
     * @param ar The archive to deserialize from.
     * @param m The matrix to populate with deserialized data.
     */
    template<class Archive, class Derived>
    inline void load(Archive &ar, Eigen::PlainObjectBase<Derived> &m, BinaryTag) requires EigenMatrix<Derived> {
        typedef Eigen::PlainObjectBase<Derived> ArrT;
        Eigen::Index rows = ArrT::RowsAtCompileTime, cols = ArrT::ColsAtCompileTime;
        if (rows == Eigen::Dynamic) ar(rows);
        if (cols == Eigen::Dynamic) ar(cols);
        if (ArrT::RowsAtCompileTime == Eigen::Dynamic || ArrT::ColsAtCompileTime == Eigen::Dynamic) {
            m.resize(rows, cols);
        }
        ar(binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(typename Derived::Scalar))));
    }

    /**
     * @brief Deserialization function for Eigen matrices for non-binary archives.
     *
     * @param ar The archive to deserialize from.
     * @param m The matrix to populate with deserialized data.
     */
    template<class Archive, class Derived>
    inline void load(Archive &ar, Eigen::PlainObjectBase<Derived> &m, NonBinaryTag) requires EigenMatrix<Derived> {
        size_type rows, cols;
        ar(cereal::make_nvp("rows", rows));
        ar(cereal::make_nvp("cols", cols));
        resize_matrix(m, rows, cols);
        std::vector<std::vector<typename Derived::Scalar>> data;
        ar(cereal::make_nvp("data", data));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                m(i, j) = data[i][j];
            }
        }
        // TODO: I don't think its necessary to optimise this, but if really needed, we could try make this stateful
        //       with the policy class maybe, and keep a statically sized intermediate container for serde.
    }
}

#endif // EIGEN_SERIALIZATION_POLICY_H
