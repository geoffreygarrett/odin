/**
 * @file state_types.hpp
 * @brief This file contains various state classes.
 * @note Currently, state classes are designed for stack allocation, thus they do not use the EIGEN_MAKE_ALIGNED_OPERATOR_NEW macro.
 *       This macro ensures proper memory alignment for dynamic allocations of Eigen types, enabling efficient SIMD optimizations.
 *       However, it's only necessary if these classes are dynamically allocated with the new operator.
 *
 *       The use of the EIGEN_MAKE_ALIGNED_OPERATOR_NEW macro doesn't affect current usage, as instances are stack-allocated.
 *       Nevertheless, if dynamic allocation is introduced in the future, each class with fixed-size Eigen types will need
 *       to include the EIGEN_MAKE_ALIGNED_OPERATOR_NEW macro to ensure proper alignment and optimal performance.
 *
 * @todo Should dynamic allocation be required in the future, include EIGEN_MAKE_ALIGNED_OPERATOR_NEW in relevant classes.
 */
#ifndef STATE_TYPES_HPP
#define STATE_TYPES_HPP

#include <Eigen/Core>
#include <utility>

template<typename ID, int Dimension, typename Scalar>
struct VecType {
    using type = Eigen::Vector<Scalar, Dimension>;
};

template<typename ID>
struct VecType<ID, 2, double> {
    using type = Eigen::Vector2d;
};

template<typename ID>
struct VecType<ID, 3, double> {
    using type = Eigen::Vector3d;
};

template<typename ID>
struct VecType<ID, 2, float> {
    using type = Eigen::Vector2f;
};

template<typename ID>
struct VecType<ID, 3, float> {
    using type = Eigen::Vector3f;
};


/**
 * @brief Base state class trait.
 */

template<typename T>
concept StateConcept = requires {
    { T::is_scalar } -> std::convertible_to<bool>;
    { T::space_dimension } -> std::convertible_to<int>;
    { T::state_dimension } -> std::convertible_to<int>;
    { T::derivative_order } -> std::convertible_to<int>;
    { T::is_translational } -> std::convertible_to<bool>;
    { T::is_rotational } -> std::convertible_to<bool>;
    { T::is_dynamic_property } -> std::convertible_to<bool>;

    // Space Dimension check
    requires (T::is_scalar || (T::space_dimension != -1));

    // State Dimension check
    requires (T::state_dimension != -1);

    // Derivative Order check
    requires (T::derivative_order >= 0 && T::derivative_order <= 1);

    // Type Check
    requires (T::is_translational || T::is_rotational || T::is_dynamic_property);
};

// Base template
template<typename T>
struct dimension_trait;

// Specialization for common scalar types
template<>
struct dimension_trait<float> {
    static constexpr int value = 1;
};

template<>
struct dimension_trait<double> {
    static constexpr int value = 1;
};

template<>
struct dimension_trait<int> {
    static constexpr int value = 1;
};

template<>
struct dimension_trait<std::complex<float>> {
    static constexpr int value = 1;
};

template<>
struct dimension_trait<std::complex<double>> {
    static constexpr int value = 1;
};

// Specialization for Eigen::Matrix types
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct dimension_trait<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
    static constexpr int value = Rows * Cols;  // The total size of the matrix
};

template<typename T>
struct is_scalar {
    static constexpr bool value = false;
};

template<typename T>
struct space_dimension {
    static constexpr int value = -1;
};

template<typename T>
struct derivative_order {
    static constexpr int value = -1;
};

template<typename T>
struct is_translational {
    static constexpr bool value = false;
};

template<typename T>
struct is_rotational {
    static constexpr bool value = false;
};

template<typename T>
struct is_dynamic_property {
    static constexpr bool value = false;
};

//template<typename ID, typename Scalar>
//struct is_scalar<Mass<ID, Scalar>> {
//static constexpr bool value = true;
//};


template<typename Derived>
struct StateBase {
//    using id_type = typename Derived::id_type;
//    using type = typename Derived::type;
//    using scalar_type = typename Derived::scalar_type;
//    static constexpr int state_dimension = dimension_trait<type>::value;
};

//template<typename Derived>
//struct StateBase {
//    // ...
//};

template<typename ID, typename Scalar = DEFAULT_SCALAR>
struct Mass : StateBase<Mass<ID, Scalar>> {
    using scalar_type = Scalar;
    using type = Scalar;
    type value;

    static constexpr bool is_scalar = true;
    static constexpr int space_dimension = -1;  // Not applicable for scalar types
    static constexpr int state_dimension = dimension_trait<type>::value;

    static constexpr int derivative_order = 0;
    static constexpr bool is_translational = false;
    static constexpr bool is_rotational = false;
    static constexpr bool is_dynamic_property = true;

    explicit Mass(const Scalar &mass) : value(mass) {
    }
};


template<typename ID, typename Scalar = DEFAULT_SCALAR>
struct MassRate : StateBase<MassRate<ID, Scalar>> {
    using scalar_type = Scalar;
    using type = Scalar;
    type value;

    static constexpr bool is_scalar = true;
    static constexpr int space_dimension = -1;  // Not applicable for scalar types
    static constexpr int state_dimension = dimension_trait<type>::value;
    static constexpr int derivative_order = 1;
    static constexpr bool is_translational = false;
    static constexpr bool is_rotational = false;
    static constexpr bool is_dynamic_property = true;

    explicit MassRate(const Scalar &mass_rate) : value(mass_rate) {
    }
};


template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct Inertia : StateBase<Inertia<ID, Dimension, Scalar>> {
    static_assert(Dimension == 2 || Dimension == 3, "Dimension must be 2 or 3 for Inertia.");
    using scalar_type = Scalar;
    using type = typename std::conditional<Dimension == 2, Scalar, Eigen::Matrix<Scalar, 3, 3>>::type;
    type value;

    static constexpr bool is_scalar = false;
    static constexpr int space_dimension = Dimension;
    static constexpr int state_dimension = dimension_trait<type>::value;
    static constexpr int derivative_order = 0;
    static constexpr bool is_translational = false;
    static constexpr bool is_rotational = true;
    static constexpr bool is_dynamic_property = true;

    Inertia(const std::initializer_list<Scalar> &init_values) {
        if constexpr (Dimension == 3) {
            static_assert(init_values.size() == 9, "Initializer list size must be 9 for 3D Inertia Matrix.");
            Eigen::Map<const Eigen::Matrix<Scalar, 3, 3>> map(init_values.begin());
            value = map;
        } else if constexpr (Dimension == 2) {
            static_assert(init_values.size() == 1, "Initializer list size must be 1 for 2D Inertia Scalar.");
            value = *init_values.begin();
        }
    }

    [[nodiscard]] auto to_vector() const {
        if constexpr (Dimension == 2) {
            return Eigen::Matrix<Scalar, 1, 1>::Constant(value);
        } else if constexpr (Dimension == 3) {
            const Eigen::Map<const Eigen::Matrix<Scalar, 3, 3, Eigen::RowMajor>> map(value.data());
            return map.template cast<Scalar>().reshaped().
                    template transpose<1, 0>().reshaped();
        }
    }

    void from_vector(const Eigen::Matrix<Scalar, (Dimension == 2) ? 1 : 9, 1> &vec) {
        if constexpr (Dimension == 2) {
            value = vec(0);
        } else if constexpr (Dimension == 3) {
            Eigen::Map<const Eigen::Matrix<Scalar, 3, 3, Eigen::RowMajor>> map(vec.data());
            value = map;
        }
    }
};

template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct InertiaRate : StateBase<InertiaRate<ID, Dimension, Scalar>> {
    static_assert(Dimension == 2 || Dimension == 3, "Dimension must be 2 or 3 for InertiaRate.");
    using scalar_type = Scalar;
    using type = typename std::conditional<Dimension == 2, Scalar, Eigen::Matrix<Scalar, 3, 3>>::type;
    type value;

    static constexpr bool is_scalar = false;
    static constexpr int space_dimension = Dimension;
//    static constexpr int state_dimension = (Dimension == 2) ? 1 : 9;
    static constexpr int derivative_order = 1;
    static constexpr bool is_translational = false;
    static constexpr bool is_rotational = true;
    static constexpr bool is_dynamic_property = true;

    InertiaRate(const std::initializer_list<Scalar> &init_values) {
        if constexpr (Dimension == 3) {
            static_assert(init_values.size() == 9, "Initializer list size must be 9 for 3D InertiaRate Matrix.");
            Eigen::Map<const Eigen::Matrix<Scalar, 3, 3>> map(init_values.begin());
            value = map;
        } else if constexpr (Dimension == 2) {
            static_assert(init_values.size() == 1, "Initializer list size must be 1 for 2D InertiaRate Scalar.");
            value = *init_values.begin();
        }
    }
};

/**
 * @section Translational State Types
 */

template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct Position : StateBase<Position<ID, Dimension, Scalar>> {
    static_assert(Dimension >= 2 && Dimension <= 3, "Dimension must be either 2 or 3.");
    using scalar_type = Scalar;
    using type = typename VecType<ID, Dimension, Scalar>::type;
    type value;

    static constexpr bool is_scalar = false;
    static constexpr int space_dimension = Dimension;
    static constexpr int state_dimension = dimension_trait<type>::value;
    static constexpr int derivative_order = 0;
    static constexpr bool is_translational = true;
    static constexpr bool is_rotational = false;
    static constexpr bool is_dynamic_property = false;

    template<typename... Args>
    explicit Position(Args... args) {
        static_assert(sizeof...(args) == Dimension, "Number of arguments must match Dimension.");
        value = type{args...};
    }
};

template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct Velocity : StateBase<Velocity<ID, Dimension, Scalar>> {
    static_assert(Dimension >= 2 && Dimension <= 3, "Dimension must be either 2 or 3.");
    using scalar_type = Scalar;
    using type = typename VecType<ID, Dimension, Scalar>::type;
    type value;

    static constexpr bool is_scalar = false;
    static constexpr int space_dimension = Dimension;
    static constexpr int state_dimension = dimension_trait<type>::value;
    static constexpr int derivative_order = 0;
    static constexpr bool is_translational = true;
    static constexpr bool is_rotational = false;
    static constexpr bool is_dynamic_property = false;

    template<typename... Args>
    explicit Velocity(Args... args) {
        static_assert(sizeof...(args) == Dimension, "Number of arguments must match Dimension.");
        value = type{args...};
    }
};

/**
 * @section Rotational State Types
 */


template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct Rotation : StateBase<Rotation<ID, Dimension, Scalar>> {
    static_assert(Dimension == 2 || Dimension == 3, "Dimension must be 2 or 3 for Rotation.");
    using scalar_type = Scalar;
    using type = typename std::conditional<Dimension == 2, std::complex<Scalar>, Eigen::Quaternion<Scalar>>::type;
    type value;

    static constexpr bool is_scalar = false;
    static constexpr int space_dimension = Dimension;
    static constexpr int state_dimension = dimension_trait<type>::value;
    static constexpr int derivative_order = 0;
    static constexpr bool is_translational = false;
    static constexpr bool is_rotational = true;
    static constexpr bool is_dynamic_property = false;

    Rotation(const std::initializer_list<Scalar> &init_values) {
        if constexpr (Dimension == 3) {
            assert(init_values.size() == 4 && "Initializer list size must be 4 for Quaternion.");
            std::array<Scalar, 4> temp_array;
            std::copy(init_values.begin(), init_values.end(), temp_array.begin());
            value = Eigen::Quaternion<Scalar>(temp_array[0], temp_array[1], temp_array[2], temp_array[3]);
        } else if constexpr (Dimension == 2) {
            assert(init_values.size() == 2 && "Initializer list size must be 2 for Complex.");
            value = std::complex<Scalar>{init_values.begin()[0], init_values.begin()[1]};
        }
    }
};

template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct AngularVelocity : StateBase<AngularVelocity<ID, Dimension, Scalar>> {
    static_assert(Dimension == 2 || Dimension == 3, "Dimension must be 2 or 3 for AngularVelocity.");
    using scalar_type = Scalar;
    using type = typename std::conditional<Dimension == 2, Scalar, Eigen::Matrix<Scalar, 3, 1>>::type;
    type value;

    static constexpr bool is_scalar = false;
    static constexpr int space_dimension = Dimension;
    static constexpr int state_dimension = dimension_trait<type>::value;
    static constexpr int derivative_order = 0;
    static constexpr bool is_translational = false;
    static constexpr bool is_rotational = true;
    static constexpr bool is_dynamic_property = false;

    AngularVelocity(const std::initializer_list<Scalar> &init_values) {
        if constexpr (Dimension == 3) {
            assert(init_values.size() == 3 && "Initializer list size must be 3 for 3D angular velocity.");
            value = Eigen::Matrix<Scalar, 3, 1>{init_values.begin()[0], init_values.begin()[1], init_values.begin()[2]};
        } else if constexpr (Dimension == 2) {
            assert(init_values.size() == 1 && "Initializer list size must be 1 for 2D angular velocity.");
            value = *init_values.begin();
        }
    }
};


struct TEST_STATE {
};
//static_assert(StateConcept<Mass<TEST_STATE>> == true, "Mass does not comply with StateConcept.");
//static_assert(StateConcept<MassRate<TEST_STATE>> == true, "Mass does not comply with StateConcept.");
//static_assert(StateConcept<Inertia<TEST_STATE>> == true, "Inertia does not comply with StateConcept.");
//static_assert(StateConcept<InertiaRate<TEST_STATE>> == true, "InertiaRate does not comply with StateConcept.");
//static_assert(StateConcept<Position<TEST_STATE>> == true, "Position does not comply with StateConcept.");
//static_assert(StateConcept<Velocity<TEST_STATE>> == true, "Velocity does not comply with StateConcept.");
//static_assert(StateConcept<Rotation<TEST_STATE>> == true, "Rotation does not comply with StateConcept.");
//static_assert(StateConcept<AngularVelocity<TEST_STATE>> == true, "AngularVelocity does not comply with StateConcept.");
#endif // STATE_TYPES_HPP