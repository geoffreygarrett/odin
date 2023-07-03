// Axis namespace
namespace frame {
    namespace handedness {
        struct right {
        };
        struct left {
        };
    }

    namespace axis {

        struct x {
        };
        struct y {
        };
        struct z {
        };
        namespace neg {
            struct x {
            };
            struct y {
            };
            struct z {
            };
        }
    };

    template<std::size_t I>
    struct arg {
        static constexpr std::size_t value = I;
    };

    // Alignment structure
    template<typename AxisType, std::size_t I>
    struct align {
        using axis = AxisType;
        using argument = arg<I>;
    };

}

template<typename T>
struct print_compile_time_type;


// Alignment Traits
template<typename T, std::size_t Dimension>
struct align_traits;

// Specialization for 2D
template<typename T>
struct align_traits<T, 2> {
    template<typename AlignType, typename... Args>
    static constexpr Eigen::Matrix<T, 2, 1> compute_vec(Args... args) {
        constexpr std::size_t index = AlignType::argument::value;
        static_assert(index <= sizeof...(args), "Arg index is out of range");

        Eigen::Matrix<T, 2, 1> vec = std::get<index>(std::tie(args...));
        if constexpr (std::is_same_v<typename AlignType::axis, frame::axis::neg::x> ||
                      std::is_same_v<typename AlignType::axis, frame::axis::neg::y>)
            vec = -vec;

        return vec.normalized();
    }
};

// Specialization for 3D
template<typename T>
struct align_traits<T, 3> {
    template<typename AlignType, typename... Args>
    static constexpr Eigen::Matrix<T, 3, 1> compute_vec(Args... args) {
        constexpr std::size_t index = AlignType::argument::value;
        static_assert(index <= sizeof...(args), "Arg index is out of range");

        Eigen::Matrix<T, 3, 1> vec = std::get<index>(std::tie(args...));
        if constexpr (std::is_same_v<typename AlignType::axis, frame::axis::neg::x> ||
                      std::is_same_v<typename AlignType::axis, frame::axis::neg::y> ||
                      std::is_same_v<typename AlignType::axis, frame::axis::neg::z>)
            vec = -vec;

        return vec.normalized();
    }
};

// Function to check that all types in a parameter pack are the same
template<auto... values>
constexpr bool all_same_value() {
    return ((values == ...));
}


template<typename T>
struct is_axis_x {
    static const bool value = std::is_same_v<T, frame::axis::x> || std::is_same_v<T, frame::axis::neg::x>;
};

template<typename T>
struct is_axis_y {
    static const bool value = std::is_same_v<T, frame::axis::y> || std::is_same_v<T, frame::axis::neg::y>;
};

template<typename T>
struct is_axis_z {
    static const bool value = std::is_same_v<T, frame::axis::z> || std::is_same_v<T, frame::axis::neg::z>;
};


template<typename T>
constexpr bool is_valid_axis() {
    return is_axis_x<T>::value || is_axis_y<T>::value || is_axis_z<T>::value;
}

template<typename A, typename B, std::size_t Dimension>
constexpr bool is_valid_3D_combination() {
    if constexpr (Dimension == 3)
        return !(std::is_same_v<typename A::axis, typename B::axis>);
    else
        return true;
}


template<typename Handedness, typename A, typename B, typename Scalar>
void populate_rotation_matrix(
        Eigen::Matrix<Scalar, 2, 2> &R,
        const Eigen::Matrix<Scalar, 2, 1> &a,
        const Eigen::Matrix<Scalar, 2, 1> &b) {
    if constexpr (std::is_same_v<Handedness, frame::handedness::right>) {
        if constexpr (is_axis_x<typename A::axis>::value)
            R << a, b;
        else if constexpr (is_axis_y<typename A::axis>::value)
            R << b, a;
    } else {
        if constexpr (is_axis_x<typename A::axis>::value)
            R << b, a;
        else if constexpr (is_axis_y<typename A::axis>::value)
            R << a, b;
    }
}

template<typename Handedness, typename A, typename B, typename Scalar>
void populate_rotation_matrix(
        Eigen::Matrix<Scalar, 3, 3> &R,
        const Eigen::Matrix<Scalar, 3, 1> &a,
        const Eigen::Matrix<Scalar, 3, 1> &b) {
    // Compute the third axis using the cross product
    Eigen::Matrix<Scalar, 3, 1> c;
    if constexpr (std::is_same_v<Handedness, frame::handedness::right>)
        c = a.cross(b);
    else
        c = b.cross(a);

    if constexpr (is_axis_x<typename A::axis>::value)
        R << a, b, c;
    else if constexpr (is_axis_y<typename A::axis>::value)
        R << c, a, b;
    else if constexpr (is_axis_z<typename A::axis>::value)
        R << b, c, a;
}


/**
 * @brief Generates a lambda function to compute a rotation matrix.
 *
 * This function template uses template parameters to generate a lambda function
 * that computes a rotation matrix. The rotation matrix aligns a chosen
 * transformation axis with a primary vector direction in the original frame.
 * The secondary axis is then derived from a second input vector through
 * orthogonal projection of the primary vector.
 *
 * @tparam A The primary axis in the original frame. Must be one of `frame::axis::x`,
 * `frame::axis::y`, or `frame::axis::z`.
 *
 * @tparam B The secondary axis in the original frame. Must be one of `frame::axis::x`,
 * `frame::axis::y`, or `frame::axis::z`. The default is `frame::axis::x`.
 *
 * @tparam Handedness Specifies the handedness of the coordinate system.
 * Must be one of `frame::handedness::right` or `frame::handedness::left`.
 * The default is `frame::handedness::right`.
 *
 * For 2D vectors, only one axis can be specified explicitly. The secondary axis
 * is automatically chosen as the remaining axis out of `x` and `y`. The z-direction
 * (inward or outward from the page) can only be set via the `Handedness` template parameter.
 *
 * @return A lambda function that takes in a variable number of arguments and computes
 * the corresponding rotation matrix.
 *
 * @note All input arguments must have the same dimension (i.e., they must all be either 2D or 3D vectors).
 */
template<typename A, typename B = frame::align<frame::axis::x, 0>, typename Handedness = frame::handedness::right>
auto make_rot_mat_lambda() {
    return [](auto... args) constexpr {
        // Ensure that all arguments have the same dimension
        static_assert(
                all_same_value<std::decay_t<decltype(args)>::RowsAtCompileTime...>(),
                "All input arguments must have the same dimension (2D or 3D vectors)");

        // Infer the Matrix type and dimension from the first argument, now that we know they're all the same.
        using ArgType = std::decay_t<decltype(std::get<0>(std::tie(args...)))>;
        using Scalar = typename ArgType::Scalar;
        constexpr std::size_t Dimension = ArgType::RowsAtCompileTime;

        // Static assertions
        static_assert(is_valid_axis<typename A::axis>(), "Unsupported primary axis type");
        static_assert(is_valid_axis<typename B::axis>(), "Unsupported secondary axis type");

        if constexpr (Dimension == 2) {
            static_assert(is_valid_3D_combination<A, B, Dimension>(),
                          "Primary and secondary axes cannot be the same in 3D space");
        }

        // Compute the primary vector
        auto a = align_traits<Scalar, Dimension>::template compute_vec<A>(args...);

        // Compute the secondary vector
        auto b_unproj = align_traits<Scalar, Dimension>::template compute_vec<B>(args...);
        auto b = (b_unproj - a * a.dot(b_unproj)).normalized();

        // Assemble the rotation matrix
        Eigen::Matrix<Scalar, Dimension, Dimension> R;
        populate_rotation_matrix<Handedness, A, B, Scalar>(R, a, b);
        return R;

    };
}
