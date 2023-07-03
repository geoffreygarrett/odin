
Ok, I have a state handler, which I want to be so general and well interfaced so that I can pass it as a "state" to a function, which can then retrieve the state through static types. I face a problem however, I want the generic template classes where the state is being used, to be completely agnostic of this "state handler" and its api. For example, I want my "BaseFilter" (used for estimation filtering) to be easily applied to simple eigen types, AS well as to my complex state handler state. Please provide an aadvanced metaprogramming technique I could use, here are the relevant classes.

BASE STATE
```C++
template<typename Derived, typename ID>
struct StateBase {
    using id_type = ID;
    using is_translational = std::false_type;
    using is_rotational = std::false_type;
    using is_rigid_body = std::false_type;
    using is_point_mass = std::false_type;
    using is_celestial_body = std::false_type;
    using is_dynamic = std::false_type;
    using dim = std::integral_constant<int, -1>;
};
```

POSITION STATE
```C++
template<typename ID, int Dimension = DEFAULT_DIMENSION, typename Scalar = DEFAULT_SCALAR>
struct Position : StateBase<Position<ID, Dimension, Scalar>, ID> {
    static_assert(Dimension >= 2 && Dimension <= 3, "Dimension must be either 2 or 3.");
    using scalar_type = Scalar;
    using is_translational = std::true_type;
    using dim = std::integral_constant<int, Dimension>;
    using type = typename VecType<ID, Dimension, Scalar>::type;
    type value;

    template<typename... Args>
    explicit Position(Args... args) {
        static_assert(sizeof...(args) == Dimension, "Number of arguments must match Dimension.");
        std::array<scalar_type, Dimension> temp_array{args...};
        value = Eigen::Map<const type>(temp_array.data());
    }
};
```

STATE HANDLER
```C++
template<typename... States>
class StateHandler {
public:
    using state_tuple = std::tuple<States...>;  ///< Tuple type for storing states.

    template<typename T>
    [[nodiscard]] static auto &get(state_tuple &states) {
        return std::get<T>(states);
    }


    template<typename T>
    [[nodiscard]] static const auto &get(const state_tuple &states) {
        return std::get<T>(states);
    }


    template<typename ID, typename Tuple, std::size_t... I>
    auto get_states_of_entity_impl(Tuple &&t, std::index_sequence<I...>) {
        return std::tuple_cat(
                (std::is_same_v<ID, typename std::tuple_element_t<I, std::decay_t<Tuple>>::id_type> ?
                 std::make_tuple(std::get<I>(t)) : std::tuple<>())...);
    }

    template<typename ID, typename Tuple>
    auto get_states_of_entity(Tuple &&t) {
        return get_states_of_entity_impl<ID>(
                std::forward<Tuple>(t),
                std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{});
    }

    [[nodiscard]] static auto split(state_tuple &states) {
        return std::apply([](auto &... state) {
            return std::make_tuple(std::addressof(state)...);
        }, states);
    }

    static void combine(state_tuple &old_states, state_tuple &&new_states) {
        old_states = std::move(new_states);
    }

    template<typename... SpecificStates>
    static void recombine(state_tuple &old_states, std::tuple<SpecificStates...> &&new_states) {
        std::apply([&old_states](auto &&... specific_state) {
            ((std::get<std::remove_reference_t<decltype(specific_state)>>(
                    old_states) = std::forward<decltype(specific_state)>(specific_state)), ...);
        }, std::move(new_states));
    }

```

BASE FILTER

The member functions should be easily usable with either Eigen matrices, or a specific override that takes the global state handler for easy compile time selection of the global state.
```C++
template<typename D, typename X, typename M, typename P>
class Filter {
public:
);

    /**
     * @brief Initialize the state and covariance matrix of the filter.
     */
    void initialize(const X &initial_state, const P &initial_covariance) {
        static_cast<D *>(this)->initialize_impl(initial_state, initial_covariance);
    }

    /**
     * @brief Perform the predict step of the filter.
     */
    template<typename... Args>
    void predict(Args &&... args) {
        static_cast<D *>(this)->predict_impl(std::forward<Args>(args)...);
    }

    /**
     * @brief Perform the update step of the filter.
     */
    template<typename... Args>
    void update(Args &&... args) {
        static_cast<D *>(this)->update_impl(std::forward<Args>(args)...);
    }
    }
```

BASE MEASUREMENT MODEL
```C++
template<typename Derived, typename StateType, typename MeasurementType, typename JacobianType>
struct MeasurementModel {

//    static constexpr int state_dim = StateType::RowsAtCompileTime;
//    static constexpr int meas_dim = MeasurementType::RowsAtCompileTime;

    constexpr MeasurementType compute_measurement(const StateType &state) const {
        return static_cast<const Derived *>(this)->compute_measurement_impl(state);
    }

    template<typename T = Derived>
    constexpr auto compute_jacobian(const StateType &state) const
    -> std::enable_if_t<has_compute_analytical_jacobian<T, StateType, JacobianType>, JacobianType> {
        return static_cast<const T *>(this)->compute_analytical_jacobian(state);
    }

    template<typename T = Derived>
    constexpr auto compute_jacobian(const StateType &state) const
    -> std::enable_if_t<!has_compute_analytical_jacobian<T, StateType, JacobianType>, JacobianType> {
        return compute_autodiff_jacobian(state);
    }


    constexpr JacobianType compute_autodiff_jacobian(const StateType &state) const {
        using namespace autodiff;
        using Scalar = typename MeasurementType::Scalar;

        using ADStateType = Eigen::Matrix<autodiff::real, state_dim, 1>;
        ADStateType state_ad = state.template cast<autodiff::real>();

        // Wrapper lambda function to pass into autodiff::jacobian
        auto f = [&](const ADStateType &state) {
            return static_cast<const Derived *>(this)->compute_measurement_impl(state);
        };

        MeasurementType result;   // the output vector F = f(state, params) evaluated together with Jacobian below

        // Compute the Jacobian
        JacobianType jacobian = autodiff::jacobian(f, wrt(state_ad), at(state_ad), result);

        return jacobian;
    }
};
```

OUR POTENTIAL SOLUTION
```C++
template<typename... States>
struct StateVariantHandler {
    using StateVariant = std::variant<States...>;

    /**
      * @brief Applies a function to the state variant in a type-safe manner.
      *
      * @details
      * This function uses std::visit to apply a function to the state variant.
      * The actual function called will depend on the active member of the variant at runtime.
      *
      * @param function The function to apply.
      * @param state The state variant.
      *
      * @return The result of applying the function to the state variant.
      */
    template<typename Function>
    auto visit(Function &&function, const StateVariant &state) const {
        return std::visit(std::forward<Function>(function), state);
    }
};

```