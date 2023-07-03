#include <tuple>

template<typename... Parameters>
class ParameterHandler {
public:
    using parameter_tuple = std::tuple<Parameters...>;  ///< Tuple type for storing parameters.

    template<typename T>
    [[nodiscard]] static auto &get(parameter_tuple &parameters) {
        return std::get<T>(parameters);
    }

    template<typename T>
    [[nodiscard]] static const auto &get(const parameter_tuple &parameters) {
        return std::get<T>(parameters);
    }

    template<typename ID, typename Tuple, std::size_t... I>
    auto get_parameters_of_entity_impl(Tuple &&t, std::index_sequence<I...>) {
        return std::tuple_cat(
                (std::is_same_v<ID, typename std::tuple_element_t<I, std::decay_t<Tuple>>::id_type> ?
                 std::make_tuple(std::get<I>(t)) : std::tuple<>())...);
    }

    template<typename ID, typename Tuple>
    auto get_parameters_of_entity(Tuple &&t) {
        return get_parameters_of_entity_impl<ID>(
                std::forward<Tuple>(t),
                std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{});
    }

    [[nodiscard]] static auto split(parameter_tuple &parameters) {
        return std::apply([](auto &... parameter) {
            return std::make_tuple(std::addressof(parameter)...);
        }, parameters);
    }

    static void combine(parameter_tuple &old_parameters, parameter_tuple &&new_parameters) {
        old_parameters = std::move(new_parameters);
    }

    template<typename... SpecificParameters>
    static void recombine(parameter_tuple &old_parameters, std::tuple<SpecificParameters...> &&new_parameters) {
        std::apply([&old_parameters](auto &&... specific_parameter) {
            ((std::get<std::remove_reference_t<decltype(specific_parameter)>>(
                    old_parameters) = std::forward<decltype(specific_parameter)>(specific_parameter)), ...);
        }, std::move(new_parameters));
    }
};

#define DEFAULT_SCALAR double


template<typename ID, typename Scalar = DEFAULT_SCALAR, bool IsEstimable = false>
struct Parameter {
    using id_type = ID;
    using scalar_type = Scalar;
    static constexpr bool is_estimable = IsEstimable;

    // Additional functionality as needed
};

template<typename Derived, typename ID, typename Scalar = DEFAULT_SCALAR>
struct MeasurementErrorModel {
    using id_type = ID;
    using scalar_type = Scalar;

    // The implementation that works with a ParameterHandler instance
    template<typename... Parameters>
    scalar_type compute_error(ParameterHandler<Parameters...> &handler) const {
        // We use decltype to get the exact type of the parameter that the model expects
        // This would typically be based on the type of Derived
        using expected_parameter_types = typename Derived::expected_parameter_types;
        auto params = handler.template get_parameters_of_entity<expected_parameter_types>(handler.parameters);
        return static_cast<const Derived *>(this)->compute_error_impl(params);
    }

};

template<typename ID, typename Scalar = DEFAULT_SCALAR>
struct ClockInstabilityErrorModel : MeasurementErrorModel<ClockInstabilityErrorModel<ID, Scalar>, ID, Scalar> {
    // This struct may need more parameters to model the error correctly.
    // You should add them as needed.
    Scalar compute_error_impl(const Parameter<ID, Scalar> &param) const override {
        // Implementation for generic clock instability error goes here
        // This should be replaced or removed when specializing this struct
        return Scalar{};
    }
};

template<typename Scalar = DEFAULT_SCALAR>
struct FrequencyOffset : Parameter<Scalar, Scalar, false> {
    using id_type = Scalar;
    Scalar value;

    explicit FrequencyOffset(const Scalar &value) : value(value) {}
    // Add any functionality specific to FrequencyOffset here
};

template<typename Scalar = DEFAULT_SCALAR>
struct NominalFrequency : Parameter<Scalar, Scalar, false> {
    using id_type = Scalar;
    Scalar value;

    explicit NominalFrequency(const Scalar &value) : value(value) {}
    // Add any functionality specific to NominalFrequency here
};

template<typename Scalar = DEFAULT_SCALAR>
struct CountTime : Parameter<Scalar, Scalar, true> {
    using id_type = Scalar;
    Scalar value;

    explicit CountTime(const Scalar &value) : value(value) {}
    // Add any functionality specific to CountTime here
};

template<typename Scalar = DEFAULT_SCALAR>
struct RTLT : Parameter<Scalar, Scalar, true> {
    using id_type = Scalar;
    Scalar value;

    explicit RTLT(const Scalar &value) : value(value) {}
    // Add any functionality specific to RTLT here
};

template<typename Scalar = DEFAULT_SCALAR>
struct ClockEpochError : Parameter<Scalar, Scalar, true> {
    using id_type = Scalar;
    Scalar value;

    explicit ClockEpochError(const Scalar &value) : value(value) {}
    // Add any functionality specific to ClockEpochError here
};

template<typename Scalar = DEFAULT_SCALAR>
struct ClockRateOffset : Parameter<Scalar, Scalar, true> {
    using id_type = Scalar;
    Scalar value;

    explicit ClockRateOffset(const Scalar &value) : value(value) {}
    // Add any functionality specific to ClockRateOffset here
};

template<typename Scalar = DEFAULT_SCALAR>
struct FractionalFrequencyDeviation : Parameter<Scalar, Scalar, true> {
    using id_type = Scalar;
    Scalar y_k, y_k_next;

    FractionalFrequencyDeviation(const Scalar &y_k, const Scalar &y_k_next)
            : y_k(y_k), y_k_next(y_k_next) {}
    // Add any functionality specific to FractionalFrequencyDeviation here
};
