//static_assert(
//        static_cast<int>(X::ColsAtCompileTime) == 1,
//        "X (`state`) must be a column vector"
//);
//
//static_assert(
//        static_cast<int>(M::ColsAtCompileTime) == 1,
//        "M (`measurement`) must be a column vector"
//);
//
//static_assert(
//        static_cast<int>(X::RowsAtCompileTime) == static_cast<int>(P::RowsAtCompileTime) &&
//        static_cast<int>(X::RowsAtCompileTime) == static_cast<int>(P::ColsAtCompileTime) &&
//        static_cast<int>(P::RowsAtCompileTime) == static_cast<int>(P::ColsAtCompileTime),
//        "P (`state_covariance` matrix) must be a square matrix with dimensions matching X (`state` column vector)"
//);


/**
 * @brief Base filter class for various types of filters, designed using Curiously Recurring Template Pattern (CRTP).
 *
 * This class template acts as a base class for filters. It enforces a consistent API across all derived filter types,
 * while allowing the derived classes to implement the specific logic of each filter type. The CRTP idiom is used to
 * achieve static (compile-time) polymorphism, which provides the benefits of polymorphism without the overhead of
 * dynamic dispatch (i.e., virtual functions).
 *
 * The base class asserts that the dimensions of the following types are compatible with each other.
 *
 * @tparam D The derived class, used for static polymorphism via CRTP.
 * @tparam X The type representing the state.
 * @tparam M The type representing the measurements.
 * @tparam P The type representing the state covariance.
 */
template<typename D, typename X, typename M, typename P>
class Filter {
public:

    /**
     * @brief Initialize the state and covariance matrix of the filter.
     *
     * @param initial_state Initial state vector.
     * @param initial_covariance Initial covariance matrix.
     */
    void initialize(const X &initial_state, const P &initial_covariance) {
        static_cast<D *>(this)->initialize_impl(initial_state, initial_covariance);
    }

    /**
     * @brief Perform the predict step of the filter.
     *
     * This method takes a variable number of arguments and forwards them to the predict_impl() method of the derived class.
     * This allows for a high degree of flexibility in the interface.
     *
     * @tparam Args The types of the arguments for the predict_impl() method in the derived class.
     * @param args The arguments for the predict_impl() method in the derived class.
     */
    template<typename... Args>
    void predict(Args &&... args) {
        static_cast<D *>(this)->predict_impl(std::forward<Args>(args)...);
    }

    /**
     * @brief Perform the update step of the filter.
     *
     * This method takes a variable number of arguments and forwards them to the update_impl() method of the derived class.
     * This allows for a high degree of flexibility in the interface.
     *
     * @tparam Args The types of the arguments for the update_impl() method in the derived class.
     * @param args The arguments for the update_impl() method in the derived class.
     */
    template<typename... Args>
    void update(Args &&... args) {
        static_cast<D *>(this)->update_impl(std::forward<Args>(args)...);
    }

    /**
     * @brief Get the current state estimate.
     *
     * @return The state estimate.
     */
    X get_state_estimate() const {
        return static_cast<const D *>(this)->get_state_estimate_impl();
    }

    /**
     * @brief Get the current covariance estimate.
     *
     * @return The covariance estimate.
     */
    P get_covariance_estimate() const {
        return static_cast<const D *>(this)->get_covariance_estimate_impl();
    }
};
