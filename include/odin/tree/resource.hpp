#ifndef ODIN_RESOURCE_HPP
#define ODIN_RESOURCE_HPP

#include <chrono>


namespace odin {

// Concept for constraints
template<typename T>
concept ConstraintConcept = requires(const T &constraint) {
    { constraint.should_terminate() } -> std::convertible_to<bool>;
};


template<typename Derived>
struct Constraint {
    [[nodiscard]] bool should_terminate() const {
        return static_cast<const Derived *>(this)->should_terminate_impl();
    }
};

// TimeConstraint derived from Constraint
struct TimeConstraint : Constraint<TimeConstraint> {
    std::chrono::duration<double>                  max_duration;
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

    explicit TimeConstraint(double max_duration) : max_duration(max_duration) {}

    [[nodiscard]] bool should_terminate_impl() const {
        return std::chrono::high_resolution_clock::now() >= start_time + max_duration;
    }
};

// IterationConstraint derived from Constraint
struct IterationConstraint : Constraint<IterationConstraint> {
    size_t max_iterations;
    size_t iterations = 0;

    explicit IterationConstraint(size_t max_iterations) : max_iterations(max_iterations) {}

    [[nodiscard]] bool should_terminate_impl() const {
        return iterations >= max_iterations;
    }
};


// VirtualConstraint for further inheritance
struct VirtualConstraint : Constraint<VirtualConstraint> {
    // Pure virtual function
    [[nodiscard]] virtual bool should_terminate_impl() const = 0;
};

}// namespace odin

#endif//ODIN_RESOURCE_HPP