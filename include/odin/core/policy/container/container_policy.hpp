#ifndef CONTAINER_POLICY_H
#define CONTAINER_POLICY_H

#include <array>
#include <map>
#include <vector>
#include <type_traits>
#include <cstddef>

template<typename Derived, typename Container>
class ContainerPolicyBase {
public:
    template<typename... Args>
    void add(Container &container, Args &&... args) {
        return static_cast<Derived *>(this)->add_impl(container, std::forward<Args>(args)...);
    }

    void clear(Container &container) {
        return static_cast<Derived *>(this)->clear_impl(container);
    }

    /**
     * @brief Returns the size in bytes of the provided container
     *
     * This method delegates the calculation of the size in bytes to the
     * derived class using the Curiously Recurring Template Pattern (CRTP).
     * In the initial implementation, a static_cast was used to upcast the
     * this pointer to the Derived type, which led to a compilation error
     * ("casts away qualifiers") when the method was invoked on a const
     * object.
     *
     * The issue was solved by using static_cast to cast this to a const
     * Derived*. This preserves the const-ness of the object, allowing the
     * method to be invoked on const objects without raising the compilation
     * error.
     *
     * @param container The container whose size in bytes is to be computed
     * @return The size in bytes of the provided container
     */
    [[nodiscard]] size_t size_in_bytes(const Container &container) const {
        return static_cast<const Derived *>(this)->size_in_bytes_impl(container);
    }

    template<typename Archive>
    void save(Archive & ar) const {
        static_cast<Derived const*>(this)->save_impl(ar);
    }

    template<typename Archive>
    void load(Archive & ar) {
        static_cast<Derived*>(this)->load_impl(ar);
    }
};

template<typename Container>
struct ContainerPolicy;

#endif // CONTAINER_POLICY_H
