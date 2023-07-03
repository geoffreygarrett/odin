#ifndef MEMORY_POLICY_H
#define MEMORY_POLICY_H

#include <string>
#include <vector>

template<typename Derived, typename Data>
class MemoryPolicy {
public:
    template<typename... Args>
    void store(const Data &data, Args &&... args) {
        return static_cast<Derived *>(this)->store_impl(data, std::forward<Args>(args)...);
    }

    template<typename... Args>
    Data retrieve(Args &&... args) {
        return static_cast<Derived *>(this)->retrieve_impl(std::forward<Args>(args)...);
    }

    template<typename... Args>
    void delete_data(Args &&... args) {
        return static_cast<Derived *>(this)->delete_data_impl(std::forward<Args>(args)...);
    }
};

#endif
