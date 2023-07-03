#ifndef DEFAULT_MEMORY_POLICY_H
#define DEFAULT_MEMORY_POLICY_H

#include "memory_policy.hpp"
#include <memory>

template<typename T, typename Container>
class DefaultMemoryPolicy : public MemoryPolicy<DefaultMemoryPolicy<T, Container>, Container> {
public:
    DefaultMemoryPolicy() : stored_data_(std::make_shared<Container>()) {}

    void store_impl(const Container &container) {
        *stored_data_ = container;
    }

    Container retrieve_impl() {
        return *stored_data_;
    }

    void delete_data_impl() {
        stored_data_->clear();
    }

private:
    std::shared_ptr<Container> stored_data_;
};

#endif
