#include "core/policy/container/container_policy.hpp"
#include "core/policy/memory/memory_policy.hpp"
#include "core/policy/memory/default_memory_policy.hpp"
#include <memory>
#include <type_traits>
#include <utility>


template<typename T, typename Container, typename ContainerPolicyType = ContainerPolicy<Container>, typename MemoryPolicy = DefaultMemoryPolicy<T, Container>, typename Allocator = std::allocator<T>>
class History {
public:
    using MemoryPolicyPtr = std::shared_ptr<MemoryPolicy>;
    using AllocatorType = typename std::allocator_traits<Allocator>::template rebind_alloc<T>;
    AllocatorType allocator_;

    explicit History(MemoryPolicyPtr memory_policy = nullptr)
            : memory_policy_(memory_policy ? std::move(memory_policy) : std::make_shared<MemoryPolicy>()) {}

    template<typename PolicyArg>
    explicit History(PolicyArg &&arg)
            : container_policy_(std::forward<PolicyArg>(arg)) {}

    virtual ~History() = default;

    template<typename... Args>
    void update(Args &&...args) {
        container_policy_.add(container_, std::forward<Args>(args)...);
        memory_policy_->store(container_);
    }

    void clear() {
        container_policy_.clear(container_);
        memory_policy_->clear(container_);
    }

    [[nodiscard]] Container get_container() const {
        if constexpr (std::is_array_v<Container>) {
            return container_policy_.get_ordered_container(container_);
        } else {
            return container_;
        }
    }

    [[nodiscard]] size_t size_in_bytes() const {
        size_t size = container_policy_.size_in_bytes(container_);  // Size of the container
        return size;
    }

    [[nodiscard]] Container move_container() {
        if constexpr (std::is_array_v<Container>) {
            return container_policy_.get_ordered_container(std::move(container_));
        } else {
            return std::move(container_);
        }
    }

protected:
    ContainerPolicyType container_policy_;   ///< The policy for manipulating the container.
    Container container_;                    ///< The container for keeping a history of all values.
    MemoryPolicyPtr memory_policy_;          ///< The memory policy for storing the container.
};
