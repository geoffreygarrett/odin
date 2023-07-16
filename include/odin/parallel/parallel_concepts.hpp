#ifndef ODIN_PARALLEL_CONCEPTS_HPP
#define ODIN_PARALLEL_CONCEPTS_HPP

#include <type_traits>

namespace odin {

    // This concept checks if an object T provides a method `thread_local_copy()`
    // that can be used to create a thread-safe copy of the object for parallel execution.
    // The copied object should be of the same type as the original object.
    template<typename T>
    concept ThreadLocalCopy = requires(T a) {
        { a.thread_local_copy() } -> std::same_as<T>;
    };

    // This concept checks if an object T provides a method `thread_local_reference()`
    // that returns a reference to an object of the same type, intended to be used as a thread-safe reference.
    // The method should return a reference to a thread-local instance of the object.
    template<typename T>
    concept ThreadLocalReference = requires(T a) {
        { a.thread_local_reference() } -> std::same_as<T &>;
    };

    // This concept checks if an object T provides a method `thread_local_pointer()`
    // that returns a raw pointer to an object of the same type, intended to be used as a thread-safe pointer.
    // The method should return a non-owning pointer to a thread-local instance of the object.
    template<typename T>
    concept ThreadLocalPointer = requires(T a) {
        { a.thread_local_pointer() } -> std::same_as<T *>;
    };

    // An object T is considered ThreadSafe if it can provide a thread-safe copy, reference, or pointer of itself.
    // This is useful in parallel computing where each thread needs its own instance of an object to prevent data races.
    template<typename T>
    concept ThreadSafe = ThreadLocalCopy<T> || ThreadLocalReference<T> || ThreadLocalPointer<T>;

    // This function attempts to obtain a thread-safe instance of the model object
    // for use in parallel computations. This is crucial in preventing data races
    // in multithreaded code.
    //
    // There are three methods by which this function tries to create a thread-safe
    // instance, in this order:
    //
    // 1. `thread_local_copy()`: This method creates a new copy of the object for each thread.
    //    This is the preferred method, as it ensures that each thread works on an independent
    //    instance of the object, eliminating the possibility of data races.
    //    This method should be implemented if the cost of copying the object is not high.
    //
    // 2. `thread_local_reference()`: This method returns a reference to an instance of the object.
    //    This is the second option, preferred when creating a copy of the object is expensive.
    //    However, this requires the referenced object to be inherently thread-safe or not modified
    //    in a way that would cause data races.
    //
    // 3. `thread_local_pointer()`: This is the last resort. It returns a raw pointer to an instance
    //    of the object. This is used when the object is large, and both copying and referencing are not suitable.
    //    Like with references, this method requires careful handling of the pointed object to prevent data races.
    //
    // Note: This function requires that the model object conforms to at least one of the three thread-safety
    // concepts defined above: `ThreadLocalCopy`, `ThreadLocalReference`, or `ThreadLocalPointer`.
    //    template<ThreadSafe T>
    //    T get_thread_local_instance(T &model) {
    //        thread_local T tl_instance;
    //        if constexpr (ThreadLocalCopy<T>) {
    //            tl_instance = model.thread_local_copy();
    //        } else if constexpr (ThreadLocalReference<T>) {
    //            tl_instance = model.thread_local_reference();
    //        } else {// ThreadLocalPointer<T>
    //            tl_instance = *model.thread_local_pointer();
    //        }
    //        return tl_instance;
    //    }

    template<ThreadLocalCopy T>
    T *get_thread_local_instance(T &model) {
        thread_local T tl_instance = model.thread_local_copy();
        return &tl_instance;
    }

    template<ThreadLocalReference T>
    T *get_thread_local_instance(T &model) {
        thread_local T &tl_instance = model.thread_local_reference();
        return &tl_instance;
    }

    template<ThreadLocalPointer T>
    T *get_thread_local_instance(T &model) {
        return model.thread_local_pointer();
    }


}// namespace odin

#endif//ODIN_PARALLEL_CONCEPTS_HPP
