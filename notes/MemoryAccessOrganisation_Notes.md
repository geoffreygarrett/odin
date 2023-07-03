Designing a data structure to organize memory based on access patterns at compile time can be challenging. In your case, where access patterns depend on the types and sequence of states given to the `CompositeState`, you're already setting yourself up for memory to be arranged in a way that reflects access patterns to a degree. The states that are together in the `CompositeState` will have their data next to each other in memory.

However, if you want to further optimize this, you would need to profile your application to understand the access patterns at runtime. This is because the order in which the states are accessed can depend on several factors including the specific computations being performed, the algorithm used, and the order of instructions issued by the CPU, which is generally not something you can predict accurately at compile time.

There is, however, a technique in template metaprogramming known as "sorting types". You can actually define a compile-time comparison function for types, then create a sorted typelist, and finally instantiate your `CompositeState` using this sorted typelist. This way, you are organizing the states based on their "importance" or "frequency of access" or any criteria you set. Here's an example:

```cpp
template <typename T, typename U>
struct compare_types {
    static constexpr bool value = T::state_dimension < U::state_dimension;
};

using SortedStates = sort_types<compare_types, Position, Velocity, Mass, MassRate>::type;

CompositeState<SortedStates...> state(position, velocity, mass, massRate);
```

In this example, `compare_types` is comparing the `state_dimension` of the types. So the `CompositeState` will be instantiated with the types sorted by their `state_dimension`.

Keep in mind that this technique requires a robust implementation of `sort_types`, and that your sorting criteria should reflect your application's specific access patterns for optimal performance. Note that reordering the types will not change the padding issue previously mentioned. And of course, runtime profiling is still necessary to validate any improvements you're hoping to achieve.

If your access patterns can't be known until runtime, then compile-time techniques won't be able to help, and you'll need to consider runtime techniques instead, such as rearranging data based on runtime profiling, using caches, or other data-oriented design strategies.



template<typename... States>
CompositeState::CompositeState(States... states)
: state_tuple(std::move(states)...) {

    // Calculate total state dimension and allocate global state vector
    global_state_vector = Eigen::VectorXd(state_dimension);

    // Helper lambda for setting up state views
    auto setup_state_view = [&](auto& state, auto offset) {
        new (&state.data) Eigen::Map<Eigen::VectorXd>(global_state_vector.data() + offset, state.data.size());
    };

    // Iterate over the state tuple and setup views into global state vector
    size_t offset = 0;
    std::apply([&](auto&... state) {
        (setup_state_view(state, offset), offset += state.data.size(), ...);
    }, state_tuple);
}
