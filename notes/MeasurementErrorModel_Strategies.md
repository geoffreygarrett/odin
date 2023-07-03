# MeasurementErrorModel_Strategies.md

In C++ generic programming, various strategies can be used to customize the behavior of a base class template. This note presents three different approaches, using a `MeasurementErrorModel` class as an example.

## Single Parameter Instance

In this strategy, we make a method in the base class that accepts a single Parameter instance:

```cpp
virtual scalar_type compute_error(const Parameter<ID, Scalar>& param) const {
    return static_cast<const Derived*>(this)->compute_error_impl(param);
}
```

This is the simplest and most straightforward approach. It provides strong type checking and simplicity of design. However, it lacks flexibility as each MeasurementErrorModel is locked to one parameter type. It's best suited when you're sure that each MeasurementErrorModel will always deal with only one parameter type.

## ParameterHandler with Specific Parameters

The second strategy is more flexible. It takes a `ParameterHandler` instance and a variadic list of parameters:

```cpp
template<typename... Params>
scalar_type compute_error(ParameterHandler<Params...>& handler, const Params&... parameters) const {
    return static_cast<const Derived*>(this)->compute_error_impl(handler.template get<parameters>()...);
}
```

This approach allows for more flexibility while maintaining explicit type checking. This method can deal with different parameter types in different contexts, which makes it more adaptable to evolving requirements.

## ParameterHandler with Derived Class Specifying Needs

The third strategy is the most flexible but also the most complex. The derived class specifies the parameters it needs:

```cpp
template<typename... Parameters>
scalar_type compute_error(ParameterHandler<Parameters...>& handler) const {
    using expected_parameter_types = typename Derived::expected_parameter_types;
    auto params = handler.template get_parameters_of_entity<expected_parameter_types>(handler.parameters);
    return static_cast<const Derived*>(this)->compute_error_impl(params);
}
```

In this case, the derived class dictates its requirements. This method is useful when you want to be able to add new MeasurementErrorModel subclasses that operate on new types of parameters without modifying the base class or the `ParameterHandler` class. This method reduces the amount of boilerplate code in each derived class but can be harder to understand due to advanced template techniques.

## Conclusion

These strategies demonstrate different ways to approach generic programming in C++. The choice of strategy often depends on the specific needs of the project, the complexity of the requirements, and the expected future changes to those requirements. Each of these strategies provides different trade-offs between simplicity, flexibility, and complexity. It's always important to consider these factors when choosing the most appropriate strategy for your project.
