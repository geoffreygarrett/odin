1. **Static Assertions**: These are assertions that are evaluated at compile time. If the condition does not hold, it
   causes a compile-time error.

   Example:
    ```cpp
    static_assert(std::is_final<Derived>::value, "Derived must be a final class!");
    ```

2. **Using `std::integral_constant` for Compile-Time Constants**: This is a way to represent a compile-time constant of
   a specified type.

   Example:
    ```cpp
    using derivative_order = std::integral_constant<int, Order>;
    ```

3. **Compile-Time Arithmetic Operations**: This refers to performing arithmetic operations during compile-time. The
   result of these operations can be used in other compile-time computations or checks.

   Example:
    ```cpp
    static constexpr int count_true_types() {
        return ((Derived::is_translational::value ? 1 : 0) +
                (Derived::is_rotational::value ? 1 : 0) +
                (Derived::is_inertia::value ? 1 : 0) +
                (Derived::is_mass::value ? 1 : 0));
    }
    ```

4. **Checking a Condition for a Set of Types with `std::disjunction`**: This is used to check if a particular condition
   holds for at least one type in a set of types. `std::disjunction` takes a set of `std::integral_constant` types and
   yields `std::true_type` if at least one of them is `std::true_type`.

   Example:
    ```cpp
    static_assert(std::disjunction<is_translational, is_rotational, is_inertia, is_mass>::value,
                  "At least one of the motion types must be true!");
    ```

5. **Checking a Type Trait for a Set of Types with `std::conjunction`**: This is used to check if a particular type
   trait holds for all types in a set of types. `std::conjunction` takes a set of `std::integral_constant` types and
   yields `std::true_type` if all of them are `std::true_type`.

   Example:
    ```cpp
    static_assert(std::conjunction<std::is_integral<T1>, std::is_integral<T2>, std::is_integral<T3>>::value,
                  "All types must be integral!");
    ```

6. **Checking if a Class is Derived from a Base Class with `std::is_base_of`**: This checks if a class is derived from a
   particular base class.

   Example:
    ```cpp
    static_assert(std::is_base_of<BaseClass, DerivedClass>::value, "DerivedClass must be derived from BaseClass!");
    ```

Sure, let's continue with more checks you can use in your C++ code.

7. **Checking Type Equality with `std::is_same`**: This can be used to check if two types are exactly the same.

   Example:
    ```cpp
    static_assert(std::is_same<T1, T2>::value, "T1 and T2 must be the same type!");
    ```

8. **Checking if a Type is a Particular Category with `std::is_` traits**: There are many `std::is_` traits that can be
   used to check if a type is a specific category. Some examples
   include `std::is_integral`, `std::is_floating_point`, `std::is_array`, `std::is_enum`, `std::is_pointer`, etc.

   Example:
    ```cpp
    static_assert(std::is_integral<T>::value, "T must be an integral type!");
    ```

9. **Checking if a Type has a Specific Property with `std::has_` traits**: There are many `std::has_` traits that can be
   used to check if a type has a specific property. Some examples
   include `std::has_virtual_destructor`, `std::has_unique_object_representations`.

   Example:
    ```cpp
    static_assert(std::has_virtual_destructor<T>::value, "T must have a virtual destructor!");
    ```

10. **Checking Type Conversion Possibilities with `std::is_convertible`**: This can be used to check if a type can be
    converted to another type.

    Example:
    ```cpp
    static_assert(std::is_convertible<T1, T2>::value, "T1 must be convertible to T2!");
    ```

11. **Checking if a Type Satisfies a Combination of Traits with `std::conditional`**: This can be used to select one of
    two types based on a condition.

    Example:
    ```cpp
    using type = std::conditional<std::is_integral<T>::value, int, float>::type;
    ```

12. **Checking for Const-Qualification with `std::is_const`**: This can be used to check if a type is const-qualified.

    Example:
    ```cpp
    static_assert(std::is_const<T>::value, "T must be const-qualified!");
    ```

13. **Checking for Specific Type Modifiers with other `std::is_` traits**: This includes traits
    like `std::is_volatile`, `std::is_signed`, `std::is_unsigned`.

    Example:
    ```cpp
    static_assert(std::is_volatile<T>::value, "T must be volatile-qualified!");
    ```

14. **Checking Size of a Type with `std::is_empty` and `sizeof`**: This can be used to check if a type is an empty
    class (i.e., has no non-static data members) or to compare the sizes of types.

    Example:
    ```cpp
    static_assert(std::is_empty<T>::value, "T must be an empty class!");
    static_assert(sizeof(T1) == sizeof(T2), "T1 and T2 must be the same size!");
    ```

15. **Checking Alignment of a Type with `alignof`**: This can be used to check if a type has a specific alignment.

    Example:
    ```cpp
    static_assert(alignof(T) == 8, "T must have an alignment of 8!");
    ```

These are not all the checks possible, but it should cover a large number of the common cases that you might encounter
when working with templates and type checking in C++.