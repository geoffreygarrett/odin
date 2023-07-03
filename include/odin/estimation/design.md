Ok paint a picture for me. Demonstrate a consistent API across the entire script. I'm envisaging two this as the start of the script to set the initial state and what states are being tracked in the ProcessStateHandler

```C++
using process_state = StateHandler<
        Position<EARTH>, 
        Velocity<EARTH>, 
        Position<DELFI>, 
        Velocity<DELFI>>;

Position<EARTH> earth_pos{10.0, 20.0, 30.0};
Position<DELFI> delfi_pos{5.0, -2.0, 7.0};
Velocity<EARTH> earth_vel{1.0, 2.0, 3.0};
Velocity<DELFI> delfi_vel{-1.0, 0.5, -0.3};
process_state::state_tuple global_states = std::make_tuple(
        earth_pos,
        earth_vel,
        delfi_pos, 
        delfi_vel);

// Example: Access and modify state using StateHandler
 ProcessState::get<Position<EARTH >>(global_states).value.x() += 10.0;
```

Then I imagine we set up the ProcessModel itself, something like this:

```C++
using process = ProcessHandler<process_state, process_state::state_tuple>;

constexpr double EARTH_SH_FILE = "earth_gravity_model.txt";
constexpr double SUN_MU = 1.32712440018e20;
constexpr double MOON_MU = 4.9028e12;

SphericalHarmonicsGravitation<EARTH, config::USE_GPU> earth_point_mass_gravity{EARTH_SH_FILE};
PointMassGravitation<SUN> sun_point_mass_gravity{SUN_MU};
PointMassGravitation<MOON> moon_point_mass_gravity{MOON_MU};

// I'm not certain what the most generic possible way is to do this,
// but we then need to map the model to the states/entity that it affects.
// ...

// Then we can set up the process model
process::model_tuple global_models = std::make_tuple(
        earth_point_mass_gravity,
        sun_point_mass_gravity,
        moon_point_mass_gravity);

// Example: Print the model equations
for (auto &model : process::model_tuple) {
    std::cout << model.equation << std::endl;
}
/*
Possible Output:
    "F = ......" // Equation for the SphericalHarmonicsGravitation model
    "F = (MU: parameter) / ((r: {relative_state}) ^ 2)" // Equation for the PointMassGravitation model (SUN)
    "F = (MU: parameter) / ((r: {relative_state}) ^ 2)" // Equation for the PointMassGravitation model (MOON)
*/

// Example: Accessing the parameters of the model
for (auto &parameter : process::parameter_tuple) {
    std::cout << parameter.name << std::endl;
}

/*
Possible Output:
    "EARTH_SH_FILE_C_0"
    "EARTH_SH_FILE_S_0"
    "EARTH_SH_FILE_C_1"
    "EARTH_SH_FILE_S_1"
    //...
    "SUN_MU"
    "MOON_MU"
*/


```

Then we can set up the measurement model, something like this:

```C++
using measurement_model = MeasurementHandler<
        RangeMeasurement<DELFI, EARTH>,
        RangeMeasurement<DELFI, MOON>>;

constexpr double NOISE_STD = 1.0;

RangeMeasurement<DELFI, EARTH> delfi_to_earth_range{NOISE_STD};
RangeMeasurement<DELFI, MOON> delfi_to_moon_range{NOISE_STD};

measurement_model::measurement_tuple global_measurements = std::make_tuple(
        delfi_to_earth_measurement,
        delfi_to_moon_measurement);
```
Estimated state and covariance are then updated as follows:

```C++
using estimated_state = StateHandler<
        Position<DELFI>,
        Velocity<DELFI>>;

using estimated_covariance = CovarianceHandler<
        Covariance<Position<DELFI>>,
        Covariance<Velocity<DELFI>>>;

Position<DELFI> delfi_estimated_pos{5.0, -2.0, 7.0};
Velocity<DELFI> delfi_estimated_vel{-1.0, 0.5, -0.3};

estimated_state::state_tuple estimated_states = std::make_tuple(
        delfi_estimated_pos,
        delfi_estimated_vel);

// Another initializer will be provided to fully specify all cross-covariances,
// if desired later on.
Covariance<Position<DELFI>> delfi_pos_covar{0.1, 0.1, 0.1};
Covariance<Velocity<DELFI>> delfi_vel_covar{0.01, 0.01, 0.01};

estimated_covariance::covariance_tuple estimated_covariances = std::make_tuple(
        delfi_pos_covar,
        delfi_vel_covar);
```

```C++
/**
 * \typedef estimator
 * 
 * \brief Typedef for an Estimator object.
 * 
 * The Estimator is defined with three template parameters:
 * - The estimated state tuple: estimated_state::state_tuple,
 * - The estimated covariance tuple: estimated_covariance::covariance_tuple,
 * - The augmentation option: augment::NONE (as placeholder; to be implemented as needed).
 * 
 * Several examples of possible augmentation methods are given below (note these are purely illustrative):
 * - augment::STATE_WITH_NOISE: Implements non-additive noise augmentation.
 * - augment::STATE_DELAYED<2>: Implements time-delayed state augmentation with a delay of 2 time steps.
 * - augment::PARAMETER<PointMassGravitation<SUN>>: Augments the state with the uncertain parameter of the PointMassGravitation<SUN> model.
 * - augment::MULTI_MODEL<estimator::extended_kalman_filter<1>,
 *                        estimator::extended_kalman_filter<2>,
 *                        estimator::extended_kalman_filter<3>>: For systems where multiple models might be suitable.
 * - augment::AUXILIARY_VARIABLE<...>: Augments the state with auxiliary variables that could simplify system representation.
 */
using estimator = Estimator<
    estimated_state::state_tuple,
    estimated_covariance::covariance_tuple,
    augment::NONE
>;

estimator estimationHandler;

estimationHandler.stateHandler = estimated_states;
estimationHandler.covarianceHandler = estimated_covariances;

// Initialize the estimator with the initial estimated state and covariance
estimationHandler.initialize(estimated_states, estimated_covariances);
```


And a very rough idea of what the final loop might look like:

```C++

using control_state = ControlHandler<
    ControlCommand<DELFI>
>;

// Define the control command with thrust magnitude and direction
ControlCommand<DELFI> delfi_thrust{10.0, {0.0, 1.0, 0.0}};

control_state::control_tuple global_controls = std::make_tuple(
    delfi_thrust
);

// Update control command using ControlHandler
control_state::get<ControlCommand<DELFI>>(global_controls).thrust = 20.0;

int main() {
    // Define the simulation time step (in seconds)
    const double timeStep = 1.0;

    // Define the total simulation duration (in seconds)
    const double duration = 100.0;

    for (double time = 0; time < duration; time += timeStep) {
        // Compute the new state using the process model
        global_states = process::update(global_states, global_models, global_controls, timeStep);

        // Generate new measurements
        global_measurements = measurement_model::measure(global_states);

        // Update the estimated state and covariance
        std::tie(estimated_states, estimated_covariances) = 
            estimationHandler.estimate(global_states, global_measurements);

        // Print the estimated position and velocity of DELFI
        std::cout << "Estimated position of DELFI at time " << time << ": " 
                  << estimated_state::get<Position<DELFI>>(estimated_states).value << std::endl;
        std::cout << "Estimated velocity of DELFI at time " << time << ": " 
                  << estimated_state::get<Velocity<DELFI>>(estimated_states).value << std::endl;

        // TODO: Implement a control strategy for modifying the global_controls based on the estimated state
        // and your mission objectives.
    }

    return 0;
}
```



```angular2html
#include <Eigen/Dense>
#include <type_traits>

template<typename Scalar, int x_dim, int v_dim, int n_dim>
struct AugmentedState {
    Eigen::Matrix<Scalar, x_dim, 1> x;
    Eigen::Matrix<Scalar, v_dim, 1> v;
    Eigen::Matrix<Scalar, n_dim, 1> n;
};

template<typename Scalar, int x_dim, int v_dim, int n_dim>
struct AugmentedCovariance {
    Eigen::Matrix<Scalar, x_dim, x_dim> P0;
    Eigen::Matrix<Scalar, v_dim, v_dim> RV;
    Eigen::Matrix<Scalar, n_dim, n_dim> RN;
};

// Base template for state transition
template<typename Scalar, int x_dim, int v_dim, int n_dim>
AugmentedState<Scalar, x_dim, v_dim, n_dim>
state_transition(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x, double dt) {
    AugmentedState<Scalar, x_dim, v_dim, n_dim> new_state;
    new_state.x = rk4.step(x.x, dt);  
    new_state.v = process_noise(x.x, dt);
    new_state.n = measurement_noise(x.x);
    return new_state;
}

// Function for process_covariance if RV and RN are variable
template<typename Scalar, int x_dim, int v_dim, int n_dim>
typename std::enable_if<(v_dim != Eigen::Dynamic && n_dim != Eigen::Dynamic), 
    AugmentedCovariance<Scalar, x_dim, v_dim, n_dim>>::type
process_covariance(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x, double dt) {
    AugmentedCovariance<Scalar, x_dim, v_dim, n_dim> Q;
    Q.P0 = original_process_covariance(x.x, dt);  
    Q.RV = calculate_process_noise_covariance(x.x, dt);
    Q.RN = calculate_measurement_noise_covariance(x.x, dt);
    return Q;
}

// Function for process_covariance if RV and RN are constant
template<typename Scalar, int x_dim, int v_dim, int n_dim>
typename std::enable_if<(v_dim == Eigen::Dynamic && n_dim == Eigen::Dynamic), 
    AugmentedCovariance<Scalar, x_dim, v_dim, n_dim>>::type
process_covariance(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x, double dt) {
    AugmentedCovariance<Scalar, x_dim, v_dim, n_dim> Q;
    Q.P0 = original_process_covariance(x.x, dt);  
    Q.RV = Eigen::Matrix<Scalar, v_dim, v_dim>::Identity();  // Or any other constant value
    Q.RN = Eigen::Matrix<Scalar, n_dim, n_dim>::Identity();  // Or any other constant value
    return Q;
}

template<typename Scalar, int x_dim, int v_dim, int n_dim>
Eigen::Matrix<Scalar, n_dim, 1>
measurement_model(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x) {
    return original_measurement_model(x.x);  
}

```


```C++
template<typename Scalar, int Dim>
struct Augmentation {
    using State = Eigen::Matrix<Scalar, Dim, 1>;
    using Covariance = Eigen::Matrix<Scalar, Dim, Dim>;

    std::pair<State, Covariance> data;
};

// Noise augmentation strategy
template<typename Scalar, int x_dim, int v_dim, int n_dim>
struct NoiseAugmentation {
    using StateAugment = Augmentation<Scalar, x_dim>;
    using VNoiseAugment = Augmentation<Scalar, v_dim>;
    using NNoiseAugment = Augmentation<Scalar, n_dim>;

    StateAugment state_augmentation;
    VNoiseAugment v_noise_augmentation;
    NNoiseAugment n_noise_augmentation;
};

// Parameter augmentation strategy
template<typename Scalar, int param_dim>
struct ParameterAugmentation {
    using ParameterAugment = Augmentation<Scalar, param_dim>;
    ParameterAugment parameter_augmentation;
};

template<typename Scalar, typename... Augmentations>
class AugmentationHandler;

// Base case for AugmentationHandler
template<typename Scalar>
class AugmentationHandler<Scalar> {
    // base behavior here
};

template<typename Scalar, int x_dim, int v_dim, int n_dim>
class AugmentationHandler {
public:
    AugmentationHandler() = default;
    
    AugmentedState<Scalar, x_dim, v_dim, n_dim>
    state_transition(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x, double dt) {
        // Implement the state_transition logic
    }
    
    AugmentedCovariance<Scalar, x_dim, v_dim, n_dim>
    process_covariance(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x, double dt) {
        // Implement the process_covariance logic
    }

    Eigen::Matrix<Scalar, n_dim, 1>
    measurement_model(const AugmentedState<Scalar, x_dim, v_dim, n_dim>& x) {
        // Implement the measurement_model logic
    }

    AugmentedState<Scalar, x_dim, v_dim, n_dim>
    construct_augmented_state() {
        return AugmentedState<Scalar, x_dim, v_dim, n_dim>();
    }

    AugmentedCovariance<Scalar, x_dim, v_dim, n_dim>
    construct_augmented_covariance() {
        return AugmentedCovariance<Scalar, x_dim, v_dim, n_dim>();
    }
};
```

Great, it's coming along nicely. Now, let's consider next steps and potential improvements.

1. **Generalize Noise and Parameter Augmentation**: Rather than hardcoding the `NoiseAugmentation`
   and `ParameterAugmentation` classes, we could create a more general `AugmentationStrategy` class that takes in a list
   of dimensions for the state and covariance. This way, we can create different augmentation strategies just by
   providing different dimension values.

2. **Improve Method Overloading**: The methods `state_transition`, `process_covariance`, and `measurement_model` in
   the `AugmentationHandler` class could be improved. Right now, they're hardcoded for specific dimensions. We could use
   template specialization or function overloading to create versions of these methods for different types
   of `AugmentedState` and `AugmentedCovariance`.

3. **Consider Dynamic Size Matrices**: Currently, we're using static-sized Eigen matrices. If the dimensions aren't
   known at compile time, we could switch to dynamic-sized matrices by using `Eigen::Dynamic`.

4. **Refactor Construct Methods**: The methods `construct_augmented_state` and `construct_augmented_covariance` are
   returning default constructed objects. We should consider refactoring these methods to initialize the objects with
   meaningful values.

5. **Introduce Compile-Time Control Flow**: Depending on the nature of our augmentation strategy, we may want to have
   compile-time control flow using `constexpr`, SFINAE, or template specialisation.

6. **Introduce Error Handling**: The code currently has no form of error handling or validation. We might want to add
   some checks to ensure the code fails gracefully when given invalid input.

7. **Add More Augmentation Strategies**: In the future, we could add more augmentation strategies such as time-delayed
   state augmentation, multi-model systems augmentation, and auxiliary variables augmentation.

8. **Test and Document**: Finally, once the major features are implemented, we'll want to write thorough unit tests to
   ensure everything is working as expected. Properly commenting the code and documenting its functionality and usage
   would also be important.

Note: Always remember to use the coding style guide (snake_case, camelCase, etc.) that is preferred in your development
environment.


Sure, we can refactor the base `Augmentation` struct to follow the Curiously Recurring Template Pattern (CRTP) paradigm, and use static polymorphism to avoid the overhead of dynamic polymorphism. This approach also allows us to keep the logic of the state transition, process covariance, and measurement model constant for all augmentations. Here is the refactored code:

```cpp
// CRTP base class
template<typename Derived, typename Scalar, int Dim>
struct Augmentation {
    using State = Eigen::Matrix<Scalar, Dim, 1>;
    using Covariance = Eigen::Matrix<Scalar, Dim, Dim>;

    std::pair<State, Covariance> data;

    // Call the state_transition function on the derived class
    State state_transition(const State& x, double dt) {
        return static_cast<Derived*>(this)->state_transition_impl(x, dt);
    }

    // Call the process_covariance function on the derived class
    Covariance process_covariance(const State& x, double dt) {
        return static_cast<Derived*>(this)->process_covariance_impl(x, dt);
    }

    // Call the measurement_model function on the derived class
    State measurement_model(const State& x) {
        return static_cast<Derived*>(this)->measurement_model_impl(x);
    }
};

template<typename Scalar, int x_dim, int v_dim, int n_dim>
struct NoiseAugmentation : Augmentation<NoiseAugmentation<Scalar, x_dim, v_dim, n_dim>, Scalar, x_dim> {
    using Base = Augmentation<NoiseAugmentation<Scalar, x_dim, v_dim, n_dim>, Scalar, x_dim>;
    using State = typename Base::State;
    using Covariance = typename Base::Covariance;

    State state_transition_impl(const State& x, double dt) {
        // Implement the state transition logic here for NoiseAugmentation
    }

    Covariance process_covariance_impl(const State& x, double dt) {
        // Implement the process covariance logic here for NoiseAugmentation
    }

    State measurement_model_impl(const State& x) {
        // Implement the measurement model logic here for NoiseAugmentation
    }
};

template<typename Scalar, int param_dim>
struct ParameterAugmentation : Augmentation<ParameterAugmentation<Scalar, param_dim>, Scalar, param_dim> {
    using Base = Augmentation<ParameterAugmentation<Scalar, param_dim>, Scalar, param_dim>;
    using State = typename Base::State;
    using Covariance = typename Base::Covariance;

    State state_transition_impl(const State& x, double dt) {
        // Implement the state transition logic here for ParameterAugmentation
    }

    Covariance process_covariance_impl(const State& x, double dt) {
        // Implement the process covariance logic here for ParameterAugmentation
    }

    State measurement_model_impl(const State& x) {
        // Implement the measurement model logic here for ParameterAugmentation
    }
};
```
In this version of the code, each augmentation strategy provides its own implementations for `state_transition_impl`, `process_covariance_impl`, and `measurement_model_impl`, which are called from the base `Augmentation` class. This allows us to avoid virtual function overhead while still providing different behaviour for each augmentation strategy.



Given the new definition for the `Augmentation` using CRTP, we can update the `AugmentationHandler` accordingly. The handler needs to know about the different augmentations and to be able to deal with multiple different ones at a time, hence the variadic template.

```cpp
template<typename Scalar, typename... Augmentations>
class AugmentationHandler;

// Base case for AugmentationHandler
template<typename Scalar>
class AugmentationHandler<Scalar> {
    // base behavior here
};

template<typename Scalar, typename... Augmentations>
class AugmentationHandler {
public:
    AugmentationHandler() = default;
    
    std::tuple<typename Augmentations::State...> 
    state_transition(const std::tuple<typename Augmentations::State...>& x, double dt) {
        return std::make_tuple((Augmentations().state_transition(std::get<typename Augmentations::State>(x), dt))...);
    }

    std::tuple<typename Augmentations::Covariance...> 
    process_covariance(const std::tuple<typename Augmentations::State...>& x, double dt) {
        return std::make_tuple((Augmentations().process_covariance(std::get<typename Augmentations::State>(x), dt))...);
    }

    std::tuple<typename Augmentations::State...> 
    measurement_model(const std::tuple<typename Augmentations::State...>& x) {
        return std::make_tuple((Augmentations().measurement_model(std::get<typename Augmentations::State>(x)))...);
    }

    std::tuple<typename Augmentations::State...> 
    construct_augmented_state() {
        return std::make_tuple(Augmentations().data.first...);
    }

    std::tuple<typename Augmentations::Covariance...> 
    construct_augmented_covariance() {
        return std::make_tuple(Augmentations().data.second...);
    }
};
```

Here we are using tuples to handle the variadic nature of our `AugmentationHandler`. Each augmentation type has its own
corresponding state and covariance, so these are stored in tuples of their respective types.

Each method in the handler uses variadic template expansion to perform the necessary operations for each augmentation
type. The state transition, process covariance, and measurement model methods are defined in terms of the corresponding
methods in the `Augmentation` base class. The `construct_augmented_state` and `construct_augmented_covariance` methods
construct a tuple of states and covariances, respectively, for each augmentation type.

Please note that you would need to provide correct initialization for each augmentation type. The above example uses
default constructor of augmentations to demonstrate how methods can be used. In a real use-case scenario, you might need
to initialize each augmentation object properly before using them.



-----

```C++
// CRTP base class
template<typename Derived, typename Scalar, int Dim>
struct Augmentation {
using State = Eigen::Matrix<Scalar, Dim, 1>;
using Covariance = Eigen::Matrix<Scalar, Dim, Dim>;

    std::pair<State, Covariance> data;

    // Call the state_transition function on the derived class
    State state_transition(const State& x, double dt) {
        return static_cast<Derived*>(this)->state_transition_impl(x, dt);
    }

    // Call the process_covariance function on the derived class
    Covariance process_covariance(const State& x, double dt) {
        return static_cast<Derived*>(this)->process_covariance_impl(x, dt);
    }

    // Call the measurement_model function on the derived class
    State measurement_model(const State& x) {
        return static_cast<Derived*>(this)->measurement_model_impl(x);
    }
};

template<typename Scalar, int x_dim, int v_dim, int n_dim>
struct NoiseAugmentation : Augmentation<NoiseAugmentation<Scalar, x_dim, v_dim, n_dim>, Scalar, x_dim> {
using Base = Augmentation<NoiseAugmentation<Scalar, x_dim, v_dim, n_dim>, Scalar, x_dim>;
using State = typename Base::State;
using Covariance = typename Base::Covariance;

    State state_transition_impl(const State& x, double dt) {
        // Implement the state transition logic here for NoiseAugmentation
    }

    Covariance process_covariance_impl(const State& x, double dt) {
        // Implement the process covariance logic here for NoiseAugmentation
    }

    State measurement_model_impl(const State& x) {
        // Implement the measurement model logic here for NoiseAugmentation
    }
};

template<typename Scalar, int param_dim>
struct ParameterAugmentation : Augmentation<ParameterAugmentation<Scalar, param_dim>, Scalar, param_dim> {
using Base = Augmentation<ParameterAugmentation<Scalar, param_dim>, Scalar, param_dim>;
using State = typename Base::State;
using Covariance = typename Base::Covariance;

    State state_transition_impl(const State& x, double dt) {
        // Implement the state transition logic here for ParameterAugmentation
    }

    Covariance process_covariance_impl(const State& x, double dt) {
        // Implement the process covariance logic here for ParameterAugmentation
    }

    State measurement_model_impl(const State& x) {
        // Implement the measurement model logic here for ParameterAugmentation
    }
};


```

PERSONAL NOTE:
> We say that a matrix is stored in row-major order if it is stored row by row. The entire first row is stored first, followed by the entire second row, and so on. Consider for example the matrix
> A=⎡⎣⎢893215244945⎤⎦⎥.
If this matrix is stored in row-major order, then the entries are laid out in memory as follows:
8 2 2 9 9 1 4 4 3 5 4 5
On the other hand, a matrix is stored in column-major order if it is stored column by column, starting with the entire first column, followed by the entire second column, and so on. If the above matrix is stored in column-major order, it is laid out as follows:
> https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
> 
> 
> 


```angular2html

    /// --- experimental ---
    /// This method directly calculates the total state dimension size by folding
    /// the state_dimension values of each State in the variadic template parameter pack States.
    /// The folding is performed by summing up the values (i.e., (+ ...)).
    /// This operation is performed at compile-time.
    /// Syntax is more concise and doesn't require a helper structure.
    static constexpr std::size_t state_dimension = (... + States::state_dimension);

    /// Get total size of the state vector
    /// This method achieves the same goal as the above method but using an auxiliary structure `StateSize`.
    /// The purpose of `StateSize` is to calculate the sum of the state_dimension of each State
    /// in the variadic template parameter pack States.
    /// It also performs this operation at compile-time.
    /// This approach makes use of an additional helper structure, which 
    /// might make it more reusable and modular,
    /// but it slightly increases the complexity of the code.
    static constexpr size_t total_state_dim() {
        return (... + States::state_dimension);
    }

    template<typename... States>
    struct StateSize {
        static constexpr int value = (... + States::state_dimension);
    };

    /// --- experimental ---
```