# ODIN: Optimal Dynamics Integration and Numerics

ODIN is a high-performance, header-only library engineered for efficient
numerical computation in dynamics, simulation,
and estimation. With the power of modern C++17 and C++20 features, it
provides a flexible and easily extensible
interface while delivering lightning-fast computation. Its generic design
ensures adaptability to a vast array of
numerical problems, reaching far beyond its origin in astrodynamics and
ordinary differential equations (ODE)
simulations.

Crafted with a compile-time-first philosophy, ODIN optimizes for maximum
performance benefits. Recognizing the inherent
balance between performance and generality, it strikes a fine equilibrium
that preserves its capability to address a
broad spectrum of computational challenges. Furthermore, with a view
towards future compatibility, including planned
Python bindings via pybind11, ODIN incorporates optional dynamic runtime
features. These are carefully designed to offer
a flexible interface, meeting diverse user needs while keeping runtime
overhead to the bare minimum.

ODIN also features optional generic serialization using Cereal, ensuring no
computational result is lost due to system
failures or interruptions. It offers cached checkpoints that support
binary, portable binary, xml, and json formats for
all Eigen types and ODIN types.

## Key Features

1. **Performance-Driven**: Designed with performance at its core, ODIN aims
   to deliver swift and efficient computational
   solutions.
2. **Generality**: With versatile algorithms, ODIN can tackle a wide array
   of numerical problems.
3. **Modern C++ Design**: Harnessing the capabilities of modern C++17 and
   C++20 features, ODIN adheres to contemporary
   development standards, promoting clean, efficient, and maintainable
   code.
4. **Flexibility**: With its adaptable architecture, ODIN can be seamlessly
   extended and modified, serving as a dynamic
   tool for numerical simulations and estimations.
5. **Resilience**: ODIN employs optional generic serialization with Cereal
   to provide cached checkpoints, ensuring the
   safety and preservation of computational results.

## Roadmap

In its quest to offer the most efficient computational solutions, ODIN is
taking significant strides towards embracing
parallelization. A key focus area is the modularization of select
computations for offloading to the GPU, a feature that
can greatly accelerate processes like the calculation of gravitational
influence in complex models. Both Threading
Building Blocks (TBB) and CUDA are part of this ambitious roadmap.

## Getting Started

### Prerequisites

Ensure you have the following prerequisites installed on your system:

- A C++ compiler that supports C++17 and C++20. This project has been
  tested with GCC 9.3.0, Clang 11.0.0, and MSVC
  19.28.
- CMake 3.25 or newer.

### Building the Project

First, clone the repository:

```bash
git clone https://github.com/geoffreygarrett/odin.git
cd odin
````

To build the project, run the following commands:

```bash
mkdir build
cd build
cmake ..
make
```

We invite you to explore, use, and contribute to the ODIN library. Your
insights and contributions will help shape its
future, cementing its position as an essential tool in the realm of
high-performance numerical computation.

### CMake Options

| Option                   | Description                                      | Default Value |
|--------------------------|--------------------------------------------------|---------------|
| `ODIN_USE_SERIALIZATION` | Enable serialization with cereal                 | `ON`          |
| `ODIN_BUILD_TESTS      ` | Build test programs                              | `ON`          |
| `ODIN_BUILD_EXAMPLES   ` | Build example programs                           | `ON`          |
| `ODIN_BUILD_BENCHMARKS ` | Build benchmark programs                         | `OFF`         |
| `ODIN_USE_AUTODIFF     ` | Enable autodiff                                  | `ON`          |
| `ODIN_USE_GLOG         ` | Enable Google's logging library glog             | `ON`          |
| `ODIN_USE_FFTW         ` | Enable the Fastest Fourier Transform in the West | `ON`          |
| `ODIN_USE_SHTNS        ` | Enable the Spherical Harmonic Transform library  | `ON`          |
| `ODIN_USE_TBB          ` | Enable Intel's Threading Building Blocks         | `ON`          |
| `ODIN_USE_GSL          ` | Enable the GNU Scientific Library                | `ON`          |

## Autodiff Options

If `ODIN_USE_AUTODIFF` is `ON`, you can use the following options:

| Option                    | Description                                       | Default Value |
|---------------------------|:--------------------------------------------------|---------------|
| `AUTODIFF_BUILD_TESTS   ` | Enable the compilation of the test files          | `OFF`         |
| `AUTODIFF_BUILD_PYTHON  ` | Enable the compilation of the python bindings     | `OFF`         |
| `AUTODIFF_BUILD_EXAMPLES` | Enable the compilation of the example files       | `OFF`         |
| `AUTODIFF_BUILD_DOCS    ` | Enable the build of the documentation and website | `OFF`         |

## Cereal Options

If `ODIN_USE_SERIALIZATION` is `ON`, you can use the following options:

| Option                               | Description                                  | Default Value |
|--------------------------------------|----------------------------------------------|---------------|
| `CEREAL_JUST_INSTALL_CEREAL        ` | Only install cereal, do nothing else         | `ON`          |
| `CEREAL_SKIP_PORTABILITY_TEST      ` | Skip portability tests                       | `ON`          |
| `CEREAL_THREAD_SAFE                ` | Enable thread safety                         | `OFF`         |
| `CEREAL_BUILD_DOC                  ` | Build documentation                          | `ON`          |
| `CEREAL_BUILD_SANDBOX              ` | Build sandbox examples                       | `ON`          |
| `CEREAL_SKIP_PERFORMANCE_COMPARISON` | Skip building performance sandbox comparison | `OFF`         |
| `CEREAL_WITH_WERROR                ` | Compile with '-Werror' C++ compiler flag     | `ON`          |
| `CEREAL_CLANG_USE_LIBCPP           ` | Use libc++ for clang compilation             | `OFF`         |
| `CEREAL_CEREAL_INSTALL             ` | Generate the install target                  | `ON`          |
| `CEREAL_BUILD_TESTS                ` | Build tests                                  | `OFF`         |

### Dependencies

If you build anything dependent on pagmo, you'll need to install gfortran.

```bash
sudo apt-get install gfortran
```