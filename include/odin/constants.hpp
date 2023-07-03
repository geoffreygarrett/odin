#include <cmath>

namespace constants {


    template<typename Scalar>
    constexpr Scalar PI = static_cast<Scalar>(M_PIl);

    template<>
    constexpr float PI<float> = M_PIf;

    template<>
    constexpr double PI<double> = M_PI;

    template<>
    constexpr long double PI<long double> = M_PIl;

    template<typename Scalar>
    constexpr Scalar PI_2 = PI<Scalar> / static_cast<Scalar>(2.0);

    template<typename Scalar>
    constexpr Scalar PI_4 = PI<Scalar> / static_cast<Scalar>(4.0);

    template<>
    constexpr float PI_2<float> = M_PI_2f;

    template<>
    constexpr float PI_4<float> = M_PI_4f;

    template<>
    constexpr double PI_2<double> = M_PI_2;

    template<>
    constexpr double PI_4<double> = M_PI_4;


    template<typename Scalar>
    constexpr static Scalar G = static_cast<Scalar>(6.67430e-11);// Gravitational constant

    template<typename Scalar>
    constexpr static Scalar c = static_cast<Scalar>(2.998e8);// Speed of light in m/s

    template<typename Scalar>
    constexpr static Scalar k = static_cast<Scalar>(1.380649e-23);// Boltzmann constant in J/K

    template<typename Scalar>
    constexpr static Scalar h = static_cast<Scalar>(6.62607015e-34);// Planck constant in J.s

    template<typename Scalar>
    constexpr static Scalar e = static_cast<Scalar>(1.602176634e-19);// Elementary charge in C

    template<typename Scalar>
    constexpr static Scalar me = static_cast<Scalar>(9.10938356e-31);// Electron mass in kg

    template<typename Scalar>
    constexpr static Scalar mp = static_cast<Scalar>(1.672621898e-27);// Proton mass in kg

    template<typename Scalar>
    constexpr static Scalar mu_0 = static_cast<Scalar>(4.0 * M_PI * 1.0e-7);// Vacuum permeability in H/m

    template<typename Scalar>
    constexpr static Scalar epsilon_0 = static_cast<Scalar>(1.0 / (mu_0<Scalar> * c<Scalar> * c<Scalar>) );// Vacuum permittivity in F/m

    template<typename Scalar>
    constexpr static Scalar M_Earth = static_cast<Scalar>(5.97219e24);// Mass of Earth in kg

    template<typename Scalar>
    constexpr static Scalar R_Earth = static_cast<Scalar>(6.3781e6);// Radius of Earth in m

    template<typename Scalar>
    constexpr static Scalar AU = static_cast<Scalar>(1.495978707e11);// Astronomical Unit in m

    // Conversion factors
    template<typename Scalar>
    constexpr static Scalar deg_to_rad = static_cast<Scalar>(M_PI / 180.0);// Degree to Radian

    template<typename Scalar>
    constexpr static Scalar rad_to_deg = static_cast<Scalar>(180.0 / M_PI);// Radian to Degree

    template<typename Scalar>
    constexpr static Scalar E = static_cast<Scalar>(M_E);// Euler's number

    // Note: M_PI and M_E are defined in math.h header file in C++

}// namespace constants