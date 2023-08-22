#ifndef ODIN_CONSTANTS_H
#define ODIN_CONSTANTS_H

#include <cmath>

namespace odin::consts {

template<typename Scalar = double>
inline constexpr Scalar PI = static_cast<Scalar>(M_PIl);

template<typename Scalar = double>
inline constexpr Scalar PI_2 = PI<Scalar> / 2.0;

template<typename Scalar = double>
inline constexpr Scalar PI_4 = PI<Scalar> / 4.0;

template<typename Scalar = double>
inline constexpr Scalar G = 6.67430e-11;// Gravitational constant

template<typename Scalar = double>
inline constexpr Scalar c = 2.998e8;// Speed of light in m/s

template<typename Scalar = double>
inline constexpr Scalar k = 1.380649e-23;// Boltzmann constant in J/K

template<typename Scalar = double>
inline constexpr Scalar h = 6.62607015e-34;// Planck constant in J.s

template<typename Scalar = double>
inline constexpr Scalar e = 1.602176634e-19;// Elementary charge in C

template<typename Scalar = double>
inline constexpr Scalar me = 9.10938356e-31;// Electron mass in kg

template<typename Scalar = double>
inline constexpr Scalar mp = 1.672621898e-27;// Proton mass in kg

template<typename Scalar = double>
inline constexpr Scalar mu_0 = 4.0 * PI<Scalar> * 1.0e-7;// Vacuum permeability in H/m

template<typename Scalar = double>
inline constexpr Scalar epsilon_0
        = 1.0 / (mu_0<Scalar> * c<Scalar> * c<Scalar>);// Vacuum permittivity in F/m

template<typename Scalar = double>
inline constexpr Scalar M_Earth = 5.97219e24;// Mass of Earth in kg

template<typename Scalar = double>
inline constexpr Scalar R_Earth = 6.3781e6;// Radius of Earth in m

template<typename Scalar = double>
inline constexpr Scalar AU = 1.495978707e11;// Astronomical Unit in m

// Conversion factors
template<typename Scalar = double>
inline constexpr Scalar DEG2RAD = PI<Scalar> / 180.0;// Degree to Radian

template<typename Scalar = double>
inline constexpr Scalar RAD2DEG = 180.0 / PI<Scalar>;// Radian to Degree

template<typename Scalar = double>
inline constexpr Scalar E = static_cast<Scalar>(M_E);// Euler's number


}// namespace odin::consts

#endif// ODIN_CONSTANTS_H