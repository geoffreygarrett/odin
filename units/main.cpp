#include "conversion.hpp"
#include "math.h"
#include <iostream>
#include <ratio>
#include <type_traits>

//template<typename Expr>
//class Expression {
//public:
//    const Expr &expr() const {
//        return static_cast<const Expr &>(*this);
//    }
//};
//
//template<typename T, typename U>
//class Quantity : public Expression<Quantity<T, U>> {
//public:
//    explicit Quantity(T value) : m_value(value) {}
//
//    T value() const { return m_value; }
//
//    template<typename NewUnit>
//    Quantity<T, NewUnit> convert_to() const {
//        double conversion = conversion_factor<U, NewUnit>();
//        return Quantity<T, NewUnit>(m_value * conversion);
//    }
//
//private:
//    T m_value;
//
//    // Conversion factor function template (replace this with actual conversion factors)
//    template<typename From, typename To>
//    static double conversion_factor() {
//        return 1.0;// Default to no conversion
//    }
//};
//
//template<typename E1, typename E2>
//class Multiply : public Expression<Multiply<E1, E2>> {
//public:
//    Multiply(const E1 &e1, const E2 &e2) : e1(e1), e2(e2) {}
//
//    auto value() const {
//        return e1.value() * e2.value();
//    }
//
//private:
//    const E1 &e1;
//    const E2 &e2;
//};
//
//template<typename E1, typename E2, std::enable_if_t<std::is_base_of_v<Expression<E1>, E1> || std::is_base_of_v<Expression<E2>, E2>, int> = 0>
//auto operator*(const E1 &e1, const E2 &e2) {
//    return Multiply<E1, E2>(e1, e2);
//}
//
//
//// Definition of compound units
//template<typename... T>
//struct compound_unit;
//
//template<typename T>
//struct compound_unit<T> {
//    using type = T;
//    static constexpr char symbol[] = T::symbol;
//};
//
//template<typename T, typename... Rest>
//struct compound_unit<T, Rest...> {
//    using type = compound_unit<T, Rest...>;
//    static constexpr char symbol[] = T::symbol + "•" + compound_unit<Rest...>::symbol;
//};


// Base units
namespace unit {
    struct meter {};   // unit of length
    struct kilogram {};// unit of mass
    struct second {};  // unit of time
    struct ampere {};  // unit of electric current
    struct kelvin {};  // unit of temperature
    struct mole {};    // unit of amount of substance
    struct candela {}; // unit of luminous intensity
}// namespace unit

// The Unit struct: defines a unit in terms of the exponents of the SI base units
template<int M, int Kg, int S = 0, int A = 0, int K = 0, int Mol = 0, int Cd = 0, typename Factor = std::ratio<1>>
struct Unit {
    static constexpr int M_exp   = M;
    static constexpr int Kg_exp  = Kg;
    static constexpr int S_exp   = S;
    static constexpr int A_exp   = A;
    static constexpr int K_exp   = K;
    static constexpr int Mol_exp = Mol;
    static constexpr int Cd_exp  = Cd;
    using scale_factor           = Factor;
};

template<typename T, typename Scalar = std::ratio<1>>
class Quantity {
public:
    explicit Quantity(double value) : m_value(value * T::scale_factor::num / T::scale_factor::den * Scalar::den / Scalar::num) {}
    double value() const { return m_value * T::scale_factor::den / T::scale_factor::num * Scalar::num / Scalar::den; }// Return the value in the unit defined by T

    // Multiply this Quantity with another Quantity
    template<typename U>
    auto operator*(const Quantity<U> &other) const {
        return Quantity<Unit<T::M_exp + U::M_exp, T::Kg_exp + U::Kg_exp, T::S_exp + U::S_exp, T::A_exp + U::A_exp, T::K_exp + U::K_exp, T::Mol_exp + U::Mol_exp, T::Cd_exp + U::Cd_exp>>(value() * other.value());
    }

    // Multiply this Quantity with a scalar
    template<typename Scalar2>
    auto operator*(const Scalar2 &scalar) const {
        return Quantity<T, Scalar2>(value() * scalar);
    }


private:
    double m_value;
};

// Scalar is a physical quantity with magnitude but without direction.
// Here, we define a Scalar as a special case of Quantity with no physical dimension.
template<int N = 1, int D = 1>
using Scalar = Quantity<Unit<0, 0, 0, 0, 0, 0, 0, std::ratio<N, D>>>;

// Fraction is used to represent constant ratios, it's a special case of Scalar.
template<int N, int D = 1>
using Fraction = Scalar<N, D>;

// Operator overloads to handle multiplication between a Quantity and a Fraction.
template<typename T, typename Scalar, int N, int D>
auto operator*(const Quantity<T, Scalar> &quantity, const Fraction<N, D> &fraction) {
    // The resultant Quantity has its Scalar modified according to the Fraction.
    return Quantity<T, std::ratio_multiply<Scalar, std::ratio<N, D>>>(quantity.value());
}

template<typename T, typename Scalar, int N, int D>
auto operator*(const Fraction<N, D> &fraction, const Quantity<T, Scalar> &quantity) {
    // The resultant Quantity has its Scalar modified according to the Fraction.
    return Quantity<T, std::ratio_multiply<Scalar, std::ratio<N, D>>>(quantity.value());
}

// Define units. Here we represent each unit as a special case of Quantity with specific dimensions.
using Meter         = Unit<1, 0, 0, 0, 0, 0, 0>;                       // The base unit of length.
using Kilometer     = Unit<1, 0, 0, 0, 0, 0, 0, std::kilo>;            // 1000 meters.
using Foot          = Unit<1, 0, 0, 0, 0, 0, 0, std::ratio<381, 1250>>;// Approximately 0.3048 meters.
using GravConstant  = Unit<3, -1, -2>;                                 // The gravitational constant, approximately 6.674 × 10⁻¹¹ m³ kg⁻¹ s⁻².
using GravParameter = Unit<3, 0, -2>;                                  // The gravitational parameter, usually represented by μ.
using Density       = Unit<-3, 1>;                                     // Mass per unit volume.

// Define types for specific quantities
using LengthMeter       = Quantity<Meter>;
using LengthKilometer   = Quantity<Kilometer>;
using DensityUnit       = Quantity<Density>;
using GravConstantUnit  = Quantity<GravConstant>;
using GravParameterUnit = Quantity<GravParameter>;


int main() {
    // Define your quantities
    LengthMeter      a(300.0);          // in meters
    LengthMeter      b(200.0);          // in meters
    LengthMeter      c(100.0);          // in meters
    DensityUnit      rho(2.0 * 1000.0); // in kg/m³
    GravConstantUnit G(6.67408 * 1e-11);// in m³/(kg·s²)

    //    // Calculate gravitational parameter
    GravParameterUnit mu = G * a * b * c * Fraction<4, 3>(1) * Scalar<1>(M_PI) * rho;
    std::cout << "mu in m³/s²: " << mu.value() << std::endl;

    // Show values
    std::cout << "a in meters: " << a.value() << std::endl;
    std::cout << "b in meters: " << b.value() << std::endl;
    std::cout << "c in meters: " << c.value() << std::endl;
    std::cout << "rho in kg/m³: " << rho.value() << std::endl;
    std::cout << "G in m³/(kg·s²): " << G.value() << std::endl;
    std::cout << "mu in m³/s²: " << mu.value() << std::endl;

    //    std::cout << "Step size in meters: " << step_size.value() << std::endl;
    //    std::cout << "Grid range in meters: " << grid_range.value() << std::endl;
    //    std::cout << "Number of steps: " << num_steps << std::endl;
    //
    //    // Demonstrate unit conversion
    //    Quantity<double, LengthUnit> a_km = a.convert_to<Kilometer>();
    //    std::cout << "a in kilometers: " << a_km.value() << std::endl;

    return 0;
}
