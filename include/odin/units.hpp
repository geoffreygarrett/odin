#include <cmath>
#include <iostream>
#include <ratio>

// Base dimensions
template<int M, int Kg, int S, int A, int K, int Mol, int Cd>
struct Dimension {
    enum { Meter    = M,
           Kilogram = Kg,
           Second   = S,
           Ampere   = A,
           Kelvin   = K,
           Mole     = Mol,
           Candela  = Cd };
};

// Scalar dimension (for pure numbers)
using ScalarDimension = Dimension<0, 0, 0, 0, 0, 0, 0>;
using Length          = Dimension<1, 0, 0, 0, 0, 0, 0>;
using Mass            = Dimension<0, 1, 0, 0, 0, 0, 0>;
using Time            = Dimension<0, 0, 1, 0, 0, 0, 0>;
using Current         = Dimension<0, 0, 0, 1, 0, 0, 0>;
using Temperature     = Dimension<0, 0, 0, 0, 1, 0, 0>;
using Amount          = Dimension<0, 0, 0, 0, 0, 1, 0>;
using Luminosity      = Dimension<0, 0, 0, 0, 0, 0, 1>;
using Angle           = Dimension<0, 0, 0, 0, 0, 0, 0>;

template<typename D, typename T = double>
struct Quantity;

template<typename Op1, typename Op2, typename Operation>
struct Expression {
    const Op1 &op1;
    const Op2 &op2;
    Operation  operation;

    Expression(const Op1 &op1, const Op2 &op2, Operation operation)
        : op1(op1), op2(op2), operation(operation) {}

    [[nodiscard]] auto value() const {
        return operation(op1.value, op2.value);
    }

    template<typename D, typename T = double>
    explicit operator Quantity<D, T>() const {
        return Quantity<D, T>(value());
    }
};

template<typename T>
concept QuantityType = requires(T a, T b) {
    { a.value } -> std::convertible_to<double>;
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
    // Define other operations you expect a Quantity to have
};

template<QuantityType Q>
Q operator+(const Q &lhs, const Q &rhs) {
    return Q(lhs.value + rhs.value);
}

template<typename T>
concept AngleType = requires(T a) {
    { a.to_radians() } -> std::convertible_to<double>;
};


template<typename D, typename T>
struct Quantity {
    T value;

    explicit Quantity(T val) : value(val) {}

    template<typename Op1, typename Op2, typename Operation>
    explicit Quantity(const Expression<Op1, Op2, Operation> &expr) : value(expr.value()) {}

    Quantity &operator=(const Quantity &rhs) {
        value = rhs.value;
        return *this;
    }

    Quantity operator+(const Quantity &rhs) const {
        return Quantity(value + rhs.value);
    }

    Quantity operator-(const Quantity &rhs) const {
        return Quantity(value - rhs.value);
    }

    template<int exponent>
    [[nodiscard]] Quantity<
            Dimension<
                    D::Meter * exponent,
                    D::Kilogram * exponent,
                    D::Second * exponent,
                    D::Ampere * exponent,
                    D::Kelvin * exponent,
                    D::Mole * exponent,
                    D::Candela * exponent>,
            T>
    pow() const {
        return Quantity<
                Dimension<
                        D::Meter * exponent,
                        D::Kilogram * exponent,
                        D::Second * exponent,
                        D::Ampere * exponent,
                        D::Kelvin * exponent,
                        D::Mole * exponent,
                        D::Candela * exponent>,
                T>(std::pow(value, exponent));
    }

    [[nodiscard]] T to_radians() const {
        static_assert(std::is_same<D, Angle>::value, "Can only convert angle quantities to radians");
        return value;
    }

    [[nodiscard]] T to_degrees() const {
        static_assert(std::is_same<D, Angle>::value, "Can only convert angle quantities to degrees");
        return value * 180.0 / M_PI;
    }

    Quantity<Dimension<0, 0, 0, 0, 0, 0, 0>, T> as_rad() const {
        static_assert(std::is_same<D, Angle>::value, "Can only convert angle quantities to radians");
        return Quantity<Dimension<0, 0, 0, 0, 0, 0, 0>, T>(to_radians());
    }

    Quantity<Dimension<0, 0, 0, 0, 0, 0, 0>, T> as_deg() const {
        static_assert(std::is_same<D, Angle>::value, "Can only convert angle quantities to degrees");
        return Quantity<Dimension<0, 0, 0, 0, 0, 0, 0>, T>(to_degrees());
    }
};

Quantity<Dimension<1, 0, 0, 0, 0, 0, 0>> operator"" _m(long double value) {
    return Quantity<Dimension<1, 0, 0, 0, 0, 0, 0>>(value);
}
Quantity<Dimension<0, 1, 0, 0, 0, 0, 0>> operator"" _kg(long double value) {
    return Quantity<Dimension<0, 1, 0, 0, 0, 0, 0>>(value);
}
Quantity<Dimension<0, 0, 1, 0, 0, 0, 0>> operator"" _s(long double value) {
    return Quantity<Dimension<0, 0, 1, 0, 0, 0, 0>>(value);
}


// Angle dimension TODO: Move to separate file
using Angle = Dimension<0, 0, 0, 0, 0, 0, 0>;// This is just a scalar, so all dimensions are 0

// Literals for angle units
Quantity<Angle> operator"" _rad(long double value) { return Quantity<Angle>(value); }
Quantity<Angle> operator"" _deg(long double value) { return Quantity<Angle>(value * M_PI / 180.0); }// Convert to radians


// Arithmetic operations for quantities
template<typename D1, typename T1, typename D2, typename T2>
auto operator/(const Quantity<D1, T1> &lhs, const Quantity<D2, T2> &rhs) {
    return Expression(lhs, rhs, [](auto l, auto r) { return l / r; });
}

template<typename D1, typename T1, typename D2, typename T2>
auto operator*(const Quantity<D1, T1> &lhs, const Quantity<D2, T2> &rhs) {
    return Expression(lhs, rhs, [](auto l, auto r) { return l * r; });
}


namespace odin {

template<typename T>
T cos(const Quantity<Angle, T> &angle) {
    return std::cos(angle.to_radians());
}

template<typename T>
T sin(const Quantity<Angle, T> &angle) {
    return std::sin(angle.to_radians());
}

template<typename T>
T tan(const Quantity<Angle, T> &angle) {
    return std::tan(angle.to_radians());
}

template<typename T>
T acos(const Quantity<Angle, T> &angle) {
    return std::acos(angle.to_radians());
}

template<typename T>
T asin(const Quantity<Angle, T> &angle) {
    return std::asin(angle.to_radians());
}

template<typename T>
T atan(const Quantity<Angle, T> &angle) {
    return std::atan(angle.to_radians());
}

template<typename T>
T atan2(const Quantity<Angle, T> &y, const Quantity<Angle, T> &x) {
    return std::atan2(y.to_radians(), x.to_radians());
}

template<typename T>
T sinh(const Quantity<Angle, T> &angle) {
    return std::sinh(angle.to_radians());
}

template<typename T>
T cosh(const Quantity<Angle, T> &angle) {
    return std::cosh(angle.to_radians());
}

template<typename T>
T tanh(const Quantity<Angle, T> &angle) {
    return std::tanh(angle.to_radians());
}

template<typename T>
T asinh(const Quantity<Angle, T> &angle) {
    return std::asinh(angle.to_radians());
}

template<typename T>
T acosh(const Quantity<Angle, T> &angle) {
    return std::acosh(angle.to_radians());
}

template<typename T>
T atanh(const Quantity<Angle, T> &angle) {
    return std::atanh(angle.to_radians());
}

}// namespace odin


using Meter      = Dimension<1, 0, 0, 0, 0, 0, 0>;
using Kilometer  = Dimension<1, 0, 0, 0, 0, 0, 0>;// Keeping the same dimension for simplicity
using Millimeter = Dimension<1, 0, 0, 0, 0, 0, 0>;// Keeping the same dimension for simplicity


namespace test {
// Conversion factor
constexpr double MetersPerFoot = 0.3048;

template<int MeterScale, int FootScale = 0, typename T = double>
struct Length {
    T value;

    explicit Length(T val) : value(val) {}

    [[nodiscard]] T to_base_units() const {
        return value * std::pow(10, MeterScale) / std::pow(MetersPerFoot, FootScale);
    }

    template<int NewMeterScale, int NewFootScale = 0>
    [[nodiscard]] Length<NewMeterScale, NewFootScale, T> convert() const {
        T newValue = to_base_units() * std::pow(10, -NewMeterScale) * std::pow(MetersPerFoot, NewFootScale);
        return Length<NewMeterScale, NewFootScale, T>(newValue);
    }

    // Conversion methods
    Length<0>                  as_m() const { return convert<0>(); }
    Length<3>                  as_km() const { return convert<3>(); }
    Length<-3>                 as_mm() const { return convert<-3>(); }
    [[nodiscard]] Length<0, 1> as_ft() const { return convert<0, 1>(); }
};
namespace literals {

// Literal operators
Length<0> operator"" _m(long double value) {
    return Length<0>(value);
}
Length<3> operator"" _km(long double value) {
    return Length<3>(value);
}
Length<-3> operator"" _mm(long double value) {
    return Length<-3>(value);
}
Length<0, 1> operator"" _ft(long double value) {
    return Length<0, 1>(value);
}

// force
using Newton = Dimension<1, 1, -2, 0, 0, 0, 0>;

Quantity<Newton> operator"" _N(long double value) {
    return Quantity<Newton>(value);
}

}// namespace literals

//Quantity<Meter> operator"" _m(long double value) {
//    return Quantity<Meter>(value);
//}
//
//Quantity<Kilometer> operator"" _km(long double value) {
//    return Quantity<Kilometer>(value);
//}
//
//Quantity<Millimeter> operator"" _mm(long double value) {
//    return Quantity<Millimeter>(value);
//}

//// Define base unit (meters)
//Length<0> operator"" _m(long double value) {
//    return Length<0>(value);
//}
//
//// Define other units (kilometers and millimeters)
//Length<3> operator"" _km(long double value) {
//    return Length<3>(value);
//}
//
//Length<-3> operator"" _mm(long double value) {
//    return Length<-3>(value);
//}

}// namespace test

int main() {
    // repeat after me: clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    //                  clarity and expressiveness - clarity and expressiveness - clarity and expressiveness
    auto dist   = 100.0_m;
    auto time   = 10.0_s;
    auto speed  = dist / time;
    auto speed2 = 100.0_m / 10.0_s;
    auto speed3 = 10.0_m / 1.0_s;

    std::cout << "Speed: " << speed.value() << " m/s" << std::endl;

    auto volume = dist.pow<3>();
    std::cout << "Volume: " << volume.value << " m^3" << std::endl;

    auto angle_deg = 90.0_deg;
    std::cout << "Angle in degrees: " << angle_deg.to_degrees() << " deg" << std::endl;
    std::cout << "Angle in radians: " << angle_deg.to_radians() << " rad" << std::endl;

    auto angle_rad = (90.0_deg).as_rad();
    std::cout << "Angle in degrees: " << angle_rad.to_degrees() << " deg" << std::endl;
    std::cout << "Angle in radians: " << angle_rad.to_radians() << " rad" << std::endl;

    //    using namespace odin;

    auto x1 = odin::acosh(1.0012301203120930129309120391023901293_rad);

    //    auto angle_deg = 90.0_deg;
    //    auto angle_rad = (90.0_deg).to_radians();
    //
    //    std::cout << "cos(90 deg): " << cos(angle_deg) << '\n';
    //    std::cout << "sin(90 deg): " << sin(angle_deg) << '\n';
    //    std::cout << "tan(90 deg): " << tan(angle_deg) << '\n';
    //    std::cout << "acos(0.0): " << acos(quantity<angle>(0.0)).to_degrees() << '\n';
    //    std::cout << "asin(1.0): " << asin(quantity<angle>(1.0)).to_degrees() << '\n';
    //    std::cout << "atan(1.0): " << atan(quantity<angle>(1.0)).to_degrees() << '\n';
    //    std::cout << "atan2(1.0, 0.0): " << atan2(quantity<angle>(1.0), quantity<angle>(0.0)).to_degrees() << '\n';
    //    std::cout << "sinh(90 deg): " << sinh(angle_deg) << '\n';
    //    std::cout << "cosh(90 deg): " << cosh(angle_deg) << '\n';
    //    std::cout << "tanh(90 deg): " << tanh(angle_deg) << '\n';
    //    std::cout << "asinh(1.0): " << asinh(quantity<angle>(1.0)).to_degrees() << '\n';
    //    std::cout << "acosh(1.0): " << acosh(quantity<angle>(1.0)).to_degrees() << '\n';
    //    std::cout << "atanh(0.5): " << atanh(quantity<angle>(0.5)).to_degrees() << '\n';

    using namespace test;
    auto length_mm = 1000.0_mm;
    std::cout << "Length in millimeters: " << length_mm.value << " mm" << std::endl;
    std::cout << "Length in meters: " << length_mm.convert<0>().value << " m" << std::endl;
    std::cout << "Length in kilometers: " << length_mm.convert<3>().value << " km" << std::endl;

    auto length_m = length_mm.convert<0>();
    std::cout << "Length in millimeters: " << length_m.convert<-3>().value << " mm" << std::endl;
    std::cout << "Length in meters: " << length_m.value << " m" << std::endl;
    std::cout << "Length in kilometers: " << length_m.convert<3>().value << " km" << std::endl;

    std::cout << "Length in millimeters: " << length_mm.value << " mm" << std::endl;
    std::cout << "Length in meters: " << length_mm.as_m().value << " m" << std::endl;
    std::cout << "Length in kilometers: " << length_mm.as_km().value << " km" << std::endl;
    std::cout << "Length in feet: " << length_mm.as_ft().value << " ft" << std::endl;
    using namespace test::literals;

    using namespace test;

    auto work_done = 100.0_N * 10.0_mm;

    auto my_function(
            (1e6_m).as_sma(),
            0.5_ecc,
            0.0_inc,
            0.0_raan,
            0.0_argp,
            0.0_ta,
            0.0_mu);


    return 0;
}