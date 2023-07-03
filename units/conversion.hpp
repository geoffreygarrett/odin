//// units
//namespace unit {
//    struct meter {
//        static constexpr char symbol[] = "m";
//    };
//    struct kilometer {
//        static constexpr char symbol[] = "km";
//    };
//    struct kilogram {
//        static constexpr char symbol[] = "kg";
//    };
//    struct second {
//        static constexpr char symbol[] = "s";
//    };
//    struct ampere {
//        static constexpr char symbol[] = "A";
//    };
//    struct kelvin {
//        static constexpr char symbol[] = "K";
//    };
//    struct mole {
//        static constexpr char symbol[] = "mol";
//    };
//    struct candela {
//        static constexpr char symbol[] = "cd";
//    };
//}// namespace unit
//
//// Conversion factor function template
//template<typename From, typename To>
//double conversion_factor();
//
//// Conversion factor for meters to meters
//template<>
//double conversion_factor<unit::meter, unit::meter>() {
//    return 1.0;
//}
//
//// Conversion factor for meters to kilometers
//template<>
//double conversion_factor<unit::meter, unit::kilometer>() {
//    return 1e-3;
//}
//
//// Conversion factor for meters to feet
//template<>
//double conversion_factor<unit::meter, unit::foot>() {
//    return 3.28084;
//}
//
//// Conversion factor for kilometers to meters
//template<>
//double conversion_factor<unit::kilometer, unit::meter>() {
//    return 1e3;
//}
//
//// Conversion factor for kilometers to kilometers
//template<>
//double conversion_factor<unit::kilometer, unit::kilometer>() {
//    return 1.0;
//}
//
//// Conversion factor for kilometers to feet
//template<>
//double conversion_factor<unit::kilometer, unit::foot>() {
//    return 3.28084 * 1e3;
//}
//
//// Conversion factor for feet to meters
//template<>
//double conversion_factor<unit::foot, unit::meter>() {
//    return 1 / 3.28084;
//}
//
//// Conversion factor for feet to kilometers
//template<>
//double conversion_factor<unit::foot, unit::kilometer>() {
//    return (1 / 3.28084) * 1e-3;
//}
//
//// Conversion factor for feet to feet
//template<>
//double conversion_factor<unit::foot, unit::foot>() {
//    return 1.0;
//}