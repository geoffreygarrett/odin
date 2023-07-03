

//#include <Eigen/Core>
//template<typename Derived, typename Float, int N>
//class Body {
//public:
//    using Scalar      = Float;
//    using Vector      = Vector3<Scalar>;
//    using VectorBatch = Vector3Batch<Scalar, N>;
//    using Matrix      = Matrix3<Scalar>;
//    using MatrixBatch = Matrix3Batch<Scalar, N>;
//
//    CelestialBody() = default;
//
//    CelestialBody(const Vector &position, const Vector &velocity, Scalar mass)
//        : position(position), velocity(velocity), mass(mass) {}
//
//    CelestialBody(const VectorBatch &position, const VectorBatch &velocity, const std::array<Scalar, N> &mass)
//        : position(position), velocity(velocity), mass(mass) {}
//
//    CelestialBody(const VectorBatch &position, const VectorBatch &velocity, const std::vector<Scalar> &mass)
//        : position(position), velocity(velocity), mass(mass) {}
//
//    CelestialBody(const CelestialBody &other) = default;
//
//    CelestialBody(CelestialBody &&other) noexcept = default;
//
//    CelestialBody &operator=(const CelestialBody &other) = default;
//
//    CelestialBody &operator=(CelestialBody &&other) noexcept = default;
//
//    ~CelestialBody() = default;
//
//    Vector position;
//    Vector velocity;
//    Scalar mass;
//};