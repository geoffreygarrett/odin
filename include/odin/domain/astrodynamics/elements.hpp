#include <Eigen/Core>
#include <array>
#include <vector>

template<typename Scalar>
using Vector3 = Eigen::Vector<Scalar, 3>;


template<typename Scalar>
struct ClassicalElements {
    Scalar p;
    Scalar e;
    Scalar i;
    Scalar raan;
    Scalar argp;
    Scalar nu;
};

template<typename Scalar, int N>
struct ClassicalElementsBatch {
    std::array<Scalar, N> p;
    std::array<Scalar, N> e;
    std::array<Scalar, N> i;
    std::array<Scalar, N> raan;
    std::array<Scalar, N> argp;
    std::array<Scalar, N> nu;
};

template<typename Scalar>
struct ClassicalElementsBatch<Scalar, -1> {
    std::vector<Scalar> p;
    std::vector<Scalar> a;
    std::vector<Scalar> e;
    std::vector<Scalar> i;
    std::vector<Scalar> raan;
    std::vector<Scalar> argp;
    std::vector<Scalar> nu;
};

template<typename Scalar>
struct CartesianElements {
    Vector3<Scalar> r;
    Vector3<Scalar> v;
};

template<typename Scalar, int N>
struct CartesianElementsBatch {
    std::array<Vector3<Scalar>, N> r;
    std::array<Vector3<Scalar>, N> v;
};

template<typename Scalar>
struct CartesianElementsBatch<Scalar, -1> {
    std::vector<Vector3<Scalar>> r;
    std::vector<Vector3<Scalar>> v;
};

template<typename Scalar>
struct EquinoctialElements {
    Scalar p;
    Scalar f;
    Scalar g;
    Scalar h;
    Scalar k;
    Scalar L;
};

template<typename Scalar, int N>
struct EquinoctialElementsBatch {
    std::array<Scalar, N> p;
    std::array<Scalar, N> f;
    std::array<Scalar, N> g;
    std::array<Scalar, N> h;
    std::array<Scalar, N> k;
    std::array<Scalar, N> L;
};

template<typename Scalar>
struct EquinoctialElementsBatch<Scalar, -1> {
    std::vector<Scalar> p;
    std::vector<Scalar> f;
    std::vector<Scalar> g;
    std::vector<Scalar> h;
    std::vector<Scalar> k;
    std::vector<Scalar> L;
};
