#ifndef ODIN_ELLIPSOID_PLANE_INTERSECTION_H
#define ODIN_ELLIPSOID_PLANE_INTERSECTION_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace odin::algorithms {

template<typename T>
bool safe_to_divide(T _, T den) {
    return std::abs(den) > std::numeric_limits<T>::epsilon();
}

template<typename U>
using vec3_type = Eigen::Vector3<U>;

template<typename U>
vec3_type<U> calculate_polar_plane_pole(
        vec3_type<U> normal, U offset, vec3_type<U> semi_axes) {
    auto n = normal.normalized();
    if (safe_to_divide(1., offset)) {
        return (semi_axes.array().pow(2).matrix().cwiseProduct(n) / offset)
                .eval();
    } else {
        return (semi_axes.array().pow(2).matrix().cwiseProduct(n)).eval();
    }
}
template<typename U>
vec3_type<U> calculate_polar_plane_pole(
        vec3_type<U> normal, U f, U a, U b, U c) {
    return calculate_polar_plane_pole(normal, f, vec3_type<U>(a, b, c));
}

template<typename U>
vec3_type<U> calculate_polar_plane_pole(U nx, U ny, U nz, U f, U a, U b, U c) {
    return calculate_polar_plane_pole(
            vec3_type<U>(nx, ny, nz), f, vec3_type<U>(a, b, c));
}

template<typename U>
std::tuple<U, U, U, U> polar_plane_from_point(
        const Eigen::Vector3d &p, U a, U b, U c) {
    // TODO: check if p is on / inside the ellipsoid
    U a2 = a * a;
    U b2 = b * b;
    U c2 = c * c;
    U f  = 1.0// clang-format off
             / std::sqrt(std::pow(p[0] / a2, 2)
                       + std::pow(p[1] / b2, 2)
                       + std::pow(p[2] / c2, 2));// clang-format on

    U nx = p[0] * f / a / a;
    U ny = p[1] * f / b / b;
    U nz = p[2] * f / c / c;
    return {nx, ny, nz, f};
}

template<typename U>
vec3_type<U> generate_orthogonal_vector(const vec3_type<U> &vi) {
    vec3_type<U> vo;
    for (int i = 0; i < 3; ++i) {
        vo = vec3_type<U>::Unit(i).cross(vi);
        if (vo.norm() > std::numeric_limits<U>::epsilon()) { break; }
    }
    return vo;
}

template<typename U>
std::tuple<vec3_type<U>, vec3_type<U>, vec3_type<U>, U, U>
ellipse_of_intersection(vec3_type<U> normal, U offset, vec3_type<U> semi_axes) {
    // Nomenclature:
    // =================
    // n : unit normal vector of plane
    // q : arbitrary point on plane, interior to the ellipse
    // r, s : arbitrary vectors perpendicular to n
    // D1 : diagonal matrix with 1/a, 1/b, 1/c
    // t : first parameter of the 2D parametric equation of the ellipse
    // u : second parameter of the 2D parametric equation of the ellipse
    // t0 : first parameter of the 2D parametric equation of the ellipse,
    //      center of ellipse
    // u0 : second parameter of the 2D parametric equation of the ellipse,
    //      center of ellipse
    // A : semi-major axis of ellipse
    // B : semi-minor axis of ellipse
    // k : distance of the plane from the origin

    // Ensure that the normal vector is normalized
    auto n = normal.normalized();

    // Define D1 as a diagonal matrix with 1/a, 1/b, 1/c
    Eigen::DiagonalMatrix<U, 3> D1(
            1 / semi_axes[0], 1 / semi_axes[1], 1 / semi_axes[2]);

    // Define q, any point on plane, interior to the ellipse
    auto q = n * offset;

    // Define r, s, arbitrary vectors perpendicular to n, satisfying
    // dot(r, r) = dot(s, s) = 1, dot(n, r) = dot(n, s) = 0, dot(r, s) = 0
    Eigen::Vector3<U> r, s;

    // iterate through the unit vectors until we find one that is not parallel to n
    r = generate_orthogonal_vector(n);
    s = n.cross(r);

    // In case (7) is not fulfilled for the chosen r and s, dot(D1 @ r, D1 @ s) != 0, we transform:
    auto D1r = (D1 * r).eval();
    auto D1s = (D1 * s).eval();

    if (D1r.dot(D1s) != 0) {
        U omega;
        if (D1r.dot(D1r) - D1s.dot(D1s) == 0) {
            omega = M_PI_4;
        } else {
            omega = 0.5
                  * std::atan(2 * D1r.dot(D1s) / (D1r.dot(D1r) - D1s.dot(D1s)));
        }
        r = (std::cos(omega) * r + std::sin(omega) * s).normalized();
        s = n.cross(r).normalized();

        // Recalculate D1r and D1s
        D1r = D1 * r;
        D1s = D1 * s;
    }

    // Calculate beta_1 and beta_2
    auto beta_1 = D1r.dot(D1r);
    auto beta_2 = D1s.dot(D1s);

    // Expression for d, k (Eq. 11)
    auto k = q.dot(n);
    auto d = k * k
           / (semi_axes.array().pow(2).matrix().dot(n.array().pow(2).matrix()));

    // Ellipse center, m (Eq. 39)
    vec3_type<U> m;
    auto         D1n = D1 * n;
    if (k == 0) {
        m = vec3_type<U>::Zero();
    } else {
        m = k * (n - D1n.dot(D1r) / beta_1 * r - D1n.dot(D1s) / beta_2 * s);
    }

    // Semi-axes, A and B (Eq. 10)
    auto A = std::sqrt((1 - d) / beta_1);
    auto B = std::sqrt((1 - d) / beta_2);

    return std::make_tuple(m, r, s, A, B);
}

template<typename U>
std::tuple<vec3_type<U>, vec3_type<U>, vec3_type<U>, U, U>
ellipse_of_intersection(U nx, U ny, U nz, U f, U a, U b, U c) {
    return ellipse_of_intersection(
            vec3_type<U>(nx, ny, nz), f, vec3_type<U>(a, b, c));
}

template<typename U>
std::tuple<vec3_type<U>, vec3_type<U>, vec3_type<U>, U, U>
ellipse_of_intersection(vec3_type<U> n, U f, U a, U b, U c) {
    return ellipse_of_intersection(n, f, vec3_type<U>(a, b, c));
}

//template<typename U>
//__global__
//        void ellipse_of_intersection(vec3_type<U> normal, vec3_type<U> semi_axes, U offset,
//                                vec3_type<U>* m, vec3_type<U>* r, vec3_type<U>* s, U* A, U* B) {
//
//    auto n = normal.normalized();
//    Eigen::DiagonalMatrix<U, 3> D1(1 / semi_axes[0], 1 / semi_axes[1], 1 / semi_axes[2]);
//    auto q = n * offset;
//
//    for (auto i = 0; i < 3; ++i) {
//        *r = Eigen::Vector3<U>::Unit(i).cross(n);
//        if (r->norm() > std::numeric_limits<U>::epsilon()) { break; }
//    }
//
//    *s = n.cross(*r);
//
//    auto D1r = D1 * *r;
//    auto D1s = D1 * *s;
//
//    if (D1r.dot(D1s) != 0) {
//        U omega;
//        if (D1r.dot(D1r) - D1s.dot(D1s) == 0) {
//            omega = M_PI_4;
//        } else {
//            omega = 0.5
//                  * std::atan(2 * D1r.dot(D1s) / (D1r.dot(D1r) - D1s.dot(D1s)));
//        }
//        *r = (std::cos(omega) * *r + std::sin(omega) * *s).normalized();
//        *s = n.cross(*r).normalized();
//
//        D1r = D1 * *r;
//        D1s = D1 * *s;
//    }
//
//    auto beta_1 = D1r.dot(D1r);
//    auto beta_2 = D1s.dot(D1s);
//    auto k = q.dot(n);
//    auto d = k * k / (semi_axes.pow(2).dot(n.pow(2)));
//
//    if (k == 0) {
//        *m = vec3_type<U>::Zero();
//    } else {
//        *m = k * (n - D1r.dot(n) / beta_1 * *r - D1s.dot(n) / beta_2 * *s);
//    }
//
//    *A = std::sqrt((1 - d) / beta_1);
//    *B = std::sqrt((1 - d) / beta_2);
//}


}// namespace odin::algorithms

#endif//ODIN_ELLIPSOID_PLANE_INTERSECTION_H