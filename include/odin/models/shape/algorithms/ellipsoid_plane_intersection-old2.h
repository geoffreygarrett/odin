#ifndef ODIN_ELLIPSOID_PLANE_INTERSECTION_H
#define ODIN_ELLIPSOID_PLANE_INTERSECTION_H

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace odin::algorithms {


template<typename U, std::size_t dim>
class Vector {
    std::array<U, dim> values;

public:
    Vector() : values{} {}

    explicit Vector(const std::array<U, dim> &arr) : values(arr) {}

    U       &operator[](std::size_t i) { return values[i]; }
    const U &operator[](std::size_t i) const { return values[i]; }

    U norm() const {
        return std::sqrt(std::inner_product(values.begin(), values.end(),
                                            values.begin(), U(0)));
    }

    Vector(std::initializer_list<U> init_list) {
        assert(init_list.size() == dim);
        std::copy(init_list.begin(), init_list.end(), values.begin());
    }

    Vector operator/(U scalar) const {
        std::array<U, dim> result;
        std::transform(values.begin(), values.end(), result.begin(),
                       [scalar](U val) { return val / scalar; });
        return Vector(result);
    }

    std::array<std::size_t, dim> argsort() const {
        std::array<std::size_t, dim> indices;
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
                  [this](std::size_t i1, std::size_t i2) {
                      return std::abs(values[i1]) < std::abs(values[i2]);
                  });
        return indices;
    }

    Vector operator*(const std::array<std::array<U, dim>, dim> &P) const {
        std::array<U, dim> result{};
        for (std::size_t i = 0; i < dim; i++) {
            for (std::size_t j = 0; j < dim; j++) {
                result[i] += P[i][j] * values[j];
            }
        }
        return Vector(result);
    }
};

template<typename U, std::size_t dim>
class Matrix {

public:
    Matrix() : values{} {}
    std::array<std::array<U, dim>, dim> values;

    Matrix(std::initializer_list<std::initializer_list<U>> list) {
        assert(list.size() <= dim);// Include <cassert> for this line
        std::size_t i = 0;
        for (const auto &inner_list: list) {
            assert(inner_list.size() <= dim);// Include <cassert> for this line
            std::copy(inner_list.begin(), inner_list.end(), values[i].begin());
            ++i;
        }
    }


    explicit Matrix(const std::array<std::array<U, dim>, dim> &arr)
        : values(arr) {}

    explicit Matrix(std::array<std::size_t, dim> indices) {
        // iterate through all, initialize, then set the 1s
        for (std::size_t i = 0; i < dim; i++) {
            std::fill(values[i].begin(), values[i].end(), 0);
        }
        for (std::size_t i = 0; i < dim; i++) {
            this->values[i][indices[i]] = 1;
        }
    }

    Matrix transpose() const {
        std::array<std::array<U, dim>, dim> transposed{};
        for (std::size_t i = 0; i < dim; i++) {
            for (std::size_t j = 0; j < dim; j++) {
                transposed[i][j] = this->values[j][i];
            }
        }
        return Matrix(transposed);
    }

    Vector<U, dim> operator*(const Vector<U, dim> &v) const {
        std::array<U, dim> result{};
        for (std::size_t i = 0; i < dim; i++) {
            for (std::size_t j = 0; j < dim; j++) {
                result[i] += this->values[i][j] * v[j];
            }
        }
        return Vector<U, dim>(result);
    }

    Matrix operator*(const Matrix<U, dim> &m) const {
        Matrix<U, dim> result;
        for (std::size_t i = 0; i < dim; i++) {
            for (std::size_t j = 0; j < dim; j++) {
                for (std::size_t k = 0; k < dim; k++) {
                    result.values[i][j] += this->values[i][k] * m.values[k][j];
                }
            }
        }
        return result;
    }
};

template<typename U, std::size_t dim>
Matrix<U, dim> permutation_matrix(const std::array<std::size_t, dim> &indices) {
    return Matrix<U, dim>(indices);
}

template<typename U>
std::tuple<U, U, U, U, U, U>
calculate_ellipse_coefficients(U h, U nx, U ny, U nz, U sx, U sy, U sz) {
    // GARRY J. TEE, "SURFACE AREA OF ELLIPSOID SEGMENT", 2005 (Eq. 51)
    // h : normalized distance from the center of the ellipsoid to the plane
    // nx, ny, nz : normalized normal vector of the plane, paper uses λ, μ, ν
    // sx, sy, sz : normalized semi-axes of the ellipsoid, paper uses a, b, c
    U A = nz * nz * sz * sz / sx / sx + nx * nx;
    U B = 2 * nx * ny;
    U C = nz * nz * sz * sz / sy / sy + ny * ny;
    U D = -2 * h * nx;
    U E = -2 * h * ny;
    U G = h * h - nz * nz * sz * sz;
    return std::make_tuple(A, B, C, D, E, G);
}

template<typename U>
std::tuple<U, U, Eigen::Vector3<U>>
compute_ellipse_properties(U A, U B, U C, U D, U E, U G, U h,
                           Eigen::Vector3<U> n_pqr) {
    // Eigenvalues of the ellipse
    U common   = std::sqrt((A - C) * (A - C) + B * B);
    U eig_val1 = (A + C + common) / 2;
    U eig_val2 = (A + C - common) / 2;

    // Semi-axes of PROJECTED x-y ellipse
    U f_const = -G
              + (A * std::pow((B * E - 2 * C * D), 2)
                 + B * (B * E - 2 * C * D) * (B * D - 2 * A * E)
                 + C * std::pow((B * D - 2 * A * E), 2))
                        / std::pow((4 * A * C - B * B), 2);

    // eig_val1 > eig_val2, so d1(eig_val2) > d2(eig_val1)
    U d1 = std::sqrt(f_const / eig_val2);
    U d2 = std::sqrt(f_const / eig_val1);

    // Center of ellipse
    U                 p0_pqr = (B * E - 2 * C * D) / (4 * A * C - B * B);
    U                 q0_pqr = (B * D - 2 * A * E) / (4 * A * C - B * B);
    Eigen::Vector3<U> center_pqr;
    center_pqr << p0_pqr, q0_pqr,
            (h - n_pqr[0] * p0_pqr - n_pqr[1] * q0_pqr) / n_pqr[2];


    return std::make_tuple(d1, d2, center_pqr);
}

template<typename U>
U calculate_theta(U A, U B, U C, U D, U f) {
    // https://math.stackexchange.com/questions/3622888/area-and-center-location-of-an-ellipse-generated-by-the-intersection-of-an-ellip
    if (B != 0) {
        return std::atan2((C - A - std::sqrt(B * B + (A - C) * (A - C))), B);
    } else {
        if (A < C) {
            return 0.;
        } else if (A > C) {
            return M_PI / 2;
        } else if (A == C) {
            if (std::abs(D) >= std::abs(f)) {
                return 0.;
            } else {
                return M_PI / 2;
            }
        }
    }
}

// rotation matrix
template<typename U, std::size_t dim>
Matrix<U, dim> rotation_matrix(std::array<std::array<U, dim>, dim> values) {
    return Matrix<U, dim>(values);
}

template<typename U>
std::tuple<U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U, U,
           Eigen::Matrix3<int>, Eigen::Matrix3<U>>
ellipse_of_intersection(U l, U m, U n, U f, U a, U b, U c) {
    // GARRY J. TEE, "SURFACE AREA OF ELLIPSOID SEGMENT", 2005 (Section 3.1)

    // Establish normal vector
    Eigen::Vector3<U> normal(l, m, n);

    // Compute the norm
    U k = normal.norm();

    // Establish semi-axes of ellipsoid
    Eigen::Vector3<U> s_xyz(a, b, c);

    // Normalize normal vector, given sign of plane offset
    Eigen::Vector3<U> n_xyz
            = (f >= 0) ? (normal / k).eval() : (-normal / k).eval();
    U h = (f >= 0) ? f / k : -f / k;

    // Find index of the maximum absolute value
    Eigen::Vector3<U>  abs_n_xyz = n_xyz.cwiseAbs();
    std::array<int, 3> indices{0, 1, 2};// indices of {0, 1, 2}
    int                max_idx;
    abs_n_xyz.maxCoeff(&max_idx);

    // Swap to make max coefficient last
    std::swap(indices[2], indices[max_idx]);

    // Frame transformation matrix
    Eigen::PermutationMatrix<3, 3> P;
    P.indices() = Eigen::Vector3i(indices[0], indices[1], indices[2]);


    std::cout << "indices[0] = " << indices[0] << std::endl;
    std::cout << "indices[1] = " << indices[1] << std::endl;
    std::cout << "indices[2] = " << indices[2] << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // Transform to new coordinate system (swap axes)
    ////////////////////////////////////////////////////////////////////////////
    // Reorder normal vector and semi-axes for the largest component to be v
    Eigen::Vector3<U> n_pqr{n_xyz[indices[0]], n_xyz[indices[1]],
                            n_xyz[indices[2]]};
    Eigen::Vector3<U> s_pqr{s_xyz[indices[0]], s_xyz[indices[1]],
                            s_xyz[indices[2]]};


    // Calculate ellipse coefficients (3.1, Eq. 51, 52)
    auto [A, B, C, D, E, G] = calculate_ellipse_coefficients(
            h, n_pqr[0], n_pqr[1], n_pqr[2], s_pqr[0], s_pqr[1], s_pqr[2]);

    // Eigenvalues of the ellipse
    auto [d1, d2, center_pqr]
            = compute_ellipse_properties(A, B, C, D, E, G, h, n_pqr);

    // Calculate intersection ellipse coefficients
    // Appendix: Circumference of Ellipse of Intersection (Eq. 125, 126)
    auto A_h                        = A;
    auto B_h                        = B * n_pqr[2];
    auto C_h                        = C * n_pqr[2] * n_pqr[2];
    auto D_h                        = D;
    auto E_h                        = E * n_pqr[2];
    auto G_h                        = G;
    auto [d1_h, d2_h, center_pqr_h] = compute_ellipse_properties(
            A_h, B_h, C_h, D_h, E_h, G_h, h, n_pqr);

    auto theta_i = calculate_theta(A_h, B_h, C_h, D_h, f);
    auto theta_0 = calculate_theta(A, B, C, D, f);

    ////////////////////////////////////////////////////////////////////////////
    // Transform back to original frame
    ////////////////////////////////////////////////////////////////////////////
    auto center_xyz = P.transpose() * center_pqr;
    //    auto center_xyz_h = P.transpose() * center_pqr_h;

    // calculate the pole of the polar plane
    U px, py, pz;
    if (std::abs(f) < 1e-6) {
        px = s_xyz[0] * s_xyz[0] * n_xyz[0];
        py = s_xyz[1] * s_xyz[1] * n_xyz[1];
        pz = s_xyz[2] * s_xyz[2] * n_xyz[2];
    } else {
        px = s_xyz[0] * s_xyz[0] * n_xyz[0] / f;
        py = s_xyz[1] * s_xyz[1] * n_xyz[1] / f;
        pz = s_xyz[2] * s_xyz[2] * n_xyz[2] / f;
    }

    // calculate transformation matrix from xyz to XYZ
    // qx, qy, qz are the coordinates of the closest point from origin to plane
    auto qx = l * f / k;
    auto qy = m * f / k;
    auto qz = n * f / k;

    //    auto omega = std::atan2(qy, qx);
    //    auto phi   = std::atan2(std::sqrt(qx * qx + qy * qy), qz);

    //    auto c_phi = std::cos(phi);
    //    auto s_phi = std::sin(phi);
    //    auto c_omg = std::cos(omega);
    //    auto s_omg = std::sin(omega);
    //    auto R     = Matrix<U, 3>({
    //            { c_phi * c_omg,  s_phi * c_omg, -s_phi * c_omg},
    //            {-c_phi * s_omg, -s_phi * s_omg,  s_phi * s_omg},
    //            {         s_phi,              0,          c_phi}
    //    });
    //    P          = P * R;
    // define origin
    //    U x = 0, y = 0, z = h / nu;

    // First create an initial estimate for eta_axis
    Eigen::Vector3<U> eta_axis = P * Eigen::Vector3<U>::UnitX();

    // If eta_axis and n_xyz are not orthogonal, adjust eta_axis
    if (eta_axis.dot(n_xyz) != 0) {
        Eigen::Vector3<U> orthogonal_component = n_xyz - n_xyz.dot(eta_axis) * eta_axis;
        orthogonal_component.normalize();
        eta_axis = orthogonal_component;
    }

    // Now create zeta_axis orthogonal to eta_axis and n_xyz
    Eigen::Vector3<U> zeta_axis = n_xyz.cross(eta_axis);
    zeta_axis.normalize();

    // Construct the η-axis, which is orthogonal to the zeta axis and the plane's normal vector
    eta_axis = zeta_axis.cross(n_xyz);
    eta_axis.normalize();

    // Form the rotation matrix
    Eigen::Matrix3<U> R;
    R.col(0) = eta_axis;
    R.col(1) = zeta_axis;
    R.col(2) = n_xyz;
    return std::make_tuple(A, B, C, D, E, G, d1, d2, d1_h, d2_h, theta_0,
                           theta_i, center_xyz[0], center_xyz[1], center_xyz[2],
                           px, py, pz, P, R);
}


}// namespace odin::algorithms

#endif//ODIN_ELLIPSOID_PLANE_INTERSECTION_H