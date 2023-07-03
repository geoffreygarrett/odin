#include <Eigen/Dense>
#include <cmath>

/**
 * @brief A class representing the Doppler shift model
 *
 * This model is used to calculate the frequency ratio
 * from a transmitting spacecraft to a receiving station,
 * taking into account the velocities and positions of both,
 * as well as the direction of signal propagation.
 *
 * @tparam ScalarType The scalar numeric type to be used in calculations
 * @tparam Dim The dimension of the space
 */
template<typename ScalarType, int Dim>
class DopplerShiftModel {
    using VectorType = Eigen::Matrix<ScalarType, Dim, 1>;
public:
    /**
     * @brief Construct a new Doppler Shift Model object
     *
     * @param c Speed of light
     */
    DopplerShiftModel(ScalarType c) : c_(c) {}

    /**
     * @brief Compute the frequency ratio as per the model
     *
     * @param v_t Velocity vector of the transmitter
     * @param v_r Velocity vector of the receiver
     * @param e Unit vector in the direction of the signal propagation
     * @param U_t Newtonian potential at the transmitter
     * @param U_r Newtonian potential at the receiver
     *
     * @return ScalarType Computed frequency ratio
     */
    ScalarType compute_frequency_ratio(const VectorType &v_t,
                                       const VectorType &v_r,
                                       const VectorType &e,
                                       ScalarType U_t,
                                       ScalarType U_r) {
        ScalarType v_t_dot_e = v_t.dot(e);
        ScalarType v_r_dot_e = v_r.dot(e);
        ScalarType v_t_sq = v_t.squaredNorm();
        ScalarType v_r_sq = v_r.squaredNorm();

        ScalarType ratio = (1 - v_r_dot_e / c_ + U_r / (c_ * c_) + v_r_sq / (2 * c_ * c_)) /
                           (1 - v_t_dot_e / c_ + U_t / (c_ * c_) + v_t_sq / (2 * c_ * c_));
        return ratio;
    }

    /**
     * @brief Compute the Jacobian of the model
     *
     * @param v_t Velocity vector of the transmitter
     * @param v_r Velocity vector of the receiver
     * @param e Unit vector in the direction of the signal propagation
     * @param U_t Newtonian potential at the transmitter
     * @param U_r Newtonian potential at the receiver
     *
     * @return VectorType Jacobian vector
     */
    vector_type compute_jacobian(const vector_type &v_t, const vector_type &v_r,
                                 const vector_type &e, scalar_type U_t, scalar_type U_r) {
        // As the calculation of the Jacobian for this specific function might be
        // complex due to its non-linearity, an analytical solution is not presented here.
        // Normally, you would calculate the partial derivatives of the function
        // with respect to each of its inputs and form a vector (or matrix, if the output is a vector).
        // For this, symbolic math software or automatic differentiation might be useful.
        // This is a placeholder.
        return vector_type::Zero();
    }

private:
    scalar_type c_;  // Speed of light
};

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

using namespace autodiff;

template<typename ScalarType>
class FrequencyRatioModel {
    using scalar_type = ScalarType;
    using vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using matrix_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
    using adiff_type = forward::real<scalar_type>;
    using adiff_vector_type = Eigen::Matrix<adiff_type, Eigen::Dynamic, 1>;

public:
    scalar_type compute_frequency_ratio(const vector_type &v_t,
                                        const vector_type &v_r,
                                        const vector_type &e,
                                        scalar_type U_t,
                                        scalar_type U_r) {
        scalar_type v_t_dot_e = v_t.dot(e);
        scalar_type v_r_dot_e = v_r.dot(e);
        scalar_type v_t_sq = v_t.squaredNorm();
        scalar_type v_r_sq = v_r.squaredNorm();
        scalar_type ratio = (1 - v_r_dot_e / c_ + U_r / (c_ * c_) + v_r_sq / (2 * c_ * c_)) /
                            (1 - v_t_dot_e / c_ + U_t / (c_ * c_) + v_t_sq / (2 * c_ * c_));
        return ratio;
    }

    vector_type compute_jacobian(const vector_type &v_t,
                                 const vector_type &v_r,
                                 const vector_type &e,
                                 scalar_type U_t,
                                 scalar_type U_r) {

        adiff_vector_type v_t_ad = v_t.template cast<adiff_type>();
        adiff_vector_type v_r_ad = v_r.template cast<adiff_type>();
        adiff_vector_type e_ad = e.template cast<adiff_type>();
        adiff_type U_t_ad = U_t;
        adiff_type U_r_ad = U_r;

        v_t_ad = forward::make_variable(v_t_ad);
        v_r_ad = forward::make_variable(v_r_ad);
        e_ad = forward::make_variable(e_ad);
        U_t_ad = forward::make_variable(U_t_ad);
        U_r_ad = forward::make_variable(U_r_ad);

        adiff_type ratio_ad = compute_frequency_ratio_autodiff(v_t_ad, v_r_ad, e_ad, U_t_ad, U_r_ad, c_);

        vector_type derivatives(v_t.size() + v_r.size() + e.size() + 2);
        int i = 0;

        for (int j = 0; j < v_t.size(); ++j, ++i)
            derivatives(i) = v_t_ad(j).derivative();

        for (int j = 0; j < v_r.size(); ++j, ++i)
            derivatives(i) = v_r_ad(j).derivative();

        for (int j = 0; j < e.size(); ++j, ++i)
            derivatives(i) = e_ad(j).derivative();

        derivatives(i++) = U_t_ad.derivative();
        derivatives(i++) = U_r_ad.derivative();

        return derivatives;
    }

private:
    adiff_type compute_frequency_ratio_autodiff(const adiff_vector_type &v_t,
                                                const adiff_vector_type &v_r,
                                                const adiff_vector_type &e,
                                                adiff_type U_t,
                                                adiff_type U_r,
                                                adiff_type c_) {
        adiff_type v_t_dot_e = v_t.dot(e);
        adiff_type v_r_dot_e = v_r.dot(e);
        adiff_type v_t_sq = v_t.squaredNorm();
        adiff_type v_r_sq = v_r.squaredNorm();
        adiff_type ratio = (1 - v_r_dot_e / c_ + U_r / (c_ * c_) + v_r_sq / (2 * c_ * c_)) /
                           (1 - v_t_dot_e / c_ + U_t / (c_ * c_) + v_t_sq / (2 * c_ * c_));
        return ratio;
    }

    scalar_type c_;
};
