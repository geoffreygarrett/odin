

#include <cmath>
#include <string>

template<typename Derived, typename Float = double>
class BaseCelestialBody {
protected:
    Float       mu_central_body;
    Float       mu_self;
    Float       radius;
    Float       safe_radius;
    std::string name;

public:
    BaseCelestialBody(Float mu_cb, Float mu_s, Float r, Float sr, std::string n)
        : mu_central_body(mu_cb), mu_self(mu_s), radius(r), safe_radius(sr), name(std::move(n)) {}

    std::tuple<Float, Float> eph(Float when) {
        return static_cast<Derived *>(this)->eph(when);
    }

    Float compute_period() {
        return 2 * M_PI * std::sqrt(std::pow(a, 3) / mu_self);
    }

    std::string human_readable_extra() {
        return static_cast<Derived *>(this)->human_readable_extra();
    }

    std::string get_name() {
        return name;
    }

    Float get_mu_central_body() {
        return mu_central_body;
    }
    Float get_mu_self() {
        return mu_self;
    }
    Float get_radius() {
        return radius;
    }
    Float get_safe_radius() {
        return safe_radius;
    }
};
