template<typename derived, typename scalar>
class NoiseModel {
public:
    NoiseModel() : noise_level_(0.0) {}

    void update() {
        static_cast<derived *>(this)->update_impl();
    }

    scalar get_noise_value() const {
        return noise_level_;
    }

private:
    scalar noise_level_;
};

#include <random>


template<typename scalar>
class GaussianNoise : public NoiseModel<GaussianNoise<scalar>, scalar> {
public:
    // The Gaussian noise model is defined by mean and standard deviation.
    GaussianNoise(scalar mean, scalar std_dev)
        : distribution_(mean, std_dev), gen_(rd_()) {}

    void update_impl() {
        // Sample a new noise value from the Gaussian distribution.
        this->noise_level_ = distribution_(gen_);
    }


private:
    std::random_device rd_;
    std::mt19937 gen_;
    std::normal_distribution<double> distribution_;
};
