#include "Eigen/Dense"
#include <pybind11/pybind11.h>

class CustomFilter : public Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>
{
public:
    using Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>::Filter;

    void initialize_impl(const Eigen::VectorXd &initial_state, const Eigen::MatrixXd &initial_covariance) override {
        PYBIND11_OVERRIDE(
                void,
                Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>,
                initialize_impl,
                initial_state, initial_covariance
        );
    }

    void predict_impl() override {
        PYBIND11_OVERRIDE_PURE(
                void,
                Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>,
                predict_impl
        );
    }

    void update_impl() override {
        PYBIND11_OVERRIDE_PURE(
                void,
                Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>,
                update_impl
        );
    }

    Eigen::VectorXd get_state_estimate_impl() const override {
        PYBIND11_OVERRIDE(
                Eigen::VectorXd,
                Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>,
                get_state_estimate_impl
        );
    }

    Eigen::MatrixXd get_covariance_estimate_impl() const override {
        PYBIND11_OVERRIDE(
                Eigen::MatrixXd,
                Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>,
                get_covariance_estimate_impl
        );
    }
};

// python wrapper
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "CustomFilter.h"

namespace py = pybind11;

PYBIND11_MODULE(your_module_name, m) {
py::class_<CustomFilter, Filter<CustomFilter, Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd>>(m, "CustomFilter")
.def(py::init<>())
.def("initialize", &CustomFilter::initialize)
.def("predict", &CustomFilter::predict_impl)
.def("update", &CustomFilter::update_impl)
.def("get_state_estimate", &CustomFilter::get_state_estimate_impl)
.def("get_covariance_estimate", &CustomFilter::get_covariance_estimate_impl);
}

// Custom python class
//from your_module_name import CustomFilter
//
//class MyFilter(CustomFilter):
//def predict(self, *args, **kwargs):
//# Your custom implementation here
//pass
//
//        def update(self, *args, **kwargs):
//# Your custom implementation here
//pass
