#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <odin/parallel/series.hpp>

namespace {

    // Mock objects
    struct MockModel {
        int value;

        explicit MockModel(int v) : value(v) {}

        // thread local copy
        [[nodiscard]] MockModel thread_local_copy() const {
            return MockModel(value);
        }
    };

    using MockIMatrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;
    using MockOMatrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

    struct BaseCalculateFn {
        virtual int operator()(MockModel &model, const Eigen::RowVectorXi &var_set) const {
            return 0;
        }
        virtual ~BaseCalculateFn() = default;
    };

    struct MockCalculateFn : BaseCalculateFn {
        int operator()(MockModel &model, const Eigen::RowVectorXi &var_set) const override {
            return model.value * var_set.sum();
        }
    };

    struct MockCalculateFnSideEffect : BaseCalculateFn {
        mutable std::atomic<int> increment = 0;

        int operator()(MockModel &model, const Eigen::RowVectorXi &var_set) const override {
            return (model.value + increment.fetch_add(1, std::memory_order_relaxed)) * var_set.sum();
        }
    };

    class SeriesCalculatorTest : public ::testing::TestWithParam<std::tuple<int, int, int>> {
    protected:
        MockModel   model;
        MockIMatrix input;
        MockOMatrix output;

        SeriesCalculatorTest() : model(std::get<2>(GetParam())) {}

        void SetUp() override {
            int rows = std::get<0>(GetParam());
            int cols = std::get<1>(GetParam());

            input  = MockIMatrix::Random(rows, cols);
            output = MockOMatrix::Zero(rows, cols == 1 ? 1 : cols);
        }
    };

    TEST_P(SeriesCalculatorTest, CalculateInSeriesMockCalculateFn) {
        MockCalculateFn calculation_func;
        odin::calculate_in_series(model, input, calculation_func, output);
        for (int i = 0; i < output.rows(); ++i) {
            for (int j = 0; j < output.cols(); ++j) {
                ASSERT_EQ(output(i, j), calculation_func(model, input.row(i)));
            }
        }
    }

//    TEST_P(SeriesCalculatorTest, CalculateInSeriesMockCalculateFnSideEffect) {
//        MockCalculateFnSideEffect calculation_func;
//        odin::calculate_in_series(model, input, calculation_func, output);
//        for (int i = 0; i < output.rows(); ++i) {
//            for (int j = 0; j < output.cols(); ++j) {
//                ASSERT_EQ(output(i, j), calculation_func(model, input.row(i)));
//            }
//        }
//    }

    INSTANTIATE_TEST_SUITE_P(
            SeriesCalculatorTestParams,
            SeriesCalculatorTest,
            ::testing::Combine(
                    ::testing::Values(1, 2, 10, 100),// rows
                    ::testing::Values(1, 2, 10, 100),// cols
                    ::testing::Values(0, 1, -1, 10)  // model values
                    ));

}// namespace