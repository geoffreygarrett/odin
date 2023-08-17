#include <Eigen/Dense>
#include <iostream>
#include <odin/logging.hpp>
#include <odin/optimization/tree/mcts.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/unordered_map.hpp>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_ellint.h>

#include <oneapi/tbb/version.h>
#include <odin/models/shape/algorithms.h>

using Eigen::MatrixXd;

// C++ includes
#include <iostream>

// autodiff include
#include <autodiff/forward/dual.hpp>
using namespace autodiff;
//
//// The single-variable function for which derivatives are needed
dual f(dual x) {
    return 1 + x + x * x + 1 / x + log(x);
}

void main2() {
    dual x = 2.0; // the input variable x
    dual u = f(x);// the output variable u

    double dudx = derivative(f, wrt(x), at(x));// evaluate the derivative du/dx

    std::cout << "u = " << u << std::endl;       // print the evaluated output u
    std::cout << "du/dx = " << dudx << std::endl;// print the evaluated derivative du/dx
}
#include <chrono>
#include <iostream>
#include <string>
#include <thread>

// Constants for the loading bar
const int         BAR_WIDTH  = 50;
const std::string BAR_CHAR   = "=";
const std::string SPACE_CHAR = " ";

//// Function to draw a loading bar
//void DrawProgressBar(float progress) {
//    // Calculate the number of characters in the bar
//    int pos = BAR_WIDTH * progress;
//
//    // Initialize the bar with spaces
//    std::string bar(BAR_WIDTH, SPACE_CHAR[0]);
//
//    // Replace the spaces with the bar character up to the current position
//    for (int i = 0; i < pos; ++i) {
//        bar[i] = BAR_CHAR[0];
//    }
//
//    // Colorize the bar
//    bar = Colorize(bar, ANSI_COLOR_GREEN);
//
//    // Print the progress bar
//    std::cout << "[" << bar << "] " << int(progress * 100.0) << " %\r";
//    std::cout << std::flush;
//}

int main(int argc, char *argv[]) {
    //    INIT_ODIN_LOGGING(argv[0], "./log/odin.log");

    std::cout << std::endl;
    std::cout << "Hello from oneTBB "
              << TBB_VERSION_MAJOR << "."
              << TBB_VERSION_MINOR << "."
              << TBB_VERSION_PATCH
              << "!" << std::endl;

    main2();
    std::cout << "Hello, Worlds!" << std::endl;
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << m << std::endl;

    ODIN_LOG_INFO << "Hello, this is information!";
    ODIN_LOG_WARNING << "Hello, this is a warning!";
    ODIN_LOG_ERROR << "Hello, this is an error!";

    for (int i = 0; i < 10; ++i) {
        ODIN_LOG_INFO << "Hello, this is information!";
        ODIN_LOG_WARNING << "Hello, this is a warning!";
        ODIN_LOG_ERROR << "Hello, this is an error!";
    }
    //    for (int i = 0; i <= 100; ++i) {
    //        DrawProgressBar(i/100.0);
    //        std::this_thread::sleep_for(std::chrono::milliseconds(100));  // Simulate work
    //    }
    std::cout << std::endl;

    std::cout << "\033[34m"
              << "This text is blue"
              << "\033[0m" << std::endl;

    std::cout << sizeof(long double) << std::endl;
    return 0;
}
