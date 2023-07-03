
#include <chrono>
#include <include/gtest/gtest.h>
#include <random>
#include <tbb/parallel_for.h>
#include <tbb/task.h>

std::vector<int> square_vector_serial(const std::vector<int> &input) {
    std::vector<int> result(input.size());
    for (size_t i = 0; i < input.size(); ++i)
        result[i] = input[i] * input[i];
    return result;
}

std::vector<int> square_vector_parallel(const std::vector<int> &input) {
    std::vector<int> result(input.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                      [&](const tbb::blocked_range<size_t> &r) {
                          for (size_t i = r.begin(); i != r.end(); ++i)
                              result[i] = input[i] * input[i];
                      });
    return result;
}

std::vector<int64_t> factorial_vector_serial(const std::vector<int> &input) {
    std::vector<int64_t> result(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        int64_t factorial = 1;
        for (int j = 1; j <= input[i]; ++j)
            factorial *= j;
        result[i] = factorial;
    }
    return result;
}

std::vector<int64_t> factorial_vector_parallel(const std::vector<int> &input) {
    std::vector<int64_t> result(input.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                      [&](const tbb::blocked_range<size_t> &r) {
                          for (size_t i = r.begin(); i != r.end(); ++i) {
                              int64_t factorial = 1;
                              for (int j = 1; j <= input[i]; ++j)
                                  factorial *= j;
                              result[i] = factorial;
                          }
                      });
    return result;
}

std::vector<int> generate_random_vector(size_t size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 100);

    std::vector<int> result(size);
    for (size_t i = 0; i < size; ++i)
        result[i] = dis(gen);
    return result;
}


TEST(TBBTest, CompareSerialAndParallel) {
    // Generate a large vector
    const size_t size = 10000000;
    std::vector<int> input = generate_random_vector(size);

    // Time and run the serial version
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> serial_result = square_vector_serial(input);
    auto stop = std::chrono::high_resolution_clock::now();
    auto serial_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cout << "Serial time: " << serial_duration << " ms\n";

    // Time and run the parallel version
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> parallel_result = square_vector_parallel(input);
    stop = std::chrono::high_resolution_clock::now();
    auto parallel_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cout << "Parallel time: " << parallel_duration << " ms\n";

    // Check that the results are the same
    EXPECT_EQ(serial_result, parallel_result);
}


TEST(TBBTest, CompareSerialAndParallelFactorial) {
    // Generate a large vector
    const size_t size = 1000000;// Be aware that a large size or large values could lead to integer overflow
    std::vector<int> input = generate_random_vector(size);

    // Time and run the serial version
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int64_t> serial_result = factorial_vector_serial(input);
    auto stop = std::chrono::high_resolution_clock::now();
    auto serial_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cout << "Serial time: " << serial_duration << " ms\n";

    // Time and run the parallel version
    start = std::chrono::high_resolution_clock::now();
    std::vector<int64_t> parallel_result = factorial_vector_parallel(input);
    stop = std::chrono::high_resolution_clock::now();
    auto parallel_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cout << "Parallel time: " << parallel_duration << " ms\n";

    // Check that the results are the same
    EXPECT_EQ(serial_result, parallel_result);
}


class FibonacciTask {
public:
    int n;           // input
    long long result;// output

    FibonacciTask(int n_) : n(n_), result(0) {}

    void run() {
        result = parallel_fibonacci(n);
    }

private:
    static long long parallel_fibonacci(int n) {
        if (n < 2) {
            return n;
        } else {
            long long x, y;
            tbb::task_group g;
            g.run([&] { x = parallel_fibonacci(n - 1); });// spawn a task
            g.run([&] { y = parallel_fibonacci(n - 2); });// spawn another task
            g.wait();                                     // wait for both tasks to complete
            return x + y;
        }
    }
};

long long serial_fibonacci(int n) {
    if (n < 2) {
        return n;
    } else {
        return serial_fibonacci(n - 1) + serial_fibonacci(n - 2);
    }
}

long long parallel_fibonacci(int n) {
    if (n < 2) {
        return n;
    } else {
        long long x, y;
        tbb::task_group g;
        g.run([&] { x = parallel_fibonacci(n - 1); });// spawn a task
        g.run([&] { y = parallel_fibonacci(n - 2); });// spawn another task
        g.wait();                                     // wait for both tasks to complete
        return x + y;
    }
}


TEST(TBBTest, ParallelFibonacci) {
    int n = 30;

    // Time and run the serial version
    auto start = std::chrono::high_resolution_clock::now();
    long long serial_result = serial_fibonacci(n);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Serial time: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms\n";

    // Time and run the parallel version
    start = std::chrono::high_resolution_clock::now();
    long long parallel_result = parallel_fibonacci(n);
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Parallel time: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms\n";

    EXPECT_EQ(serial_result, parallel_result);
}