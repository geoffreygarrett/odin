#include <iostream>
#include <odin/eigen.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>

int main() {
    // Generate 32M random numbers serially.
    thrust::default_random_engine         rng(1337);
    thrust::uniform_int_distribution<int> dist;
    thrust::host_vector<int>              h_vec(32 << 20);

    // Using Thrust to generate numbers
    thrust::generate(h_vec.begin(), h_vec.end(), [&] { return dist(rng); });

    // Transfer data to the device.
    thrust::device_vector<int> d_vec = h_vec;

    // Sort data on the device (GPU).
    thrust::sort(d_vec.begin(), d_vec.end());

    // Transfer data back to host.
    thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());

    // Use Thrust to compute the sum of the numbers
    int sum = thrust::reduce(h_vec.begin(), h_vec.end());

    std::cout << "Sum: " << sum << std::endl;
}