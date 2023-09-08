//
//  main.cpp
//
//  Created by Srimukh Sripada on 03.12.21.
//  Edited by Geoffrey Garrett on 28.08.23.
//      - Bazel build system
//
#include <vector>// <--- Add this line
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <chrono>// For time measurement
#include <iostream>
#include <vector>

#include "metal_adder.hpp"

void add_arrays(const float* inA, const float* inB, float* result, int length) {
    for (int index = 0; index < length; index++) { result[index] = inA[index] + inB[index]; }
}

int main(int argc, const char* argv[]) {
    NS::AutoreleasePool* p_pool = NS::AutoreleasePool::alloc()->init();
    MTL::Device*         device = MTL::CreateSystemDefaultDevice();
    auto*                adder  = new metal_adder();
    adder->init_with_device(device);

//    const int array_length = 1 << 30;

    std::vector<float> dataA(array_length, 1.0f);
    std::vector<float> dataB(array_length, 2.0f);
    std::vector<float> result_cpu(array_length);

    // Initialize vectors
    for (int index = 0; index < array_length; index++) {
        dataA[index] = static_cast<float>(index);
        dataB[index] = static_cast<float>(index);
    }

    // Time the CPU-based array addition
    auto start_cpu = std::chrono::high_resolution_clock::now();
    add_arrays(dataA.data(), dataB.data(), result_cpu.data(), array_length);
    auto                                      end_cpu = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed_cpu = end_cpu - start_cpu;

    // Time the GPU-based array addition
    adder->prepare_data_with_vectors(dataA, dataB);
    auto start_gpu = std::chrono::high_resolution_clock::now();
    adder->send_compute_command();
    auto                                      end_gpu = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed_gpu = end_gpu - start_gpu;

    // Print time comparison
    std::cout << "Time taken by CPU-based addition: " << elapsed_cpu.count() << " ms\n";
    std::cout << "Time taken by GPU-based addition: " << elapsed_gpu.count() << " ms\n";

    auto result_ptr = (float*) adder->m_buffer_result->contents();

    // Verify results
    for (int index = 0; index < array_length; index++) {
        if (result_ptr[index] != result_cpu[index]) {
            std::cout << "Error: GPU result does not match CPU result at index " << index << "\n";
            std::cout << "GPU result: " << result_ptr[index] << "\n";
            std::cout << "CPU result: " << result_cpu[index] << "\n";
            return 1;
        }
    }
    // You can uncomment the next lines to print the results if you wish
//    /*
//    for (int index = 0; index < array_length; index++) {
//        std::cout << result_ptr[index] << " ";
//    }
//    */

    std::cout << "Execution finished.\n";
    p_pool->release();

    delete adder;

    return 0;
}