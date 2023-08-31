//
//  MetalAdder.cpp
//
//  Created by Srimukh Sripada on 04.12.21.
//

#include "metal_adder.hpp"

#include <iostream>

const unsigned int array_length = 1 << 5;
const unsigned int buffer_size  = array_length * sizeof(float);


void metal_adder::init_with_device(MTL::Device* device) {
    m_device = device;
    NS::Error* error;
    auto       default_library = m_device->newDefaultLibrary();

    if (!default_library) {
        std::cerr << "Failed to load default library.";
        // NOTES:
        // https://stackoverflow.com/questions/36204360/metal-default-library-not-found
        // Metal only looks in the app’s main bundle for default.metallib when you call newDefaultLibrary.
        std::exit(-1);
    }

    auto function_name = NS::String::string("add_arrays", NS::ASCIIStringEncoding);
    auto add_function  = default_library->newFunction(function_name);

    if (!add_function) { std::cerr << "Failed to find the adder function."; }

    m_add_function_pso = m_device->newComputePipelineState(add_function, &error);
    m_command_queue    = m_device->newCommandQueue();
};


void metal_adder::prepare_data() {
    m_buffer_A      = m_device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
    m_buffer_B      = m_device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
    m_buffer_result = m_device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);

    generate_random_float_data(m_buffer_A);
    generate_random_float_data(m_buffer_B);
}


void metal_adder::generate_random_float_data(MTL::Buffer* buffer) {
    auto* data_ptr = (float*) buffer->contents();
    for (unsigned long index = 0; index < array_length; index++) {
        data_ptr[index] = (float) rand() / (float) (RAND_MAX);
    }
}

void metal_adder::send_compute_command() {
    MTL::CommandBuffer* command_buffer = m_command_queue->commandBuffer();
    //    assert(command_buffer != nullptr);
    MTL::ComputeCommandEncoder* compute_encoder = command_buffer->computeCommandEncoder();
    encode_add_command(compute_encoder);
    compute_encoder->endEncoding();
    //    MTL::CommandBufferStatus status = command_buffer->status();
    //    std::cout << status << std::endl;
    command_buffer->commit();
    command_buffer->waitUntilCompleted();


    verify_results();
}

void metal_adder::encode_add_command(MTL::ComputeCommandEncoder* compute_encoder) const {
    compute_encoder->setComputePipelineState(m_add_function_pso);
    compute_encoder->setBuffer(m_buffer_A, 0, 0);
    compute_encoder->setBuffer(m_buffer_B, 0, 1);
    compute_encoder->setBuffer(m_buffer_result, 0, 2);

    MTL::Size grid_size = MTL::Size(array_length, 1, 1);

    NS::UInteger thread_group_size_ = m_add_function_pso->maxTotalThreadsPerThreadgroup();
    if (thread_group_size_ > array_length) { thread_group_size_ = array_length; }

    MTL::Size thread_group_size = MTL::Size(thread_group_size_, 1, 1);

    compute_encoder->dispatchThreads(grid_size, thread_group_size);
}


void metal_adder::prepare_data_with_vectors(const std::vector<float> &dataA, const std::vector<float> &dataB) {
    // Ensure the Metal buffers have been created
    m_buffer_A      = m_device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
    m_buffer_B      = m_device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
    m_buffer_result = m_device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);

    // Copy data from std::vector to Metal buffer
    std::memcpy(m_buffer_A->contents(), dataA.data(), buffer_size);
    std::memcpy(m_buffer_B->contents(), dataB.data(), buffer_size);
}

void metal_adder::verify_results() const {
    auto a      = (float*) m_buffer_A->contents();
    auto b      = (float*) m_buffer_B->contents();
    auto result = (float*) m_buffer_result->contents();

    for (unsigned long index = 0; index < array_length; index++) {
        if (result[index] != (a[index] + b[index])) {
            std::cout << "Compute ERROR: index=" << index << "result=" << result[index] << "vs "
                      << a[index] + b[index] << "=a+b\n";
            assert(result[index] == (a[index] + b[index]));
        }
    }
    std::cout << "Compute results as expected\n";
}