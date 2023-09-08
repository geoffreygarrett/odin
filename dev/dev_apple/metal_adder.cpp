//
//  MetalAdder.cpp
//
//  Created by Srimukh Sripada on 04.12.21.
//

#include "metal_adder.hpp"

#include <iostream>


void metal_adder::init_with_device(MTL::Device* device) {
    m_device         = device;
    NS::Error* error = nullptr;

    auto default_library = m_device->newDefaultLibrary();
    //
    //    NS::String* default_path    = NS::String::string("default", NS::UTF8StringEncoding);
    //    auto        default_library = m_device->newLibrary(default_path, &error);

    //    NSURL *baseURL = [NSURL fileURLWithPath:@"file:///path/to/user/"];
    //    NSURL *URL = [NSURL URLWithString:@"folder/file.html" relativeToURL:baseURL];
    //    NSLog(@"absoluteURL = %@", [URL absoluteURL]);
    //
    //    auto base_url = NS::URL::URLWithString(NS::String::string("file:///path/to/user/", NS::UTF8StringEncoding));
    //    auto url = NS::URL::URLWithString(NS::String::string("folder/file.metallib", NS::UTF8StringEncoding), base_url);
    //    auto default_library = m_device->newLibrary(url, &error);

    // get cwd
    //    auto cwd = NS::CW

    //    auto default_path = NS::String::string("./default.metallib", NS::UTF8StringEncoding);
    //    auto default_path = NS::String::string("file:///Users/geoffreygarrett/CLionProjects/bor/bazel-bin/freyr/examples/example_metal_gpgpu.runfiles/bor/freyr/assets/shaders/metal/default.metallib", NS::UTF8StringEncoding);
    //    auto default_path = NS::String::string("file://./default.metallib", NS::UTF8StringEncoding);
    //    auto url = NS::URL::URLWithString(default_path);
    //    auto url = NS::URL::fileURLWithPath(default_path);
    //    auto default_library = m_device->newLibrary(url, &error);

    if (error) {
        std::cerr << "Failed to load default library. Detailed error information:" << std::endl;
        std::cerr << "  - Description: " << error->localizedDescription() << std::endl;
        std::cerr << "  - Recovery Options: " << error->localizedRecoveryOptions() << std::endl;
        std::cerr << "  - Recovery Suggestion: " << error->localizedRecoverySuggestion()
                  << std::endl;
        std::cerr << "  - Failure Reason: " << error->localizedFailureReason() << std::endl;
        std::cerr << "  - Error Code: " << error->code() << std::endl;
        std::cerr << "  - Error Domain: " << error->domain() << std::endl;
        std::cerr << "  - User Info: " << error->userInfo() << std::endl;

        // Additional notes or links for debugging
        std::cerr << "Notes:" << std::endl;
        std::cerr << "  - Metal only looks in the app’s main bundle for default.metallib when you "
                     "call newDefaultLibrary."
                  << std::endl;
        std::cerr << "  - Further Info: "
                     "https://stackoverflow.com/questions/36204360/metal-default-library-not-found"
                  << std::endl;

        std::exit(-1);
    } else if (!default_library) {
        std::cerr << "Failed to load default library. No error information available." << std::endl;
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
    assert(command_buffer != nullptr);
    MTL::ComputeCommandEncoder* compute_encoder = command_buffer->computeCommandEncoder();
    encode_add_command(compute_encoder);
    compute_encoder->endEncoding();
    MTL::CommandBufferStatus status = command_buffer->status();
    std::cout << status << std::endl;
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

    std::cout << "thread_group_size_ " << thread_group_size_ << std::endl;

    MTL::Size thread_group_size = MTL::Size(thread_group_size_, 1, 1);

    compute_encoder->dispatchThreads(grid_size, thread_group_size);
}

void metal_adder::prepare_data_with_vectors(
        const std::vector<float>& dataA, const std::vector<float>& dataB) {
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