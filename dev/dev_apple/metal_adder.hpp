//
//  MetalAdder.hpp
//  metalnet
//
//  Created by Srimukh Sripada on 04.12.21.
//

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

const unsigned int array_length = 1 << 25;
const unsigned int buffer_size  = array_length * sizeof(float);


class metal_adder {
public:
    MTL::Device               *m_device;
    MTL::ComputePipelineState *m_add_function_pso;
    MTL::CommandQueue         *m_command_queue;

    MTL::Buffer *m_buffer_A;
    MTL::Buffer *m_buffer_B;
    MTL::Buffer *m_buffer_result;

    void init_with_device(MTL::Device *);
    void prepare_data();
    void send_compute_command();
    void prepare_data_with_vectors(const std::vector<float> &dataA, const std::vector<float> &dataB);


private:
    static void generate_random_float_data(MTL::Buffer *buffer);
    void encode_add_command(MTL::ComputeCommandEncoder *compute_encoder) const;
    void verify_results() const;
};