// File: ellipsoid_plane_intersection.metal
#include <metal_stdlib>
using namespace metal;

// No template in MSL, assume float3 input
float3 generate_orthogonal_vector(float3 vi) {
    float3 vo;
    float3 unit[3];
    unit[0] = float3(1.0f, 0.0f, 0.0f);
    unit[1] = float3(0.0f, 1.0f, 0.0f);
    unit[2] = float3(0.0f, 0.0f, 1.0f);

    for (int i = 0; i < 3; ++i) {
        vo = cross(unit[i], vi);
        if (length(vo) > FLT_EPSILON) { break; }
    }
    return vo;
}

// No template in MSL, assume float3 input
float3 calculate_polar_plane_pole(
        float3 normal, float offset, float3 semi_axes) {
    float3 n = normalize(normal);
    if (1.0f / offset > FLT_EPSILON) {
        return (semi_axes * semi_axes * n) / offset;
    } else {
        return semi_axes * semi_axes * n;
    }
}