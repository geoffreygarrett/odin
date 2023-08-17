#ifndef ODIN_MARCHING_CUBES_H
#define ODIN_MARCHING_CUBES_H

#include <odin/core/point.h>

#include <array>
#include <map>
#include <valarray>
#include <vector>

#include "../shape_concepts.h"
#include "face_vertex.h"
#include "face_vertex_impl.h"

w
namespace odin {

// TODO: Refactor odin::shape namespace for point.h
using namespace odin::shape;

struct Point3D {
    float x, y, z;
};

struct Triangle {
    Point3D vertices[3];
};

template<Point P>
struct Voxel {
    P     position{};
    bool  is_active{};
    float value{};
};

// http://paulbourke.net/geometry/polygonise/
constexpr std::array<std::array<int, 2>, 12> edge_vertex_indices = {
  // clang-format off
        {{{0, 1}},
         {{1, 2}},
         {{2, 3}},
         {{3, 0}},
         {{4, 5}},
         {{5, 6}},
         {{6, 7}},
         {{7, 4}},
         {{0, 4}},
         {{1, 5}},
         {{2, 6}},
         {{3, 7}}}  // clang-format on
};

// Marching Cubes lookup tables https://gist.github.com/dwilliamson/c041e3454a713e58baf6e4f8e5fffecd
constexpr std::array<int, 256> edge_table = {
        // clang-format off
        0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905,
        0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99, 0x393, 0x29a,
        0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93,
        0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
        0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9,
        0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6,
        0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f,
        0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc, 0xdfc, 0xcf5,
        0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a,
        0x256, 0x35f, 0x55, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
        0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
        0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9,
        0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc, 0x1c5, 0x2cf, 0x3c6,
        0x4ca, 0x5c3, 0x6c9, 0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f,
        0xf55, 0xe5c, 0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5,
        0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a,
        0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663,
        0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
        0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0, 0xd30, 0xc39,
        0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636,
        0x13a, 0x33, 0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f,
        0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605,
        0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0// clang-format on
};

constexpr std::array<std::array<int, 16>, 256> tri_table = {
        {// clang-format off
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  8,  3,  9,  8,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  3,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 9,  2, 10,  0,  2,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 2,  8,  3,  2, 10,  8, 10,  9,  8, -1, -1, -1, -1, -1, -1, -1},
        { 3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0, 11,  2,  8, 11,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  9,  0,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1, 11,  2,  1,  9, 11,  9,  8, 11, -1, -1, -1, -1, -1, -1, -1},
        { 3, 10,  1, 11, 10,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0, 10,  1,  0,  8, 10,  8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        { 3,  9,  0,  3, 11,  9, 11, 10,  9, -1, -1, -1, -1, -1, -1, -1},
        { 9,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  3,  0,  7,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  1,  9,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  1,  9,  4,  7,  1,  7,  3,  1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 10,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 3,  4,  7,  3,  0,  4,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1},
        { 9,  2, 10,  9,  0,  2,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1},
        { 2, 10,  9,  2,  9,  7,  2,  7,  3,  7,  9,  4, -1, -1, -1, -1},
        { 8,  4,  7,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11,  4,  7, 11,  2,  4,  2,  0,  4, -1, -1, -1, -1, -1, -1, -1},
        { 9,  0,  1,  8,  4,  7,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1},
        { 4,  7, 11,  9,  4, 11,  9, 11,  2,  9,  2,  1, -1, -1, -1, -1},
        { 3, 10,  1,  3, 11, 10,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1},
        { 1, 11, 10,  1,  4, 11,  1,  0,  4,  7, 11,  4, -1, -1, -1, -1},
        { 4,  7,  8,  9,  0, 11,  9, 11, 10, 11,  0,  3, -1, -1, -1, -1},
        { 4,  7, 11,  4, 11,  9,  9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        { 9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 9,  5,  4,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  5,  4,  1,  5,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 8,  5,  4,  8,  3,  5,  3,  1,  5, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 10,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 3,  0,  8,  1,  2, 10,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1},
        { 5,  2, 10,  5,  4,  2,  4,  0,  2, -1, -1, -1, -1, -1, -1, -1},
        { 2, 10,  5,  3,  2,  5,  3,  5,  4,  3,  4,  8, -1, -1, -1, -1},
        { 9,  5,  4,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0, 11,  2,  0,  8, 11,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1},
        { 0,  5,  4,  0,  1,  5,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1},
        { 2,  1,  5,  2,  5,  8,  2,  8, 11,  4,  8,  5, -1, -1, -1, -1},
        {10,  3, 11, 10,  1,  3,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1},
        { 4,  9,  5,  0,  8,  1,  8, 10,  1,  8, 11, 10, -1, -1, -1, -1},
        { 5,  4,  0,  5,  0, 11,  5, 11, 10, 11,  0,  3, -1, -1, -1, -1},
        { 5,  4,  8,  5,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1},
        { 9,  7,  8,  5,  7,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 9,  3,  0,  9,  5,  3,  5,  7,  3, -1, -1, -1, -1, -1, -1, -1},
        { 0,  7,  8,  0,  1,  7,  1,  5,  7, -1, -1, -1, -1, -1, -1, -1},
        { 1,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 9,  7,  8,  9,  5,  7, 10,  1,  2, -1, -1, -1, -1, -1, -1, -1},
        {10,  1,  2,  9,  5,  0,  5,  3,  0,  5,  7,  3, -1, -1, -1, -1},
        { 8,  0,  2,  8,  2,  5,  8,  5,  7, 10,  5,  2, -1, -1, -1, -1},
        { 2, 10,  5,  2,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1},
        { 7,  9,  5,  7,  8,  9,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1},
        { 9,  5,  7,  9,  7,  2,  9,  2,  0,  2,  7, 11, -1, -1, -1, -1},
        { 2,  3, 11,  0,  1,  8,  1,  7,  8,  1,  5,  7, -1, -1, -1, -1},
        {11,  2,  1, 11,  1,  7,  7,  1,  5, -1, -1, -1, -1, -1, -1, -1},
        { 9,  5,  8,  8,  5,  7, 10,  1,  3, 10,  3, 11, -1, -1, -1, -1},
        { 5,  7,  0,  5,  0,  9,  7, 11,  0,  1,  0, 10, 11, 10,  0, -1},
        {11, 10,  0, 11,  0,  3, 10,  5,  0,  8,  0,  7,  5,  7,  0, -1},
        {11, 10,  5,  7, 11,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  3,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 9,  0,  1,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  8,  3,  1,  9,  8,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1},
        { 1,  6,  5,  2,  6,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  6,  5,  1,  2,  6,  3,  0,  8, -1, -1, -1, -1, -1, -1, -1},
        { 9,  6,  5,  9,  0,  6,  0,  2,  6, -1, -1, -1, -1, -1, -1, -1},
        { 5,  9,  8,  5,  8,  2,  5,  2,  6,  3,  2,  8, -1, -1, -1, -1},
        { 2,  3, 11, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11,  0,  8, 11,  2,  0, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1},
        { 0,  1,  9,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1},
        { 5, 10,  6,  1,  9,  2,  9, 11,  2,  9,  8, 11, -1, -1, -1, -1},
        { 6,  3, 11,  6,  5,  3,  5,  1,  3, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8, 11,  0, 11,  5,  0,  5,  1,  5, 11,  6, -1, -1, -1, -1},
        { 3, 11,  6,  0,  3,  6,  0,  6,  5,  0,  5,  9, -1, -1, -1, -1},
        { 6,  5,  9,  6,  9, 11, 11,  9,  8, -1, -1, -1, -1, -1, -1, -1},
        { 5, 10,  6,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  3,  0,  4,  7,  3,  6,  5, 10, -1, -1, -1, -1, -1, -1, -1},
        { 1,  9,  0,  5, 10,  6,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1},
        {10,  6,  5,  1,  9,  7,  1,  7,  3,  7,  9,  4, -1, -1, -1, -1},
        { 6,  1,  2,  6,  5,  1,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2,  5,  5,  2,  6,  3,  0,  4,  3,  4,  7, -1, -1, -1, -1},
        { 8,  4,  7,  9,  0,  5,  0,  6,  5,  0,  2,  6, -1, -1, -1, -1},
        { 7,  3,  9,  7,  9,  4,  3,  2,  9,  5,  9,  6,  2,  6,  9, -1},
        { 3, 11,  2,  7,  8,  4, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1},
        { 5, 10,  6,  4,  7,  2,  4,  2,  0,  2,  7, 11, -1, -1, -1, -1},
        { 0,  1,  9,  4,  7,  8,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1},
        { 9,  2,  1,  9, 11,  2,  9,  4, 11,  7, 11,  4,  5, 10,  6, -1},
        { 8,  4,  7,  3, 11,  5,  3,  5,  1,  5, 11,  6, -1, -1, -1, -1},
        { 5,  1, 11,  5, 11,  6,  1,  0, 11,  7, 11,  4,  0,  4, 11, -1},
        { 0,  5,  9,  0,  6,  5,  0,  3,  6, 11,  6,  3,  8,  4,  7, -1},
        { 6,  5,  9,  6,  9, 11,  4,  7,  9,  7, 11,  9, -1, -1, -1, -1},
        {10,  4,  9,  6,  4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4, 10,  6,  4,  9, 10,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1},
        {10,  0,  1, 10,  6,  0,  6,  4,  0, -1, -1, -1, -1, -1, -1, -1},
        { 8,  3,  1,  8,  1,  6,  8,  6,  4,  6,  1, 10, -1, -1, -1, -1},
        { 1,  4,  9,  1,  2,  4,  2,  6,  4, -1, -1, -1, -1, -1, -1, -1},
        { 3,  0,  8,  1,  2,  9,  2,  4,  9,  2,  6,  4, -1, -1, -1, -1},
        { 0,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 8,  3,  2,  8,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1},
        {10,  4,  9, 10,  6,  4, 11,  2,  3, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  2,  2,  8, 11,  4,  9, 10,  4, 10,  6, -1, -1, -1, -1},
        { 3, 11,  2,  0,  1,  6,  0,  6,  4,  6,  1, 10, -1, -1, -1, -1},
        { 6,  4,  1,  6,  1, 10,  4,  8,  1,  2,  1, 11,  8, 11,  1, -1},
        { 9,  6,  4,  9,  3,  6,  9,  1,  3, 11,  6,  3, -1, -1, -1, -1},
        { 8, 11,  1,  8,  1,  0, 11,  6,  1,  9,  1,  4,  6,  4,  1, -1},
        { 3, 11,  6,  3,  6,  0,  0,  6,  4, -1, -1, -1, -1, -1, -1, -1},
        { 6,  4,  8, 11,  6,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 7, 10,  6,  7,  8, 10,  8,  9, 10, -1, -1, -1, -1, -1, -1, -1},
        { 0,  7,  3,  0, 10,  7,  0,  9, 10,  6,  7, 10, -1, -1, -1, -1},
        {10,  6,  7,  1, 10,  7,  1,  7,  8,  1,  8,  0, -1, -1, -1, -1},
        {10,  6,  7, 10,  7,  1,  1,  7,  3, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2,  6,  1,  6,  8,  1,  8,  9,  8,  6,  7, -1, -1, -1, -1},
        { 2,  6,  9,  2,  9,  1,  6,  7,  9,  0,  9,  3,  7,  3,  9, -1},
        { 7,  8,  0,  7,  0,  6,  6,  0,  2, -1, -1, -1, -1, -1, -1, -1},
        { 7,  3,  2,  6,  7,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 2,  3, 11, 10,  6,  8, 10,  8,  9,  8,  6,  7, -1, -1, -1, -1},
        { 2,  0,  7,  2,  7, 11,  0,  9,  7,  6,  7, 10,  9, 10,  7, -1},
        { 1,  8,  0,  1,  7,  8,  1, 10,  7,  6,  7, 10,  2,  3, 11, -1},
        {11,  2,  1, 11,  1,  7, 10,  6,  1,  6,  7,  1, -1, -1, -1, -1},
        { 8,  9,  6,  8,  6,  7,  9,  1,  6, 11,  6,  3,  1,  3,  6, -1},
        { 0,  9,  1, 11,  6,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 7,  8,  0,  7,  0,  6,  3, 11,  0, 11,  6,  0, -1, -1, -1, -1},
        { 7, 11,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 3,  0,  8, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  1,  9, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 8,  1,  9,  8,  3,  1, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1},
        {10,  1,  2,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 10,  3,  0,  8,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1},
        { 2,  9,  0,  2, 10,  9,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1},
        { 6, 11,  7,  2, 10,  3, 10,  8,  3, 10,  9,  8, -1, -1, -1, -1},
        { 7,  2,  3,  6,  2,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 7,  0,  8,  7,  6,  0,  6,  2,  0, -1, -1, -1, -1, -1, -1, -1},
        { 2,  7,  6,  2,  3,  7,  0,  1,  9, -1, -1, -1, -1, -1, -1, -1},
        { 1,  6,  2,  1,  8,  6,  1,  9,  8,  8,  7,  6, -1, -1, -1, -1},
        {10,  7,  6, 10,  1,  7,  1,  3,  7, -1, -1, -1, -1, -1, -1, -1},
        {10,  7,  6,  1,  7, 10,  1,  8,  7,  1,  0,  8, -1, -1, -1, -1},
        { 0,  3,  7,  0,  7, 10,  0, 10,  9,  6, 10,  7, -1, -1, -1, -1},
        { 7,  6, 10,  7, 10,  8,  8, 10,  9, -1, -1, -1, -1, -1, -1, -1},
        { 6,  8,  4, 11,  8,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 3,  6, 11,  3,  0,  6,  0,  4,  6, -1, -1, -1, -1, -1, -1, -1},
        { 8,  6, 11,  8,  4,  6,  9,  0,  1, -1, -1, -1, -1, -1, -1, -1},
        { 9,  4,  6,  9,  6,  3,  9,  3,  1, 11,  3,  6, -1, -1, -1, -1},
        { 6,  8,  4,  6, 11,  8,  2, 10,  1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 10,  3,  0, 11,  0,  6, 11,  0,  4,  6, -1, -1, -1, -1},
        { 4, 11,  8,  4,  6, 11,  0,  2,  9,  2, 10,  9, -1, -1, -1, -1},
        {10,  9,  3, 10,  3,  2,  9,  4,  3, 11,  3,  6,  4,  6,  3, -1},
        { 8,  2,  3,  8,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1},
        { 0,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  9,  0,  2,  3,  4,  2,  4,  6,  4,  3,  8, -1, -1, -1, -1},
        { 1,  9,  4,  1,  4,  2,  2,  4,  6, -1, -1, -1, -1, -1, -1, -1},
        { 8,  1,  3,  8,  6,  1,  8,  4,  6,  6, 10,  1, -1, -1, -1, -1},
        {10,  1,  0, 10,  0,  6,  6,  0,  4, -1, -1, -1, -1, -1, -1, -1},
        { 4,  6,  3,  4,  3,  8,  6, 10,  3,  0,  3,  9, 10,  9,  3, -1},
        {10,  9,  4,  6, 10,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  9,  5,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  3,  4,  9,  5, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1},
        { 5,  0,  1,  5,  4,  0,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11,  7,  6,  8,  3,  4,  3,  5,  4,  3,  1,  5, -1, -1, -1, -1},
        { 9,  5,  4, 10,  1,  2,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1},
        { 6, 11,  7,  1,  2, 10,  0,  8,  3,  4,  9,  5, -1, -1, -1, -1},
        { 7,  6, 11,  5,  4, 10,  4,  2, 10,  4,  0,  2, -1, -1, -1, -1},
        { 3,  4,  8,  3,  5,  4,  3,  2,  5, 10,  5,  2, 11,  7,  6, -1},
        { 7,  2,  3,  7,  6,  2,  5,  4,  9, -1, -1, -1, -1, -1, -1, -1},
        { 9,  5,  4,  0,  8,  6,  0,  6,  2,  6,  8,  7, -1, -1, -1, -1},
        { 3,  6,  2,  3,  7,  6,  1,  5,  0,  5,  4,  0, -1, -1, -1, -1},
        { 6,  2,  8,  6,  8,  7,  2,  1,  8,  4,  8,  5,  1,  5,  8, -1},
        { 9,  5,  4, 10,  1,  6,  1,  7,  6,  1,  3,  7, -1, -1, -1, -1},
        { 1,  6, 10,  1,  7,  6,  1,  0,  7,  8,  7,  0,  9,  5,  4, -1},
        { 4,  0, 10,  4, 10,  5,  0,  3, 10,  6, 10,  7,  3,  7, 10, -1},
        { 7,  6, 10,  7, 10,  8,  5,  4, 10,  4,  8, 10, -1, -1, -1, -1},
        { 6,  9,  5,  6, 11,  9, 11,  8,  9, -1, -1, -1, -1, -1, -1, -1},
        { 3,  6, 11,  0,  6,  3,  0,  5,  6,  0,  9,  5, -1, -1, -1, -1},
        { 0, 11,  8,  0,  5, 11,  0,  1,  5,  5,  6, 11, -1, -1, -1, -1},
        { 6, 11,  3,  6,  3,  5,  5,  3,  1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 10,  9,  5, 11,  9, 11,  8, 11,  5,  6, -1, -1, -1, -1},
        { 0, 11,  3,  0,  6, 11,  0,  9,  6,  5,  6,  9,  1,  2, 10, -1},
        {11,  8,  5, 11,  5,  6,  8,  0,  5, 10,  5,  2,  0,  2,  5, -1},
        { 6, 11,  3,  6,  3,  5,  2, 10,  3, 10,  5,  3, -1, -1, -1, -1},
        { 5,  8,  9,  5,  2,  8,  5,  6,  2,  3,  8,  2, -1, -1, -1, -1},
        { 9,  5,  6,  9,  6,  0,  0,  6,  2, -1, -1, -1, -1, -1, -1, -1},
        { 1,  5,  8,  1,  8,  0,  5,  6,  8,  3,  8,  2,  6,  2,  8, -1},
        { 1,  5,  6,  2,  1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  3,  6,  1,  6, 10,  3,  8,  6,  5,  6,  9,  8,  9,  6, -1},
        {10,  1,  0, 10,  0,  6,  9,  5,  0,  5,  6,  0, -1, -1, -1, -1},
        { 0,  3,  8,  5,  6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10,  5,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11,  5, 10,  7,  5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11,  5, 10, 11,  7,  5,  8,  3,  0, -1, -1, -1, -1, -1, -1, -1},
        { 5, 11,  7,  5, 10, 11,  1,  9,  0, -1, -1, -1, -1, -1, -1, -1},
        {10,  7,  5, 10, 11,  7,  9,  8,  1,  8,  3,  1, -1, -1, -1, -1},
        {11,  1,  2, 11,  7,  1,  7,  5,  1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  3,  1,  2,  7,  1,  7,  5,  7,  2, 11, -1, -1, -1, -1},
        { 9,  7,  5,  9,  2,  7,  9,  0,  2,  2, 11,  7, -1, -1, -1, -1},
        { 7,  5,  2,  7,  2, 11,  5,  9,  2,  3,  2,  8,  9,  8,  2, -1},
        { 2,  5, 10,  2,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1},
        { 8,  2,  0,  8,  5,  2,  8,  7,  5, 10,  2,  5, -1, -1, -1, -1},
        { 9,  0,  1,  5, 10,  3,  5,  3,  7,  3, 10,  2, -1, -1, -1, -1},
        { 9,  8,  2,  9,  2,  1,  8,  7,  2, 10,  2,  5,  7,  5,  2, -1},
        { 1,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  7,  0,  7,  1,  1,  7,  5, -1, -1, -1, -1, -1, -1, -1},
        { 9,  0,  3,  9,  3,  5,  5,  3,  7, -1, -1, -1, -1, -1, -1, -1},
        { 9,  8,  7,  5,  9,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 5,  8,  4,  5, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1},
        { 5,  0,  4,  5, 11,  0,  5, 10, 11, 11,  3,  0, -1, -1, -1, -1},
        { 0,  1,  9,  8,  4, 10,  8, 10, 11, 10,  4,  5, -1, -1, -1, -1},
        {10, 11,  4, 10,  4,  5, 11,  3,  4,  9,  4,  1,  3,  1,  4, -1},
        { 2,  5,  1,  2,  8,  5,  2, 11,  8,  4,  5,  8, -1, -1, -1, -1},
        { 0,  4, 11,  0, 11,  3,  4,  5, 11,  2, 11,  1,  5,  1, 11, -1},
        { 0,  2,  5,  0,  5,  9,  2, 11,  5,  4,  5,  8, 11,  8,  5, -1},
        { 9,  4,  5,  2, 11,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 2,  5, 10,  3,  5,  2,  3,  4,  5,  3,  8,  4, -1, -1, -1, -1},
        { 5, 10,  2,  5,  2,  4,  4,  2,  0, -1, -1, -1, -1, -1, -1, -1},
        { 3, 10,  2,  3,  5, 10,  3,  8,  5,  4,  5,  8,  0,  1,  9, -1},
        { 5, 10,  2,  5,  2,  4,  1,  9,  2,  9,  4,  2, -1, -1, -1, -1},
        { 8,  4,  5,  8,  5,  3,  3,  5,  1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  4,  5,  1,  0,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 8,  4,  5,  8,  5,  3,  9,  0,  5,  0,  3,  5, -1, -1, -1, -1},
        { 9,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4, 11,  7,  4,  9, 11,  9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        { 0,  8,  3,  4,  9,  7,  9, 11,  7,  9, 10, 11, -1, -1, -1, -1},
        { 1, 10, 11,  1, 11,  4,  1,  4,  0,  7,  4, 11, -1, -1, -1, -1},
        { 3,  1,  4,  3,  4,  8,  1, 10,  4,  7,  4, 11, 10, 11,  4, -1},
        { 4, 11,  7,  9, 11,  4,  9,  2, 11,  9,  1,  2, -1, -1, -1, -1},
        { 9,  7,  4,  9, 11,  7,  9,  1, 11,  2, 11,  1,  0,  8,  3, -1},
        {11,  7,  4, 11,  4,  2,  2,  4,  0, -1, -1, -1, -1, -1, -1, -1},
        {11,  7,  4, 11,  4,  2,  8,  3,  4,  3,  2,  4, -1, -1, -1, -1},
        { 2,  9, 10,  2,  7,  9,  2,  3,  7,  7,  4,  9, -1, -1, -1, -1},
        { 9, 10,  7,  9,  7,  4, 10,  2,  7,  8,  7,  0,  2,  0,  7, -1},
        { 3,  7, 10,  3, 10,  2,  7,  4, 10,  1, 10,  0,  4,  0, 10, -1},
        { 1, 10,  2,  8,  7,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  9,  1,  4,  1,  7,  7,  1,  3, -1, -1, -1, -1, -1, -1, -1},
        { 4,  9,  1,  4,  1,  7,  0,  8,  1,  8,  7,  1, -1, -1, -1, -1},
        { 4,  0,  3,  7,  4,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 4,  8,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 9, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 3,  0,  9,  3,  9, 11, 11,  9, 10, -1, -1, -1, -1, -1, -1, -1},
        { 0,  1, 10,  0, 10,  8,  8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        { 3,  1, 10, 11,  3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  2, 11,  1, 11,  9,  9, 11,  8, -1, -1, -1, -1, -1, -1, -1},
        { 3,  0,  9,  3,  9, 11,  1,  2,  9,  2, 11,  9, -1, -1, -1, -1},
        { 0,  2, 11,  8,  0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 3,  2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 2,  3,  8,  2,  8, 10, 10,  8,  9, -1, -1, -1, -1, -1, -1, -1},
        { 9, 10,  2,  0,  9,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 2,  3,  8,  2,  8, 10,  0,  1,  8,  1, 10,  8, -1, -1, -1, -1},
        { 1, 10,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 1,  3,  8,  9,  1,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  9,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        { 0,  3,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
}// clang-format on
};

//template<Point P>
//P interpolate_edge(
//        float iso_level, const Voxel<P> &voxel1, const Voxel<P> &voxel2) {
//    if (std::abs(iso_level - voxel1.value) < 0.00001) return voxel1.position;
//    if (std::abs(iso_level - voxel2.value) < 0.00001) return voxel2.position;
//    if (std::abs(voxel1.value - voxel2.value) < 0.00001) return voxel1.position;
//
//    float mu = (iso_level - voxel1.value) / (voxel2.value - voxel1.value);
//    P     point;
//    point.x = voxel1.position.x + mu * (voxel2.position.x - voxel1.position.x);
//    point.y = voxel1.position.y + mu * (voxel2.position.y - voxel1.position.y);
//    point.z = voxel1.position.z + mu * (voxel2.position.z - voxel1.position.z);
//
//    return point;
//}

//template< typename U, int dim>
//Eigen::Matrix<U, dim, 1> interpolate_edge(U iso_level,
//        const Eigen::Matrix<U, dim, 1> &pos1, U val1,
//        const Eigen::Matrix<U, dim, 1> &pos2, U val2) {
//    if (std::abs(iso_level - val1) < 0.00001) return pos1;
//    if (std::abs(iso_level - val2) < 0.00001) return pos2;
//    if (std::abs(val1 - val2) < 0.00001) return pos1;
//
//    U                        mu    = (iso_level - val1) / (val2 - val1);
//    Eigen::Matrix<U, dim, 1> point = pos1 + mu * (pos2 - pos1);
//
//    return point;
//}

//template<typename U, int dim>
//Eigen::Matrix<U, dim, 1> interpolate_edge(U iso_level,
//        const Eigen::Matrix<U, dim, 1>     &pos1,
//        U                                   val1,
//        const Eigen::Matrix<U, dim, 1>     &pos2,
//        U                                   val2) {
//
//    // Copy input vectors to allow modifications
//    Eigen::Matrix<U, dim, 1> v1 = pos1;
//    Eigen::Matrix<U, dim, 1> v2 = pos2;
//
//    // Ensure v1, val1 always represents the point with smaller coordinates
//    if (v2[0] < v1[0] || v2[1] < v1[1] || (dim == 3 && v2[2] < v1[2])) {
//        std::swap(v1, v2);
//        std::swap(val1, val2);
//    }
//
//    // Perform linear interpolation
//    Eigen::Matrix<U, dim, 1> point;
//    if (std::abs(val1 - val2) > 0.00001) {
//        U mu  = (iso_level - val1) / (val2 - val1);
//        point = v1 + mu * (v2 - v1);
//    } else {
//        point = v1;
//    }
//
//    return point;
//}
//template<typename U, int dim>
//Eigen::Matrix<U, dim, 1> interpolate_edge(U iso_level,
//        const Eigen::Matrix<U, dim, 1>     &pos1,
//        U                                   val1,
//        const Eigen::Matrix<U, dim, 1>     &pos2,
//        U                                   val2) {
//
//    // Check if the values are close enough to consider them equal
//    bool values_are_close = std::abs(val1 - val2) < 0.00001;
//
//    // If the iso_level is close to val1 or val2, return the corresponding position
//    if (std::abs(iso_level - val1) < 0.00001) { return pos1; }
//    if (std::abs(iso_level - val2) < 0.00001) { return pos2; }
//
//    // If the values are close, return pos1
//    if (values_are_close) { return pos1; }
//
//    // Perform linear interpolation
//    U mu = (iso_level - val1) / (val2 - val1 + 1e-10);
//    return pos1 + mu * (pos2 - pos1);
//}


template<typename U, int dim>
Eigen::Matrix<U, dim, 1> interpolate_edge(U iso_level,
        const Eigen::Matrix<U, dim, 1>     &pos1,
        U                                   val1,
        const Eigen::Matrix<U, dim, 1>     &pos2,
        U                                   val2,
        U                                   eps = 0.00001) {

    if (std::abs(iso_level - val1) < eps) { return pos1; }
    if (std::abs(iso_level - val2) < eps) { return pos2; }
    if (std::abs(val1 - val2) < eps) { return pos1; }
    U mu = (iso_level - val1) / (val2 - val1);
    return pos1 + mu * (pos2 - pos1);
}


std::array<std::array<int, 3>, 8> vertex_offsets = {
        {{0, 0, 0},
         {1, 0, 0},
         {1, 1, 0},
         {0, 1, 0},
         {0, 0, 1},
         {1, 0, 1},
         {1, 1, 1},
         {0, 1, 1}}
};


template<typename U, int dim>
std::pair<Eigen::Matrix<U, dim, 1>, U> get_vertex(
        const std::array<Eigen::Tensor<U, 3>, 3> &grid_positions,
        const Eigen::Tensor<U, 3>                &grid_values,
        const Eigen::array<size_t, 3>            &indices,
        int                                       vertex_index) {

    Eigen::array<size_t, 3> vertex_indices = indices;
    for (int j = 0; j < 3; j++) {
        vertex_indices[j] += vertex_offsets[vertex_index][j];
    }

    return std::make_pair(
            Eigen::Matrix<U, dim, 1>(grid_positions[0](vertex_indices[0],
                                             vertex_indices[1],
                                             vertex_indices[2]),
                    grid_positions[1](vertex_indices[0],
                            vertex_indices[1],
                            vertex_indices[2]),
                    grid_positions[2](vertex_indices[0],
                            vertex_indices[1],
                            vertex_indices[2])),
            grid_values(
                    vertex_indices[0], vertex_indices[1], vertex_indices[2]));
}


#include <spdlog/spdlog.h>

template<typename U, int dim = 3>
mesh::face_vertex<U, Eigen::Matrix<U, dim, 1>> marching_cubes(
        const std::array<Eigen::Tensor<U, 3>, 3> &grid_positions,
        const Eigen::Tensor<U, 3>                &grid_values,
        U                                         iso_level) {
    spdlog::info("Started marching_cubes method");
    mesh::face_vertex<U, Eigen::Matrix<U, dim, 1>> mesh;
    using facet_type =
            typename mesh::face_vertex<U, Eigen::Matrix<U, dim, 1>>::facet_type;

    constexpr int     num_vertices_per_voxel = 8;
    const std::size_t x_dim                  = grid_values.dimension(0);
    const std::size_t y_dim                  = grid_values.dimension(1);
    const std::size_t z_dim                  = grid_values.dimension(2);

    for (std::size_t x = 0; x < x_dim-1; ++x) {
        for (std::size_t y = 0; y < y_dim-1; ++y) {
            for (std::size_t z = 0; z < z_dim-1; ++z) {
                spdlog::debug("Processing voxel ({}, {}, {})", x, y, z);
                Eigen::array<size_t, 3> indices{x, y, z};
                std::vector<std::pair<Eigen::Matrix<U, dim, 1>, U>> vertices(
                        num_vertices_per_voxel);

                for (int i = 0; i < num_vertices_per_voxel; i++) {
                    vertices[i] = get_vertex<U, dim>(
                            grid_positions, grid_values, indices, i);
                }

                int cube_index = 0;
                for (int i = 0; i < 8; ++i) {
                    if (vertices[i].second < iso_level) cube_index |= 1 << i;
                }

                std::array<Eigen::Matrix<U, dim, 1>, 12> edge_intersections;
                for (int i = 0; i < 12; ++i) {
                    // check if its on the edge of the grid
                    if (edge_table[cube_index] & (1 << i)) {
                        edge_intersections[i] = interpolate_edge<U, dim>(
                                iso_level,
                                vertices[edge_vertex_indices[i][0]].first,
                                vertices[edge_vertex_indices[i][0]].second,
                                vertices[edge_vertex_indices[i][1]].first,
                                vertices[edge_vertex_indices[i][1]].second);
                    }
                }

                std::map<std::tuple<int, int, int, int>, int> edge_to_vertex;
                auto get_vertex_index = [&](int edge_index) -> int {
                    std::tuple<int, int, int, int> edge_key
                            = std::make_tuple(x, y, z, edge_index);
                    auto it = edge_to_vertex.find(edge_key);
                    if (it != edge_to_vertex.end()) {
                        return it->second;
                    } else {
                        int vertex_index = mesh.vertices.size();
                        mesh.vertices.push_back(edge_intersections[edge_index]);
                        edge_to_vertex[edge_key] = vertex_index;
                        return vertex_index;
                    }
                };

                for (int i = 0; tri_table[cube_index][i] != -1; i += 3) {
                    facet_type face;
                    for (int j = 0; j < 3; ++j) {
                        face[j] = get_vertex_index(
                                tri_table[cube_index][i + j]);
                    }
                    mesh.faces.push_back(face);
                }
            }
        }
    }

    spdlog::info("Finished marching_cubes method");
    return mesh;
}
//
//emplate<typename U, int dim = 3>
//        mesh::face_vertex<U, Eigen::Matrix<U, dim, 1>> marching_cubes(
//                const std::array<Eigen::Tensor<U, 3>, 3> &grid_positions,
//                const Eigen::Tensor<U, 3>                &grid_values,
//                U                                         iso_level) {
//    spdlog::info("Started marching_cubes method");
//    mesh::face_vertex<U, Eigen::Matrix<U, dim, 1>> mesh;
//    using facet_type =
//            typename mesh::face_vertex<U, Eigen::Matrix<U, dim, 1>>::facet_type;
//
//    constexpr int     num_vertices_per_voxel = 8;
//    const std::size_t x_dim                  = grid_values.dimension(0);
//    const std::size_t y_dim                  = grid_values.dimension(1);
//    const std::size_t z_dim                  = grid_values.dimension(2);
//
//    for (std::size_t x = 0; x < x_dim-1; ++x) {
//        for (std::size_t y = 0; y < y_dim-1; ++y) {
//            for (std::size_t z = 0; z < z_dim-1; ++z) {
//                spdlog::debug("Processing voxel ({}, {}, {})", x, y, z);
//                Eigen::array<size_t, 3> indices{x, y, z};
//                std::vector<std::pair<Eigen::Matrix<U, dim, 1>, U>> vertices(
//                        num_vertices_per_voxel);
//
//                for (int i = 0; i < num_vertices_per_voxel; i++) {
//                    vertices[i] = get_vertex<U, dim>(
//                            grid_positions, grid_values, indices, i);
//                }
//
//                int cube_index = 0;
//                for (int i = 0; i < 8; ++i) {
//                    if (vertices[i].second < iso_level) cube_index |= 1 << i;
//                }
//
//                std::array<Eigen::Matrix<U, dim, 1>, 12> edge_intersections;
//                for (int i = 0; i < 12; ++i) {
//                    // check if its on the edge of the grid
//                    if (edge_table[cube_index] & (1 << i)) {
//                        edge_intersections[i] = interpolate_edge<U, dim>(
//                                iso_level,
//                                vertices[edge_vertex_indices[i][0]].first,
//                                vertices[edge_vertex_indices[i][0]].second,
//                                vertices[edge_vertex_indices[i][1]].first,
//                                vertices[edge_vertex_indices[i][1]].second);
//                    }
//                }
//
//                std::map<std::tuple<int, int, int, int>, int> edge_to_vertex;
//                auto get_vertex_index = [&](int edge_index) -> int {
//                    std::tuple<int, int, int, int> edge_key
//                            = std::make_tuple(x, y, z, edge_index);
//                    auto it = edge_to_vertex.find(edge_key);
//                    if (it != edge_to_vertex.end()) {
//                        return it->second;
//                    } else {
//                        int vertex_index = mesh.vertices.size();
//                        mesh.vertices.push_back(edge_intersections[edge_index]);
//                        edge_to_vertex[edge_key] = vertex_index;
//                        return vertex_index;
//                    }
//                };
//
//                for (int i = 0; tri_table[cube_index][i] != -1; i += 3) {
//                    facet_type face;
//                    for (int j = 0; j < 3; ++j) {
//                        face[j] = get_vertex_index(
//                                tri_table[cube_index][i + j]);
//                    }
//                    mesh.faces.push_back(face);
//                }
//            }
//        }
//    }
//
//    spdlog::info("Finished marching_cubes method");
//    return mesh;
//}

}// namespace odin

#endif//ODIN_MARCHING_CUBES_H