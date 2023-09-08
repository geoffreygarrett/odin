#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// Sample code, assumes you have defined Vector3 and Voxel types
struct Vector3 {
    float x, y, z;
};

struct Voxel {
    Vector3 minCorner;
    float   size;
    float   density;// An example field, to be determined by your actual data
};

class OctreeNode {
public:
    Voxel                       voxel;
    std::unique_ptr<OctreeNode> children[8];// 8 children for octree

    OctreeNode(const Voxel& v)
        : voxel(v) {}

    // Determine if this node is at the boundary (you can change this condition)
    bool isBoundary() const {
        return std::abs(voxel.density) < 0.1;// Dummy condition
    }

    // Subdivide this node into 8 children
    void subdivide() {
        if (!isBoundary()) { return; }

        float newSize = voxel.size / 2.0;
        for (int x = 0; x < 2; ++x) {
            for (int y = 0; y < 2; ++y) {
                for (int z = 0; z < 2; ++z) {
                    Vector3 newMinCorner = {voxel.minCorner.x + x * newSize,
                            voxel.minCorner.y + y * newSize,
                            voxel.minCorner.z + z * newSize};
                    Voxel newVoxel = {newMinCorner, newSize, 0.0f /* fill density appropriately */};
                    children[x * 4 + y * 2 + z] = std::make_unique<OctreeNode>(newVoxel);
                    children[x * 4 + y * 2 + z]->subdivide();
                }
            }
        }
    }
};

// A placeholder for your Marching Cubes function
void runMarchingCubesOnLeaf(const Voxel& voxel) {
    // Perform Marching Cubes on this voxel
    std::cout << "Running Marching Cubes on voxel at (" << voxel.minCorner.x << ", "
              << voxel.minCorner.y << ", " << voxel.minCorner.z << ") with size " << voxel.size
              << "\n";
}

void processOctreeNode(const OctreeNode& node) {
    if (node.isBoundary()) {
        if (node.children[0] == nullptr) {
            runMarchingCubesOnLeaf(node.voxel);
        } else {
            for (const auto& child: node.children) {
                if (child != nullptr) { processOctreeNode(*child); }
            }
        }
    }
}

int main() {
    Voxel rootVoxel = {
            {0, 0, 0},
            1.0, 0.0
    };// Starting voxel
    OctreeNode rootNode(rootVoxel);
    rootNode.subdivide();

    processOctreeNode(rootNode);

    return 0;
}
