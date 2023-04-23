#include <cmath>

#include "point.h"
#include "voxel.h"

using namespace Symmetry;



// CONSTRUCTORS
// Default constructor without parameters.
Voxel::Voxel() :
    x(0),
    y(0),
    z(0),
    interesting(false),
    material(false),
    checked(false),
    inSymmetry(false)
{}

// Constructor with all three coordinates.
Voxel::Voxel(uint x, uint y, uint z) :
    x(x),
    y(y),
    z(z),
    interesting(false),
    material(false),
    checked(false),
    inSymmetry(false)
{}

// Constructor with all three coordinates and booleans for a super-interesting and an interesting voxel.
Voxel::Voxel(const uint x, const uint y, const uint z, const bool si, const bool i) :
    x(x),
    y(y),
    z(z),
    interesting(si),
    material(i),
    checked(false),
    inSymmetry(false)
{}

// Constructor with all three coordinates and a boolean for in-symmetry.
Voxel::Voxel(const uint x, const uint y, const uint z, const bool is) :
    x(x),
    y(y),
    z(z),
    interesting(false),
    material(false),
    checked(false),
    inSymmetry(is)
{}

// Copy constructor with Point<> and a boolean for an interesting voxel.
Voxel::Voxel(const IntPoint& p, const bool i) :
    x(p.x),
    y(p.y),
    z(p.z),
    interesting(i),
    material(i),
    checked(false),
    inSymmetry(false)
{}

// Copy constructor with Point.
Voxel::Voxel(const Point& p) :
    x(static_cast<uint>(p.x)),
    y(static_cast<uint>(p.y)),
    z(static_cast<uint>(p.z)),
    interesting(false),
    material(false),
    checked(false),
    inSymmetry(false)
{}


// OVERLOADED OPERATORS
// Checking whether the two voxels have the same coordinates.
bool Voxel::operator == (const Voxel& v) const {
    if (x == v.x && y == v.y && z == v.z) {
        return true;
    }

    return false;
}

// Operator < serves for the sort method as the comparator.
bool Voxel::operator < (const Voxel& v) const {
    return (
        (x < v.x) ||
        (x == v.x && y < v.y) ||
        (x == v.x && y == v.y && z < v.z)
    );
}

// Operator * serves for multiplication.
Voxel Voxel::operator * (const uint factor) const {
    return Voxel(x * factor, y * factor, z * factor);
}


// OBJECT METHODS
// Getting the coordinate from the voxel center point.
Point Voxel::centerCoordinate(const VoxelMesh& voxelMesh) const {
    const double pointX = x + (0.5 * voxelMesh.voxelSideSize);
    const double pointY = y + (0.5 * voxelMesh.voxelSideSize);
    const double pointZ = z + (0.5 * voxelMesh.voxelSideSize);

    return Point(pointX, pointY, pointZ);
}

// Getting all voxel coordinates in an array.
std::array<uint, 3> Voxel::getNormalizedCoordinates(const VoxelMesh& voxelMesh) const {
    const uint xCoordinate = static_cast<uint>(x / voxelMesh.voxelSideSize);
    const uint yCoordinate = static_cast<uint>(y / voxelMesh.voxelSideSize);
    const uint zCoordinate = static_cast<uint>(z / voxelMesh.voxelSideSize);
    
    return std::array<uint, 3>({ xCoordinate, yCoordinate, zCoordinate });
}

// Getting voxel coordinates in an array.
std::array<uint, 3> Voxel::getCoordinates() const {
    return std::array<uint, 3>({ x, y, z });
}

// Getting a normalized voxel.
Voxel Voxel::getNormalizedVoxel(const uint factor) const {
    const uint xCoordinate = static_cast<uint>(x / factor);
    const uint yCoordinate = static_cast<uint>(y / factor);
    const uint zCoordinate = static_cast<uint>(z / factor);

    return Voxel(xCoordinate, yCoordinate, zCoordinate);
}

// Resetting all voxel properties.
void Voxel::resetAllProperties() {
    material = interesting = inSymmetry = checked = false;
}



// VOXEL FUNCTIONS NAMESPACE
// Getting the data whether the voxel should be painted brightly.
bool VoxelFunctions::isBrightNeighbor(const uint x, const uint y, const uint z) {
    // First bright option.
    if (z % 2 == 0 && y % 2 == 0 && x % 2 == 0) {
        return true;
    }

    // Second bright option.
    if (z % 2 == 0 && y % 2 == 1 && x % 2 == 1) {
        return true;
    }

    // Third bright option.
    if (z % 2 == 1 && y % 2 == 1 && x % 2 == 0) {
        return true;
    }

    // Fourth bright option.
    if (z % 2 == 1 && y % 2 == 0 && x % 2 == 1) {
        return true;
    }

    return false;
}

// Getting the layer of the voxel mesh according to the Z coordinate.
uint VoxelFunctions::getLayerInVoxelMesh(const double z, const VoxelMesh& vm) {
    const uint layer = static_cast<uint>((z - vm.minZ) / vm.voxelSideSize);
    return layer;
}

// Returning interesting points from voxel vector.
std::vector<Point> VoxelFunctions::getInterestingPointsFromVoxels(
    const VoxelVector& voxels,
    const double voxelEdge,
    const bool mustBeSuperInteresting
)
{
    std::vector<Point> points;

    // Creating a new point from a voxel.
    for (uint i = 0; i < voxels.size(); i++) {
        for (uint j = 0; j < voxels[i].size(); j++) {
            for (uint k = 0; k < voxels[i][j].size(); k++) {
                if (!mustBeSuperInteresting || voxels[i][j][k].interesting) {
                    points.push_back(
                        Point(
                            voxels[i][j][k].x + voxelEdge / 2,
                            voxels[i][j][k].y + voxelEdge / 2,
                            voxels[i][j][k].z + voxelEdge / 2
                        )
                    );
                }
            }
        }
    }

    return points;
}

// Clearing all variables in a voxel vector.
void VoxelFunctions::clearVoxelVector(VoxelVector& voxelVector) {
    for (auto& vector2D : voxelVector) {
        for (auto& vector1D : vector2D) {
            for (Voxel& voxel : vector1D) {
                voxel.resetAllProperties();
            }
        }
    }
}

// Checking whether the voxel is in bounds of the voxel mesh.
bool VoxelFunctions::isVoxelInVoxelMesh(const Voxel& voxel, const VoxelMesh& vm) {
    if (voxel.x < 0 || voxel.x > vm.maxX || voxel.y < 0 || voxel.y > vm.maxY || voxel.z < 0 || voxel.z > vm.maxZ) {
        return false;
    }

    return true;
}