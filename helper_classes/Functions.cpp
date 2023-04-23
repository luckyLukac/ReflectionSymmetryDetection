#include "Functions.hpp"

using namespace Symmetry;



// POINT FUNCTIONS
// Calculation of the Euclidean distance between two points.
double PointFunctions::distance(const Point& p1, const Point& p2) {
    return (
        std::sqrt(
            std::pow(p1.x - p2.x, 2) +
            std::pow(p1.y - p2.y, 2) +
            std::pow(p1.z - p2.z, 2)
        )
        );
}

// Calculation of the midpoint of a line segment.
Point PointFunctions::calculateMidPoint(const Point& p1, const Point& p2) {
    // Midpoint can be calculated as the "average point" between two points: (p1 + p2) / 2.
    Point midpoint((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0, (p1.z + p2.z) / 2.0);

    return midpoint;
}

// Calculation of a barycenter point (centroid).
Point PointFunctions::calculateBarycenterPoint(const std::vector<Point>& points) {
    Point p;
    for (const Point& point : points) {
        p += point;
    }
    p /= static_cast<double>(points.size());

    return p;
}

// Calculation of the voxel position according to voxel mesh.
Voxel PointFunctions::voxelFromPoint(const Point& point, const VoxelMesh& voxelMesh) {
    const uint x = static_cast<uint>(point.x / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;
    const uint y = static_cast<uint>(point.y / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;
    const uint z = static_cast<uint>(point.z / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;

    return Voxel(x, y, z);
}




// VECTOR FUNCTIONS
// Cross product of two vectors.
Vector3d VectorFunctions::crossProduct(const Vector3d& lhs, const Vector3d& rhs) {
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x
    };
}

// Getting the angle between two vectors (represented as Points).
double VectorFunctions::angle(const Vector3d& a, const Vector3d& b) {
    const double dot = a.x * b.x + a.y * b.y;
    const double det = a.x * b.y - a.y * b.x;
    const double angle = std::atan2(det, dot);

    // Angle lies inside the interval [0, 2 * PI].
    return angle >= 0 ? angle : angle + 2 * PI;
}