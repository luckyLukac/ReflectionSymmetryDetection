#pragma once

#include "point.h"
#include "voxel.h"


namespace Symmetry {
    // Namespace for Point related functions.
    namespace PointFunctions {
        double distance(const Point& p1, const Point& p2);                     // Calculation of the Euclidean distance between two points.
        Point calculateMidPoint(const Point& p1, const Point& p2);             // Calculation of the midpoint of a line segment.
        Point calculateBarycenterPoint(const std::vector<Point>& points);      // Calculation of a barycenter point (centroid).
        Voxel voxelFromPoint(const Point& point, const VoxelMesh& voxelMesh);  // Calculation of the voxel position according to voxel mesh.
    };

    // Namespace for Vector related functions.
    namespace VectorFunctions {
        Vector3d crossProduct(const Vector3d& lhs, const Vector3d& rhs);  // Cross product of two vectors.
        double angle(const Vector3d& a, const Vector3d& b);               // Getting the angle between two vectors (represented as Points).
    };
}