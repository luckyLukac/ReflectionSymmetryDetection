#ifndef VOXEL_H
#define VOXEL_H

#include <array>
#include <limits>
#include <vector>

#include "constants.h"
#include "point.h"


namespace Symmetry {
    // Structure for dividing the voxels into a mesh.
    struct VoxelMesh {
        uint voxelCount = 0;                               // Number of voxels.
        uint voxelSideSize = 0;                            // Length of a voxel edge.
        uint voxelX = 0;                                   // Number of voxels by X.
        uint voxelY = 0;                                   // Number of voxels by Y.
        uint voxelZ = 0;                                   // Number of voxels by Z.
        double minX = std::numeric_limits<double>::max();  // Minimum X.
        double maxX = std::numeric_limits<double>::min();  // Maximum X.
        double minY = std::numeric_limits<double>::max();  // Minimum Y.
        double maxY = std::numeric_limits<double>::min();  // Maximum Y.
        double minZ = std::numeric_limits<double>::max();  // Minimum Z.
        double maxZ = std::numeric_limits<double>::min();  // Maximum Z.
        double deltaX = 0.0;                               // Distance between the maximum and the minimum X coordinate.
        double deltaY = 0.0;                               // Distance between the maximum and the minimum Y coordinate.
        double deltaZ = 0.0;                               // Distance between the maximum and the minimum Z coordinate.
    };


    // Voxel for discretization in the 3D space.
    struct Voxel {
        // OBJECT VARIABLES
        uint x;                              // X coordinate.
        uint y;                              // Y coordinate.
        uint z;                              // Z coordinate.
        bool interesting;                    // True if the voxel contains at least one point and is not surrounded with other interesting voxels.
        bool material;                       // True if the voxel contains at least one point and is not surrounded with other interesting voxels.
        bool checked;                        // True if the voxel has been checked during the geometry search.
        bool inSymmetry;                     // True if the voxel is a part of a symmetry.
        
        // CONSTRUCTORS
        Voxel();                                                                       // Default constructor without parameters.
        Voxel(const uint x, const uint y, const uint z);                               // Constructor with all three coordinates.
        Voxel(const uint x, const uint y, const uint z, const bool si, const bool i);  // Constructor with all three coordinates and booleans for a super-interesting and an interesting voxel.
        Voxel(const uint x, const uint y, const uint z, const bool is);                // Constructor with all three coordinates and a boolean for in-symmetry.
        Voxel(const IntPoint& p, const bool i);                                        // Copy constructor with Point and a boolean for an interesting voxel.
        Voxel(const Point& p);                                                         // Copy constructor with Point.

        // OVERLOADED OPERATORS
        bool operator == (const Voxel& v) const;  // Checking whether the two voxels have the same coordinates.
        bool operator < (const Voxel& v) const;   // Operator < serves for the sort method as the comparator.
        Voxel operator * (const uint factor) const;   // Operator * serves for multiplication.
        
        // OBJECT METHODS
        Point centerCoordinate(const VoxelMesh& voxelMesh) const;                         // Getting the coordinate from the voxel center point.
        std::array<uint, 3> getNormalizedCoordinates(const VoxelMesh& voxelMesh) const;   // Getting normalized voxel coordinates in an array.
        std::array<uint, 3> getCoordinates() const;                                       // Getting voxel coordinates in an array.
        Voxel getNormalizedVoxel(const uint factor) const;                                // Getting a normalized voxel.
        void resetAllProperties();                                                        // Resetting all voxel properties.
    };


    // ALIASES
    using VoxelVector = std::vector<std::vector<std::vector<Voxel>>>;


    // VOXEL FUNCTION NAMESPACE
    // Namespace for voxel functions.
    namespace VoxelFunctions {
        bool isBrightNeighbor(const uint x, const uint y, const uint z);                                                                                 // Getting the data whether the voxel should be painted brightly.
        uint getLayerInVoxelMesh(const double z, const VoxelMesh& vm);                                                                                   // Getting the layer of the voxel mesh according to the Z coordinate.
        std::vector<Point> getInterestingPointsFromVoxels(const VoxelVector& voxels, const double voxelEdge, const bool mustBeSuperInteresting = true);  // Returning interesting points from voxel vector.
        void clearVoxelVector(VoxelVector& voxelVector);                                                                                                 // Clearing all variables in a voxel vector.
        bool isVoxelInVoxelMesh(const Voxel& voxel, const VoxelMesh& vm);                                                                                // Checking whether the voxel is in bounds of the voxel mesh.
    }
};


#endif // VOXEL_H
