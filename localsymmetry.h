#ifndef LOCALSYMMETRY_H
#define LOCALSYMMETRY_H

#include <set>
#include <stack>
#include <string>
#include <vector>

#include "helper_classes/linesegment.h"
#include "helper_classes/point.h"
#include "helper_classes/voxel.h"


namespace Symmetry {
    // Class for local symmetries that contains methods for
    // both reflection symmetries and rotational symmetries.
    class LocalSymmetry {
    protected:
        // CLASS VARIABLES
        std::vector<Point> m_Points;             // Vector of LAS points.
        VoxelMesh m_VoxelMesh;                   // Object with the data about the voxel mesh.
        VoxelVector m_Voxels;                    // 3D vector of voxels.
        std::vector<Voxel> m_InterestingVoxels;  // Vector of interesting voxels.
        std::vector<Voxel> m_MaterialVoxels;     // Vector of material voxels that are not interesting.

    public:
        // GETTERS AND SETTERS
        void setPoints(std::vector<Point>& points);   // Setting points.
        const std::vector<Point>& getPoints() const;  // Point vector getter.
        VoxelMesh& getVoxelMesh();                    // Voxel mesh object getter.
        VoxelVector& getVoxels();                     // 3D voxel vector getter.
        std::vector<Voxel>& getInterestingVoxels();   // Interesting voxel getters.
        std::vector<Voxel>& getMaterialVoxels();      // Material voxel getters.
    };


    // SYMMETRY FUNCTIONS NAMESPACE
    namespace SymmetryFunctions {
        void calculateVoxelMeshByVoxelSideLength(VoxelVector& voxels, VoxelMesh& voxelMesh, const int userInput);                                                                                                     // Voxel mesh calculation by voxel side length.
        void calculateVoxelMeshByMaximumVoxelCount(VoxelVector& voxels, VoxelMesh& voxelMesh, const uint userInput);                                                                                                  // Voxel mesh calculation by maximum voxel count.
        void buildVoxelVector(VoxelVector& voxels, const VoxelMesh& voxelMesh);                                                                                                                                       // Building a 3D voxel vector from a voxel mesh.
        void findInterestingVoxels(const std::vector<Point>& points, std::vector<Voxel>& interesting, std::vector<Voxel>& material, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minClusterSize);      // Searching super-interesting voxels according to points in the 3D voxel vector.
        SplitLineSegments splitLineSegmentVectorIntoParts(const std::vector<LineSegment>& lineSegments, const double tolerance);                                                                                      // Splitting a line segment vector into parts, where the lenghts of line segments is below the tolerance.
        VoxelVector getNormalisedVoxelVector(const VoxelVector& voxelVector, const VoxelMesh& voxelMesh);                                                                                                             // Getting a normalized 3D voxel vector (voxel edge size equals 1).
        void findInterestingClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z);                          // Clustering step.
        std::vector<std::set<Voxel>> findInterestingClusters(VoxelVector& normalizedVoxels, const VoxelMesh& voxelMesh);                                                                                              // Finding clusters of symmetry voxels in symmetry.
        void removeSmallInterestingClusters(std::vector<Voxel>& interestingVoxels, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minimumClusterSize);                                                   // Removing small clusters of interesting voxels.
        std::vector<LineSegment> calculateLineSegmentsBetweenPoints(const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh, const double tolerance, const double minDistance);  // Getting line segments between all pairs of points (fully connected graph).
        double calculateInterestingPercentInLineSegment(const LineSegment& lineSegment, const VoxelVector& voxels, const VoxelMesh& voxelMesh);                                                                       // Calculating the percent of interesting line segment sections.
    }
};

#endif // LOCALSYMMETRY_H
