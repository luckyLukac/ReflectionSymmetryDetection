#ifndef REFLECTIONSYMMETRY_H
#define REFLECTIONSYMMETRY_H

#include <set>
#include <vector>

#include "helper_classes/linesegment.h"
#include "helper_classes/plane.h"
#include "helper_classes/positionfromplane.h"
#include "localsymmetry.h"


namespace Symmetry {
    // Class for reflection symmetries.
    template <typename T>
    class ReflectionSymmetryT : public LocalSymmetry {
    private:
        // OBJECT VARIABLES
        uint m_Index;
        Plane m_Plane;
        std::vector<LineSegment> m_LineSegments;
        std::vector<Voxel> m_SymmetryVoxels;

    public:
        // OBJECT METHODS
        void clusterRecursiveStep(std::set<Voxel>& cluster, VoxelVector& voxelVector, const uint x, const uint y, const uint z) const;  // Clustering recursive step.
        std::vector<std::set<Voxel>> findClusters(VoxelVector voxels, const VoxelMesh& voxelMesh) const;                                // Finding clusters of symmetry voxels in symmetry.

        // CONSTRUCTORS
        ReflectionSymmetryT();                                                      // Empty constructor.
        ReflectionSymmetryT(Plane plane, std::vector<LineSegment> lineSegments);    // Default constructor.
        ReflectionSymmetryT(const Plane& plane, const std::vector<Voxel>& voxels);  // Constructor with the plane and the voxels.

        // OVERLOADED OPERATORS
        bool operator == (ReflectionSymmetryT& s);                 // Operator == is used for checking whether the two symmetries are equal.
        ReflectionSymmetryT operator &= (ReflectionSymmetryT& s);  // Operator &= is used for merging the current symmetry with another.

        // GETTERS AND SETTERS
        uint getIndex() const;                                // Index getter.
        void setIndex(const uint index);                      // Index setter.
        Plane getPlane() const;                               // Symmetry plane getter.
        std::vector<LineSegment> getLineSegments() const;     // Line segments getter.
        u128 getLineSegmentCount() const;                     // Getting number of line segments.
        void addLineSegment(const LineSegment& lineSegment);  // Adding a new line segment.
        std::vector<Voxel>& getVoxels();                      // Voxel vector getter.
        u128 getVoxelCount() const;                           // Voxel count getter.
        void addSymmetryVoxel(const Voxel& voxel);            // Symmetry voxel adder.

        // OBJECT METHODS
        std::vector<PositionFromPlane> calculatePositionsFromPlane(const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh, const Plane& plane, const bool validPlane, u128& pointCount, const bool normalize);  // Calculating the point positions according to a plane.
    };


    // ALIASES
    using ReflectionSymmetry = ReflectionSymmetryT<double>;


    // REFLECTION SYMMETRY FUNCTIONS
    namespace ReflectionSymmetryFunctions {
        std::vector<ReflectionSymmetry> findSimpleSymmetries(const std::vector<LineSegment>& lineSegments, const double tolerance);                                                                                                                                                                  // Finding simple symmetries.
        bool calculateSimpleReflectionSymmetry(ReflectionSymmetry& symmetry, const LineSegment& l1, const LineSegment& l2, const double tolerance);                                                                                                                                                  // Checking whether two line segments represent a reflection symmetry.
        std::vector<ReflectionSymmetry> mergeSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const double tolerance);                                                                                                                                                              // Merging reflection symmetries.
        void getVoxelsInSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelMesh& voxelMesh);                                                                                                                                                                               // Getting voxels for each reflection symmetry.
        void addMaterialVoxelsToSymmetry(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelVector& voxels, const VoxelMesh& voxelMesh, std::vector<Voxel>& materialVoxels);                                                                                                                                           // Adding interesting voxels back to symmetry (if in symmetry).
        void findInterestingClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z);                                                                                                         // Clustering step.
        void removeSmallClusters(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minimumClusterSize);                                                                                                                       // Removing small clusters of symmetries.
        void processPointsOnThePlaneAndVoxelEdge(std::vector<ReflectionSymmetry>& reflectionSymmetries, const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh);                                                                                               // Processing possible points that lie on the voxel edge and the symmetry plane simultaneously.
        void removeVoxelsWithoutPairs(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelMesh& voxelMesh);                                                                                                                                                                            // Removing voxels without their pairs in reflection symmetries.
        std::vector<ReflectionSymmetry> findPartialSymmetries(const std::vector<ReflectionSymmetry>& reflectionSymmetries, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minSymmetrySize);                                                                                             // Finding partial symmetries from all reflection symmetries.
        std::vector<ReflectionSymmetry> calculateReflectionSymmetries(const std::vector<Point>& points, const VoxelVector& voxels, const VoxelMesh& voxelMesh, std::vector<Voxel>& materialVoxels, const double tolerance, const uint minimumClusterSize, const double minimumLineSegmentDistance);  // Calculating reflection symmetries.
        std::vector<ReflectionSymmetry> calculateNearestSymmetriesByPlane(const std::vector<ReflectionSymmetry>& symmetries, const Plane& desiredPlane, const VoxelMesh& voxelMesh);                                                                                                                                                    // Calculation of nearest symmetries according to the desired symmetry plane.
    };
};

#endif // REFLECTIONSYMMETRY_H
