#pragma once
#include <algorithm>
#include <cmath>
#include <set>

#include "helper_classes/tolerance.h"
#include "reflectionsymmetry.h"

using namespace Symmetry;



// OBJECT METHODS
// Clustering recursive step.
template <typename T>
void ReflectionSymmetryT<T>::clusterRecursiveStep(std::set<Voxel>& cluster, VoxelVector& voxelVector, const uint x, const uint y, const uint z) const {
    // If the voxel coordinates lie outside of the voxel mesh or the
    // voxel has been already checked, the recursion unfolds.
    if (x < 0 || x >= m_VoxelMesh.voxelX ||
        y < 0 || y >= m_VoxelMesh.voxelY ||
        z < 0 || z >= m_VoxelMesh.voxelZ ||
        voxelVector[z][y][x].checked
    )
    {
        return;
    }

    // Setting the current voxel to checked.
    voxelVector[z][y][x].checked = true;

    // If the current voxel is not in symmetry, the recursion starts to unfold.
    if (!voxelVector[z][y][x].inSymmetry) {
        return;
    }

    // Inserting the current voxel to the cluster.
    cluster.insert(Voxel(x, y, z, true));

    // 26 recursive steps to for depth-first search of the neighborhood.
    // 9 in the upper layer, 8 (all but current) in the current layer, 9 in the lower layer.
    clusterRecursiveStep(cluster, voxelVector, x - 1, y - 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y - 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y - 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 0, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 0, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 0, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 1, z + 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y - 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y - 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y - 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 0, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 0, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 0, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 1, z + 0);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y - 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y - 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y - 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 0, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 0, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 0, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x - 1, y + 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 0, y + 1, z - 1);
    clusterRecursiveStep(cluster, voxelVector, x + 1, y + 1, z - 1);
}

// Finding clusters of symmetry voxels in symmetry.
template <typename T>
std::vector<std::set<Voxel>> ReflectionSymmetryT<T>::findClusters(VoxelVector voxels, const VoxelMesh& voxelMesh) const {
    std::vector<std::set<Voxel>> clusters;

    // Cleaning the checked property.
    for (uint z = 0; z < voxels.size(); z++) {
        for (uint y = 0; y < voxels[z].size(); y++) {
            for (uint x = 0; x < voxels[z][y].size(); x++) {
                voxels[z][y][x].checked = false;
            }
        }
    }

    // Getting the 3D voxel vector.
    for (const Voxel& v : m_SymmetryVoxels) {
        // Setting all voxels in symmetry to true.
        const auto [x, y, z] = v.getNormalizedCoordinates(voxelMesh);
        voxels[z][y][x].inSymmetry = true;
    }

    // Searching for clusters in each voxel in symmetry (depth-first-search).
    for (const Voxel& v : m_SymmetryVoxels) {
        // Getting the coordinates of the voxel.
        const auto [x, y, z] = v.getNormalizedCoordinates(voxelMesh);

        // If the current voxel has not been checked and is in symmetry, a new cluster has been found.
        if (!voxels[z][y][x].checked) {
            std::set<Voxel> cluster;         // Creating a new cluster.
            std::stack<Voxel> stack({ v });  // Creating a stack for flood fill.

            // While stack is not empty, we push the new items on it.
            while (!stack.empty()) {
                // Getting the top item from the stack.
                const Voxel top = stack.top();
                stack.pop();

                // Casting the coordinates of the top voxel.
                const auto [topX, topY, topZ] = top.getNormalizedCoordinates(voxelMesh);

                // 26 steps to for depth-first search of the neighborhood.
                // 9 in the upper layer, 8 (all but current) in the current layer, 9 in the lower layer.
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY - 1, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY - 1, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY - 1, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY + 0, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY + 0, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY + 0, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY + 1, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY + 1, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY + 1, topZ + 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY - 1, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY - 1, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY - 1, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY + 0, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY + 0, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY + 0, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY + 1, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY + 1, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY + 1, topZ + 0);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY - 1, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY - 1, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY - 1, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY + 0, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY + 0, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY + 0, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX - 1, topY + 1, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 0, topY + 1, topZ - 1);
                ReflectionSymmetryFunctions::findInterestingClusterItem(cluster, voxels, voxelMesh, stack, topX + 1, topY + 1, topZ - 1);
            }

            clusters.push_back(cluster);  // Adding a new cluster to the list.
        }
    }

    return clusters;
}



// CONSTRUCTORS
// Empty constructor.
template <typename T>
ReflectionSymmetryT<T>::ReflectionSymmetryT() :
    m_Plane(Plane()),
    m_LineSegments(std::vector<LineSegment>())
{}

// Default constructor.
template <typename T>
ReflectionSymmetryT<T>::ReflectionSymmetryT(
    Plane plane,
    std::vector<LineSegment> lineSegments
) :
    m_Plane(plane),
    m_LineSegments(lineSegments)
{}

// Constructor with the plane and the voxels.
template <typename T>
ReflectionSymmetryT<T>::ReflectionSymmetryT(const Plane& plane, const std::vector<Voxel>& voxels) :
    m_Plane(plane),
    m_SymmetryVoxels(voxels)
{}



// OVERLOADED OPERATORS
// Operator == is used for checking whether the two symmetries are equal.
template <typename T>
bool ReflectionSymmetryT<T>::operator == (ReflectionSymmetryT& s) {
    // First condition: same planes.
    if (m_Plane == s.m_Plane) {
        // Second condition: same voxel count.
        if (m_SymmetryVoxels.size() == s.m_SymmetryVoxels.size()) {
            // Third condition: same voxels.
            for (uint i = 0; i < m_SymmetryVoxels.size(); i++) {
                if (m_SymmetryVoxels[i] == s.m_SymmetryVoxels[i]) {
                    continue;
                }
                
                return false;
            }

            return true;
        }
    }

    return false;
}

// Operator &= is used for merging the current symmetry with another.
template <typename T>
ReflectionSymmetryT<T> ReflectionSymmetryT<T>::operator &= (ReflectionSymmetryT& s) {
    // Adding the line segments of the second symmetry.
    std::vector<LineSegment> lineSegments = s.getLineSegments();
    for (uint i = 0; i < lineSegments.size(); i++) {
        this->addLineSegment(lineSegments[i]);
    }

    return *this;
}



// GETTERS AND SETTERS
// Index getter.
template <typename T>
uint ReflectionSymmetryT<T>::getIndex() const {
    return m_Index;
}

// Index setter.
template <typename T>
void ReflectionSymmetryT<T>::setIndex(const uint index) {
    m_Index = index;
}

// Symmetry plain vector getter.
template <typename T>
Plane ReflectionSymmetryT<T>::getPlane() const {
    return m_Plane;
}

// Line segments getter.
template <typename T>
std::vector<LineSegment> ReflectionSymmetryT<T>::getLineSegments() const {
    return m_LineSegments;
}

// Getting number of line segments.
template <typename T>
u128 ReflectionSymmetryT<T>::getLineSegmentCount() const {
    return m_LineSegments.size();
}

// Adding a new line segment.
template <typename T>
void ReflectionSymmetryT<T>::addLineSegment(const LineSegment& lineSegment) {
    m_LineSegments.push_back(lineSegment);
}

// Voxel vector getter.
template <typename T>
std::vector<Voxel>& ReflectionSymmetryT<T>::getVoxels() {
    return m_SymmetryVoxels;
}

// Voxel count getter.
template <typename T>
u128 ReflectionSymmetryT<T>::getVoxelCount() const {
    return m_SymmetryVoxels.size();
}

// Symmetry voxel adder.
template <typename T>
void ReflectionSymmetryT<T>::addSymmetryVoxel(const Voxel& voxel) {
    m_SymmetryVoxels.push_back(voxel);
}



// OBJECT METHODS
// Calculating the point positions according to a plane.
template <typename T>
std::vector<PositionFromPlane> ReflectionSymmetryT<T>::calculatePositionsFromPlane(const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh, const Plane& plane, const bool validPlane, u128& pointCount, const bool normalize) {
    const u128 numberOfPoints = points.size();  // Getting the point vector size.
    const uint normalizeFactor = normalize ? voxelMesh.voxelSideSize : 1;

    // If there is no plane, all points are undefined.
    if (!validPlane) {
        return std::vector<PositionFromPlane>(numberOfPoints, PositionFromPlane::undefined);
    }

    // Creating a new vector with positions.
    std::vector<PositionFromPlane> positions(numberOfPoints, PositionFromPlane::undefined);

    for (u128 i = 0; i < numberOfPoints; i++) {
        // Getting each single point and its voxel coordinates.
        const Point& point = points[i];
        const Voxel voxel = PointFunctions::voxelFromPoint(point, voxelMesh);
        const auto [x, y, z] = voxel.getNormalizedCoordinates(voxelMesh);

        // If the voxel is not at least material, its position is set to undefined.
        if (!voxelVector[z][y][x].material) {
            positions[i] = PositionFromPlane::undefined;
            continue;
        }

        // If the current point lies outside of an interesting voxel, its position is set to undefined.
        if (std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxel.getNormalizedVoxel(normalizeFactor)) == m_SymmetryVoxels.end()) {
            if (x == 0 ||
                y == 0 ||
                (
                    plane.a == 1 &&
                    (
                        !Tolerance::isInTolerance(fmod(point.y / static_cast<double>(voxelMesh.voxelSideSize), voxelMesh.voxelSideSize), 0.0, 0.0001) ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y - 1][x - 1].getNormalizedVoxel(normalizeFactor)) == m_SymmetryVoxels.end() ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y - 1][x].getNormalizedVoxel(normalizeFactor)) == m_SymmetryVoxels.end()
                    )
                )
            )
            {
                positions[i] = PositionFromPlane::undefined;
                continue;
            }
            else if (
                x == 0 ||
                y == 0 ||
                (
                    plane.a == 0 &&
                    (
                        !Tolerance::isInTolerance(fmod(point.x / static_cast<double>(voxelMesh.voxelSideSize), voxelMesh.voxelSideSize), 0.0, 0.0001) ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y - 1][x - 1].getNormalizedVoxel(normalizeFactor)) == m_SymmetryVoxels.end() ||
                        std::find(m_SymmetryVoxels.begin(), m_SymmetryVoxels.end(), voxelVector[z][y][x - 1].getNormalizedVoxel(normalizeFactor)) == m_SymmetryVoxels.end())
                    )
                )
            {
                positions[i] = PositionFromPlane::undefined;
                continue;
            }
            else if (plane.a != 0 && plane.a != 1) {
                positions[i] = PositionFromPlane::undefined;
                continue;
            }
        }

        pointCount++;

        // If the point and the projection point have the same coordinates, the
        // plane goes through the point. Therefore, its position is set to center.
        const Point pointProjection = plane.calculateProjectionPoint(point, voxelMesh);
        if (point == pointProjection) {
            positions[i] = PositionFromPlane::center;
            continue;
        }

        // Calculating the center coordinates of the voxel where the point is located.
        const Point voxelCenter = voxel.centerCoordinate(voxelMesh);

        // Calculating the projection point to the symmetry plane.
        const Point projection = plane.calculateProjectionPoint(voxelCenter, voxelMesh);
        const auto [voxelProjectionX, voxelProjectionY, voxelProjectionZ] = PointFunctions::voxelFromPoint(projection, voxelMesh).getNormalizedCoordinates(voxelMesh);

        // If the projection point lies within the same voxel as the point, the position is set to center.
        // Note: voxels that are only touched by the symmetry plane on one edge are NOT center voxels.
        if (x == voxelProjectionX && !Tolerance::isInTolerance(voxel.x, projection.x, 0.001) &&
            y == voxelProjectionY && !Tolerance::isInTolerance(voxel.y, projection.y, 0.001) &&
            z == voxelProjectionZ && !Tolerance::isInTolerance(voxel.z, projection.z, 0.001)
        )
        {
            positions[i] = PositionFromPlane::center;
            continue;
        }

        // If the point is on the left side of the plane, left
        // position is added, right otherwise.
        if (plane.isPointOnTheLeftSide(point, voxelMesh)) {
            positions[i] = PositionFromPlane::left;
        }
        else {
            positions[i] = PositionFromPlane::right;
        }
    }

    return positions;
}



// REFLECTION SYMMETRY FUNCTIONS
// Finding simple symmetries.
std::vector<ReflectionSymmetry> ReflectionSymmetryFunctions::findSimpleSymmetries(const std::vector<LineSegment>& lineSegments, const double tolerance) {
    std::vector<ReflectionSymmetry> symmetries;  // Symmetry vector.

    // Splitting line segments into parts.
    auto splitLineSegments = SymmetryFunctions::splitLineSegmentVectorIntoParts(lineSegments, tolerance);

    // Iterating through all line segment vector parts.
    ReflectionSymmetry symmetry;
    for (std::vector<LineSegment>& part : splitLineSegments) {
        // Checking every possible combination.
        //#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(part.size() - 1); i++) {
            for (uint j = i + 1; j < part.size(); j++) {
                const bool found = ReflectionSymmetryFunctions::calculateSimpleReflectionSymmetry(symmetry, part[i], part[j], tolerance);
                if (found) {
                    symmetries.push_back(symmetry);
                }
            }
        }
    }

    return symmetries;
}

// Checking whether two line segments represent a reflection symmetry.
bool ReflectionSymmetryFunctions::calculateSimpleReflectionSymmetry(ReflectionSymmetry& symmetry, const LineSegment& l1, const LineSegment& l2, const double tolerance) {
    // At first, we have to connect the two line segments.
    LineSegment conn1(l1.p1, l2.p1);  // Making a first connection between line segments.
    LineSegment conn2(l1.p2, l2.p2);  // Making a second connection between line segments.

    // As we only search symmetries that lie in the same layer,
    // NULL pointer is returned if the two Z coordinates are not equal.
    if (!Tolerance::isInTolerance(l1.p1.z, l2.p1.z, tolerance)) {
        return false;
    }

    // If the connections intersect with each other, it is necessary
    // to reconect the both line segments with other two points.
    if (LineSegmentFunctions::doLineSegmentsIntersect(conn1, conn2, tolerance, true)) {
        conn1.p1 = l1.p1;
        conn1.p2 = l2.p2;
        conn2.p1 = l1.p2;
        conn2.p2 = l2.p1;
    }

    // Calculating center points of the both line segments.
    const Point center1 = conn1.getCenterPoint();
    const Point center2 = conn2.getCenterPoint();

    // Calculating vectors.
    Vector3d axis = (center1 - center2).normalize();
    Vector3d V1 = (conn1.p1 - conn1.p2).normalize();
    Vector3d V2 = (conn2.p1 - conn2.p2).normalize();

    // If axis cannot be calculated, there is no symmetry.
    if (std::isnan(axis.x) || std::isnan(axis.y) || std::isnan(axis.z)) {
        return false;
    }

    // Calculating the angles between the center line segment and the connections.
    const double angle1 = std::acos(axis.dot(V1));
    const double angle2 = std::acos(axis.dot(V2));

    // If both angles are equal (with tolerance) to PI/2, a symmetry is found.
    if (
        (std::isnan(angle1) || Tolerance::isInTolerance(angle1, PI / 2, tolerance)) &&
        (std::isnan(angle2) || Tolerance::isInTolerance(angle2, PI / 2, tolerance))
    )
    {
        // Creating a plane and a vector of line segments.
        Vector3d point(center1.x, center1.y, center1.z);
        if (point.x == 0 && point.y == 0 && point.z == 0) {
            point.x = center2.x;
            point.y = center2.z;
            point.z = center2.z;
        }
        const Plane plane(point, axis, Vector3d(0, 0, 1));
        std::vector<LineSegment> lineSegments({ l1, l2 });

        symmetry = ReflectionSymmetry(plane, lineSegments);

        return true;
    }

    return false;
}

// Merging reflection symmetries.
std::vector<ReflectionSymmetry> ReflectionSymmetryFunctions::mergeSymmetries(
    std::vector<ReflectionSymmetry>& reflectionSymmetries,
    const double tolerance
)
{
    // Vector of vectors of symmetries with certain full quotients:
    // { [0, 0.25], (0.25, 0.50], (0.50, 0.75], (0.75, 1.00] }.
    std::vector<ReflectionSymmetry> mergedSymmetries;

    // Merging the found symmetries by Z axis and center point.
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
        // Searching for a potential symmetry in the list, where X and
        // Y coordinates are the same as in the current symmetry, the
        // symmetry level is also the same.
        auto itr = std::find_if(
            mergedSymmetries.begin(),
            mergedSymmetries.end(),
            [&symmetry, &tolerance](ReflectionSymmetry& curSymmetry) {
                const Plane& p1 = symmetry.getPlane();
                const Plane& p2 = curSymmetry.getPlane();

                return (
                    (
                        Tolerance::isInTolerance(p1.a, p2.a, tolerance) &&
                        Tolerance::isInTolerance(p1.b, p2.b, tolerance) &&
                        Tolerance::isInTolerance(p1.c, p2.c, tolerance) &&
                        Tolerance::isInTolerance(p1.d, p2.d, tolerance)
                    )
                    ||
                    (
                        Tolerance::isInTolerance(p1.a, -p2.a, tolerance) &&
                        Tolerance::isInTolerance(p1.b, -p2.b, tolerance) &&
                        Tolerance::isInTolerance(p1.c, -p2.c, tolerance) &&
                        Tolerance::isInTolerance(p1.d, -p2.d, tolerance)
                    )
                );
            }
        );

        // Adding the symmetry to the list if not exists.
        if (itr == mergedSymmetries.end()) {
            mergedSymmetries.push_back(symmetry);
        }
        // Updating the Symmetry object with new limits.
        else {
            const u128 i = std::distance(mergedSymmetries.begin(), itr);
            mergedSymmetries[i] &= symmetry;
        }
    }

    return mergedSymmetries;
}

// Getting voxels for each reflection symmetry.
void ReflectionSymmetryFunctions::getVoxelsInSymmetries(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelMesh& voxelMesh) {
    // Getting voxels in each reflection symmetry.
    for (ReflectionSymmetry& reflectionSymmetry : reflectionSymmetries) {
    //#pragma omp parallel for
    //for (int i = 0; i < reflectionSymmetries.size(); i++) {
        //ReflectionSymmetry reflectionSymmetry = reflectionSymmetries[i];
        const std::vector<LineSegment> lineSegments = reflectionSymmetry.getLineSegments();
        std::vector<Voxel>& voxels = reflectionSymmetry.getVoxels();

        voxels.clear();  // Clearing the previous possible artefacts in the list. Better safe than sorry.

        // Adding voxels according to reflection symmetry line segments.
        for (const LineSegment& ls : lineSegments) {
            // First line segment point.
            {
                // Getting the line segment point and calculating the voxel coordinates.
                const Point p = ls.p1;
                const uint x = static_cast<uint>(p.x / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;
                const uint y = static_cast<uint>(p.y / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;
                const uint z = static_cast<uint>(p.z / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;

                // Creating a new voxel.
                Voxel v(x, y, z, true, true);
                v.inSymmetry = true;

                // If the voxel is not already present in the symmetry, it is added to the list.
                if (std::find(voxels.begin(), voxels.end(), v) == voxels.end()) {
                    voxels.push_back(Voxel(x, y, z, true, true));
                }
            }
            // Second line segment point.
            {
                // Getting the line segment point and calculating the voxel coordinates.
                Point p = ls.p2;
                const uint x = static_cast<uint>(p.x / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;
                const uint y = static_cast<uint>(p.y / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;
                const uint z = static_cast<uint>(p.z / voxelMesh.voxelSideSize) * voxelMesh.voxelSideSize;

                // Creating a new voxel.
                Voxel v(x, y, z, true, true);
                v.inSymmetry = true;

                // If the voxel is not already present in the symmetry, it is added to the list.
                if (std::find(voxels.begin(), voxels.end(), v) == voxels.end()) {
                    voxels.push_back(Voxel(x, y, z, true, true));
                }
            }
        }
    }
}

// Adding interesting voxels back to symmetry (if in symmetry).
void ReflectionSymmetryFunctions::addMaterialVoxelsToSymmetry(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelVector& voxels, const VoxelMesh& voxelMesh, std::vector<Voxel>& materialVoxels) {
    if (materialVoxels.empty()) {
        return;
    }

    // Adding interesting voxels for each reflection symmetry.
    //#pragma omp parallel for
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
    //#pragma omp parallel for
    //for (int i = 0; i < reflectionSymmetries.size(); i++) {
        //ReflectionSymmetry symmetry = reflectionSymmetries[i];
        std::vector<std::tuple<double, Point, Voxel>> distances;  // Distances from the plane to interesting voxels and plane projection points.
        std::vector<Voxel>& symmetryVoxels = symmetry.getVoxels();   // Getting voxels in symmetry.
        const Plane plane = symmetry.getPlane();                     // Getting symmetry plane.

        // Checking voxels by each layer.
        for (uint z = 0; z < voxelMesh.voxelZ; z++) {
            distances.clear();  // Clearing distances from the previous iteration.

            for (uint y = 0; y < voxelMesh.voxelY; y++) {
                for (uint x = 0; x < voxelMesh.voxelX; x++) {
                    // Getting a voxel from the mesh according to indices.
                    const Voxel& v = voxels[z][y][x];

                    // If the voxel is not interesting, we move to the next one as soon as possible.
                    if (!v.material) {
                        continue;
                    }

                    // Checking whether the voxel is already present in the list.
                    bool exists = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), v) != symmetryVoxels.end();

                    // If the voxel is (super) interesting and is not a part of the symmetry,
                    // it is checked if the voxel has its counterpart voxel on the other side
                    // of the symmetry plane.
                    //if (!exists) {
                        // Calculating voxel indices and the plane projection point.
                        const Point voxelPoint = voxels[z][y][x].centerCoordinate(voxelMesh);
                        const Point opposite = plane.calculateOppositePoint(voxelPoint, voxelMesh);

                        // Calculation projection point voxel indices.
                        const Voxel oppositeVoxel = PointFunctions::voxelFromPoint(opposite, voxelMesh);
                        const auto [xOpposite, yOpposite, zOpposite] = oppositeVoxel.getNormalizedCoordinates(voxelMesh);

                        // Checking whether the voxel is already present in the list.
                        bool oppositeExists = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), oppositeVoxel) != symmetryVoxels.end();

                        if (VoxelFunctions::isVoxelInVoxelMesh(oppositeVoxel, voxelMesh) && voxels[zOpposite][yOpposite][xOpposite].material && !oppositeExists) {
                            if (!exists) {
                                symmetryVoxels.push_back(v);
                            }
                            if (v != oppositeVoxel) {
                                symmetryVoxels.push_back(oppositeVoxel);
                            }

                            const Point oppositeLeftUpper = opposite - Vector3d(-0.5 * voxelMesh.voxelSideSize, -0.5 * voxelMesh.voxelSideSize, 0);
                            const Voxel oppositeLeftUpperVoxel = PointFunctions::voxelFromPoint(oppositeLeftUpper, voxelMesh);
                            bool oppositeLeftUpperExists = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), oppositeLeftUpperVoxel) != symmetryVoxels.end();
                            const auto [xOppositeLeftUpper, yOppositeLeftUpper, zOppositeLeftUpper] = oppositeLeftUpperVoxel.getNormalizedCoordinates(voxelMesh);
                            if (VoxelFunctions::isVoxelInVoxelMesh(oppositeLeftUpperVoxel, voxelMesh) && voxels[zOppositeLeftUpper][yOppositeLeftUpper][xOppositeLeftUpper].material && !oppositeLeftUpperExists) {
                                symmetryVoxels.push_back(oppositeLeftUpperVoxel);
                            }

                            const Point oppositeRightUpper = opposite - Vector3d(0.5 * voxelMesh.voxelSideSize, -0.5 * voxelMesh.voxelSideSize, 0);
                            const Voxel oppositeRightUpperVoxel = PointFunctions::voxelFromPoint(oppositeRightUpper, voxelMesh);
                            bool oppositeRightUpperExists = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), oppositeRightUpperVoxel) != symmetryVoxels.end();
                            const auto [xOppositeRightUpper, yOppositeRightUpper, zOppositeRightUpper] = oppositeRightUpperVoxel.getNormalizedCoordinates(voxelMesh);
                            if (VoxelFunctions::isVoxelInVoxelMesh(oppositeRightUpperVoxel, voxelMesh) && voxels[zOppositeRightUpper][yOppositeRightUpper][xOppositeRightUpper].material && !oppositeRightUpperExists) {
                                symmetryVoxels.push_back(oppositeRightUpperVoxel);
                            }

                            const Point oppositeLeftLower = opposite - Vector3d(-0.5 * voxelMesh.voxelSideSize, 0.5 * voxelMesh.voxelSideSize, 0);
                            const Voxel oppositeLeftLowerVoxel = PointFunctions::voxelFromPoint(oppositeLeftLower, voxelMesh);
                            bool oppositeLeftLowerExists = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), oppositeLeftLowerVoxel) != symmetryVoxels.end();
                            const auto [xOppositeLeftLower, yOppositeLeftLower, zOppositeLeftLower] = oppositeLeftLowerVoxel.getNormalizedCoordinates(voxelMesh);
                            if (VoxelFunctions::isVoxelInVoxelMesh(oppositeLeftLowerVoxel, voxelMesh) && voxels[zOppositeLeftLower][yOppositeLeftLower][xOppositeLeftLower].material && !oppositeLeftLowerExists) {
                                symmetryVoxels.push_back(oppositeLeftLowerVoxel);
                            }

                            const Point oppositeRightLower = opposite - Vector3d(0.5 * voxelMesh.voxelSideSize, 0.5 * voxelMesh.voxelSideSize, 0);
                            const Voxel oppositeRightLowerVoxel = PointFunctions::voxelFromPoint(oppositeRightLower, voxelMesh);
                            bool oppositeRightLowerExists = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), oppositeRightLowerVoxel) != symmetryVoxels.end();
                            const auto [xOppositeRightLower, yOppositeRightLower, zOppositeRightLower] = oppositeRightLowerVoxel.getNormalizedCoordinates(voxelMesh);
                            if (VoxelFunctions::isVoxelInVoxelMesh(oppositeRightLowerVoxel, voxelMesh) && voxels[zOppositeRightLower][yOppositeRightLower][xOppositeRightLower].material && !oppositeRightLowerExists) {
                                symmetryVoxels.push_back(oppositeRightLowerVoxel);
                            }
                        }
                    //}
                }
            }

            // If there are no distances in the list, no pairs have to be checked. Yippee ki-yay!!!
            if (distances.empty()) {
                continue;
            }
        }
    }
}

// Clustering step.
void ReflectionSymmetryFunctions::findInterestingClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z) {
    // Checking whether the point is in voxel mesh bounds.
    if (x >= 0 && x < voxelMesh.voxelX &&
        y >= 0 && y < voxelMesh.voxelY &&
        z >= 0 && z < voxelMesh.voxelZ &&
        !voxelVector[z][y][x].checked && voxelVector[z][y][x].inSymmetry
    )
    {
        voxelVector[z][y][x].checked = true;    // Setting the current voxel to checked. 
        stack.push(voxelVector[z][y][x]);       // Pushing a voxel to the stack.
        cluster.emplace(Voxel(x, y, z, true));  // Adding a voxel to the current cluster.
    }
}

// Removing small clusters of symmetries.
void ReflectionSymmetryFunctions::removeSmallClusters(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minimumClusterSize) {
    // If the minimum cluster size equals 1, no clusters have to be removed. Yaaayyy!
    if (minimumClusterSize == 1) {
        return;
    }

    // Removing clusters for each symmetry.
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
    //#pragma omp parallel for
    //for (int i = 0; i < reflectionSymmetries.size(); i++) {
        //ReflectionSymmetry symmetry = reflectionSymmetries[i];
        std::vector<Voxel>& symmetryVoxels = symmetry.getVoxels();
        std::vector<std::set<Voxel>> clusters = symmetry.findClusters(voxels, voxelMesh);  // Finding clusters for the current symmetry.

        // Removing clusters that contain too few voxels.
        for (uint i = 0; i < clusters.size(); i++) {
            const u128 clusterSize = clusters[i].size();  // Reading the current cluster size (number of voxels).

            // If the current cluster size is smaller than the minimum cluster size, the cluster is removed from the symmetry.
            if (clusterSize < minimumClusterSize) {
                for (const Voxel& v : clusters[i]) {
                    const auto [x, y, z] = v.getNormalizedCoordinates(voxelMesh);
                    std::vector<Voxel>::iterator position = std::find(symmetryVoxels.begin(), symmetryVoxels.end(), Voxel(x, y, z, true, true));
                    if (position != symmetryVoxels.end()) {
                        symmetryVoxels.erase(position);
                    }
                }
            }
        }
    }

    // Removing empty symmetries.
    for (uint i = 0; i < reflectionSymmetries.size(); i++) {
        if (reflectionSymmetries[i].getVoxels().size() <= 0) {
            reflectionSymmetries.erase(reflectionSymmetries.begin() + i);
            i--;
        }
    }

    // Sorting reflection symmetries by voxel count in each symmetry.
    std::sort(
        reflectionSymmetries.begin(),
        reflectionSymmetries.end(),
        [](const ReflectionSymmetry& s1, const ReflectionSymmetry& s2) {
            return s1.getVoxelCount() > s2.getVoxelCount();
        }
    );

    for (uint i = 0; i < reflectionSymmetries.size(); i++) {
        reflectionSymmetries[i].setIndex(i);
    }
}

// Processing possible points that lie on the voxel edge and the symmetry plane simultaneously.
void ReflectionSymmetryFunctions::processPointsOnThePlaneAndVoxelEdge(std::vector<ReflectionSymmetry>& reflectionSymmetries, const std::vector<Point>& points, const VoxelVector& voxelVector, const VoxelMesh& voxelMesh) {
    // Processing points for each reflection symmetry.
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
        const Plane& plane = symmetry.getPlane();
        
        std::vector<Voxel>& voxels = symmetry.getVoxels();
        
        // Getting the points that lie on the plane and a voxel edge.
        const std::vector<uint>& pointsOnPlaneAndVoxelEdge = symmetry.getPlane().getPointsIndicesOnPlaneAndVoxelEdge(points, voxelMesh);

        // Processing each point.
        for (const uint pointIndex : pointsOnPlaneAndVoxelEdge) {
            const Point& p = points[pointIndex];  // Retrieving the point.

            // Getting voxel coordinates.
            uint voxelX = static_cast<uint>(p.x / voxelMesh.voxelSideSize);
            uint voxelY = static_cast<uint>(p.y / voxelMesh.voxelSideSize);
            uint voxelZ = static_cast<uint>(p.z / voxelMesh.voxelSideSize);
            if (voxelX == voxelMesh.voxelX) {
                voxelX--;
            }
            if (voxelY == voxelMesh.voxelY) {
                voxelY--;
            }
            if (voxelZ == voxelMesh.voxelZ) {
                voxelX--;
            }

            if (
                Tolerance::isInTolerance(fmod(voxelX, 1.0), 0.0, 0.0001) &&
                Tolerance::isInTolerance(fmod(voxelY, 1.0), 0.0, 0.0001) &&
                plane.a == 0
            )
            {
                if (
                    !voxelVector[voxelZ][voxelY][voxelX - 1].material &&
                    !voxelVector[voxelZ][voxelY - 1][voxelX - 1].material
                )
                {
                    voxels.push_back(Voxel(voxelX - 1, voxelY - 1, voxelZ));
                    voxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
                }

                if (
                    !voxelVector[voxelZ][voxelY][voxelX].material &&
                    !voxelVector[voxelZ][voxelY - 1][voxelX].material
                )
                {
                    voxels.push_back(Voxel(voxelX, voxelY - 1, voxelZ));
                    voxels.push_back(Voxel(voxelX, voxelY, voxelZ));
                }
            }
            else if (
                Tolerance::isInTolerance(fmod(voxelX, 1.0), 0.0, 0.0001) &&
                Tolerance::isInTolerance(fmod(voxelY, 1.0), 0.0, 0.0001) &&
                plane.a == 1
            )
            {
                if (
                    !voxelVector[voxelZ][voxelY - 1][voxelX - 1].material &&
                    !voxelVector[voxelZ][voxelY - 1][voxelX].material
                )
                {
                    voxels.push_back(Voxel(voxelX - 1, voxelY - 1, voxelZ));
                    voxels.push_back(Voxel(voxelX, voxelY - 1, voxelZ));
                }

                if (
                    !voxelVector[voxelZ][voxelY][voxelX - 1].material &&
                    !voxelVector[voxelZ][voxelY][voxelX].material
                )
                {
                    voxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
                    voxels.push_back(Voxel(voxelX, voxelY, voxelZ));
                }
            }
            else if (
                Tolerance::isInTolerance(fmod(voxelX, 1.0), 0.0, 0.0001) &&
                !voxelVector[voxelZ][voxelY][voxelX].material &&
                !voxelVector[voxelZ][voxelY][voxelX - 1].material
            )
            {
                voxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
                voxels.push_back(Voxel(voxelX, voxelY, voxelZ));
            }
            else if (
                Tolerance::isInTolerance(fmod(voxelY, 1.0), 0.0, 0.0001) &&
                !voxelVector[voxelZ][voxelY][voxelX].material &&
                !voxelVector[voxelZ][voxelY - 1][voxelX].material
            )
            {
                voxels.push_back(Voxel(voxelX - 1, voxelY, voxelZ));
                voxels.push_back(Voxel(voxelX, voxelY, voxelZ));
            }
        }
    }
}

// Removing voxels without their pairs.
void ReflectionSymmetryFunctions::removeVoxelsWithoutPairs(std::vector<ReflectionSymmetry>& reflectionSymmetries, const VoxelMesh& voxelMesh) {
    // Checking each reflection symmetry.
    for (ReflectionSymmetry& symmetry : reflectionSymmetries) {
        std::vector<Voxel>& voxels = symmetry.getVoxels();

        // Checking all voxels in the symmetry.
        for (uint i = 0; i < voxels.size(); i++) {
            // Getter of the current voxel and calculation of the projection point to the symmetry plane.
            const Voxel& v = voxels[i];
            const Point voxelPoint((v.x + 0.5) * voxelMesh.voxelSideSize, (v.y + 0.5) * voxelMesh.voxelSideSize, (v.z + 0.5) * voxelMesh.voxelSideSize);
            const Point projection = symmetry.getPlane().calculateProjectionPoint(voxelPoint, voxelMesh);

            // Calculation of the opposite voxel across the symmetry plane.
            const Vector3d vector(projection.x - voxelPoint.x, projection.y - voxelPoint.y, projection.z - voxelPoint.z);
            const Vector3d vectorMultiplied = vector * 2.0;
            const Point oppositePoint = voxelPoint + vectorMultiplied;
            const uint oppositeX = static_cast<uint>(floor((oppositePoint.x / voxelMesh.voxelSideSize) - voxelMesh.minX));
            const uint oppositeY = static_cast<uint>(floor((oppositePoint.y / voxelMesh.voxelSideSize) - voxelMesh.minY));
            const uint oppositeZ = static_cast<uint>(floor((oppositePoint.z / voxelMesh.voxelSideSize) - voxelMesh.minZ));
            const Voxel opposite(oppositeX, oppositeY, oppositeZ);

            auto it = std::find(voxels.begin(), voxels.end(), opposite);  // Searching whether the opposite voxel is a part of symmetry.

            // If a voxel that should be opposite is not a part of the symmetry,
            // the current voxel is removed from the list.
            if (it == voxels.end()) {
                voxels.erase(voxels.begin() + i);
            }
        }
    }
}

// Calculating reflection symmetries.
std::vector<ReflectionSymmetry> ReflectionSymmetryFunctions::calculateReflectionSymmetries(
    const std::vector<Point>& points,
    const VoxelVector& voxels,
    const VoxelMesh& voxelMesh,
    std::vector<Voxel>& materialVoxels,
    const double tolerance,
    const uint minimumClusterSize,
    const double minimumLineSegmentDistance
)
{
    // Algorithm for finding reflection symmetries.
    std::vector<LineSegment> lineSegments = SymmetryFunctions::calculateLineSegmentsBetweenPoints(points, voxels, voxelMesh, tolerance, minimumLineSegmentDistance);  // Calculation of line segments betweeen all pairs of points.
    std::vector<ReflectionSymmetry> symmetries = ReflectionSymmetryFunctions::findSimpleSymmetries(lineSegments, tolerance);                                          // Finding simple (trivial) symmetries according to line segments.
    std::vector<ReflectionSymmetry> mergedSymmetries = ReflectionSymmetryFunctions::mergeSymmetries(symmetries, tolerance);                                           // Merging the simple (trivial) symmetries with the given tolerance.
    ReflectionSymmetryFunctions::getVoxelsInSymmetries(mergedSymmetries, voxelMesh);                                                                                  // Adding voxels to each symmetry.
    ReflectionSymmetryFunctions::addMaterialVoxelsToSymmetry(mergedSymmetries, voxels, voxelMesh, materialVoxels);                                                 // Adding additional interesting voxels to each symmetry.
    ReflectionSymmetryFunctions::removeSmallClusters(mergedSymmetries, voxels, voxelMesh, minimumClusterSize);                                                        // Removing voxels in symmetries that are in too small clusters.
    //ReflectionSymmetryFunctions::processPointsOnThePlaneAndVoxelEdge(mergedSymmetries, points, voxels, voxelMesh);                                                  // Processing points that lie on the plane and voxel edge.

    return mergedSymmetries;
}

// Finding partial symmetries from all reflection symmetries.
std::vector<ReflectionSymmetry> ReflectionSymmetryFunctions::findPartialSymmetries(const std::vector<ReflectionSymmetry>& reflectionSymmetries, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minSymmetrySize) {
    std::vector<ReflectionSymmetry> partialSymmetries;

    // Iterating through all the reflection symmetries.
    for (const ReflectionSymmetry& symmetry : reflectionSymmetries) {
        const std::vector<std::set<Voxel>> clusters = symmetry.findClusters(voxels, voxelMesh);
        const Plane plane = symmetry.getPlane();

        // Iterating through all the clusters.
        for (const std::set<Voxel>& cluster : clusters) {
            // If the cluster size is equal or bigger than the minimal
            // symmetry size, a new symmetry is added to the list.
            if (cluster.size() >= minSymmetrySize) {
                std::vector<Voxel> voxels(cluster.begin(), cluster.end());
                //for (uint i = 0; i < voxels.size(); i++) {
                //    voxels[i] = voxels[i] * voxelMesh.voxelSideSize;
                //}

                for (const Voxel& voxel : voxels) {
                    // Calculation of voxel and projection points.
                    const double voxelX = (voxel.x + 0.5) * voxelMesh.voxelSideSize + voxelMesh.minX;
                    const double voxelY = (voxel.y + 0.5) * voxelMesh.voxelSideSize + voxelMesh.minY;
                    const double voxelZ = (voxel.z + 0.5) * voxelMesh.voxelSideSize + voxelMesh.minZ;
                    const Point voxelPoint(voxelX, voxelY, voxelZ);
                    const Point projection = plane.calculateProjectionPoint(voxelPoint, voxelMesh);

                    // Calculating the projection point voxel coordinates.
                    const uint xp = static_cast<uint>(floor((projection.x - voxelMesh.minX) / voxelMesh.voxelSideSize));  // X index.
                    const uint yp = static_cast<uint>(floor((projection.y - voxelMesh.minY) / voxelMesh.voxelSideSize));  // Y index.
                    const uint zp = static_cast<uint>(floor((projection.z - voxelMesh.minZ) / voxelMesh.voxelSideSize));  // Z index.

                    // If the voxel and the projection point voxel are the same,
                    // a reflection symmetry is also a partial symmetry.
                    if (voxel.x == xp && voxel.y == yp && voxel.z && zp) {
                        ReflectionSymmetry sym(plane, voxels);
                        sym.setIndex(static_cast<uint>(partialSymmetries.size()));
                        partialSymmetries.push_back(sym);
                        break;
                    }
                }
            }
        }
    }

    std::sort(partialSymmetries.begin(), partialSymmetries.end(), [](const ReflectionSymmetry& s1, const ReflectionSymmetry& s2) { return s1.getVoxelCount() > s2.getVoxelCount(); });

    return partialSymmetries;
}


// Calculation of nearest symmetries according to the desired symmetry plane.
std::vector<ReflectionSymmetry> ReflectionSymmetryFunctions::calculateNearestSymmetriesByPlane(const std::vector<ReflectionSymmetry>& symmetries, const Plane& desiredPlane, const VoxelMesh&) {
    std::vector<ReflectionSymmetry> nearestSymmetries(symmetries.size());
    std::vector<std::pair<double, uint>> areas(symmetries.size());

    for (int i = 0; i < symmetries.size(); i++) {
        const Plane plane = symmetries[i].getPlane();
        
        double angle = 180 * VectorFunctions::angle(plane.calculateParallelVector(), desiredPlane.calculateParallelVector()) / PI;
        if (angle > 180) {
            angle = 180 - Difference::difference(90, angle);
        }
        else if (angle > 90) {
            angle = 90 - Difference::difference(90, angle);
        }
        double diff = Difference::difference(plane.d, desiredPlane.d);
        areas[i] = std::make_pair(angle + diff, i);
    }

    std::sort(
        areas.begin(),
        areas.end(),
        [](auto area1, auto area2) {
            auto [distance1, index1] = area1;
            auto [distance2, index2] = area2;

            return distance1 < distance2;
        }
    );

    //for (auto area : areas) {
    for (int i = 0; i < areas.size(); i++) {
        const auto [distance, index] = areas[i];

        nearestSymmetries[i] = symmetries[index];
    }

    return nearestSymmetries;
}



// TEMP!!!
template class ReflectionSymmetryT<double>;