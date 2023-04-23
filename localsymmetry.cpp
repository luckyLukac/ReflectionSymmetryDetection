#include <algorithm>
#include <set>
#include <stack>

#include "helper_classes/tolerance.h"
#include "localsymmetry.h"

using namespace Symmetry;



// GETTERS
// Point vector getter.
const std::vector<Point>& LocalSymmetry::getPoints() const {
    return m_Points;
}

// Setting points.
void LocalSymmetry::setPoints(std::vector<Point>& points) {
    m_Points = points;
    m_VoxelMesh = {};

    for (Point& lasPoint : LocalSymmetry::m_Points) {
        // Getting point coordinates.
        const double x = lasPoint.x;
        const double y = lasPoint.y;
        const double z = lasPoint.z;

        // Finding minimum and maximum coordinates
        // of points to create the bounding box.
        if (x < m_VoxelMesh.minX)
            m_VoxelMesh.minX = x;
        if (x > m_VoxelMesh.maxX)
            m_VoxelMesh.maxX = x;
        if (y < m_VoxelMesh.minY)
            m_VoxelMesh.minY = y;
        if (y > m_VoxelMesh.maxY)
            m_VoxelMesh.maxY = y;
        if (z < m_VoxelMesh.minZ)
            m_VoxelMesh.minZ = z;
        if (z > m_VoxelMesh.maxZ)
            m_VoxelMesh.maxZ = z;
    }

    // Calculation of the difference between minimum
    // and maximum point coordinates.
    m_VoxelMesh.deltaX = m_VoxelMesh.maxX - m_VoxelMesh.minX;
    m_VoxelMesh.deltaY = m_VoxelMesh.maxY - m_VoxelMesh.minY;
    m_VoxelMesh.deltaZ = m_VoxelMesh.maxZ - m_VoxelMesh.minZ;

    for (u128 i = 0; i < points.size(); i++) {
        m_Points[i].x = m_Points[i].x - m_VoxelMesh.minX;
        m_Points[i].y = m_Points[i].y - m_VoxelMesh.minY;
        m_Points[i].z = m_Points[i].z - m_VoxelMesh.minZ;
    }
    m_VoxelMesh.minX = m_VoxelMesh.minY = m_VoxelMesh.minZ = 0;
    m_VoxelMesh.maxX = m_VoxelMesh.deltaX;
    m_VoxelMesh.maxY = m_VoxelMesh.deltaY;
    m_VoxelMesh.maxZ = m_VoxelMesh.deltaZ;
}

// Voxel mesh object getter.
VoxelMesh& LocalSymmetry::getVoxelMesh() {
    return m_VoxelMesh;
}

// 3D voxel vector getter.
VoxelVector& LocalSymmetry::getVoxels() {
    return m_Voxels;
}

// Interesting voxel getters.
std::vector<Voxel>& LocalSymmetry::getInterestingVoxels() {
    return m_InterestingVoxels;
}

// Material voxel getters.
std::vector<Voxel>& LocalSymmetry::getMaterialVoxels() {
    return m_MaterialVoxels;
}



// SYMMETRY FUNCTIONS
// Voxel mesh calculation by voxel side length.
void SymmetryFunctions::calculateVoxelMeshByVoxelSideLength(VoxelVector& voxels, VoxelMesh& voxelMesh, const int userInput) {
    // Calculating X, Y, Z and total voxel count.
    const int x = static_cast<int>(floor(voxelMesh.deltaX / userInput)) + 1;  // X coordinate voxel count.
    const int y = static_cast<int>(floor(voxelMesh.deltaY / userInput)) + 1;  // Y coordinate voxel count.
    const int z = static_cast<int>(floor(voxelMesh.deltaZ / userInput)) + 1;  // Z coordinate voxel count.
    const int count = x * y * z;                                              // Total voxel count.

    // Saving the optimal values to the voxel mesh.
    voxelMesh.voxelCount = count;
    voxelMesh.voxelSideSize = userInput;
    voxelMesh.voxelX = x;
    voxelMesh.voxelY = y;
    voxelMesh.voxelZ = z;

    // Building a 3D voxel vector.
    SymmetryFunctions::buildVoxelVector(voxels, voxelMesh);
}

// Maximum voxel count calculation.
void SymmetryFunctions::calculateVoxelMeshByMaximumVoxelCount(VoxelVector& voxels, VoxelMesh& voxelMesh, const uint userInput) {
    uint bestVoxelCount = 0;
    uint bestVoxelSideSize = 0;
    uint bestVoxelX = 0;
    uint bestVoxelY = 0;
    uint bestVoxelZ = 0;

    // Fitting the voxel mesh by X.
    for (int i = 1; i < static_cast<int>(voxelMesh.deltaX); i++) {
        voxelMesh.voxelSideSize = static_cast<int>(ceil(voxelMesh.deltaX / i));
        voxelMesh.voxelX = std::min(i, static_cast<int>(ceil(voxelMesh.deltaX / voxelMesh.voxelSideSize))) + 1;
        voxelMesh.voxelY = static_cast<int>(ceil(voxelMesh.deltaY / voxelMesh.voxelSideSize)) + 1;
        voxelMesh.voxelZ = static_cast<int>(ceil(voxelMesh.deltaZ / voxelMesh.voxelSideSize)) + 1;
        voxelMesh.voxelCount = voxelMesh.voxelX * voxelMesh.voxelY * voxelMesh.voxelZ;

        // If the voxel count is larger than the
        // maximum input value, the search is over.
        if (voxelMesh.voxelCount > userInput) {
            break;
        }

        // If a better option is found, it is stored in the voxel mesh.
        if (voxelMesh.voxelCount > bestVoxelCount) {
            bestVoxelCount = voxelMesh.voxelCount;
            bestVoxelSideSize = voxelMesh.voxelSideSize;
            bestVoxelX = voxelMesh.voxelX;
            bestVoxelY = voxelMesh.voxelY;
            bestVoxelZ = voxelMesh.voxelZ;

            // Voxel edge size is minimum 1, therefore the search is over.
            if (voxelMesh.voxelSideSize == 1) {
                break;
            }
        }
    }

    // Fitting the voxel mesh by Y.
    for (int i = 1; i < static_cast<int>(voxelMesh.deltaY); i++) {
        voxelMesh.voxelSideSize = static_cast<int>(ceil(voxelMesh.deltaY / i));
        voxelMesh.voxelX = static_cast<int>(ceil(voxelMesh.deltaX / voxelMesh.voxelSideSize)) + 1;
        voxelMesh.voxelY = std::min(i, static_cast<int>(ceil(voxelMesh.deltaY / voxelMesh.voxelSideSize))) + 1;
        voxelMesh.voxelZ = static_cast<int>(ceil(voxelMesh.deltaZ / voxelMesh.voxelSideSize)) + 1;
        voxelMesh.voxelCount = voxelMesh.voxelX * voxelMesh.voxelY * voxelMesh.voxelZ;

        // If the voxel count is larger than the
        // maximum input value, the search is over.
        if (voxelMesh.voxelCount > userInput) {
            break;
        }

        // If a better option is found, it is stored in the voxel mesh.
        if (voxelMesh.voxelCount > bestVoxelCount) {
            bestVoxelCount = voxelMesh.voxelCount;
            bestVoxelSideSize = voxelMesh.voxelSideSize;
            bestVoxelX = voxelMesh.voxelX;
            bestVoxelY = voxelMesh.voxelY;
            bestVoxelZ = voxelMesh.voxelZ;

            // Voxel edge size is minimum 1, therefore the search is over.
            if (voxelMesh.voxelSideSize == 1) {
                break;
            }
        }
    }

    // Fitting the voxel mesh by Z.
    for (int i = 1; i < static_cast<int>(voxelMesh.deltaZ); i++) {
        voxelMesh.voxelSideSize = static_cast<int>(ceil(voxelMesh.deltaZ / i));
        voxelMesh.voxelX = static_cast<int>(ceil(voxelMesh.deltaX / voxelMesh.voxelSideSize)) + 1;
        voxelMesh.voxelY = static_cast<int>(ceil(voxelMesh.deltaY / voxelMesh.voxelSideSize)) + 1;
        voxelMesh.voxelZ = std::min(i, static_cast<int>(ceil(voxelMesh.deltaZ / voxelMesh.voxelSideSize))) + 1;
        voxelMesh.voxelCount = voxelMesh.voxelX * voxelMesh.voxelY * voxelMesh.voxelZ;

        // Voxel edge size is minimum 1, therefore the search is over.
        if (voxelMesh.voxelCount > userInput) {
            break;
        }

        // If a better option is found, it is stored in the voxel mesh.
        if (voxelMesh.voxelCount > bestVoxelCount) {
            bestVoxelCount = voxelMesh.voxelCount;
            bestVoxelSideSize = voxelMesh.voxelSideSize;
            bestVoxelX = voxelMesh.voxelX;
            bestVoxelY = voxelMesh.voxelY;
            bestVoxelZ = voxelMesh.voxelZ;

            // Voxel edge size is minimum 1, therefore the search is over.
            if (voxelMesh.voxelSideSize == 1) {
                break;
            }
        }
    }

    // If no solution is found, the whole scene is represented as 1 voxel.
    // Not nice for a user but nice for the algorithm :)
    if (bestVoxelCount == 0) {
        voxelMesh.voxelCount = 1;
        voxelMesh.voxelSideSize = static_cast<int>(std::max({ voxelMesh.deltaX, voxelMesh.deltaY, voxelMesh.deltaZ }));
        voxelMesh.voxelX = 1;
        voxelMesh.voxelY = 1;
        voxelMesh.voxelZ = 1;
    }
    else {
        // Saving the optimal values to the voxel mesh.
        voxelMesh.voxelCount = bestVoxelCount;
        voxelMesh.voxelSideSize = bestVoxelSideSize;
        voxelMesh.voxelX = bestVoxelX;
        voxelMesh.voxelY = bestVoxelY;
        voxelMesh.voxelZ = bestVoxelZ;
    }

    // Building a 3D voxel vector.
    SymmetryFunctions::buildVoxelVector(voxels, voxelMesh);
}

// Building a 3D voxel vector from a voxel mesh.
void SymmetryFunctions::buildVoxelVector(VoxelVector& voxels, const VoxelMesh& voxelMesh) {
    voxels = VoxelVector(voxelMesh.voxelZ, std::vector<std::vector<Voxel>>(voxelMesh.voxelY, std::vector<Voxel>(voxelMesh.voxelX, Voxel())));  // Clearing a previous obsolete voxel vector.

    // Setting voxel coordinates.
    for (uint z = 0; z < voxelMesh.voxelZ; z++) {
        for (uint y = 0; y < voxelMesh.voxelY; y++) {
            for (uint x = 0; x < voxelMesh.voxelX; x++) {
                // Calculation of voxel coordinates.
                const uint voxelX = static_cast<uint>(floor(x * voxelMesh.voxelSideSize));  // X coordinate.
                const uint voxelY = static_cast<uint>(floor(y * voxelMesh.voxelSideSize));  // Y coordinate.
                const uint voxelZ = static_cast<uint>(floor(z * voxelMesh.voxelSideSize));  // Z coordinate.

                // Setting voxel coordinates.
                voxels[z][y][x].x = voxelX;
                voxels[z][y][x].y = voxelY;
                voxels[z][y][x].z = voxelZ;
            }
        }
    }
}

// Searching interesting voxels according to points in the 3D voxel vector.
void SymmetryFunctions::findInterestingVoxels(const std::vector<Point>& points, std::vector<Voxel>& interesting, std::vector<Voxel>& material, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minClusterSize) {
    interesting.clear();                       // Reset of interesting voxels.
    material.clear();                          // Reset of material voxels.
    VoxelFunctions::clearVoxelVector(voxels);  // Clearing voxel vector.

    // Adding all interesting voxels.
    for (const Point& point : points) {
        const uint x = static_cast<uint>(floor((point.x - voxelMesh.minX) / voxelMesh.voxelSideSize));
        const uint y = static_cast<uint>(floor((point.y - voxelMesh.minY) / voxelMesh.voxelSideSize));
        const uint z = static_cast<uint>(floor((point.z - voxelMesh.minZ) / voxelMesh.voxelSideSize));

        // If the current point lies on the voxel edge, it is ignored.
        if (
            Tolerance::isInTolerance(x, (point.x - voxelMesh.minX) / voxelMesh.voxelSideSize, 0.0001) ||
            Tolerance::isInTolerance(y, (point.y - voxelMesh.minY) / voxelMesh.voxelSideSize, 0.0001)
        )
        {
            continue;
        }

        // Searching the voxel in the vector and setting
        // super-interesting (temporary) and interesting.
        if (!voxels[z][y][x].interesting) {
            voxels[z][y][x].interesting = true;
            voxels[z][y][x].material = true;

            interesting.push_back(Voxel(x * voxelMesh.voxelSideSize, y * voxelMesh.voxelSideSize, z * voxelMesh.voxelSideSize));
        }
    }

    // Removing too small clusters.
    SymmetryFunctions::removeSmallInterestingClusters(interesting, voxels, voxelMesh, minClusterSize);

    material = interesting;

    // Removing interesting voxels in the middle,
    // so that the whole planes are not searched.
    VoxelVector tempVoxels(voxels);
    for (uint z = 0; z < voxels.size(); z++) {
        for (uint y = 1; y < voxels[z].size() - 1; y++) {
            for (uint x = 1; x < voxels[z][y].size() - 1; x++) {
                if (voxels[z][y].size() >= 3 &&
                    voxels[z].size() >= 3 &&
                    tempVoxels[z][y][x].material &&
                    tempVoxels[z][y][x - 1].material &&
                    tempVoxels[z][y + 1][x].material &&
                    tempVoxels[z][y + 1][x + 1].material &&
                    tempVoxels[z][y - 1][x].material &&
                    voxels[z][y][x].interesting
                )
                {
                    voxels[z][y][x].interesting = false;

                    // Moving the voxel from interesting to material voxels.
                    interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
                }
            }
        }
    }

    // Removing super-interesting voxels in the middle,
    // so that the whole planes are not searched.
    for (uint z = 1; z < voxels.size() - 1; z++) {
        for (uint y = 0; y < voxels[z].size(); y++) {
            for (uint x = 1; x < voxels[z][y].size() - 1; x++) {
                if (
                    tempVoxels[z][y][x].material &&
                    tempVoxels[z][y][x - 1].material &&
                    tempVoxels[z + 1][y][x].material &&
                    tempVoxels[z][y][x + 1].material &&
                    tempVoxels[z - 1][y][x].material &&
                    voxels[z][y][x].interesting
                )
                {
                    voxels[z][y][x].interesting = false;

                    // Moving the voxel from interesting to material voxels.
                    interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
                }
            }
        }
    }

    // Removing super-interesting voxels in the middle,
    // so that the whole planes are not searched.
    for (uint z = 1; z < voxels.size() - 1; z++) {
        for (uint y = 1; y < voxels[z].size() - 1; y++) {
            for (uint x = 0; x < voxels[z][y].size(); x++) {
                if (
                    tempVoxels[z][y][x].material &&
                    tempVoxels[z][y - 1][x].material &&
                    tempVoxels[z + 1][y][x].material &&
                    tempVoxels[z][y + 1][x].material &&
                    tempVoxels[z - 1][y][x].material &&
                    voxels[z][y][x].interesting

                )
                {
                    voxels[z][y][x].interesting = false;

                    // Moving the voxel from interesting to material voxels.
                    interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
                }
            }
        }
    }

    // Removing super-interesting voxels in the middle,
    // so that the whole planes are not searched.
    for (uint z = 1; z < voxels.size() - 1; z++) {
        for (uint y = 1; y < voxels[z].size() - 1; y++) {
            for (uint x = 1; x < voxels[z][y].size() - 1; x++) {
                if (
                    (tempVoxels[z][y][x].material &&
                    tempVoxels[z][y][x - 1].material &&
                    tempVoxels[z][y][x + 1].material &&
                    tempVoxels[z + 1][y + 1][x].material &&
                    tempVoxels[z - 1][y - 1][x].material)
                    ||
                    (tempVoxels[z][y][x].material &&
                    tempVoxels[z][y][x - 1].material &&
                    tempVoxels[z][y][x + 1].material &&
                    tempVoxels[z - 1][y + 1][x].material &&
                    tempVoxels[z + 1][y - 1][x].material)
                    ||
                    (tempVoxels[z][y][x].material &&
                    tempVoxels[z][y - 1][x].material &&
                    tempVoxels[z - 1][y][x + 1].material &&
                    tempVoxels[z][y - 1][x].material &&
                    tempVoxels[z + 1][y][x - 1].material)
                    ||
                    (tempVoxels[z][y][x].material &&
                    tempVoxels[z][y - 1][x].material &&
                    tempVoxels[z + 1][y][x + 1].material &&
                    tempVoxels[z][y - 1][x].material &&
                    tempVoxels[z - 1][y][x - 1].material)
                    ||
                    (tempVoxels[z][y][x].material &&
                    tempVoxels[z + 1][y][x].material &&
                    tempVoxels[z - 1][y][x].material &&
                    tempVoxels[z][y + 1][x - 1].material &&
                    tempVoxels[z][y - 1][x + 1].material)
                    ||
                    (tempVoxels[z][y][x].material &&
                    tempVoxels[z + 1][y][x].material &&
                    tempVoxels[z - 1][y][x].material &&
                    tempVoxels[z][y - 1][x - 1].material &&
                    tempVoxels[z][y + 1][x + 1].material)
                    &&
                    voxels[z][y][x].interesting
                )
                {
                    voxels[z][y][x].interesting = false;

                    // Moving the voxel from interesting to material voxels.
                    interesting.erase(std::remove(interesting.begin(), interesting.end(), voxels[z][y][x]), interesting.end());
                }
            }
        }
    }
}

// Splitting a line segment vector into parts, where the lenghts of line segments is below the tolerance.
std::vector<std::vector<LineSegment>> SymmetryFunctions::splitLineSegmentVectorIntoParts(
    const std::vector<LineSegment>& lineSegments,
    const double tolerance
)
{
    // Creating a vector and pushing a first part.
    std::vector<std::vector<LineSegment>> splitLineSegments;
    splitLineSegments.push_back(std::vector<LineSegment>());

    // Iterating through all line segments.
    for (uint i = 0; i < lineSegments.size(); i++) {
        // If no elements in vector, a first line segment is added.
        if (splitLineSegments.back().size() == 0) {
            splitLineSegments.back().push_back(lineSegments[i]);
        }
        // If a line segment length is inside of the allowed tolerance,
        // it is added to the same part as the previous one.
        else if (
            Tolerance::isInTolerance(
                splitLineSegments.back()[0].length(),
                lineSegments[i].length(),
                tolerance
            )
        )
        {
            splitLineSegments.back().push_back(lineSegments[i]);
        }
        // If a line segment length is outside of the allowed tolerance,
        // a new part is created, line segment is added to the new part.
        else {
            // If the last part contains only one line
            // segment, there will be no symmetry.
            if (splitLineSegments.back().size() == 1) {
                splitLineSegments.pop_back();
            }

            splitLineSegments.push_back(std::vector<LineSegment>());
            splitLineSegments.back().push_back(lineSegments[i]);
        }
    }

    // If the last part contains only one line
    // segment, there will be no symmetry.
    if (splitLineSegments.back().size() == 1) {
        splitLineSegments.pop_back();
    }

    return splitLineSegments;
}

// Getting a normalized 3D voxel vector (voxel edge size equals 1).
VoxelVector SymmetryFunctions::getNormalisedVoxelVector(const VoxelVector& voxelVector, const VoxelMesh& voxelMesh) {
    // Vector for storing normalized voxels.
    VoxelVector normalisedVoxelVector(voxelMesh.voxelZ, std::vector<std::vector<Voxel>>(voxelMesh.voxelY, std::vector<Voxel>(voxelMesh.voxelX, Voxel())));

    for (uint z = 0; z < voxelMesh.voxelZ; z++) {
        for (uint y = 0; y < voxelMesh.voxelY; y++) {
            for (uint x = 0; x < voxelMesh.voxelX; x++) {
                normalisedVoxelVector[z][y][x] = Voxel(x, y, z, voxelVector[z][y][x].interesting, voxelVector[z][y][x].material);
            }
        }
    }

    return normalisedVoxelVector;
}

// Clustering step.
void SymmetryFunctions::findInterestingClusterItem(std::set<Voxel>& cluster, VoxelVector& voxelVector, const VoxelMesh& voxelMesh, std::stack<Voxel>& stack, const uint x, const uint y, const uint z) {
    // Checking whether the point is in voxel mesh bounds.
    if (x >= 0 && x < voxelMesh.voxelX &&
        y >= 0 && y < voxelMesh.voxelY &&
        z >= 0 && z < voxelMesh.voxelZ &&
        !voxelVector[z][y][x].checked && voxelVector[z][y][x].material
    )
    {
        voxelVector[z][y][x].checked = true;    // Setting the current voxel to checked. 
        stack.push(voxelVector[z][y][x]);       // Pushing a voxel to the stack.
        cluster.emplace(Voxel(x, y, z, true));  // Adding a voxel to the current cluster.
    }
}

// Finding clusters of symmetry voxels in symmetry.
std::vector<std::set<Voxel>> SymmetryFunctions::findInterestingClusters(VoxelVector& normalizedVoxels, const VoxelMesh& voxelMesh) {
    std::vector<std::set<Voxel>> clusters;

    // Searching for clusters in each voxel in symmetry (depth-first-search).
    for (uint z = 0; z < normalizedVoxels.size(); z++) {
        for (uint y = 0; y < normalizedVoxels[z].size(); y++) {
            for (uint x = 0; x < normalizedVoxels[z][y].size(); x++) {
                // If the current voxel has not been checked and is in symmetry, a new cluster has been found.
                if (!normalizedVoxels[z][y][x].checked && normalizedVoxels[z][y][x].material) {
                    std::set<Voxel> cluster;                                 // Creating a new cluster.
                    std::stack<Voxel> stack({ normalizedVoxels[z][y][x] });  // Creating a stack for flood fill.

                    // While stack is not empty, we push the new items on it.
                    while (!stack.empty()) {
                        // Getting the top item from the stack.
                        const Voxel top = stack.top();
                        stack.pop();

                        // Casting the coordinates of the top voxel.
                        const uint topX = static_cast<uint>(top.x / voxelMesh.voxelSideSize);
                        const uint topY = static_cast<uint>(top.y / voxelMesh.voxelSideSize);
                        const uint topZ = static_cast<uint>(top.z / voxelMesh.voxelSideSize);

                        // 26 steps to for depth-first search of the neighborhood.
                        // 9 in the upper layer, 8 (all but current) in the current layer, 9 in the lower layer.
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY - 1, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY - 1, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY - 1, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 0, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 0, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 0, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 1, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 1, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 1, topZ + 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY - 1, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY - 1, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY - 1, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 0, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 0, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 0, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 1, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 1, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 1, topZ + 0);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY - 1, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY - 1, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY - 1, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 0, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 0, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 0, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX - 1, topY + 1, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 0, topY + 1, topZ - 1);
                        findInterestingClusterItem(cluster, normalizedVoxels, voxelMesh, stack, topX + 1, topY + 1, topZ - 1);
                    }


                    clusters.push_back(cluster);  // Adding a new cluster to the list.
                }
            }
        }
    }

    return clusters;
}

// Removing small clusters of interesting voxels.
void SymmetryFunctions::removeSmallInterestingClusters(std::vector<Voxel>& interestingVoxels, VoxelVector& voxels, const VoxelMesh& voxelMesh, const uint minimumClusterSize) {
    // If minimum cluster size equals 1, no clusters have to be removed. Yaaayyy!
    if (minimumClusterSize == 1) {
        return;
    }

    // Searching for clusters of interesting voxels.
    const std::vector<std::set<Voxel>> clusters = SymmetryFunctions::findInterestingClusters(voxels, voxelMesh);

    // Removing clusters that contain too few voxels.
    for (uint i = 0; i < clusters.size(); i++) {
        if (clusters[i].size() < minimumClusterSize) {
            for (const Voxel& v : clusters[i]) {
                voxels[v.z][v.y][v.x].material = false;
                voxels[v.z][v.y][v.x].interesting = false;

                // Removing the obsolete interesting voxel from the vector.
                interestingVoxels.erase(std::remove(interestingVoxels.begin(), interestingVoxels.end(), voxels[v.z][v.y][v.x]), interestingVoxels.end());
            }
        }
    }
}

// Getting line segments between all pairs of points (fully connected graph).
std::vector<LineSegment> SymmetryFunctions::calculateLineSegmentsBetweenPoints(
    const std::vector<Point>& points,
    const VoxelVector& voxelVector,
    const VoxelMesh& voxelMesh,
    const double tolerance,
    const double minDistance
)
{
    std::vector<LineSegment> lineSegments;

    // Calculation of every pair of points.
    for (uint i = 0; i < points.size() - 1; i++) {
        for (uint j = i + 1; j < points.size(); j++) {
            // Cycles are not allowed in graph.
            // Line segments between points with greater Z coordinate
            // difference than tolerance are ignored.
            if (i != j && Tolerance::isInTolerance(points[i].z, points[j].z, tolerance)) {
                const LineSegment ls(points[i], points[j]);
                if (ls.length() < minDistance) {
                    continue;
                }

                // Calculation of the interesting ratio of the line segment in the voxel space.
                const double percent = SymmetryFunctions::calculateInterestingPercentInLineSegment(ls, voxelVector, voxelMesh);

                // If the ratio is greater or equal to 80 %, a new line segment is appended to the list.
                if (percent >= 0.8) {
                    lineSegments.push_back(ls);
                }
            }
        }
    }

    // Sorting all line segments by length.
    std::sort(
        lineSegments.begin(),
        lineSegments.end(),
        [](const LineSegment& ls1, const LineSegment& ls2) {
            return ls1.length() < ls2.length();
        }
    );

    return lineSegments;
}


// Calculating the percent of interesting line segment sections.
double SymmetryFunctions::calculateInterestingPercentInLineSegment(const LineSegment& lineSegment, const VoxelVector& voxels, const VoxelMesh& voxelMesh) {
    // Calculating the number of segments.
    uint interestingSegments = 0;
    const uint segments = static_cast<uint>(10 * lineSegment.length());

    // Calculating the vector from P1 to P2.
    const Point p1 = lineSegment.p1;
    const Point p2 = lineSegment.p2;
    const Vector3d difference = p2 - p1;

    // Marching on the line segment and calculating
    // interesting line segment sections.
    for (uint i = 0; i < segments; i++) {
        double factor = static_cast<double>(i) / segments;
        const Vector3d p = Vector3d(
            p1.x + factor * difference.x,
            p1.y + factor * difference.y,
            p1.z + factor * difference.z
        );

        // Calculating a voxel position.
        const uint x = static_cast<uint>(floor((p.x - voxelMesh.minX) / voxelMesh.voxelSideSize));
        const uint y = static_cast<uint>(floor((p.y - voxelMesh.minY) / voxelMesh.voxelSideSize));
        const uint z = static_cast<uint>(floor((p.z - voxelMesh.minZ) / voxelMesh.voxelSideSize));

        // If a voxel is interesting, the number of
        // interesting sections is incremented.
        if (voxels[z][y][x].material) {
            interestingSegments++;
        }
    }

    return static_cast<double>(interestingSegments) / segments;
}