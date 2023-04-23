#ifndef PLANE_H
#define PLANE_H

#include <iomanip>
#include <string>

#include "point.h"
#include "voxel.h"


namespace Symmetry {
    // Class for plane (ax + by + cz - d == 0).
    template <typename T>
    struct PlaneT {
        // OBJECT VARIABLES
        T a;
        T b;
        T c;
        T d;


        // CONSTRUCTORS
        // Default constructor with no arguments.
        PlaneT() {
            a = b = c = d = static_cast<T>(std::numeric_limits<T>::max());
        }

        // Default constructor with two vectors and a point.
        PlaneT(const Point& p, const Vector3d& v1, const Vector3d& v2) {
            // Calculating a, b, c and d coefficients.
            // According to Linear Algebra (TM), (a, b, c) represents the plane normal vector.
            const Vector3d normal = VectorFunctions::crossProduct(v1, v2);
            a = static_cast<T>(normal.x);
            b = static_cast<T>(normal.y);
            c = static_cast<T>(normal.z);
            d = static_cast<T>(p.dot(normal));

            // To make our lives easier, A coefficient should always be >=0.
            if (a < 0) {
                a = static_cast<T>(-a);
                b = static_cast<T>(-b);
                c = static_cast<T>(-c);
                d = static_cast<T>(-d);
            }
        }


        // OVERLOADED OPERATORS
        // Checking if the two planes have the same coefficients.
        bool operator == (const PlaneT<T>& p) const {
            return (
                (a == p.a || a == -p.a) &&
                (b == p.b || b == -p.b) &&

                (c == p.c || c == -p.c) &&
                (d == p.d || d == -p.d)
            );
        }


        // OBJECT METHODS
        // Converting a plane to string.
        std::string toString() const {
            std::stringstream ss;
            ss << std::setprecision(2);  // Setting the maximum of 2 decimal places.

            // A coefficient output.
            ss << a << "x";

            // B coefficient output.
            if (b >= 0) {
                ss << "+" << std::abs(b) << "y";
            }
            else {
                ss << b << "y";
            }

            // C coefficient output.
            if (c >= 0) {
                ss << "+" << std::abs(c) << "z=";
            }
            else {
                ss << c << "z=";
            }

            // D coefficient output.
            if (d >= 0) {
                ss << std::abs(d);
            }
            else {
                ss << d;
            }

            return ss.str();
        }

        // Calculating a plane normal vector;
        Vector3d calculateNormalVector() const {
            return Vector3d(a, b, c);
        }

        // Calculating a vector, parallel to a plane.
        Vector3d calculateParallelVector() const {
            Vector3d v(-b, a, 0);  // Calculating a parallel vector.

            // To make our lives easier again, the X component
            // of the vector should be >=0. Always.
            if (v.x < 0) {
                v.x = -v.x;
                v.y = -v.y;
            }

            return v;
        }

        // Calculating a vector that goes along the Z axis.
        Vector3d calculateZVector() const {
            return Vector3d(0, 0, 1);
        }

        // Calculating a point on the plane with Y==max if possible.
        bool calculateTopBoundingBoxIntersection(Point& p, const VoxelMesh& vm) const {
            const Vector3d v = calculateParallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the Y axis, it
            // will never reach the top side of the bounding box.
            if (v.y == 0 && d != 0) {
                return false;
            }

            // Calculating the coordinates of the intersection point.
            const double x = (d - b * vm.maxY) / a;
            const double y = vm.maxY;
            const double z = 0;

            // If X lies outside the bounding box, there is no intersection.
            if (x < vm.minX || x > vm.maxX) {
                return false;
            }

            // Setting the coordinates.
            p.x = x;
            p.y = y;
            p.z = z;

            return true;
        }

        // Calculating a point on the plane with Y==min if possible.
        bool calculateBottomBoundingBoxIntersection(Point& p, const VoxelMesh& vm) const {
            const Vector3d v = calculateParallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the Y axis, it
            // will never reach the bottom side of the bounding box.
            if (v.y == 0 && d != 0) {
                return false;
            }

            // Calculating the coordinates of the intersection point.
            const double x = (d - b * vm.minY) / a;
            const double y = vm.minY;
            const double z = 0;

            // If X lies outside the bounding box, there is no intersection.
            if (x < vm.minX || x > vm.maxX) {
                return false;
            }

            // Setting the coordinates.
            p.x = x;
            p.y = y;
            p.z = z;

            return true;
        }

        // Calculating a point on the plane with X==min if possible.
        bool calculateLeftBoundingBoxIntersection(Point& p, const VoxelMesh& vm) const {
            const Vector3d v = calculateParallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the X axis, it
            // will never reach the left side of the bounding box.
            if (v.x == 0 && d != 0) {
                return false;
            }

            // Calculating the coordinates of the intersection point.
            const double x = vm.minX;
            const double y = (d - a * vm.minX) / b;
            const double z = 0;

            // If Y lies outside the bounding box, there is no intersection.
            if (y < vm.minY || y > vm.maxY) {
                return false;
            }

            // Setting the coordinates.
            p.x = x;
            p.y = y;
            p.z = z;

            return true;
        }

        // Calculating a point on the plane with X==max if possible.
        bool calculateRightBoundingBoxIntersection(Point& p, const VoxelMesh& vm) const {
            // The Y component of the point is calculated as d/b.
            const Vector3d v = calculateParallelVector();  // Calculating a parallel vector.

            // If the plane is constant along the X axis, it
            // will never reach the left side of the bounding box.
            if (v.x == 0 && d != 0) {
                return false;
            }

            // Calculating the coordinates of the intersection point.
            const double x = vm.maxX;
            const double y = (d - a * vm.maxX) / b;
            const double z = 0;

            // If Y lies outside the bounding box, there is no intersection.
            if (y < vm.minY || y > vm.maxY) {
                return false;
            }

            // Setting the coordinates.
            p.x = x;
            p.y = y;
            p.z = z;

            return true;
        }

        // Returning true if a point is on the left side of the plane.
        bool isPointOnTheLeftSide(const Point& p, const VoxelMesh& vm) const {
            Point pTop;
            Point pBottom;
            Point pLeft;

            // Calculating top, bottom and leftbounding box intersections.
            // Since all planes intersect with the bounding box at least
            // twice, 3 intersection points are enough for sure.
            const bool top = calculateTopBoundingBoxIntersection(pTop, vm);           // Calculating top intersection.
            const bool bottom = calculateBottomBoundingBoxIntersection(pBottom, vm);  // Calculating bottom intersection.
            const bool left = calculateLeftBoundingBoxIntersection(pLeft, vm);        // Calculating left intersection.
            double cross = 0.0;                                                       // Cross product between two vectors.

            // Calculation of the side with a help of the top intersection point.
            const Vector3d parallel = calculateParallelVector();
            if (top) {
                const Point p1 = pTop - parallel;
                const Vector3d v1 = (Vector3d(calculateParallelVector().x, calculateParallelVector().y, 0)).normalize();
                const Vector3d v2 = (Vector3d(p.x - p1.x, p.y - p1.y, 0)).normalize();
                cross = v1.x * v2.y - v2.x * v1.y;
            }
            // Calculation of the side with a help of the bottom intersection point.
            else if (bottom) {
                const Point p1 = pBottom - parallel;
                const Vector3d v1 = (Vector3d(calculateParallelVector().x, calculateParallelVector().y, 0)).normalize();
                const Vector3d v2 = (Vector3d(p.x - p1.x, p.y - p1.y, 0)).normalize();
                cross = v1.x * v2.y - v2.x * v1.y;
            }
            // Calculation of the side with a help of the left intersection point.
            else if (left) {
                const Point p1 = pLeft - parallel;
                const Vector3d v1 = (Vector3d(calculateParallelVector().x, calculateParallelVector().y, 0)).normalize();
                const Vector3d v2 = (Vector3d(p.x - p1.x, p.y - p1.y, 0)).normalize();
                cross = v1.x * v2.y - v2.x * v1.y;
            }

            // If cross product is <0, the point lies on the left side of the plane.
            return cross < 0;
        }

        // Calculating the projection point to the plane.
        Point calculateProjectionPoint(const Point& p, const VoxelMesh& vm) const {
            Point pTop;
            Point pBottom;
            Point pLeft;
            Point pRight;
            Point projection;

            // Calculating top, bottom and leftbounding box intersections.
            const bool top = calculateTopBoundingBoxIntersection(pTop, vm);           // Calculating the top intersection.
            const bool bottom = calculateBottomBoundingBoxIntersection(pBottom, vm);  // Calculating the bottom intersection.
            const bool left = calculateLeftBoundingBoxIntersection(pLeft, vm);        // Calculating the left intersection.
            const bool right = calculateRightBoundingBoxIntersection(pRight, vm);     // Calculating the right intersection.

            // Calculation of the projection point with a help of the top intersection point.
            if (top) {
                // Base vector calculation.
                const Vector3d v1(calculateParallelVector());
                const Vector3d v2(pTop.x - p.x + vm.minX, pTop.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pTop.x - p.x + vm.minX, pTop.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = -v1.dot(vn);
                projection.x = pTop.x + dot * v1.x * length + vm.minX;
                projection.y = pTop.y + dot * v1.y * length + vm.minX;
                projection.z = p.z;
            }
            // Calculation of the projection point with a help of the bottom intersection point.
            else if (bottom) {
                // Base vector calculation.
                const Vector3d v1(calculateParallelVector());
                const Vector3d v2(pBottom.x - p.x + vm.minX, pBottom.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pBottom.x - p.x + vm.minX, pBottom.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = -v1.dot(vn);
                projection.x = pBottom.x + dot * v1.x * length + vm.minX;
                projection.y = pBottom.y + dot * v1.y * length + vm.minY;
                projection.z = p.z;
            }
            // Calculation of the projection point with a help of the left intersection point.
            else if (left) {
                // Base vector calculation.
                const Vector3d v1(calculateParallelVector());
                const Vector3d v2(pLeft.x - p.x + vm.minX, pLeft.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pLeft.x - p.x + vm.minX, pLeft.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = std::abs(v1.dot(vn));
                projection.x = pLeft.x + dot * v1.x * length + vm.minX;
                projection.y = pLeft.y + dot * v1.y * length + vm.minY;
                projection.z = p.z;
            }
            // Calculation of the projection point with a help of the right intersection point.
            else if (right) {
                // Base vector calculation.
                const Vector3d v1(calculateParallelVector());
                const Vector3d v2(pRight.x - p.x + vm.minX, pRight.y - p.y + vm.minY, 0);
                const Vector3d vn = Vector3d(pRight.x - p.x + vm.minX, pRight.y - p.y + vm.minY, 0).normalize();
                const double length = v2.length();

                // Projection calculation.
                const double dot = std::abs(v1.dot(vn));
                projection.x = pRight.x - dot * v1.x * length + vm.minX;
                projection.y = pRight.y - dot * v1.y * length + vm.minY;
                projection.z = p.z;
            }

            return projection;
        }

        // Calculating the projection point to the plane.
        Point calculateOppositePoint(const Point& p, const VoxelMesh& vm) const {
            const Point projectionPoint = calculateProjectionPoint(p, vm);  // Calculation of the projection point.
            const Vector3d vector = projectionPoint - p;                    // Calculation of a vector from the point to the projection point.
            const Point oppositePoint = projectionPoint + vector;           // Calculation of the opposite point.

            return oppositePoint;
        }

        // Calculating the Y coordinate on the plane according to the given X coordinate.
        double calculateYCoordinateAtX(const double x) const {
            // Y can be calculated pretty easily from the plane equation (ax+by+cz=d).
            const double y = (d - a * x) / b;
            return y;
        }

        // Calculating the X coordinate on the plane according to the given Y coordinate.
        double calculateXCoordinateAtY(const double y) const {
            // X can be calculated pretty easily from the plane equation (ax+by+cz=d).
            const double x = (d - b * y) / a;
            return x;
        }

        // Getting the indices of the points that lie on the symmetry plane and voxel edge simultaneously.
        std::vector<uint> getPointsIndicesOnPlaneAndVoxelEdge(const std::vector<Point>& points, const VoxelMesh& vm) const {
            std::vector<uint> indices;

            for (uint i = 0; i < points.size(); i++) {
                // If the projection point and the point do not have
                // the same coordinates, the point does not lie on the plane.
                const Point projection = calculateProjectionPoint(points[i], vm);
                if (points[i] != projection) {
                    continue;
                }

                // If the quotient between the point X coordinate
                if (!Tolerance::isInTolerance(std::fmod(points[i].x / vm.voxelSideSize, 1), 0.0, 0.0001) &&
                    !Tolerance::isInTolerance(std::fmod(points[i].y / vm.voxelSideSize, 1), 0.0, 0.0001)
                )
                {
                    continue;
                }

                indices.push_back(i);
            }

            return indices;
        }
        
        // Calculating the starting point of the plane in a symmetry.
        std::tuple<Vector3d, double, double, double> calculateStartPoint(const VoxelMesh& vm) const {
            Point pTop;
            Point pBottom;
            Point pLeft;
            Point pRight;

            // Calculating the intersections betweeen the symmetry plane and the bounding box.
            const bool top = calculateTopBoundingBoxIntersection(pTop, vm);
            const bool bottom = calculateBottomBoundingBoxIntersection(pBottom, vm);
            const bool left = calculateLeftBoundingBoxIntersection(pLeft, vm);
            const bool right = calculateRightBoundingBoxIntersection(pRight, vm);

            // Setting the position of a plane.
            Vector3d position;
            const double RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE = 0.5;
            const double planeY = calculateParallelVector().y;
            double distance = 0.0;
            double distanceX = 0.0;
            double distanceY = 0.0;

            // Positioning the plane according to all the posibilities (8).
            if (planeY < 0) {
                // Option 1: left and bottom bounding box intersections.
                if (left && bottom) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pBottom.x - pLeft.x;
                    distanceY = pLeft.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pLeft.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 2: top and bottom bounding box intersections.
                else if (top && bottom) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pBottom.x - pTop.x;
                    distanceY = pTop.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pTop.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pTop.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 3: top and right bounding box intersections.
                else if (top && right) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pTop.x;
                    distanceY = pTop.y - pRight.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pTop.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pTop.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 4: top and right bounding box intersections.
                else if (left && right) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pLeft.x;
                    distanceY = pLeft.y - pRight.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX;
                    position.y = pLeft.y + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
            }
            else {
                // Option 5: bottom and right bounding box intersections.
                if (bottom && right) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pBottom.x;
                    distanceY = pRight.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pRight.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - 2 * distanceX;
                    position.y = pRight.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY - distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 6: bottom and top bounding box intersections.
                else if (bottom && top) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pTop.x - pBottom.x;
                    distanceY = pTop.y - pBottom.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pBottom.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - distanceX;
                    position.y = pBottom.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 7: left and top bounding box intersections.
                else if (left && top) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pTop.x - pLeft.x;
                    distanceY = pTop.y - pLeft.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - distanceX;
                    position.y = pLeft.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
                // Option 8: left and right bounding box intersections.
                else if (left && right) {
                    // Calculating the distances (by X, by Y, and total).
                    distanceX = pRight.x - pLeft.x;
                    distanceY = pRight.y - pLeft.y;
                    distance = (2 * RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE + 1) * sqrt(pow(distanceX, 2) + pow(distanceY, 2));

                    // Positioning the plane.
                    position.x = pLeft.x + RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceX - distanceX;
                    position.y = pLeft.y - RELATIVE_LENGTH_OUTSIDE_BBOX_EACH_SIDE * distanceY;
                    position.z = vm.minZ - (vm.deltaZ / 2);
                }
            }

            return std::make_tuple(position, distance, distanceX, distanceY);
        }
    };


    // ALIASES
    using Plane = PlaneT<double>;
};

#endif // PLANE_H