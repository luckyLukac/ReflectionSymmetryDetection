#ifndef LINESEGMENT_H
#define LINESEGMENT_H

#include <sstream>
#include <string>
#include <vector>

#include "Functions.hpp"
#include "point.h"


namespace Symmetry {
    // Line segment with two end points.
    template <typename T>
    struct LineSegmentT {
        // CLASS VARIABLES
        Point p1;      // First end point that limits the line segment.
        Point p2;      // Second end point that limits the line segment.
       
        // CONSTRUCTORS
        // Main constructor.
        LineSegmentT() :
            p1(Point()),
            p2(Point())
        {}

        // Constructor with both end points.
        LineSegmentT(const Point& p1, const Point& p2) :
            p1(p1),
            p2(p2)
        {}


        // OBJECT METHODS
        // Line segment length getter.
        T length() const {
            return static_cast<T>(std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2) + std::pow(p1.z - p2.z, 2)));
        }

        // Converting line segment to string.
        std::string toString() const {
            std::stringstream ss;
            ss << "LineSegment(" << p1.toString() << ", " << p2.toString() << ")";
            return ss.str();
        }

        // Checking whether the two line segments are equal.
        bool operator == (const LineSegmentT& ls) const {
            return (
                (p1 == ls.p1 && p2 == ls.p2) ||
                (p1 == ls.p2 && p2 == ls.p1)
            );
        }

        // Checking if line segment is smaller than the other.
        //bool operator < (const LineSegmentT& ls) const {
        //    return p1.x < ls.p1.x || p1.y < ls.p1.y || p1.z < ls.p1.z || p2.x < ls.p2.x || p2.y < ls.p2.y || p2.z < ls.p2.z;
        //}

        // Calculating the angle between the two line segments.
        double calculateAngle(const LineSegmentT& ls) const {
            const Vector3d v1 = p2 - p1;
            const Vector3d v2 = ls.p2 - ls.p1;

            return VectorFunctions::angle(v1, v2);
        }

        // Calculating the center point of the line segment.
        Point getCenterPoint() const {
            return (
                Point(
                    (p1.x + p2.x) / 2,
                    (p1.y + p2.y) / 2,
                    (p1.z + p2.z) / 2
                )
            );
        }

        // Calculation of the bisection.
        LineSegmentT<T> calculateBisection() const {
            const Vector3d lsVector = p2 - p1;                                                              // Calculating a line segment vector.
            const Vector3d perpendicular(-lsVector.y, lsVector.x, 0);                                       // Calculating a perpendicular vector to the line segment.
            const Point center = getCenterPoint();                                                          // Calculating the center point of the line segment.
            const LineSegmentT<T> bisection(center - perpendicular * 1000, center + perpendicular * 1000);  // Calculating the bisection.

            return bisection;
        }
    };


    // ALIASES
    using LineSegment = LineSegmentT<double>;
    using SplitLineSegments = std::vector<std::vector<LineSegment>>;


    // NAMESPACE FUNCTIONS
    namespace LineSegmentFunctions {
        // Returning true if the two line segments have an intersection point.
        template <typename T>
        inline bool doLineSegmentsIntersect(const LineSegmentT<T>& l1, const LineSegmentT<T>& l2, const double tolerance, const bool edgesIncluded = false) {
            // If the two line segments are parallel and at the same location,
            // they have infinite number of intersections. Yaaaaaayyyy!
            if (l1 == l2) {
                return true;
            }

            // Getting points from line segments.
            const Point P1 = l1.p1;
            const Point P2 = l1.p2;
            const Point P3 = l2.p1;
            const Point P4 = l2.p2;

            // Calculating vectors between points.
            const Vector3d V12 = P2 - P1;
            const Vector3d V34 = P4 - P3;
            const Vector3d V31 = P1 - P3;

            // Coefficient calculation.
            const T D = (V12.x * V34.y) - (V12.y * V34.x);
            const T A = (V34.x * V31.y) - (V34.y * V31.x);
            const T B = (V12.x * V31.y) - (V12.y * V31.x);

            // If D == 0, the line segments are parallel, therefore,
            // they do not have an intersection point.
            if (Tolerance::isInTolerance(D, 0.0, tolerance)) {
                return false;
            }

            // Normalizing coefficients.
            const T Ua = A / D;
            const T Ub = B / D;

            if (!edgesIncluded &&
                (Tolerance::isInTolerance(Ua, 0, 0.0001) || Tolerance::isInTolerance(Ua, 1, 0.0001) ||
                 Tolerance::isInTolerance(Ub, 0, 0.0001) || Tolerance::isInTolerance(Ub, 1, 0.0001)))
            {
                return false;
            }

            // If the two line segment intersect,
            // true is returned as a result.
            if (Ua > -0.0001 && Ua < 1.0001 && Ub > -0.0001 && Ub < 1.0001) {
                return true;
            }

            return false;
        }

        // Calculation of an intersetion point between the two line segments.
        template <typename T>
        inline bool calculateIntersetion(Point& point, const LineSegmentT<T>& l1, const LineSegmentT<T>& l2, const double tolerance = 0.0001) {
            // Getting points from line segments.
            const Point P1 = l1.p1;
            const Point P2 = l1.p2;
            const Point P3 = l2.p1;
            const Point P4 = l2.p2;

            // Calculating vectors between points.
            const Vector3d V12 = P2 - P1;
            const Vector3d V34 = P4 - P3;
            const Vector3d V31 = P1 - P3;

            // Coefficient calculation.
            const T D = (V12.x * V34.y) - (V12.y * V34.x);
            const T A = (V34.x * V31.y) - (V34.y * V31.x);
            const T B = (V12.x * V31.y) - (V12.y * V31.x);

            // If D == 0, the line segments are parallel, therefore,
            // they do not have an intersection point.
            if (Tolerance::isInTolerance(D, 0.0, tolerance)) {
                return false;
            }

            // Normalizing coefficients.
            const T Ua = A / D;
            const T Ub = B / D;

            // If the two line segment intersect,
            // true is returned as a result.
            if (Ua >= 0 && Ua <= 1 && Ub >= 0 && Ub <= 1) {
                const T x = P1.x + Ua * (P2.x - P1.x);
                const T y = P1.y + Ua * (P2.y - P1.y);
                const T z = P1.z + Ua * (P2.z - P1.z);

                point.x = x;
                point.y = y;
                point.z = z;

                return true;
            }

            return false;
        }
    };
};

#endif // LINESEGMENT_H
