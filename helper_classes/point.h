#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "constants.h"
#include "tolerance.h"


namespace Symmetry {
    // Point with X, Y and Z coordinate.
    template <class T>
    struct PointT {
        // OBJECT VARIABLES
        T x;  // X coordinate.
        T y;  // Y coordinate.
        T z;  // Z coordinate.
       

        // CONSTRUCTORS
        // Default constructor without parameters.
        PointT() :
            x(static_cast<T>(0)),
            y(static_cast<T>(0)),
            z(static_cast<T>(0))
        {}

        // Constructor with X, Y and Z coordinates.
        PointT(const T x, const T y, const T z) :
            x(x),
            y(y),
            z(z)
        {}

        // Copy constructor.
        PointT(const PointT<T>& p) :
            x(p.x),
            y(p.y),
            z(p.z)
        {}


        // OVERLOADED OPERATORS
        // Sum of point coordinates and a number.
        PointT<T> operator + (const double number) const {
            return PointT<T>(x + number, y + number, z + number);
        }

        // Adding a vector to the point.
        PointT<T> operator + (const PointT<T>& p) const {
            return PointT<T>(x + p.x, y + p.y, z + p.z);
        }

        // Sum of two points.
        PointT<T> operator += (const PointT<T>& p) {
            x += p.x;
            y += p.y;
            z += p.z;

            return *this;
        }

        // Difference of point coordinates and a number.
        PointT<T> operator - (const double number) const {
            return PointT<T>(x - number, y - number, z - number);
        }

        // Difference between two points is calculated
        // as the distance between the pairs of coordinates.
        PointT<T> operator - (const PointT<T>& p) const {
            return PointT<T>(x - p.x, y - p.y, z - p.z);
        }

        // Multiplying of a point and a scalar.
        PointT<T> operator * (const double factor) const {
            return PointT<T>(factor * x, factor * y, factor * z);
        }

        // Multiplying of a vector and a point.
        PointT<T> operator * (const PointT<T>& p) const {
            return PointT<T>(p.x * x, p.y * y, p.z * z);
        }

        // Division of point coordinates and a number.
        PointT<T> operator / (const double factor) const {
            return PointT<T>(x / factor, y / factor, z / factor);
        }

        // Division of point coordinates and a number.
        PointT<T> operator /= (const double factor) {
            x /= factor;
            y /= factor;
            z /= factor;

            return *this;
        }

        // Checking whether the two points are the same.
        bool operator == (const PointT<T>& p) const {
            return Tolerance::isInTolerance(x, p.x, 0.01) && Tolerance::isInTolerance(y, p.y, 0.01) && Tolerance::isInTolerance(z, p.z, 0.01);
        }

        // Checking whether the two points are not the same.
        bool operator != (const PointT<T>& p) const {
            return !Tolerance::isInTolerance(x, p.x, 0.01) || !Tolerance::isInTolerance(y, p.y, 0.01) || !Tolerance::isInTolerance(z, p.z, 0.01);
        }
        

        // OBJECT METHODS
        // Converting a point to a string.
        std::string toString() const {
            std::stringstream ss;
            ss << "Point(" << x << "," << y << "," << z << ")";
            return ss.str();
        }

        // Dot product of two vectors.
        double dot(const PointT<T>& p1) const {
            const T xDot = p1.x * x;
            const T yDot = p1.y * y;
            const T zDot = p1.z * z;

            return xDot + yDot + zDot;
        }

        // Getting length of a vector (point).
        double length() const {
            const T xDot = x * x;
            const T yDot = y * y;
            const T zDot = z * z;
            const T dot = xDot + yDot + zDot;

            return std::sqrt(dot);
        }

        // Vector normalization.
        PointT<T> normalize() {
            const double len = length();
            return PointT<T>(x / len, y / len, z / len);
        }

        // Converting the point to a certain type.
        template <typename U>
        PointT<U> convert() const {
            return PointT<U>(static_cast<U>(x), static_cast<U>(y), static_cast<U>(z));
        }
    };


    // ALIASES
    using IntPoint = PointT<int>;
    using Point = PointT<double>;
    using Vector3f = PointT<float>;
    using Vector3d = PointT<double>;
};

#endif // POINT_H
