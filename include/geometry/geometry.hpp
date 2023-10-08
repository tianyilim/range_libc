#ifndef RANGELIB_GEOM_UTILS_HPP_
#define RANGELIB_GEOM_UTILS_HPP_

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace geom {

struct Point2d {
    float x;
    float y;
};
struct Line2d {
    Point2d p1;
    Point2d p2;
};

/// @brief distance of a ray cast to a line segment.
/// @param[in] pt starting point of ray
/// @param[in] theta angle of ray in radians
/// @param[in] line segment to check for intersection
/// @param[out] result intersection point
/// @returns distance if there is an intersection, NaN if not.
double distRayCast(const Point2d& pt, const double theta, Line2d lineSeg, Point2d& res);

/// @brief distance from point to convex polygon.
float distToPolygon(const Point2d& pt, const double theta, std::vector<Point2d>& poly);

}  // namespace geom

#endif