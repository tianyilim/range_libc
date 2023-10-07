#include <math.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "RangeLib.h"
#include "rangelib/lookup_table.hpp"
#include "rangelib/range_method.hpp"
#include "rangelib/ray_casting.hpp"

// some testing thinges
#define MB (1024.0 * 1024.0)

#define GRID_STEP 10
#define GRID_RAYS 40
#define GRID_SAMPLES 1
#define RANDOM_SAMPLES 200000

struct Point2d {
    float x;
    float y;
};
struct Line2d {
    Point2d p1;
    Point2d p2;
};

double distance(Line2d ray, Line2d lineSeg, Point2d& res)
{
    const auto& p1 = ray.p1;
    const auto& p2 = ray.p2;
    const auto& p3 = lineSeg.p1;
    const auto& p4 = lineSeg.p2;

    // std::cout << "P1: " << p1.x << ", " << p1.y << " ";
    // std::cout << "P2: " << p2.x << ", " << p2.y << " ";
    // std::cout << "P3: " << p3.x << ", " << p3.y << " ";
    // std::cout << "P4: " << p4.x << ", " << p4.y << std::endl;

    double den = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);
    if (std::abs(den) < 1e-6) {
        // std::cout << "Den: " << den << " Lines do not intersect" << std::endl;
        return nanf("");
    }

    double t = ((p1.x - p3.x) * (p3.y - p4.y) - (p1.y - p3.y) * (p3.x - p4.x)) / den;
    if (t < -1e-6 || t > 1) {
        // std::cout << "t: " << t << " Not in forward projection of ray" << std::endl;
        return nanf("");
    }

    double u = ((p1.x - p3.x) * (p1.y - p2.y) - (p1.y - p3.y) * (p1.x - p2.x)) / den;
    if (u < -1e-6 || u > 1) {
        // std::cout << "u: " << u << " Not within target line segement" << std::endl;
        return nanf("");
    }

    double pxt = p1.x + t * (p2.x - p1.x);
    double pyt = p1.y + t * (p2.y - p1.y);

    double pxu = p3.x + u * (p4.x - p3.x);
    double pyu = p3.y + u * (p4.y - p3.y);

    if (fabs(pxt - pxu) > 100 * std::numeric_limits<float>::epsilon() ||
        fabs(pyt - pyu) > 100 * std::numeric_limits<float>::epsilon()) {
        std::cout << "Disagreement t: (" << pxt << ", " << pyt << ") u: (" << pxu << ", " << pyu
                  << ")" << std::endl;
        return nanf("");
    }

    res.x = pxu;
    res.y = pyu;
    float range = std::hypot(pxu - p1.x, pyu - p1.y);

    // std::cout << "Intersection at " << res.x << ", " << res.y << " with range: " << range << " t:
    // ("
    //           << pxt << ", " << pyt << ") u: (" << pxu << ", " << pyu << ")" << std::endl;
    return range;
}

float distToProjection(Line2d l, std::vector<Point2d>& poly)
{
    float min_dist = std::numeric_limits<float>::max();

    for (unsigned i = 0; i < poly.size(); ++i) {
        Line2d l2 = {poly[i], poly[(i + 1) % poly.size()]};

        Point2d intPoint;
        float dist = distance(l, l2, intPoint);

        if (isfinite(dist) && dist < min_dist) {
            min_dist = dist;
        }
    }

    return min_dist;
}

int main()
{
    const std::string filename = "../../maps/tests/symm_box.png";
    ranges::WorldValues WorldVals(1.0, 0, 0.0, -120);

    // cube coords in map frame
    std::vector<Point2d> points{{21, 21}, {99, 21}, {99, 99}, {21, 99}};

    const float maxRange = 999;

    ranges::OMap BLTestOMap(filename);
    BLTestOMap.setWorldValues(WorldVals);

    ranges::BresenhamsLine Bresenham(BLTestOMap, maxRange);
    std::cout << BLTestOMap.height() << "x" << BLTestOMap.width() << std::endl;
    ranges::RayMarching RayMarching(BLTestOMap, maxRange, 0.01);
    // ranges::GiantLUTCast GiantLUTCast(BLTestOMap, maxRange, 180.0);

    // Let's view ranges
    for (float x : {0, 22}) {
        for (float y : {0, -22, -99}) {
            for (float heading = 0.05; heading < M_PI * 1.9; heading += M_PI / 4) {
                float map_x, map_y, map_heading;
                std::tie(map_x, map_y, map_heading) = Bresenham.mapCoordinates(x, y, heading);

                // Fancy modulo
                map_heading = atan2(sin(map_heading), cos(map_heading));

                Line2d ray = {
                    {map_x, map_y},
                    {map_x + maxRange * cosf(map_heading), map_y + maxRange * sinf(map_heading)}};

                // Compute analytic range. This will be useful to see the noise of each method
                float analytic = std::min(distToProjection(ray, points), maxRange);

                float ins[3] = {x, y, heading};
                float bl[1] = {maxRange}, rm[1] = {maxRange}, glt[1] = {maxRange};
                Bresenham.numpy_calc_range(ins, bl, 1);
                RayMarching.numpy_calc_range(ins, rm, 1);
                // GiantLUTCast.numpy_calc_range(ins, glt, 1);

                if (analytic != maxRange) {
                    std::cout << "Range at (map coords) x: " << map_x << " y: " << map_y
                              << " heading: " << map_heading / M_PI * 180.0 << ": BL: " << bl[0]
                              << " RM: " << rm[0] << " GLT: " << glt[0] << " Analytic: " << analytic
                              << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    // Option for using real ROS map
    // const std::string filename = "../../maps/tests/f.png";
    // ranges::WorldValues WorldVals(0.05, 0, -11.60654, -26.520793);
}