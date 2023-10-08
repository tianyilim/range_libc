#include "geometry/geometry.hpp"

namespace geom {

double distRayCast(const Point2d& pt, const double theta, Line2d lineSeg, Point2d& res)
{
    const Point2d& p1 = pt;
    const Point2d& p2 = {p1.x + cosf(theta), p1.y + sinf(theta)};
    const Point2d& p3 = lineSeg.p1;
    const Point2d& p4 = lineSeg.p2;

    // check if point is on the line
    const auto distP1P3 = std::hypot(p1.x - p3.x, p1.y - p3.y);
    const auto distP1P4 = std::hypot(p1.x - p4.x, p1.y - p4.y);
    const auto distP3P4 = std::hypot(p3.x - p4.x, p3.y - p4.y);
    const auto dist1 = distP1P3 + distP1P4;

    // pt is on the line
    if (std::abs(dist1 - distP3P4) < 1e-6) {
        res = p1;
        return 0.0;
    }

    // std::cout << "P1: " << p1.x << ", " << p1.y << " ";
    // std::cout << "P2: " << p2.x << ", " << p2.y << " ";
    // std::cout << "P3: " << p3.x << ", " << p3.y << " ";
    // std::cout << "P4: " << p4.x << ", " << p4.y << std::endl;

    double den = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);
    if (std::abs(den) < 1e-12) {
        // std::cout << "Den: " << den << " Lines do not intersect" << std::endl;
        return nanf("");
    }

    double t = ((p1.x - p3.x) * (p3.y - p4.y) - (p1.y - p3.y) * (p3.x - p4.x)) / den;
    if (t < -1e-6) {
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
        std::cerr << "Disagreement t: (" << pxt << ", " << pyt << ") u: (" << pxu << ", " << pyu
                  << ")" << std::endl;
        return INFINITY;
    }

    res.x = pxu;
    res.y = pyu;
    float range = std::hypot(pxu - p1.x, pyu - p1.y);

    // std::cout << "Intersection at " << res.x << ", " << res.y << " with range: " << range << " t:
    // ("
    //           << pxt << ", " << pyt << ") u: (" << pxu << ", " << pyu << ")" << std::endl;
    return range;
}

float distToPolygon(const Point2d& pt, const double theta, std::vector<Point2d>& poly)
{
    float min_dist = std::numeric_limits<float>::max();

    for (unsigned i = 0; i < poly.size(); ++i) {
        Line2d l2 = {poly[i], poly[(i + 1) % poly.size()]};

        Point2d intPoint;
        float dist = distRayCast(pt, theta, l2, intPoint);

        if (std::isfinite(dist) && dist < min_dist) {
            min_dist = dist;
        }
    }

    return min_dist;
}
}  // namespace geom