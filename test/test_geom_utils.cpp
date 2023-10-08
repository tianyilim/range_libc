#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "geometry/geometry.hpp"

using namespace geom;

TEST(GeometryTest, DistanceLineTest)
{
    Point2d pt{0, 0};
    Line2d line;
    double dist;
    Point2d res;

    // Regular line
    line = {{1, 1}, {1, -1}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_FLOAT_EQ(dist, 1.0);
    EXPECT_FLOAT_EQ(res.x, 1);
    EXPECT_FLOAT_EQ(res.y, 0);

    // far away
    line = {{100, 1}, {100, -1}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_FLOAT_EQ(dist, 100.0);
    EXPECT_FLOAT_EQ(res.x, 100);
    EXPECT_FLOAT_EQ(res.y, 0);

    // behind (no intersection)
    line = {{-1, 1}, {-1, -1}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_TRUE(std::isnan(dist));

    // parallel lines
    line = {{0, 1}, {1, 1}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_TRUE(std::isnan(dist));

    // no intersection
    line = {{1, 1}, {1, 2}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_TRUE(std::isnan(dist));

    // on the line
    line = {{0, 1}, {0, -1}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_FLOAT_EQ(dist, 0.0);
    EXPECT_FLOAT_EQ(res.x, 0);
    EXPECT_FLOAT_EQ(res.y, 0);

    // 90 degrees
    pt = {1, 1};
    line = {{0, 2}, {2, 2}};
    dist = distRayCast(pt, M_PI / 2, line, res);
    EXPECT_FLOAT_EQ(dist, 1.0);
    EXPECT_FLOAT_EQ(res.x, 1);
    EXPECT_FLOAT_EQ(res.y, 2);

    // -90 degrees
    pt = {1, 1};
    line = {{0, -2}, {2, -2}};
    dist = distRayCast(pt, -M_PI / 2, line, res);
    EXPECT_FLOAT_EQ(dist, 3.0);
    EXPECT_FLOAT_EQ(res.x, 1);
    EXPECT_FLOAT_EQ(res.y, -2);

    // Origin on the line (y) but different angle
    pt = {0, 0.5};
    line = {{0, 0}, {0, 1}};
    dist = distRayCast(pt, M_PI, line, res);
    EXPECT_FLOAT_EQ(dist, 0.0);
    EXPECT_FLOAT_EQ(res.x, 0);
    EXPECT_FLOAT_EQ(res.y, 0.5);

    // Origin on the line (x) but different angle
    pt = {0.5, 0};
    line = {{0, 0}, {1, 0}};
    dist = distRayCast(pt, M_PI / 2, line, res);
    EXPECT_FLOAT_EQ(dist, 0.0);
    EXPECT_FLOAT_EQ(res.x, 0.5);
    EXPECT_FLOAT_EQ(res.y, 0);

    // Collinear
    pt = {0.5, 0};
    line = {{0, 0}, {1, 0}};
    dist = distRayCast(pt, 0.0, line, res);
    EXPECT_FLOAT_EQ(dist, 0.0);
    EXPECT_FLOAT_EQ(res.x, 0.5);
    EXPECT_FLOAT_EQ(res.y, 0);

    // Origin lies on the line
    pt = {0.5, 0};
    line = {{0, 0}, {1, 0}};
    dist = distRayCast(pt, 0.1, line, res);
    EXPECT_FLOAT_EQ(dist, 0.0);
    EXPECT_FLOAT_EQ(res.x, 0.5);
    EXPECT_FLOAT_EQ(res.y, 0);

    // Collinear but different directions
    pt = {0.5, 0};
    line = {{1, 0}, {2, 0}};
    dist = distRayCast(pt, M_PI, line, res);
    EXPECT_TRUE(std::isnan(dist));
}

TEST(GeometryTest, DistToPolygonTest)
{
    std::vector<Point2d> poly = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    Point2d pt;
    float dist;

    // Points on the polygon
    pt = {0, 0};
    dist = distToPolygon(pt, 0.0, poly);
    EXPECT_EQ(dist, 0);

    pt = {0, 0.5};
    dist = distToPolygon(pt, M_PI / 2, poly);
    EXPECT_EQ(dist, 0);

    pt = {0.5, 0};
    dist = distToPolygon(pt, M_PI, poly);
    EXPECT_EQ(dist, 0);
}