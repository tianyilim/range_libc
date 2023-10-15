#include <chrono>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <limits>

#include "RangeLib.h"
#include "geometry/geometry.hpp"
#include "rangelib/lookup_table.hpp"
#include "rangelib/range_method.hpp"
#include "rangelib/ray_casting.hpp"

// some testing thinges
#define MB (1024.0 * 1024.0)

#define GRID_STEP 10
#define GRID_RAYS 40
#define GRID_SAMPLES 1
#define RANDOM_SAMPLES 200000

int main()
{
    const std::string filename = "../../maps/tests/symm_box.png";
    ranges::WorldValues WorldVals(1.0, 0, 0.0, -120);

    // cube coords in map frame
    std::vector<geom::Point2d> points{{21, 21}, {99, 21}, {99, 99}, {21, 99}};

    const float maxRange = 50;

    ranges::OMap BLTestOMap(filename);
    BLTestOMap.setWorldValues(WorldVals);

    ranges::BresenhamsLine Bresenham(BLTestOMap, maxRange);
    std::cout << BLTestOMap.height() << "x" << BLTestOMap.width() << std::endl;
    ranges::RayMarching RayMarching(BLTestOMap, maxRange, 0.01);
    // ranges::GiantLUTCast GiantLUTCast(BLTestOMap, maxRange, 180.0);

    // TODO look into how the choice of theta discretization changes things

    // Let's view ranges
    for (float x : {0, 21}) {
        for (float y : {0, -21, -99}) {
            for (float heading = 0.05; heading < M_PI * 1.9; heading += M_PI / 4) {
                float map_x, map_y, map_heading;
                std::tie(map_x, map_y, map_heading) = Bresenham.mapCoordinates(x, y, heading);

                // Fancy modulo
                map_heading = atan2(sin(map_heading), cos(map_heading));

                // geom::Line2d ray = {
                //     {map_x, map_y},
                //     {map_x + maxRange * cosf(map_heading), map_y + maxRange *
                //     sinf(map_heading)}};

                // Compute analytic range. This will be useful to see the noise of each method
                float analytic =
                    std::min(geom::distToPolygon({map_x, map_y}, map_heading, points), maxRange);

                ranges::VectorFloat glt(1);
                ranges::PosesFloat ins(3, 1);

                ins(0, 0) = x;
                ins(0, 1) = y;
                ins(0, 2) = heading;

                ranges::VectorFloat bl = Bresenham.batchCalcRange(ins);
                ranges::VectorFloat rm = RayMarching.batchCalcRange(ins);
                // glt =  GiantLUTCast.numpy_calc_range(ins);

                if (analytic != maxRange) {
                    std::cout << "Range at (map coords) x: " << map_x << " y: " << map_y
                              << " heading: " << map_heading / M_PI * 180.0 << ": BL: " << bl(0)
                              << " RM: " << rm(0) << " GLT: " << glt(0) << " Analytic: " << analytic
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