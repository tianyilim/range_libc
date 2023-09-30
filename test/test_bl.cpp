#include <gtest/gtest.h>

#include "RangeLib.h"

class BresenhamTest : public ::testing::Test {
   protected:
    void SetUp() override { BLTestOMap = ranges::OMap(filename); }

    const std::string filename = "../../maps/tests/symm_box.png";

    const float maxRange = 25.0f;
    ranges::OMap BLTestOMap;

    ranges::WorldValues WorldVals;
};

TEST_F(BresenhamTest, SimpleConstructor)
{
    ranges::BresenhamsLine Bresenham(BLTestOMap, maxRange);

    EXPECT_NO_THROW(Bresenham.memory());
    EXPECT_NO_THROW(Bresenham.report());
    EXPECT_FLOAT_EQ(Bresenham.maxRange(), maxRange);

    // test equality of distance transform
    EXPECT_EQ(Bresenham.getMap().height(), BLTestOMap.height());
    EXPECT_EQ(Bresenham.getMap().width(), BLTestOMap.width());
    EXPECT_EQ(Bresenham.getMap().filename(), BLTestOMap.filename());
    EXPECT_EQ(Bresenham.getMap().worldValues(), BLTestOMap.worldValues());
}

TEST_F(BresenhamTest, RangeTest)
{
    ranges::BresenhamsLine Bresenham(BLTestOMap, maxRange);

    EXPECT_FLOAT_EQ(Bresenham.calc_range(0, 0, 0), maxRange);
    EXPECT_FLOAT_EQ(Bresenham.calc_range(0, 0, M_PIf / 2), maxRange);

    EXPECT_FLOAT_EQ(Bresenham.calc_range(0, 21, 0), 21.0f);             // from edge of square (x)
    EXPECT_FLOAT_EQ(Bresenham.calc_range(21, 0, M_PIf / 2), 21.0f);     // from edge of square (y)
    EXPECT_FLOAT_EQ(Bresenham.calc_range(119, 21, M_PIf), 20.0f);       // from edge of square (x)
    EXPECT_FLOAT_EQ(Bresenham.calc_range(21, 119, -M_PIf / 2), 20.0f);  // from edge of square (y)

    // angles pointing out of map
    EXPECT_FLOAT_EQ(Bresenham.calc_range(0, 0, -M_PIf / 2), maxRange);
    EXPECT_FLOAT_EQ(Bresenham.calc_range(119, 0, 0), maxRange);
    EXPECT_FLOAT_EQ(Bresenham.calc_range(0, 119, M_PIf / 2), maxRange);
    EXPECT_FLOAT_EQ(Bresenham.calc_range(119, 119, 0), maxRange);

    // out-of-range
    EXPECT_THROW(Bresenham.calc_range(120, 0, 0), std::invalid_argument);
    EXPECT_THROW(Bresenham.calc_range(0, 120, 0), std::invalid_argument);
    EXPECT_THROW(Bresenham.calc_range(120, 120, 0), std::invalid_argument);
}

TEST_F(BresenhamTest, OtherFunctionTest)
{
    ranges::BresenhamsLine Bresenham(BLTestOMap, maxRange);

    double sensorModel[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
    int tableWidth = 10;

    EXPECT_NO_THROW(Bresenham.set_sensor_model(sensorModel, tableWidth));
}

TEST_F(BresenhamTest, MapToWorldConversion)
{
    ranges::BresenhamsLine Bresenham(BLTestOMap, maxRange);
}