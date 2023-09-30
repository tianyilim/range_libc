#include <gtest/gtest.h>

#include "RangeLib.h"

class RayMarchingTest : public ::testing::Test {
   protected:
    void SetUp() override { BLTestOMap = ranges::OMap(filename); }

    const std::string filename = "../../maps/tests/symm_box.png";

    const float maxRange = 25.0f;
    ranges::OMap BLTestOMap;

    ranges::WorldValues WorldVals;
};

TEST_F(RayMarchingTest, SimpleConstructor)
{
    ranges::RayMarching RayMarching(BLTestOMap, maxRange);

    EXPECT_NO_THROW(RayMarching.memory());
    EXPECT_NO_THROW(RayMarching.report());
    EXPECT_FLOAT_EQ(RayMarching.maxRange(), maxRange);

    // test equality of distance transform
    EXPECT_EQ(RayMarching.getMap().height(), BLTestOMap.height());
    EXPECT_EQ(RayMarching.getMap().width(), BLTestOMap.width());
    EXPECT_EQ(RayMarching.getMap().filename(), BLTestOMap.filename());
    EXPECT_EQ(RayMarching.getMap().worldValues(), BLTestOMap.worldValues());
}

TEST_F(RayMarchingTest, RangeTest)
{
    ranges::RayMarching RayMarching(BLTestOMap, maxRange);

    EXPECT_FLOAT_EQ(RayMarching.calc_range(0, 0, 0), maxRange);
    EXPECT_FLOAT_EQ(RayMarching.calc_range(0, 0, M_PIf / 2), maxRange);

    EXPECT_FLOAT_EQ(RayMarching.calc_range(0, 21, 0), 21.0f);             // from edge of square (x)
    EXPECT_FLOAT_EQ(RayMarching.calc_range(21, 0, M_PIf / 2), 21.0f);     // from edge of square (y)
    EXPECT_FLOAT_EQ(RayMarching.calc_range(119, 21, M_PIf), 20.0f);       // from edge of square (x)
    EXPECT_FLOAT_EQ(RayMarching.calc_range(21, 119, -M_PIf / 2), 20.0f);  // from edge of square (y)

    // angles pointing out of map
    EXPECT_FLOAT_EQ(RayMarching.calc_range(0, 0, -M_PIf / 2), maxRange);
    EXPECT_FLOAT_EQ(RayMarching.calc_range(119, 0, 0), maxRange);
    EXPECT_FLOAT_EQ(RayMarching.calc_range(0, 119, M_PIf / 2), maxRange);
    EXPECT_FLOAT_EQ(RayMarching.calc_range(119, 119, 0), maxRange);

    // out-of-range
    EXPECT_NEAR(RayMarching.calc_range(120, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(RayMarching.calc_range(0, 120, 0), maxRange, 0.01);
    EXPECT_NEAR(RayMarching.calc_range(120, 120, 0), maxRange, 0.01);
}

TEST_F(RayMarchingTest, OtherFunctionTest)
{
    ranges::RayMarching RayMarching(BLTestOMap, maxRange);

    double sensorModel[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
    int tableWidth = 10;

    EXPECT_NO_THROW(RayMarching.set_sensor_model(sensorModel, tableWidth));
}
