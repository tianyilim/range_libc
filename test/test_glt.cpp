#include <gtest/gtest.h>

#include "RangeLib.h"

class GiantLUTTest : public ::testing::Test {
   protected:
    void SetUp() override
    {
        ranges::WorldValues wv;
        wv.worldScale = 0.5;

        BLTestOMap = ranges::OMap(filename);
        BLTestOMap.setWorldValues(wv);
    }

    const std::string filename = "../../maps/tests/symm_box.png";

    const float maxRange = 25.0f;
    const unsigned thetaDiscretization = 90;
    ranges::OMap BLTestOMap;

    ranges::WorldValues WorldVals;
};

TEST_F(GiantLUTTest, SimpleConstructor)
{
    ranges::GiantLUTCast GiantLUT(BLTestOMap, maxRange, thetaDiscretization);

    EXPECT_NO_THROW(GiantLUT.memory());
    EXPECT_NO_THROW(GiantLUT.report());
    EXPECT_FLOAT_EQ(GiantLUT.maxRange(), maxRange);

    // test equality of distance transform
    EXPECT_EQ(GiantLUT.getMap().height(), BLTestOMap.height());
    EXPECT_EQ(GiantLUT.getMap().width(), BLTestOMap.width());
    EXPECT_EQ(GiantLUT.getMap().filename(), BLTestOMap.filename());
    EXPECT_EQ(GiantLUT.getMap().worldValues(), BLTestOMap.worldValues());
}

TEST_F(GiantLUTTest, RangeTest)
{
    ranges::GiantLUTCast GiantLUT(BLTestOMap, maxRange, thetaDiscretization);

    EXPECT_NEAR(GiantLUT.calc_range(0, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(0, 0, M_PIf / 2), maxRange, 0.01);

    // From edges of the square
    EXPECT_NEAR(GiantLUT.calc_range(0, 21, 0), 21.0f / 2, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(21, 0, M_PIf / 2), 21.0f / 2, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(119, 21, M_PIf), 20.0f / 2, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(21, 119, -M_PIf / 2), 20.0f / 2, 0.01);

    // angles pointing out of map
    // THIS FAILS
    EXPECT_NEAR(GiantLUT.calc_range(0, 0, -M_PIf / 2), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(119, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(0, 119, M_PIf / 2), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(119, 119, 0), maxRange, 0.01);

    // out-of-range
    EXPECT_NEAR(GiantLUT.calc_range(120, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(0, 120, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(120, 120, 0), maxRange, 0.01);
}

TEST_F(GiantLUTTest, OtherFunctionTest)
{
    ranges::GiantLUTCast GiantLUT(BLTestOMap, maxRange, thetaDiscretization);

    double sensorModel[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
    int tableWidth = 10;

    EXPECT_NO_THROW(GiantLUT.set_sensor_model(sensorModel, tableWidth));
}
