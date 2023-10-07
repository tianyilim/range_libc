#include <gtest/gtest.h>

#include "RangeLib.h"

class GiantLUTTest : public ::testing::Test {
   protected:
    void SetUp() override { GLTTestOMap = ranges::OMap(filename); }

    inline float degToRad(float r) { return r * M_PI / 180.0; }

    const std::string filename = "../../maps/tests/symm_box.png";

    const float maxRange = 25.0f;
    const unsigned thetaDiscretization = 90;
    ranges::OMap GLTTestOMap;

    // ranges::WorldValues WorldVals;
};

TEST_F(GiantLUTTest, SimpleConstructor)
{
    ranges::GiantLUTCast GiantLUT(GLTTestOMap, maxRange, thetaDiscretization);

    EXPECT_NO_THROW(GiantLUT.memory());
    EXPECT_NO_THROW(GiantLUT.report());
    EXPECT_FLOAT_EQ(GiantLUT.maxRange(), maxRange);

    // test equality of distance transform
    EXPECT_EQ(GiantLUT.getMap().height(), GLTTestOMap.height());
    EXPECT_EQ(GiantLUT.getMap().width(), GLTTestOMap.width());
    EXPECT_EQ(GiantLUT.getMap().filename(), GLTTestOMap.filename());
    EXPECT_EQ(GiantLUT.getMap().worldValues(), GLTTestOMap.worldValues());
}

TEST_F(GiantLUTTest, IndexTest)
{
    ranges::GiantLUTCast GiantLUT(GLTTestOMap, maxRange, thetaDiscretization);

    EXPECT_EQ(GiantLUT.lutWidth(), 120);
    EXPECT_EQ(GiantLUT.lutHeight(), 120);
    EXPECT_EQ(GiantLUT.lutThetaDiscretization(), thetaDiscretization);
    EXPECT_EQ(GiantLUT.lutArraySize(), 120 * 120 * thetaDiscretization);
    EXPECT_EQ(GiantLUT.getMap().width() * GiantLUT.getMap().height() * thetaDiscretization,
              GiantLUT.lutArraySize());

    EXPECT_EQ(GiantLUT.lutHeight(), GiantLUT.getMap().height());
    EXPECT_EQ(GiantLUT.lutWidth(), GiantLUT.getMap().width());

    unsigned i1, i2;  // index counters

    i1 = GiantLUT.getLutIdx(0, 0, 0);
    i2 = GiantLUT.getLutIdx(0, 0, 1);
    EXPECT_EQ(i1 + 1, i2);

    i1 = GiantLUT.getLutIdx(0, 0, 0);
    i2 = GiantLUT.getLutIdx(0, 1, 0);
    EXPECT_EQ(i1 + thetaDiscretization, i2);

    i1 = GiantLUT.getLutIdx(0, 0, 0);
    i2 = GiantLUT.getLutIdx(1, 0, 0);
    EXPECT_EQ(i1 + 120 * thetaDiscretization, i2);

    i1 = GiantLUT.getLutIdx(119, 119, thetaDiscretization - 1);
    EXPECT_EQ(i1, GiantLUT.lutArraySize() - 1);
}

TEST_F(GiantLUTTest, ThetaDiscretizeTest)
{
    ranges::GiantLUTCast GiantLUT(GLTTestOMap, maxRange, thetaDiscretization);

    int res;
    res = GiantLUT.discretize_theta(degToRad(0.0));
    EXPECT_EQ(res, 0);
    res = GiantLUT.discretize_theta(degToRad(90.0));
    EXPECT_EQ(res, 23);
    res = GiantLUT.discretize_theta(degToRad(-270.0));
    EXPECT_EQ(res, 23);
    res = GiantLUT.discretize_theta(degToRad(270.0));
    EXPECT_EQ(res, 68);
    res = GiantLUT.discretize_theta(degToRad(-90.0));
    EXPECT_EQ(res, 68);
    res = GiantLUT.discretize_theta(degToRad(180.0));
    EXPECT_EQ(res, 45);
    res = GiantLUT.discretize_theta(degToRad(360.0));
    EXPECT_EQ(res, 0);
    res = GiantLUT.discretize_theta(degToRad(319.25));
    EXPECT_EQ(res, 80);
    res = GiantLUT.discretize_theta(degToRad(370.0));
    EXPECT_EQ(res, 3);
}

TEST_F(GiantLUTTest, RangeTest)
{
    ranges::GiantLUTCast GiantLUT(GLTTestOMap, maxRange, thetaDiscretization);

    EXPECT_NEAR(GiantLUT.calc_range(0, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(0, 0, M_PIf / 2), maxRange, 0.01);

    // From edges of the square
    EXPECT_NEAR(GiantLUT.calc_range(0, 21, 0), 21.0f, 0.01);             // from edge of square (x)
    EXPECT_NEAR(GiantLUT.calc_range(21, 0, M_PIf / 2), 21.0f, 0.01);     // from edge of square (y)
    EXPECT_NEAR(GiantLUT.calc_range(119, 21, M_PIf), 20.0f, 0.01);       // from edge of square (x)
    EXPECT_NEAR(GiantLUT.calc_range(21, 119, -M_PIf / 2), 20.0f, 0.01);  // from edge of square (y)

    // angles pointing out of map
    EXPECT_NEAR(GiantLUT.calc_range(0, 0, -M_PIf / 2), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(119, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(0, 119, M_PIf / 2), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(119, 119, 0), maxRange, 0.01);

    // out-of-range
    EXPECT_NEAR(GiantLUT.calc_range(120, 0, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(0, 120, 0), maxRange, 0.01);
    EXPECT_NEAR(GiantLUT.calc_range(120, 120, 0), maxRange, 0.01);
}

TEST_F(GiantLUTTest, SensorModelSetTest)
{
    ranges::GiantLUTCast GiantLUT(GLTTestOMap, maxRange, thetaDiscretization);

    double sensorModel[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
    int tableWidth = 10;

    EXPECT_NO_THROW(GiantLUT.set_sensor_model(sensorModel, tableWidth));
}
