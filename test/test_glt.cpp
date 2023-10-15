#include <gtest/gtest.h>

#include "RangeLib.h"
#include "rangelib/lookup_table.hpp"

class GiantLUTTest : public ::testing::Test {
   protected:
    void SetUp() override
    {
        m_GLTTestOMap.reset(new ranges::OMap(m_FILENAME));
        m_GiantLUT.reset(new ranges::GiantLUTCast(*m_GLTTestOMap, m_MAX_RANGE, m_THETA_DISC));
    }

    inline float degToRad(float r) { return r * M_PI / 180.0; }

    const std::string m_FILENAME = "../../maps/tests/symm_box.png";

    const float m_MAX_RANGE = 25.0f;
    const unsigned m_THETA_DISC = 90;

    std::unique_ptr<ranges::OMap> m_GLTTestOMap;
    std::unique_ptr<ranges::GiantLUTCast> m_GiantLUT;
};

TEST_F(GiantLUTTest, SimpleConstructor)
{
    EXPECT_NO_THROW(m_GiantLUT->memory());
    EXPECT_NO_THROW(m_GiantLUT->report());
    EXPECT_FLOAT_EQ(m_GiantLUT->maxRange(), m_MAX_RANGE);

    // test equality of distance transform
    EXPECT_EQ(m_GiantLUT->getMap().height(), m_GLTTestOMap->height());
    EXPECT_EQ(m_GiantLUT->getMap().width(), m_GLTTestOMap->width());
    EXPECT_EQ(m_GiantLUT->getMap().filename(), m_GLTTestOMap->filename());
    EXPECT_EQ(m_GiantLUT->getMap().worldValues(), m_GLTTestOMap->worldValues());
}

TEST_F(GiantLUTTest, IndexTest)
{
    EXPECT_EQ(m_GiantLUT->lutWidth(), 120);
    EXPECT_EQ(m_GiantLUT->lutHeight(), 120);
    EXPECT_EQ(m_GiantLUT->lutThetaDiscretization(), m_THETA_DISC);
    EXPECT_EQ(m_GiantLUT->lutArraySize(), 120 * 120 * m_THETA_DISC);
    EXPECT_EQ(m_GiantLUT->getMap().width() * m_GiantLUT->getMap().height() * m_THETA_DISC,
              m_GiantLUT->lutArraySize());

    EXPECT_EQ(m_GiantLUT->lutHeight(), m_GiantLUT->getMap().height());
    EXPECT_EQ(m_GiantLUT->lutWidth(), m_GiantLUT->getMap().width());

    unsigned i1, i2;  // index counters

    i1 = m_GiantLUT->getLutIdx(0, 0, 0);
    i2 = m_GiantLUT->getLutIdx(0, 0, 1);
    EXPECT_EQ(i1 + 1, i2);

    i1 = m_GiantLUT->getLutIdx(0, 0, 0);
    i2 = m_GiantLUT->getLutIdx(0, 1, 0);
    EXPECT_EQ(i1 + m_THETA_DISC, i2);

    i1 = m_GiantLUT->getLutIdx(0, 0, 0);
    i2 = m_GiantLUT->getLutIdx(1, 0, 0);
    EXPECT_EQ(i1 + 120 * m_THETA_DISC, i2);

    i1 = m_GiantLUT->getLutIdx(119, 119, m_THETA_DISC - 1);
    EXPECT_EQ(i1, m_GiantLUT->lutArraySize() - 1);
}

TEST_F(GiantLUTTest, ThetaDiscretizeTest)
{
    int res;
    res = m_GiantLUT->discretize_theta(degToRad(0.0));
    EXPECT_EQ(res, 0);
    res = m_GiantLUT->discretize_theta(degToRad(90.0));
    EXPECT_EQ(res, 23);
    res = m_GiantLUT->discretize_theta(degToRad(-270.0));
    EXPECT_EQ(res, 23);
    res = m_GiantLUT->discretize_theta(degToRad(270.0));
    EXPECT_EQ(res, 68);
    res = m_GiantLUT->discretize_theta(degToRad(-90.0));
    EXPECT_EQ(res, 68);
    res = m_GiantLUT->discretize_theta(degToRad(180.0));
    EXPECT_EQ(res, 45);
    res = m_GiantLUT->discretize_theta(degToRad(360.0));
    EXPECT_EQ(res, 0);
    res = m_GiantLUT->discretize_theta(degToRad(319.25));
    EXPECT_EQ(res, 80);
    res = m_GiantLUT->discretize_theta(degToRad(370.0));
    EXPECT_EQ(res, 3);
}

TEST_F(GiantLUTTest, RangeTest)
{
    EXPECT_NEAR(m_GiantLUT->calc_range(0, 0, 0), m_MAX_RANGE, 0.01);
    EXPECT_NEAR(m_GiantLUT->calc_range(0, 0, M_PIf / 2), m_MAX_RANGE, 0.01);

    // From edges of the square
    EXPECT_NEAR(m_GiantLUT->calc_range(0, 21, 0), 21.0f, 0.01);          // from edge of square (x)
    EXPECT_NEAR(m_GiantLUT->calc_range(21, 0, M_PIf / 2), 21.0f, 0.01);  // from edge of square (y)
    EXPECT_NEAR(m_GiantLUT->calc_range(119, 21, M_PIf), 20.0f, 0.01);    // from edge of square (x)
    EXPECT_NEAR(m_GiantLUT->calc_range(21, 119, -M_PIf / 2), 20.0f,
                0.01);  // from edge of square (y)

    // angles pointing out of map
    EXPECT_NEAR(m_GiantLUT->calc_range(0, 0, -M_PIf / 2), m_MAX_RANGE, 0.01);
    EXPECT_NEAR(m_GiantLUT->calc_range(119, 0, 0), m_MAX_RANGE, 0.01);
    EXPECT_NEAR(m_GiantLUT->calc_range(0, 119, M_PIf / 2), m_MAX_RANGE, 0.01);
    EXPECT_NEAR(m_GiantLUT->calc_range(119, 119, 0), m_MAX_RANGE, 0.01);

    // out-of-range
    EXPECT_NEAR(m_GiantLUT->calc_range(120, 0, 0), m_MAX_RANGE, 0.01);
    EXPECT_NEAR(m_GiantLUT->calc_range(0, 120, 0), m_MAX_RANGE, 0.01);
    EXPECT_NEAR(m_GiantLUT->calc_range(120, 120, 0), m_MAX_RANGE, 0.01);
}

TEST_F(GiantLUTTest, SensorModelSetTest)
{
    double sensorModel[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
    int tableWidth = 10;

    EXPECT_NO_THROW(m_GiantLUT->set_sensor_model(sensorModel, tableWidth));
}
