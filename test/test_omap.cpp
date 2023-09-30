#include <gtest/gtest.h>

#include "rangelib/omap.hpp"

class OMapTestFixture : public ::testing::Test {
   protected:
    void SetUp() override {}

    void TearDown() override {}
};

TEST_F(OMapTestFixture, SizeConstructorFixture) {}

TEST(OMapTest, SizeConstructor)
{
    ranges::OMap TestOMap;
    EXPECT_EQ(TestOMap.width(), 0) << "Empty Constructor Width not 0, was " << TestOMap.width();
    EXPECT_EQ(TestOMap.height(), 0) << "Empty Constructor Height not 0, was " << TestOMap.height();
    EXPECT_EQ(TestOMap.filename(), "")
        << "Empty Constructor Filename not empty, was " << TestOMap.filename();

    for (const auto w : {1, 5, 7, 9, 23}) {
        for (const auto h : {2, 6, 10, 99, 133}) {
            ranges::OMap TestOMap1(w, h);
            EXPECT_EQ(TestOMap1.width(), w)
                << "Sizes Constructor Width not " << w << ", was " << TestOMap1.width();
            EXPECT_EQ(TestOMap1.height(), h)
                << "Sizes Constructor Height not " << h << ", was " << TestOMap1.height();
            EXPECT_EQ(TestOMap1.filename(), "")
                << "Sizes Constructor Filename not empty, was " << TestOMap1.filename();
        }
    }
}

TEST(OMapTest, DistTransformSizeConstructor)
{
    ranges::DistanceTransform TestDT;
    EXPECT_EQ(TestDT.width(), 0) << "Empty Constructor Width not 0, was " << TestDT.width();
    EXPECT_EQ(TestDT.height(), 0) << "Empty Constructor Height not 0, was " << TestDT.height();
    EXPECT_EQ(TestDT.filename(), "")
        << "Empty Constructor Filename not empty, was " << TestDT.filename();

    EXPECT_EQ(TestDT.grid().size(), TestDT.SDF().size())
        << "Empty Constructor SDF and Occupancy Grid sizes not the same";

    for (const auto w : {1, 5, 7, 9, 23}) {
        for (const auto h : {2, 6, 10, 99, 133}) {
            ranges::DistanceTransform TestDT_1(w, h);
            EXPECT_EQ(TestDT_1.width(), w)
                << "Sizes Constructor Width not " << w << ", was " << TestDT_1.width();
            EXPECT_EQ(TestDT_1.height(), h)
                << "Sizes Constructor Height not " << h << ", was " << TestDT_1.height();
            EXPECT_EQ(TestDT_1.filename(), "")
                << "Sizes Constructor Filename not empty, was " << TestDT_1.filename();

            EXPECT_EQ(TestDT_1.grid().size(), TestDT_1.SDF().size())
                << "Sizes Constructor SDF and Occupancy Grid sizes not the same";
        }
    }
}

TEST(OMapTest, ImageConstructor)
{
    const std::string filename = "../../maps/basement_hallways_10cm.png";
    const int e_width = 600, e_height = 600;

    ranges::OMap FilepathOMap(filename);
    EXPECT_EQ(FilepathOMap.width(), e_width)
        << "Filepath Constructor Width not " << e_width << ", was " << FilepathOMap.width();
    EXPECT_EQ(FilepathOMap.height(), e_height)
        << "Filepath Constructor Height not " << e_height << ", was " << FilepathOMap.height();
    EXPECT_EQ(FilepathOMap.filename(), filename)
        << "Filepath Constructor Filename not " << filename << ", was " << FilepathOMap.filename();
}

TEST(OMapTest, DistTransformImageConstructor)
{
    const std::string filename = "../../maps/basement_hallways_10cm.png";
    const int e_width = 600, e_height = 600;

    ranges::OMap FileOMap(filename);

    ranges::DistanceTransform OMapDT(FileOMap);

    EXPECT_EQ(OMapDT.width(), e_width)
        << "Filepath Constructor Width not " << e_width << ", was " << OMapDT.width();
    EXPECT_EQ(OMapDT.height(), e_height)
        << "Filepath Constructor Height not " << e_height << ", was " << OMapDT.height();
    EXPECT_EQ(OMapDT.filename(), filename)
        << "Filepath Constructor Filename not " << filename << ", was " << OMapDT.filename();

    EXPECT_EQ(OMapDT.grid().size(), OMapDT.SDF().size())
        << "OMap Constructor SDF and Occupancy Grid sizes not the same";

    // TODO test more distanceTransform related elements here
}

// TODO consider a test fixture