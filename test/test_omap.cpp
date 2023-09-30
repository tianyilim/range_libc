#include <gtest/gtest.h>

#include "rangelib/omap.hpp"

TEST(OMapTest, SizeConstructor)
{
    ranges::OMap SizeConstructor_OMap;
    EXPECT_EQ(SizeConstructor_OMap.width(), 0)
        << "Empty Constructor Width not 0, was " << SizeConstructor_OMap.width();
    EXPECT_EQ(SizeConstructor_OMap.height(), 0)
        << "Empty Constructor Height not 0, was " << SizeConstructor_OMap.height();
    EXPECT_EQ(SizeConstructor_OMap.filename(), "")
        << "Empty Constructor Filename not empty, was " << SizeConstructor_OMap.filename();

    EXPECT_THROW(SizeConstructor_OMap.isOccupied(1, 0), std::invalid_argument)
        << "Accessing out of bound width did not throw!";

    EXPECT_THROW(SizeConstructor_OMap.isOccupied(0, 1), std::invalid_argument)
        << "Accessing out of bound height did not throw!";

    // Test some random numbers
    for (const auto w : {1, 5, 7, 9, 23}) {
        for (const auto h : {1, 2, 6, 10, 99, 133}) {
            ranges::OMap SizeConstructor_OMap_loop(w, h);
            EXPECT_EQ(SizeConstructor_OMap_loop.width(), w)
                << "Sizes Constructor Width not " << w << ", was "
                << SizeConstructor_OMap_loop.width();
            EXPECT_EQ(SizeConstructor_OMap_loop.height(), h)
                << "Sizes Constructor Height not " << h << ", was "
                << SizeConstructor_OMap_loop.height();
            EXPECT_EQ(SizeConstructor_OMap_loop.filename(), "")
                << "Sizes Constructor Filename not empty, was "
                << SizeConstructor_OMap_loop.filename();

            EXPECT_THROW(SizeConstructor_OMap_loop.isOccupied(w + 1, 0), std::invalid_argument)
                << "Accessing out of bound width did not throw!";

            EXPECT_THROW(SizeConstructor_OMap_loop.isOccupied(0, h + 1), std::invalid_argument)
                << "Accessing out of bound height did not throw!";

            EXPECT_NO_THROW(SizeConstructor_OMap_loop.isOccupied(w - 1, h - 1))
                << "Accessing in-bounds width and height did throw!";
        }
    }
}

TEST(OMapTest, DistTransformSizeConstructor)
{
    ranges::DistanceTransform DTSize_DT;
    EXPECT_EQ(DTSize_DT.width(), 0) << "Empty Constructor Width not 0, was " << DTSize_DT.width();
    EXPECT_EQ(DTSize_DT.height(), 0)
        << "Empty Constructor Height not 0, was " << DTSize_DT.height();
    EXPECT_EQ(DTSize_DT.filename(), "")
        << "Empty Constructor Filename not empty, was " << DTSize_DT.filename();

    EXPECT_EQ(DTSize_DT.grid().size(), DTSize_DT.SDF().size())
        << "Empty Constructor SDF and Occupancy Grid sizes not the same";

    for (const auto w : {1, 5, 7, 9, 23}) {
        for (const auto h : {1, 2, 6, 10, 99, 133}) {
            ranges::DistanceTransform DTSize_DT_1(w, h);
            EXPECT_EQ(DTSize_DT_1.width(), w)
                << "Sizes Constructor Width not " << w << ", was " << DTSize_DT_1.width();
            EXPECT_EQ(DTSize_DT_1.height(), h)
                << "Sizes Constructor Height not " << h << ", was " << DTSize_DT_1.height();
            EXPECT_EQ(DTSize_DT_1.filename(), "")
                << "Sizes Constructor Filename not empty, was " << DTSize_DT_1.filename();

            EXPECT_EQ(DTSize_DT_1.grid().size(), DTSize_DT_1.SDF().size())
                << "Sizes Constructor SDF and Occupancy Grid sizes not the same";

            EXPECT_THROW(DTSize_DT_1.isOccupied(w + 1, 0), std::invalid_argument)
                << "Accessing out of bound width did not throw!";
            EXPECT_THROW(DTSize_DT_1.isOccupied(0, h + 1), std::invalid_argument)
                << "Accessing out of bound height did not throw!";
            EXPECT_NO_THROW(DTSize_DT_1.isOccupied(w - 1, h - 1))
                << "Accessing in-bounds width and height did throw!";

            EXPECT_THROW(DTSize_DT_1.signedDistanceValue(w + 1, 0), std::invalid_argument)
                << "Accessing out of bound width did not throw!";
            EXPECT_THROW(DTSize_DT_1.signedDistanceValue(0, h + 1), std::invalid_argument)
                << "Accessing out of bound height did not throw!";
            EXPECT_NO_THROW(DTSize_DT_1.signedDistanceValue(w - 1, h - 1))
                << "Accessing in-bounds width and height did throw!";
        }
    }
}

TEST(OMapTest, ImageConstructor)
{
    const std::string filename = "../../maps/tests/symm_box.png";
    const int e_width = 120, e_height = 120;

    ranges::OMap FilepathOMap(filename);
    EXPECT_EQ(FilepathOMap.width(), e_width)
        << "Filepath Constructor Width not " << e_width << ", was " << FilepathOMap.width();
    EXPECT_EQ(FilepathOMap.height(), e_height)
        << "Filepath Constructor Height not " << e_height << ", was " << FilepathOMap.height();
    EXPECT_EQ(FilepathOMap.filename(), filename)
        << "Filepath Constructor Filename not " << filename << ", was " << FilepathOMap.filename();

    // Check occupancy at a few points.
    // Test object is a square from (21, 99) in x and y
    // x, y bounds are (0, 199)
    EXPECT_EQ(FilepathOMap.isOccupied(0, 0), false);
    EXPECT_EQ(FilepathOMap.isOccupied(0, 21), false);
    EXPECT_EQ(FilepathOMap.isOccupied(0, 21), FilepathOMap.isOccupied(21, 0));
    EXPECT_EQ(FilepathOMap.isOccupied(0, 50), false);
    EXPECT_EQ(FilepathOMap.isOccupied(0, 119), false);
    EXPECT_EQ(FilepathOMap.isOccupied(21, 21), true);
    EXPECT_EQ(FilepathOMap.isOccupied(21, 0), false);
    EXPECT_EQ(FilepathOMap.isOccupied(21, 50), true);
    EXPECT_EQ(FilepathOMap.isOccupied(21, 99), true);
    EXPECT_EQ(FilepathOMap.isOccupied(119, 119), false);
    EXPECT_EQ(FilepathOMap.isOccupied(110, 0), false);
}

// TODO check ranges::OMap::save() function

TEST(OMapTest, DistTransformImageConstructor)
{
    const std::string filename = "../../maps/tests/symm_box.png";
    const int e_width = 120, e_height = 120;

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

    // Check distance transform at a few points.
    // Test object is a square from (21, 99) in x and y
    // x, y bounds are (0, 199)
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(0, 0), std::hypot(21, 21));      // top left corner
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(0, 119), std::hypot(21, 20));    // top right corner
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(119, 0), std::hypot(20, 21));    // bot left corner
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(119, 119), std::hypot(20, 20));  // bot right corner

    // Corners of square
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(20, 20), std::hypot(1, 1));
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(20, 100), std::hypot(1, 1));
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(100, 100), std::hypot(1, 1));
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(100, 20), std::hypot(1, 1));

    // Various other tests
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(0, 50), 21.0);  // along the border of sq
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(21, 21), 0.0);  // top left of square
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(21, 0), 21.0);
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(21, 50), 0.0);  // within square
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(21, 99), 0.0);
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(21, 100), 1.0);
    EXPECT_FLOAT_EQ(OMapDT.signedDistanceValue(10, 0), std::hypot(11, 21));
}

TEST(OMapTest, EdgeMapTest)
{
    const std::string filename = "../../maps/tests/symm_box.png";

    ranges::OMap EdgeMap_OMap(filename);

    ranges::OMap edgemap;
    ASSERT_NO_THROW(edgemap = EdgeMap_OMap.makeEdgeMap());

    EXPECT_EQ(EdgeMap_OMap.grid().size(), edgemap.grid().size())
        << "OMap and Edge Map sizes not the same";

    // TODO check functionality of this edge map
}
