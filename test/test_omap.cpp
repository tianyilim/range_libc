#include <gtest/gtest.h>

#include "RangeLib.h"

TEST(OMapTest, SizeConstructor)
{
    ranges::OMap TestOMap;
    EXPECT_EQ(TestOMap.width(), 0) << "Empty Constructor Width not 0, was " << TestOMap.width();
    EXPECT_EQ(TestOMap.height(), 0) << "Empty Constructor Height not 0, was " << TestOMap.height();
    EXPECT_EQ(TestOMap.filename(), "") << "Empty Constructor Filename not empty, was " << TestOMap.filename();

    for (const auto w : {1, 5, 7, 9, 23}) {
        for (const auto h : {2, 6, 10, 99, 133}) {
            ranges::OMap TestOMap1(w, h);
            EXPECT_EQ(TestOMap1.width(), w) << "Sizes Constructor Width not " << w << ", was " << TestOMap.width();
            EXPECT_EQ(TestOMap1.height(), h) << "Sizes Constructor Height not " << h << ", was " << TestOMap.height();
            EXPECT_EQ(TestOMap1.filename(), "") << "Sizes Constructor Filename not empty, was " << TestOMap.filename();
        }
    }
}
