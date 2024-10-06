#include <gtest/gtest.h>
#include "Grid2D.h"

// A simple test case for Grid2D
TEST(Grid2DTest, GridInitialization) {
	Grid2D grid(10, 10, 0.01, 0.01);
	EXPECT_EQ(grid.nx, 10);  // Test if grid size is correct
	EXPECT_EQ(grid.ny, 10);
	EXPECT_DOUBLE_EQ(grid.dx, 0.01);  // Test if spacing is correct
}
