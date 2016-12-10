// Genom1Test
#include <gtest/gtest.h>
#include "sequence_test.hpp"
#include "matrixprotein_test.hpp"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
