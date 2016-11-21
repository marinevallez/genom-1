// Genom1_test.cpp 
//#include "vec.h"
#include <gtest/gtest.h>
#include "../src/sequence.hpp"
#include <iostream>

using namespace std;

TEST(SequenceTest, MotifSizeComparison)
{
	vector<char> seq1;
	vector<char> seq2;
	
	compare(seq1, seq2);
	EXPECT_EQ(seq1.size(), seq2.size());
	
}

TEST(SequenceTest, ReverseComplementarySize)
{
	
			vector<char> sequence("ACGTACGTACGT");
			vector<char> reverse(giveReverseComplementarySeq(sequence));
			EXPECT_EQ(sequence.size(),reverse.size());	
}

TEST(SequenceTest, MotifLengthFromBedFile)
{
	
}


/* TEST(MatrixTest, ) */


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
