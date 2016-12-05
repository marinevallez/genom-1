// Genom1_test.cpp 
//#include "vec.h"
#include <gtest/gtest.h>
#include "../src/sequence.hpp"
#include <iostream>

using namespace std;

const vector<char> seq1({'A','C','T','T','C','G','A','T','C'});
const vector<char> seq2({'A','T','T','T','G','A','A','C','C'});
const vector<char> seq1Complement({'T','G','A','A','G','C','T','A','G'});
const vector<char> reverseSeq1({'G','A','T','C','G','A','A','G','T'});


const vector<PosDir> goodVec({{5,1,"chrT",'+',{0.0,"CCCTTTG"}},{40,1,"chrT",'-',{0.0,"CCCTTTG"}}});

bool isEqual(const vector<char>& v1, const vector<char>& v2);
bool isEqual(const vector<PosDir>& v1, const vector<PosDir>& v2);

TEST(SequenceTest, SequenceComparison)
{
	
	ASSERT_FALSE(isEqual(seq1, seq2));
	
}

TEST(SequenceTest, MotifRecognitionInFasta)
{
	Sequence seq;
	vector<PosDir> results = seq.motifRecognition("CCCTTTG", "fasta_test.fasta");
	ASSERT_TRUE(isEqual(goodVec, results));
	
}

TEST(SequenceTest, Complementary)
{
	vector<char> compl_;
	for(size_t i(0); i < seq1.size(); ++i)
	{
		compl_.push_back(giveComplementaryBase(seq1[i]));
	}
	ASSERT_TRUE(isEqual(seq1Complement, compl_));
}

TEST(SequenceTest, ReverseComplementary)
{
	Sequence seq;
	vector<char> rev_ = seq.giveReverseComplementarySeq(seq1);
	ASSERT_TRUE(isEqual(rev_, reverseSeq1));
}


/* TEST(MatrixTest, ) */

bool isEqual(const vector<char>& v1, const vector<char>& v2)
{
	if(v1.size() != v2.size())   {return false;}
	else
	{
		for(size_t i(0); i < v1.size(); ++i)
		{
			if(v1[i] != v2[i]) {return false;}
		}
		return true;
	}
}

bool isEqual(const vector<PosDir>& v1, const vector<PosDir>& v2)
{
	if(v1.size() != v2.size())   {return false;}
	else 
	{
		for(size_t i(0); i < v1.size(); ++i)
		{
			if((v1[i].pos != v2[i].pos) 
				or (v1[i].seqNb != v2[i].seqNb)
				or (v1[i].chrNb != v2[i].chrNb) 
				or (v1[i].dir != v2[i].dir)
				or ((v1[i].pattern).bScore != (v2[i].pattern).bScore)
				or ((v1[i].pattern).site != (v2[i].pattern).site))
				{
					return false;
				} 
		}
		return true;
	}
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
