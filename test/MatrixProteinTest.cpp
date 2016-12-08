// MatrixProteinTest.cpp 
#include <gtest/gtest.h>
#include "../src/matrixprotein.hpp"
#include <iostream>

using namespace std;

const matrix absolute_PWM ({0.25, 0.5, 0.1, 0.4},
							{0.25, 0.25, 0.5, 0});
const matrix relative_PWM ({1, 0.4, 0.4, 0.2},
							{1, 0.4, 0.4, 0.2});
const matrix not_possible_matrix ({1, 2, 3, 4, 5});
							
const vector <bool> bool_absolute_PWM ({1,0});
const vector <bool> bool_relative_PWM ({1,1});

bool isEqual(const vector<char>& v1, const vector<char>& v2);



TEST(MatrixTest, CheckStatus1)
{
	vector <bool> test;
	test = matrix_status(absolute_PWM);
	ASSERT_TRUE(isEqual(test, bool_absolute_PWM)));
	
}

TEST(MatrixTest, CheckStatus2)
{
	vector <bool> test;
	test = matrix_status(relative_PWM);
	ASSERT_TRUE(isEqual(test, bool_relative_PWM)));
	
}

TEST(MatrixTest, CheckPossible)
{
	bool a(0);
	test = possible(not_possible_matrix);
	ASSERT_EQ(a,test)
	
}
	
TEST(MatrixTest, RandomlGeneratesDNA) 
{
	double length(100);
	vector <char> actg(seq_generator(100));
	ASSERT_EQ(actg.size(), length);
}


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


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
