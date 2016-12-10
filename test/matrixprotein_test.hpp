// MatrixProtein Tests
#include <gtest/gtest.h>
#include "../src/sequence.hpp"
#include "../src/matrixprotein.hpp"
#include <iostream>

using namespace std;

//Const variables

const matrix pwmrel({{0, 1, 0},{0.25, 0.25, 0.5},{1, 0, 0},{0.5, 0.25, 0.25}});
const matrix pssmrel({{-100, 2, -100},{0, 0, 1},{2, -100, -100},{1, 0, 0}});

//Prototypes

bool isEqualM(const matrix& m1, const matrix& m2);

//Tests

TEST(MatrixProteinTest, SwapToPwm)
{
	MatrixProtein mp1;
	mp1.setpssm_rel(pssmrel);
	matrix pwmtest;
	mp1.swaptopwm(mp1.getpssm_rel());
	pwmtest = mp1.getmx();
	ASSERT_TRUE(isEqualM(pwmrel, pwmtest));	
}

TEST(MatrixProteinTest, SwapToPssm)
{
	MatrixProtein mp;
	mp.setpwm_rel(pwmrel);
	matrix pssmtest;
	mp.swaptopssm(mp.getpwm_rel());
	pssmtest = mp.getmx();
	ASSERT_TRUE(isEqualM(pssmrel, pssmtest));	
}

//Definitions

bool isEqualM(const matrix& m1, const matrix& m2)
{
	if(m1.size() != m2.size())    {return false;}
	else
	{
		for (unsigned int i(0); i < m1.size() ; ++i)
		{
			if(m1[i].size() != m2[i].size())	{return false;}
			else
			{
				for (unsigned int j(0); j < m1[i].size(); ++j)
				{
					//cout << m2[i][j] << " "; 
					if (m1[i][j] != m2[i][j])
					{
						return false;
					}
				}  
			} 
		}
	}
	return true; 
}
