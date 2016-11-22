#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

//This code gives 3 functions, one that checks if a PSSM is absolute or relative, an other one
//that checks if a matrix is a PWM or a PSSM and an other one that checks if the matrix is possible.
//I did a  4th function that uses the 3 previous ones to give a matrix status, it indicates if it's
//an impossible matrix, an absolute PSSM, a relative PSSM, an absolute PWM or a relative PWM. The function
//creer_matrice is just here to test. If one value is near 1 and every other values are near 0, the 
//programm doesn't see a difference between ansolute and relative.

using namespace std;

vector<vector<double> > creer_matrice();
bool absolute(vector<vector<double> > matrice);
void PWM_to_PSSM(vector<vector<double> >& matrice);
void PWM_to_PSSM_2(vector<vector<double> >& matrice);
bool which_PWM_to_PSSM(vector<vector<double> > matrice);
bool PWM(vector<vector<double> > matrice);
bool possible(vector<vector<double> > matrice);
vector <bool> matrix_status(vector<vector<double> > matrice);
vector<vector<double> > creer_matrice();


int main()
{
	vector<vector<double> > M;
	
	M = creer_matrice();
	
	
		
	return 0;
}

//======================================================================

//1 if the matrix is absolute, 0 otherwise. This function checks every
//line of the matrix. A matrix with relative lines and absolute lines
//would be categorised as relative. The input has to be a PSSM. The reason 
//why I check every line is because the line ( 1.0 0.0 0.0 0.0 ) where
//1.0 can be at any place, is the same when it's converted to relative.

bool absolute(vector<vector<double> > matrice)
{
	assert( possible(matrice) == 1 );
	
	size_t b(0);
	for(size_t i(0); i<matrice.size(); ++i)
	{
		double a(0);
		for( size_t j(0); j<matrice[0].size(); ++j)
		{
			a += matrice[i][j];
		}
		if (a > 0.999 and a < 1.001)
		{//Since we always compute a limited amount of digits, the sum of
		 //the values of each column may be slightly different from 1
		 //that's why I didn't put a == 1.0
			++b;
		}
	}
	if ( b == matrice.size() )
	{
		return 1;
	}
	
	return 0;
}

//======================================================================

//1 if it's a PWM. The function doesn't launch if the matrix is impossible.
//There's 3 tests, a PSSM can't have a sum of its columns higher than 4
//because the maximum is (1.0 1.0 1.0 1.0) which would be (0.25 0.25 0.25 0.25)
//converted to a relative matrix. If it's an absolute matrix a.k.a
//a matrix that have no 1.0 values, then the sum of its values must be 1.
//The last test is obvious, it can't have any negative values

bool PWM(vector<vector<double> > matrice)
{
	assert ( possible(matrice) == 1 );
	
	for(size_t i(0); i<matrice.size(); ++i)
	{
		double a(matrice[i][0]);
		double b(matrice[i][1]);
		double c(matrice[i][2]);
		double d(matrice[i][3]);
		
		double sum(a + b + c + d);
		if ( sum < 0.999 or sum > 4.0 )
		{
			return 1;
		}
		if (a != 1 and b != 1 and c != 1 and d != 1 and ( sum < 0.999 or sum > 1.001))
		{
			return 1;
		}
		if (a < 0 or b < 0 or c < 0 or d < 0)
		{
			return 1;
		}
	}
	
	return 0;
}			

//======================================================================

//Just here to test the importants functions

vector<vector<double> > creer_matrice()
{
	size_t nb_lignes (0);
	
	do
	{
	  cout<<"Nombre de lignes : ";
	  cin>>nb_lignes;
    }
    while (nb_lignes<1);
    
	size_t nb_colonnes (0);
	
	do
	{
	cout<<"Nombre de colonnes : ";
	cin>>nb_colonnes;
    }
    while (nb_colonnes<1);
    
	vector <vector<double> > M (nb_lignes, vector<double> (nb_colonnes));
	for ( size_t i(0); i<nb_lignes ; ++i)
	{
		for ( size_t j(0); j<nb_colonnes; ++j)
		{
			cout<<"M["<<i+1<<","<<j+1<<"]=";
			cin>>M[i][j];
		}
	}
	return M;
}	

//======================================================================

//Used in the possible function. Sometimes, they don't multiply by 0.25 so
//I added a function that does it without to do a third test in the possible function.

void PWM_to_PSSM(vector<vector<double> >& matrice) 
{
	for (size_t i(0); i < matrice.size() ; ++i) 
	{
	    for (size_t j(0); j < matrice[i].size(); ++j) 
	    { 
			matrice[i][j] = (pow(2, matrice[i][j]))*0.25;
		}  
	}
}

void PWM_to_PSSM_2(vector<vector<double> >& matrice) 
{
	for (size_t i(0); i < matrice.size() ; ++i) 
	{
	    for (size_t j(0); j < matrice[i].size(); ++j) 
	    { 
			matrice[i][j] = (pow(2, matrice[i][j]));
		}  
	}
}

//======================================================================

//1 if it's a possible matrix. The idea is simple; we test if the matrix is a PSSM,
//if it's a PWM it may not pass the test. After that we convert the matrix to a PSSM
//and we do the test again. If it doesn't pass the test twice, the matrix is impossible.
//If it's a PSSM it would 100% succeed the first test and probably not the second. If it's a
//PWM it would probably not pass the first test but it would 100% succeed the second.

bool possible(vector<vector<double> > matrice)
{
	bool a1(1);
	bool a2(1);
	bool a3(1);
	
	for(size_t i(0); i<matrice.size(); ++i)
	{
		double a(matrice[i][0]);
		double b(matrice[i][1]);
		double c(matrice[i][2]);
		double d(matrice[i][3]);
		
		double sum(a + b + c + d);
		if ( sum < 0.999 or sum > 4.0 )
		{
			a1 = 0;
		}
		if (a != 1 and b != 1 and c != 1 and d != 1 and ( sum < 0.999 or sum > 1.001))
		{
			a1 = 0;
		}
		if (a < 0 or b < 0 or c < 0 or d < 0)
		{
			a1 = 0;
		}
		if (a > 1 or b > 1 or c > 1 or d > 1)
		{
			a1 = 0;
		}
	} 
	
	vector<vector<double> > matrice_2;
	matrice_2 = matrice;
	PWM_to_PSSM(matrice_2);
	
	for(size_t j(0); j<matrice_2.size(); ++j)
	{
		double a(matrice_2[j][0]);
		double b(matrice_2[j][1]);
		double c(matrice_2[j][2]);
		double d(matrice_2[j][3]);
		
		double sum(a + b + c + d);

		if ( sum < 0.99 or sum > 4.01 )
		{
			a2 = 0;
		}
		if ((a != 1.0) and (b != 1.0) and (c != 1.0) and (d != 1.0) and (sum < 0.999 or sum > 1.001))
		{
			a2 = 0;
		}
	}
	
	PWM_to_PSSM_2(matrice);
	
	for(size_t k(0); k<matrice.size(); ++k)
	{
		double a(matrice[k][0]);
		double b(matrice[k][1]);
		double c(matrice[k][2]);
		double d(matrice[k][3]);
		
		double sum(a + b + c + d);

		if ( sum < 0.99 or sum > 4.01 )
		{
			a3 = 0;
		}
		if ((a != 1.0) and (b != 1.0) and (c != 1.0) and (d != 1.0) and (sum < 0.999 or sum > 1.001))
		{
			a3 = 0;
		}
	}
	
	return (a1 + a2 + a3);
}

//======================================================================

vector <bool> matrix_status(vector<vector<double> > matrice)
{
	vector <bool> a(2);
	
	assert( possible(matrice) == 1 );
	    if ( PWM(matrice) == true )
	    {
			a[0] = 1;
			PWM_to_PSSM(matrice);
	    }
	    else
	    {
		    a[0] = 0;
	    }
	    if ( absolute(matrice) == true )
	    {
			a[1] = 0;
		}
		else
		{
			a[1] = 1;
		}
    
    
    return a;
}
	
//======================================================================

//A function to determine if the transformation from PSSM to PWM was done 
//with or without the use of the 0.25. 1 is with, 0 is without.

bool which_PWM_to_PSSM(vector<vector<double> > matrice)
{
	assert ( possible(matrice) == 1 );

	PWM_to_PSSM(matrice);
	
	if (possible(matrice) == 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}





