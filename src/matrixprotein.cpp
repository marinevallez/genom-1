#include "matrixprotein.hpp"
#include <vector>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cmath>
#include <cassert>

using namespace std;

// ==============================================================================================================CONSTRUCTEUR ET DESTRUCTEUR

MatrixProtein::MatrixProtein() 
{
	/*Il nous faut une méthode qui choisi loadmatrix_fromscore ou loadmatrix_fromfile
	 * pour génerer les différentes matrices*/
}

MatrixProtein::~MatrixProtein()
{
	
}

// ==============================================================================================================METHODES

void MatrixProtein::fillPattern(double const& bScore, vector<char> const& site)
{
	pattern.bScore = bScore;
	pattern.listOfSites = site;
}


void MatrixProtein::fillVectorPatterns(Pattern pattern)					//à revoir en fonction de si on pré remplit le tableau ou pas
{
    /*if (matrix.empty()) {
		matrix.push_back(pattern);
	} else {
		for(size_t i(0); i < matrix.size() ; ++i) {
			if (matrix[i].empty()) { //comment vérifier empty() pour une seule case?
				matrix[i] = pattern;
			} else {
				++i;
			}
		}
    }	*/
}


double MatrixProtein::probas(double n, double tot)
{
    return n/tot;
}



void MatrixProtein::loadmatrix_fromscore()
{
    
    double compteurA(0);
    double compteurT(0);
    double compteurG(0);
    double compteurC(0);
    double nombreTot(0);
    std::vector <double> tabA;
    std::vector <double> tabT;
    std::vector <double> tabG;
    std::vector <double> tabC;
    matrix finale;
    
    /* Motif m1 = {{'A', 'T', 'G', 'C', 'T', 'G', 'T'}, 0.8};
     sequences.push_back( m1) ;
     Motif m2 = {{'T', 'C', 'A', 'C', 'T','C','C'}, 0.2};
     sequences.push_back( m2) ;
     Motif m3 = {{'T', 'T', 'G', 'C' , 'C','A', 'C'}, 0.4};
     sequences.push_back( m3) ;
     Motif m4 = {{ 'A' , 'A', 'C', 'G', 'G','T', 'C'}, 0.7};
     sequences.push_back( m4) ;
     */
    
    
    for (size_t j(0) ; j < 7 ; ++j)
    {
        
        for ( size_t i(0) ; i < patterns.size() ; ++i)
        {
            
            nombreTot = nombreTot + patterns[i].bScore ;
            if ( patterns[i].listOfSites[j] == 'A')
            {
                compteurA = compteurA + patterns[i].bScore ;
            }
            if ( patterns[i].listOfSites[j] == 'T')
            {
                compteurT = compteurT + patterns[i].bScore ;
            }
            if ( patterns[i].listOfSites[j] == 'G')
            {
                compteurG = compteurG + patterns[i].bScore ;
            }
            if ( patterns[i].listOfSites[j] == 'C')
            {
                compteurC = compteurC + patterns[i].bScore ;
            }
        }
        
        tabA.push_back(probas(compteurA , nombreTot));
        tabT.push_back(probas(compteurT , nombreTot));
        tabG.push_back(probas(compteurG , nombreTot));
        tabC.push_back(probas(compteurC , nombreTot));
        compteurA = 0;
        compteurT = 0;
        compteurG = 0;
        compteurC = 0;
        nombreTot = 0;
    }
    
    finale.push_back(tabA);
    finale.push_back(tabT);
    finale.push_back(tabG);
    finale.push_back(tabC);
    
    mx = finale; 
    matrix_generation(); 
    
}


/*void MatrixProtein::display_PWM(matrix finale)
{
    
    for ( int i(0); i < finale.size() ; ++i)
    {
        
        for (int j(0) ; j < finale[i].size() ; ++j)
        {
            
            cout  << setprecision(5) << setw(5) << finale[i][j] << " | " ; 
        } 
        
        std::cout << std::endl; 
    }
    
}*/

void MatrixProtein::loadmatrix_fromfile(string Data){ // the function stores data from a file containing a PWM or a PSSM assuming that the bases are stored in column and in the ACGT order

    int colmn;
    int row;
    string fichier (Data);
    vector<double> temp; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open(fichier);
  
    
    if (file.fail())
    {
        throw std::runtime_error("Erreur de lecture du fichier de donnée ");
    }
    
    string var;
    while (!file.eof()) 
    {
		while (!file.eof()) 
		{
            file >> var >> ws;
            double dbl = To_double(var); // allows to get the data in double
            temp.push_back(dbl);
        }
    }
    file.close();

    int z(0);
    colmn = 4;
    row = (temp.size())/4;  // here we make the assumption that the base are stored in column
    vector<double> tmp (colmn,0.0 );
    matrix mtrix (row,tmp); // all the case are initialize at 0.0
    for (int i(0); i < row; ++i) {
        for (int j(0); j < colmn; ++j) {
            mtrix[i][j]=temp[z];
            ++z;
        }
    }
    mx = mtrix;
    matrix_generation();
}

vector<Pattern> MatrixProtein::getPatterns()
{
	return patterns; 
}


void MatrixProtein::matrix_generation()
{
	std::vector<bool> check(2);
	assert (possible(mx)==1);
	check=matrix_status(mx);
	if( (check[0] == 0) & (check[1] == 0))
	{
		pssm_abs=mx;
		swaptorelative(mx);
		pssm_rel = mx;
		swaptopwm(mx);
		pwm_rel = mx;
		mx=pssm_abs;
		swaptopwm(mx);
		pwm_abs=mx;
		
	}
	else if ( (check[0] == 0) & (check[1] == 1))
	{
		pssm_rel=mx;
		swaptoabsolute(mx);
		pssm_abs = mx;
		swaptopwm(mx);
		pwm_abs = mx;
		mx = pssm_rel;
		swaptopwm(mx);
		pwm_rel=mx;
		
	}
	else if ( (check[0] == 1) & (check[1] == 0))
	{
		pwm_abs=mx;
		swaptopssm(mx);
		pssm_abs = mx;
		swaptorelative(mx);
		pssm_rel = mx;
		swaptopwm(mx);
		pwm_rel=mx;		
	}
	else 
	{
		pwm_rel=mx;
		swaptopssm(mx);
		pssm_rel = mx;
		swaptoabsolute(mx);
		pssm_abs=mx;
		swaptopwm(mx);
		pwm_abs=mx;
		}	
};

//Since the computer doesn't like to compute things with -inf, here's a
//function that replaces these values by -100, 2^-100 will be approximated
//by 0.

void MatrixProtein::readjust_values(matrix& mtx)
{
	for (size_t i(0); i < mtx.size(); ++i)
	{
		for (size_t j(0); j < mtx[0].size(); ++j)
		{
			if ( mtx[i][j] < -100 )
			{
				mtx[i][j] = -100;
			}
		}
	}
}

//A function to determine if the transformation from PSSM to PWM was done 
//with or without the use of the 0.25. 1 is with, 0 is without.

bool MatrixProtein::which_PWM_to_PSSM(matrix matrice)
{
	assert (possible(matrice));

	PWM_to_PSSM(matrice);
	
	if (possible(matrice))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void MatrixProtein::swaptorelative(matrix& mtx)
{
	
	assert (absolute(mtx));
	assert (possible(mtx));
	for(size_t i(0); i < mtx.size(); ++i)
	{
		double max(0.0);
		
		for(size_t j(0); j<4; ++j) 
		{                                    
			if ( mtx[i][j] > max )
			{
				max = mtx[i][j];
			}
		}
		
		for(size_t k(0); k<4; ++k)
		{
			mtx[i][k] = mtx[i][k]/max;
		}
	}
}

void MatrixProtein::swaptoabsolute(matrix& mtx)
{
	
	assert (absolute(mtx)==0);
	assert (possible(mtx));
	for(size_t i(0); i < mtx.size(); ++i)
	{	
		double total(0.0);
		
		for(size_t j(0); j<4; ++j)
		{
			total += mtx[i][j];
		}
		double x(1/total);
		for(size_t k(0); k<4; ++k)
		{
			mtx[i][k] *= x;
		}
	}
		
}

void MatrixProtein::swaptopssm(matrix& mtx){
	assert (check_if_pmworpssm(mtx));
	for (unsigned int i(0); i < mtx.size() ; ++i)
	{
		for (unsigned int j(0); j < mtx[i].size(); ++j)
		{ 
			if ( mtx[i][j] == 0 )
			{
				mtx[i][j] = -100;
			}
			else
			{
			    mtx[i][j] = log2(mx[i][j]/0.25);
			}
		}  
	}

	
}
void MatrixProtein::swaptopwm(matrix& mtx){
	assert (check_if_pmworpssm(mtx) == 0);
	for (unsigned int i(0); i < mx.size() ; ++i)
	{
		for (unsigned int j(0); j < mx[i].size(); ++j)
		{
			mx[i][j] = exp2(mx[i][j])*0.25;  
		}   // we choose 0.25 as a backgroud because each aa has the same probability to appear randomly
	}	
}


vector <bool> MatrixProtein::matrix_status(matrix matrice)
{
	
	vector <bool> a(2);
	
	assert(possible(matrice));
	if (check_if_pmworpssm(matrice))
		{
			a[0] = 1;
			PWM_to_PSSM(matrice);
	    }
	    else
	    {
		    a[0] = 0;
	    }
	    if (absolute(matrice))
	    {
			a[1] = 0;
		}
		else
		{
			a[1] = 1;
		}
    
    
    return a;
}

//1 if it's a possible matrix. The idea is simple; we test if the matrix is a PSSM,
//if it's a PWM it may not pass the test. After that we convert the matrix to a PSSM
//and we do the test again. If it doesn't pass the test twice, the matrix is impossible.
//If it's a PSSM it would 100% succeed the first test and probably not the second. If it's a
//PWM it would probably not pass the first test but it would 100% succeed the second.

bool MatrixProtein::possible(matrix matrice){
	
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
	
	matrix matrice_2;
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

		if ( sum < 0.999 or sum > 4 )
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

void MatrixProtein::PWM_to_PSSM_2(matrix& matrice) 
{
	for (size_t i(0); i < matrice.size() ; ++i) 
	{
	    for (size_t j(0); j < matrice[i].size(); ++j) 
	    { 
			matrice[i][j] = (pow(2, matrice[i][j]));
		}  
	}
}

//Used in the possible function. Sometimes, they don't multiply by 0.25 so
//We added a function that does it without to do a third test in the possible function.

void MatrixProtein::PWM_to_PSSM(matrix& matrice) 
{
	for (size_t i(0); i < matrice.size() ; ++i) 
	{
	    for (size_t j(0); j < matrice[i].size(); ++j) 
	    { 
			matrice[i][j] = (pow(2, matrice[i][j]))*0.25;
		}  
	}
}

//1 if the matrix is absolute, 0 otherwise. This function checks every
//line of the matrix. A matrix with relative lines and absolute lines
//would be categorised as relative. The input has to be a PSSM. The reason 
//why I check every line is because the line ( 1.0 0.0 0.0 0.0 ) where
//1.0 can be at any place, is the same when it's converted to relative.

bool MatrixProtein::absolute(matrix matrice)
{
	assert(possible(matrice));
	
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

//Could be better in another global utilitary file
double MatrixProtein::To_double( const string& string ){ // allows a convertion from string to double

    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}


/*
Matrix::Matrix(matrix matrice)
: mx(matrice) {
	readjust_values(mx);
	matrix_generation();}
*/


//1 if it's a PWM. T~he function doesn't launch if the matrix is impossible.
//There's 3 tests, a PSSM can't have a sum of its columns higher than 4
//because the maximum is (1.0 1.0 1.0 1.0) which would be (0.25 0.25 0.25 0.25)
//converted to a relative matrix. If it's an absolute matrix a.k.a
//a matrix that have no 1.0 values, then the sum of its values must be 1.
//The last test is obvious, it can't have any negative values

bool MatrixProtein::check_if_pmworpssm(matrix matrice)
{
	assert (possible(matrice));
	
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
