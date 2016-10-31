//Valentin - in progress (
#include <iostream>
#include "matrix.hpp"
#include <vector>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cmath>

using namespace std;

//Could be better in another global utilitary file
double To_double( const string& string ){ // allows a convertion from string to double

    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}


Matrix::~Matrix(){};

Matrix::Matrix(){};


void Matrix::setrw(int value){
	
	rw = value;
	
}

//i kept our last idea of class methods
void Matrix::swaptopssm(){
	
	if (status){
		for (unsigned int i(0); i < mx.size() ; ++i){
			for (unsigned int j(0); j < mx[i].size(); ++j){ 
				mx[i][j] = log2(mx[i][j]/0.25);
			}  
		}
	}else{
		cout<<"Already PSSM";
	}
	
}
// same here
void Matrix::swaptopwm(){
	
	if (status){
		for (unsigned int i(0); i < mx.size() ; ++i){
			for (unsigned int j(0); j < mx[i].size(); ++j){
				mx[i][j] = exp2(mx[i][j])*0.25;  
			}   // we choose 0.25 as a backgroud because each aa has the same probability to appear randomly
		}	
}

void Matrix::swaptoabsolute(matrix& mtx){
	
	for(size_t i(0); i<rw; ++i)
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

void Matrix::swaptorelative(matrix& mtx){
	
	for(size_t i(0); i<rw; ++i)
	{
		double max(0.0);
		
		for(size_t j(0); j<4; ++j) 
		{                                    
			if ( matrix[i][j] > max )
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
