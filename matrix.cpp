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
	}else{
		cout<<"Already PWM";
	}
		
	
	
}

//****Function to load a matrix from a file, maybe better as a void modifing mx attribute****

matrix Matrix::loadmatrix(string Data){ // the function stores data from a file containing a PWM or a PSSM assuming that the bases are stored in column and in the ACGT order

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
    while (!file.eof()) {
        
        while (!file.eof()) {
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
    setrw(row);
    return mtrix;
}
/* (idea before monday) We need to add something to get the status of the matrix (pssm or pwm) in this function or in constructor */
