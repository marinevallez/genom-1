#ifndef MATRIX_HPP
#define MATRIX_HPP
//Valentin - in progress
#include <iostream>
#include <vector>

typedef std::vector<std::vector<double> > matrix;

class Matrix {
	
	public : 
	
	Matrix();
	~Matrix();
	void setrw(int value);
	void swaptopssm(); //previous idea
	void swaptopwm();
	matrix loadmatrix(std::string Data);
	
	private :
	
	int rw;
	matrix mx;
	bool status; //true if PWM false if PSSM  ?
	/*need to think about what Matrix's attributes could be, so i kept previous ideas*/
};
#endif


//Question to myself 
// How can we cange class of matrix from PWM to PSSM and (<-), + constructors, shall i define operator= for matrix;
