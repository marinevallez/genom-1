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
	void loadpssm(); 
	void loadpwm();
	void swaptoabsolute(matrix& mtx);
	void swaptorelative(matrix& mtx);
	bool absolute(matrix matrice);
	void PWM_to_PSSM(matrix& matrice);
	void PWM_to_PSSM_2(matrix& matrice);
	bool which_PWM_to_PSSM(matrix matrice);
	bool check_if_pmworpssm(matrix matrice);
	bool possible(matrix matrice);
	std::vector <bool> matrix_status(matrix matrice);
	
	private :
	
	int rw;
	matrix mx;
	matrix pssm_abs;
	matrix pwm_abs;
	matrix pssm_rel; 
	matrix pwm_rel;
};
#endif


//Question to myself 
// constructors, shall i define operator= for matrix;
