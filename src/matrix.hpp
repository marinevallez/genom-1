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
	bool absolute(vector<vector<double> > matrice);
	void PWM_to_PSSM(vector<vector<double> >& matrice);
	void PWM_to_PSSM_2(vector<vector<double> >& matrice);
	bool which_PWM_to_PSSM(vector<vector<double> > matrice);
	bool PWM(vector<vector<double> > matrice);
	bool possible(vector<vector<double> > matrice);
	vector <bool> matrix_status(vector<vector<double> > matrice);
	
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
