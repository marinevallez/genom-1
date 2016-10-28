// Conversion from PWM to PSSM qnd from PSSM to PWM
// Marine

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


void matrix_conversion(vector<vector<double> > Matrix, bool PMW);

/* int main()
{
	vector<vector<double> > Matrix(3, vector<double>(3,5));     // test
	bool PMW(true);   
                  
    matrix_conversion(Matrix, PMW);  
	
	return 0;
	
} */


void matrix_conversion(vector<vector<double> > Matrix, bool PMW) {
 
 // boolean has to be true if you are imputing a PMW

	if(PMW) {
		for (unsigned int i(0); i < Matrix.size() ; ++i) {
			for (unsigned int j(0); j < Matrix[i].size(); ++j) { 
				Matrix[i][j] = log2(Matrix[i][j]/0.25);
				}  
			}
			cout<< " Here is the PSSM : " <<endl;
		}
	

	else {
		for (unsigned int i(0); i < Matrix.size() ; ++i) {
			for (unsigned int j(0); j < Matrix[i].size(); ++j) {
				Matrix[i][j] = exp2(Matrix[i][j])*0.25;  
				}   // we choose 0.25 as a backgroud because each aa has the same probability to appear randomly
			}
			cout<< " Here is the PMW : " <<endl;
		}

 for (unsigned int i(0); i < Matrix.size() ; ++i) {
	    cout<< endl;
			for (unsigned int j(0); j < Matrix[i].size(); ++j) {
				cout<< Matrix[i][j] << "  ";
				}
			}	

}


