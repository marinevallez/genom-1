

//#include "matrix_reader.hpp"


#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
 #include <sstream>


using namespace std;

double To_double( const string& string ) // allows a convertion from string to double
{
    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}

vector<vector<double>> loadmatrix(string Data) // the function stores data from a file containing a PWM or a PSSM assuming that the bases are stored in column and in the ACGT order
{
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
    vector<vector<double>> Matrix (row,tmp); // all the case are initialize at 0.0
    for (int i(0); i < row; ++i) {
        for (int j(0); j < colmn; ++j) {
            Matrix[i][j]=temp[z];
            ++z;
        }
    }
    return Matrix;
}
   
   
    

/*int main ()
{
    vector<vector<double>> test (loadmatrix("/Users/oriane/Downloads/DBP_PPM\(1\).mat"));
    cout << "        A       C        G        T" << endl;
    int nbrline(2);
    cout << 1 << "   ";
    for (int i(0); i < test.size() ; ++i) {
        for (int j(0); j < test[i].size(); ++j) {
            cout << test[i][j] << " ";
        }
        cout << endl;
        cout << nbrline << "   ";
        ++nbrline;
    }
    
}*/
