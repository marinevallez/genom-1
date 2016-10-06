//
//  matrix reader.cpp
//  readingDNA
//
//  Created by Oriane Peter on 05.10.16.
//  Copyright © 2016 Oriane Peter. All rights reserved.
//

#include "matrix reader.hpp"


#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
 #include <sstream>


using namespace std;

double To_double( const string& string ) // permet une convertion des string en double
{
    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}

vector<vector<double>> loadmatrix(string Data)
{
    int colmn;
    int row;
    string fichier (Data);
    vector<double> temp; // stockage dans une matrice 1x1 pour commencer (permet de récuperer le nombre de donnée)
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
            double dbl = To_double(var); // permet d'avoir les donnée en double
            temp.push_back(dbl);
        }
    }
    
     file.close();

    int z(0);
    colmn = 4;
    row = (temp.size())/4;  // il y a 4 base (A C G T) en partant du principe que tout les fichiers mettent les bases en colonnes
    vector<double> tmp (colmn,0.0 );
    vector<vector<double>> Matrix (row,tmp); // tout les cases sont initialiser à 0.0
    for (int i(0); i < row; ++i) {
        for (int j(0); j < colmn; ++j) {
            Matrix[i][j]=temp[z];
            ++z;
        }
    }
    return Matrix;
}
   
   
    

int main ()
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
    
}
