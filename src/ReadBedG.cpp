// #include "ReadBedG.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <sstream>
#include <string>


using namespace std;

double To_double( const string& string ) // allows a convertion from string to double
{
    istringstream stream(string);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}

double To_int( const string& string ) // allows a convertion from string to int
{
    istringstream stream(string);
    int a;
    if (!(stream >> a))
        return 0;
    return a;
}


struct Coordinate {
    string chromosome;
    int start;
    int end;
    double score;
};

vector<Coordinate> ReadBedGraph(string Data) // the function stores data from a file containing a chromosome n°, 2 positions and a score
{
    
    string fichier(Data);
    
    vector<Coordinate> Coordinates;
    
    vector<string> temp; // stock in a 1x1 Matrix to calculate the number of data
    ifstream file;
    file.open(fichier);
    
    
    if (file.fail()) {
        cerr << "This file can't be openned" <<endl;
    }
    
    else {
        string var;
        
        while (!file.eof()) {
            
            while (!file.eof()) {
                file >> var >> ws;
                temp.push_back(var);
            }
        }
        
        file.close();
        
        
        unsigned int z(0);
        
        while (z < temp.size()) {
            Coordinate c;
            c.chromosome = temp[z];
            ++z;
            c.start = To_int(temp[z]);
            ++z;
            c.end = To_int(temp[z]);
            ++z;
            c.score = To_double(temp[z]);
            ++z;
            Coordinates.push_back(c);
            
        }
    }
    return Coordinates;
    
}

double addition_intervalle(int pos, vector<Coordinate> Coordinates)
{
	
	double a(0);
	size_t i(0);
	
	while (Coordinates[i].end < (pos-50))
	{
		++i;
	}
	
	size_t j(i);
	
	while (Coordinates[j].end < (pos+50))
	{
		++j;
	} 
	
	for (size_t k(i); k <= j; ++k)
	{
		int start(Coordinates[k].start);
		int end(Coordinates[k].end);
		
		if (start < (pos-50))
		{
			start = (pos-50);
		}
		
		if (end > (pos+50))
		{
			end = (pos+50);
		}
		
		a += (end-start+1)*Coordinates[k].score;
	}
	
	return a;
}

	
