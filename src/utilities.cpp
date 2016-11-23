#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "matrixprotein.hpp"
#include "sequence.hpp"
#include <vector>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cmath>
#include <map>
#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <iomanip>
#include <iostream>
using namespace std;


//Conversions methods

double to_Double(const string& str) // allows a convertion from string to double
{
    istringstream stream(str);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}

int toInt(const string& str) // allows a convertion from string to int
{
    istringstream stream(str);
    int a;
    if (!(stream >> a))
        return 0;
    return a;
}

vector<char> toVector(string str) 
{
	vector<char> vec;
	for (size_t i(0); i <= str.size(); ++i)
	{
		vec.push_back(str[i]);
	}
	return vec;
}

string toString(vector<char> vec)
{
	string str;
	for (size_t i(0); i <= vec.size(); ++i)
	{
		str += vec[i];
	}
	return str;
}

#endif
