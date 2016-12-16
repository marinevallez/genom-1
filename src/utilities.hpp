#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <string>
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

typedef vector<vector<double>> matrix;

/*!
 * Allows a conversion from string to double.
 * */
double to_Double(const string&);

/*!
 * Allows a conversion from string to int.
 * */
int toInt(const string&);

/*!
 * Allows a conversion from string to vector.
 * */
vector<char> toVector(const string& str);

/*!
 * Allows a conversion from vector to string.
 * */
string toString(const vector<char>&);

double calculateScore(matrix, vector<char>);

#endif
