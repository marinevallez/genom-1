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

double to_Double(const string&); // allows a convertion from string to double

int toInt(const string&);	// allows a convertion from string to int

vector<char> toVector(const string& str);

string toString(const vector<char>&);

double calculateScore(matrix, vector<char>);

#endif
