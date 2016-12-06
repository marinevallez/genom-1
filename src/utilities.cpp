#include "utilities.hpp"

//Conversions methods

double to_Double(const string& str) // allows a convertion from string to double
{
    istringstream stream(str);
    double dbl;
    if (!(stream >> dbl))
        return 0;
    return dbl;
}

int toInt(const string& str)  // allows a convertion from string to int
{
    istringstream stream(str);
    int a;
    if (!(stream >> a))
        return 0;
    return a;
}

vector<char> toVector(const string& str)  
{
	vector<char> vec;
	for (size_t i(0); i < str.size(); ++i)
	{
		vec.push_back(str[i]);
	}
	return vec;
}

string toString(const vector<char>& vec) 
{
	string str;
	for (size_t i(0); i < vec.size(); ++i)
	{
		str += vec[i];
	}
	return str;
} 

