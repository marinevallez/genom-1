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

double calculateScore(matrix mx, vector<char> seq_)
{
	if (mx.size() < seq_.size()) // checks that the sequence is entirly contained in the matrix
    {
        throw runtime_error("Error : The Matrix doesn't fit the sequence"); //throw exeption
    }
    if (mx.size()==0) { // checks if the Matrix is configurated
        throw runtime_error("Error : The Matrix isn't configurated "); //throw exeption
    }
    
    double score(0);
    for (unsigned int i(0); i < seq_.size(); ++i) {
        switch (seq_[i]) {
            case 'A':
                score += mx[i][0];
                break;
            case 'C':
                score += mx[i][1];
                break;
            case 'G':
                score += mx[i][2];
                break;
            case 'T':
                score += mx[i][3];
                break;
            default:
                break;
        }
    }
    return score;
}

