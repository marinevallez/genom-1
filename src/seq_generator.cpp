#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <vector>

using namespace std;

static const char nucleotides[] =		// list of nucleotides
"ATCG"
;

int stringLength = sizeof(nucleotides) - 1;

char genRandomChar()  					// Generates a random nucleotide
{

    return nucleotides[rand() % stringLength];
}


vector<char> seq_generator()   
{

	cout << "Length of the sequence ?" << endl;
	int length(0);
	cin >> length; 							// length of the sequence
	vector <char> sequence;
	
    while (cin.fail() or length < 0)		// checks the input
    {
		
		cerr << "Please select a positive number " << endl;
		cin.clear();
		cin.ignore(numeric_limits<streamsize>::max(), '\n' );
		cin >> length;
		
	}
	
				
    for(int i=0; i < length; i++)      				// creates the sequence
    {
        sequence.push_back(genRandomChar());

    }

    return sequence;
}

int main()								// test

{
	seq_generator();
	
	return 0;
}

