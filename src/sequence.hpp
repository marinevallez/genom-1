#include <vector>
#include <string>
using namespace std;

class Sequence {
	
private:
	
	vector<char> seq1;	//the two sequences we are working with
	vector<char> seq2;
	
public:

	void loadSeq();					// read two sequences from a file
	size_t searchSeq(string) const;	//find the position of subsequence within the bigger sequence
	vector<int> searchSeq_(string) const; 	//anotther version of searchSeq method
	void display();					//display the two sequences
	
};
