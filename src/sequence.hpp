#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <vector>
#include <string>
#include "Input.hpp"
using namespace std;

class Sequence : public Input {
	
private:
	
	vector<char> seq1;							//the two sequences we are working with
	vector<char> seq2;
	
public:

    Sequence() {}								//constructor + destructor
    virtual ~Sequence() {}
    virtual void loadFile(string) override; 			//a redefinition of the method from the Input class
	vector<int> searchMotif(string, int) const; //outputs the position of where the motif is found within the sequence
												//int parameter allows you to switch between sequence, 1 to search the first one, 2 to search the second
	void display();								//display the two sequences
	
	
	//From here : methods impl. to give the complementary of a sequence
	vector<char> giveComplementarySeq(vector<char>); //method giving the REVERSE complementary sequence from one of the two DNA strand
	vector<char> getSequence1(); 		 // Getters are needed for giveComplementarySeq to access the class attributes
	vector<char> getSequence2(); 
};

#endif
