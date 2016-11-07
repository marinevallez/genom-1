#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <vector>
#include <string>
#include "Input.hpp"
using namespace std;

struct PosDir 	//stands for position and direction; position of the motif on the sequence + direction of the sequence
{
	size_t pos;
	char dir;
};

class Sequence {
	
private:
	
	vector<char> seq1;							//the two sequences we are working with
	vector<char> seq2;
	
public:

    Sequence();							//constructor + destructor
    virtual ~Sequence();
    virtual void loadFile(string); 			//a redefinition of the method from the Input class
    
/*!
 * The searchMotif method looks for a particular motif within the two sequences.
 */
	vector<int> searchMotif(const vector<char>&, const string&) const; //outputs the position of where the motif is found within the sequence
												//int parameter allows you to switch between sequence, 1 to search the first one, 2 to search the second
/*!
 * The display method displays the two sequences from a .fasta file  on the terminal.
 */
	void display();								//display the two sequences
	
/*! 
 * Th outputSite method outputs a file containing sequence numbers (from a .fasta file) and the position on which a motif is found, 
 * as well as the direction of the reading frame (+ means forward, - means reverse).
 */
	void outputSites(const string&, const string&);
	
/*!
 * The quickRead method reads a given .fasta file and outputs on the terminal both sequences of the file, withough saving them somewhere for efficiency.
 */
	vector<vector<char> > quickRead(string) const;
	
	
	vector<PosDir> motifRecognition(const string&) const; 
	
	//From here : methods impl. to give the complementary of a sequence
	
/*!
 * The giveComplementarySeq method gives the reverse complementary sequence of any nucleotidic sequences, whether a motif or a genomic sequence.
 */
	vector<char> giveComplementarySeq(vector<char>); //method giving the REVERSE complementary sequence from one of the two DNA strand
	
/*!
 * This method makes accessible the first DNA sequence from a .fasta file.
 */
	vector<char> getSequence1(); 		 // Getters are needed for giveComplementarySeq to access the class attributes
	
/*!
 * This method makes accessible the second DNA sequence from a .fasta file.
 */
	vector<char> getSequence2(); 
};

#endif
