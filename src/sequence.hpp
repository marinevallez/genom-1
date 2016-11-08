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
	
public:

    Sequence();							//constructor + destructor
    virtual ~Sequence();
	
/*! 
 * Th outputSite method outputs a file containing sequence numbers (from a .fasta file) and the position on which a motif is found, 
 * as well as the direction of the reading frame (+ means forward, - means reverse).
 */
	void outputSites(const vector<PosDir>&) const;
	
/*!
 * motifRecognition opens a .fasta file and find a given motif on a sequence (and its compliment) in the .fasta
 */

	vector<PosDir> motifRecognition(const string&) const; 
	
	
/*!
 * The giveComplementarySeq method gives the reverse complementary sequence of any nucleotidic sequences, whether a motif or a genomic sequence.
 */
	vector<char> giveReverseComplementarySeq(const vector<char>&) const; //method giving the REVERSE complementary sequence from one of the two DNA strand
	
};

#endif
