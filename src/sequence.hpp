#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Input.hpp"
#include "Protein.hpp"
using namespace std;

struct PosDir 	//stands for position and direction; position of the motif on the sequence + direction of the sequence
{
	size_t pos;
	size_t seqNb;
	string chrNb;
	char dir;
	Pattern pattern;
};

class Sequence {

private:
	
public:

    Sequence();							//constructor + destructor
    virtual ~Sequence();
	
/*! 
 * The outputSite method outputs a file containing sequence numbers (from a .fasta file) and the position on which a motif is found, 
 * as well as the direction of the reading frame (+ means forward, - means reverse).
 */
	void outputSites(const vector<PosDir>&) const;
	
/*!
 *The motifRecognition method opens a .fasta file and find a given motif on a sequence (and its complment) in the .fasta file.
 */

	vector<PosDir> motifRecognition(const string&, const string& fileName) const; 
	
/*!
 * The motifRecognition method can also open a .fasta file and, with a list of possible motifs, find which one are present within the sequences of the .fasta file.
 */
	
	vector<PosDir> motifRecognition(Protein protein, const string& fileName) const;
	
/*!
 * The giveComplementarySeq method gives the reverse complementary sequence of any nucleotidic sequences, whether a motif or a genomic sequence.
 */
	vector<char> giveReverseComplementarySeq(const vector<char>&) const; //method giving the REVERSE complementary sequence from one of the two DNA strand
	
};

#endif
