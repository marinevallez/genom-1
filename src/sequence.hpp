#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Input.hpp"
#include "matrixprotein.hpp"
using namespace std;

//----Different struct/objects----
struct PosDir 			//stands for position and direction; position of the motif on the sequence + direction of the sequence
{
	size_t pos;
	size_t seqNb;
	string chrNb;
	char dir;
	Pattern pattern;
};

struct Coordinate  				  //numerical position(beginning and end) of an unknown motif from a .bed file
{
    string chromosome;
    int start;
    int end;
    double score;
};

class Sequence {

	private:
	
	public:
	
	//-----Constructor & Destructor
    Sequence();							
    virtual ~Sequence();
	
/*! 
 * The outputSite method outputs a file containing sequence numbers (from a .fasta file) and the position on which a motif is found, 
 * as well as the direction of the reading frame (+ means forward, - means reverse).
 */
	void outputSites(const vector<PosDir>&) const;
	
/*!
 *The motifRecognition method opens a .fasta file and find a given motif on a sequence (and its complement) in the .fasta file.
 */

	vector<PosDir> motifRecognition(const string&, const string& fileName) const; 
	
/*!
 * The motifRecognition method can also open a .fasta file and, with a list of possible motifs, find which one are present within the sequences of the .fasta file.
 */
	
	vector<PosDir> motifRecognition(const vector<Coordinate>& coordinates, const string& fileName) const;  //not sure yet
	
/*!
 * The giveComplementarySeq method gives the reverse complementary sequence of any nucleotidic sequences, whether a motif or a genomic sequence.
 */
	vector<char> giveReverseComplementarySeq(const vector<char>&) const; //method giving the REVERSE complementary sequence from one of the two DNA strand

/*!
 * The readBedGraph method lets us read a .bed file to find the start and end positions of a motif on a chromosome, along with a score.
 * */
 
	vector<Coordinate> readBedGraph(const string& fileName);
    double To_double(const string& string); // allows a convertion from string to double
	
	/*!
	 * The findMotifs method uses a list of coordinates from a .bedgrah file to scan each chromosome sequence in a .fasta file for the corresponding motifs.
	 * */
	
	vector<string> findMotifs(vector<Coordinate>& coordinates, const string& fileName);};

//Conversions




double To_int(const string& string); // allows a convertion from string to int


#endif
