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

struct PosDir 								//stands for position and direction; position of the motif on the sequence + direction of the sequence
{
	size_t pos;
	size_t seqNb;
	string chrNb;
	char dir;
	Pattern pattern;
};

struct Coordinate  						    //numerical position(beginning and end) and score of an unknown motif from a .bedGraph file
{
    string chromosome;
    int start;
    int end;
    double score;
};

struct BedCoordinate       				    //numerical position(beginning and end) of an unknown motif from a .bed file
{
    string chromosome;
    int start;
    int end;
};

class Sequence {

	private:
	
	public:
	
	//---Constructor & Destructor---
    Sequence();							
    virtual ~Sequence();
    
    /*//---Conversions---	
	vector<char> toVector(string str);
	string toString(vector<char> vec);
    double toDouble(const string& string);
	int toInt(const string& str);*/
	
	//---Methods---
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
 * The readBedGraph method lets us read a .bedGraph file to find the start and end positions of a motif on a chromosome, along with a score.
 * */
    
	vector<Coordinate> readBedGraph(const string& fileName);
    
    /*!
     * The ReadBed method lets us read a .bed file to find the start and end positions of a motif on a chromosome (without any score).
     * */
    
    vector<BedCoordinate> ReadBed(const string& fileName);
    

	/*!
	 * The scanFasta method uses a list of coordinates from a .bedgrah file to scan each chromosome genomic sequence in a .fasta file for regions listed in the .bedgraph.
	 * */
    
	vector<string> scanFasta(vector<Coordinate>& coordinates, const string& fileName, int nbrOfSeq);

    /*!
     * loadResultsOnFile method loads on a given file all the informations related to a given sequence/several sequences
     * */

    void loadResultsOnFile(const string& fileName, PosDir posdir, double sommeScores);

    /*!
     * loadMatrixOnFile method loads on a given file a given matrix
     * */

    void loadMatrixOnFile(const string& fileName, matrix matrice);
    
    /*!
     * find method take a matrix, loads it, goes through sequences in fasta and gets all motifs with sufficient scores with their positions and direction
     * */
    
    vector<PosDir> find(MatrixProtein&, const string&, const string&, const double&);


};

#endif
