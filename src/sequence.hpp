#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "matrixprotein.hpp"
using namespace std;

//----Different struct/objects----

/*!
 * 
 * The PosDir struct - which stands for Position/Direction - groups together all the informations related 
 * to one occurence of motif in order to display them and write them on files. PosDir contains the position 
 * of the motif on the sequence, the number of the sequence on which it was found, the reading direction of the sequence- 
 * that is the strand on which the motif is found- the motif itself and its binding score.
 * 
 * */

struct PosDir 										//stands for position and direction; position of the motif on the sequence + direction of the sequence
{
    size_t pos;
    size_t seqNb;
    string chrNb;
    char dir;
    string sequence;
    double bindingscore;
};

/*!
 * 
 * The Coordinate struct groups all the informations that are give by .bedgraph file, that is the chromosome 
 * on which the motif is located, the starting and ending positions of the motif (its 'coordinates' on the chromosome) 
 * found in genomic .fasta file (i.e. the whole chromosome) and its score.
 * 
 * */

struct Coordinate  						 		   //numerical position(beginning and end) and score of an unknown motif from a .bedGraph file
{
    string chromosome;
    int start;
    int end;
    double score;
};

/*!
 * 
 * The BedCoordinate struct is very similar to Coordinate, yet lacks the score. 
 * Its content may be given from a .bed file but also a .bedgraph file.
 * 
 * */
 
struct BedCoordinate       				    //numerical position(beginning and end) of an unknown motif from a .bed file
{
    string chromosome;
    int start;
    int end;
};

/*!
 * /brief The class Sequence represents the informations and operations that are relevant to a DNA sequence. 
 * Sequence manages the reading and scanning of genomic sequences, conversions into complementary DNA sequences as well as 
 * the organisation of the results found from scanning a sequence, whether it is looking for a specified sequence or for potential sequences
 * if a Position-Weight Matrix is given. The Sequence class also produces .fasta files and makes sure that such files are well-formatted. 
 * 
 * */
class Sequence {

	private:

	/*!
	 * The private attribute motifs4output gathers the ultimate set of motifs that are going to be given to the user, along with all the necessary informations.
	 * 
	 * */
     vector<PosDir> motifs4output;		//the motifs we are going to export into a .txt file
    
	public:
    
    //---Constructor & Destructor---
    Sequence();
    virtual ~Sequence();
    
    //---Methods---
    
/*!
* The outputSite method produces a .txt. file containing sequence numbers (from a .fasta file) and the position on which 
* a motif is found, as well as the direction of the reading frame (+ means forward, - means reverse, 
* depending on which DNA strand the chromosome was found).
*/

	 void outputSites(const vector<PosDir>&) const;
	 
/*!
* The motifRecognition method opens a .fasta file and find a given (specified by the user) 
* motif on a sequence in the .fasta file.
*/
 
	vector<PosDir> motifRecognition(const string&, const string& fileName) const;
    
/*!
* The motifRecognition method can also open a .fasta file and, with a list of generated possible motifs, 
* find which one are present within the sequences of the .fasta file.
*/
    
    vector<PosDir> motifRecognition(const vector<Coordinate>& coordinates, const string& fileName) const;  
     

/*!
* The fastaPlusMatrix method opens a .fasta file (non-genomic, a short .fasta) and 
* using a matrix fills out motifs4output with motifs from the .fasta that have a sufficient score 
* (higher than a threshold provided by the user).
*/
    
    void fastaPlusMatrix(const string&,matrix&, const double&);    
    
/*!
* The giveComplementarySeq method gives the reverse complementary sequence of any nucleotidic sequences, 
* whether a motif or a genomic sequence.
*/
    vector<char> giveReverseComplementarySeq(const vector<char>&) const;
    
/*!
* The readBedGraph method reads a .bedGraph file to find the start and end positions of a motif on a chromosome, along with a score. 
*/
    
    vector<Coordinate> readBedGraph(const string& fileName, string);
    
/*!
* The readBed method read a .bed file to find the start and end positions of a motif on a chromosome (without any score).
* */
    
    vector<Coordinate> readBed(const string& fileName, string);

/*!
* The scanFasta method uses a list of coordinates from a .bedgrah file to scan each chromosomal genomic sequence 
* in a .fasta file for regions listed in the .bedgraph. 
* */
    
    vector<string> scanFasta(vector<Coordinate>& coordinates, const string& fileName);

/*!
* The scanFasta method uses a list of coordinates from a .bedgrah file to scan each chromosomal genomic sequence 
* in a .fasta file for regions listed in the .bedgraph, looking for a specific chromosome. 
* */
    
    vector<string> scanFasta(vector<Coordinate>& coordinates, const string& fileName, int);
    
/*!
 * The makeFasta method uses the list of sequences given by scanFasta and generates a .fasta file containing all the regions that were found 
 * in .bed and .bedgraph files for easier manipulations.
 * */    
    
    void makeFasta(const vector<string>&, const vector<Coordinate>&) const;
    
/*!
* The delete_vectors_too_small method deletes sequences given by scanFasta (from .bed and .bedgraph) 
* that are to small to be used for the EM algorithm, that is, sequences smaller than the motif of interest.
*/
    
    vector<vector<char>> delete_vectors_too_small(size_t n, vector<vector<char>> target);
    
/*!
* The loadResultsOnFile method loads on a given file all the informations related to a given sequence/several sequences
* loadResultsOnFile is overloaded for the readBedGraph output and one for the readBed output.
*/
    
    void loadResultsOnFile(const string& fileName, const vector<PosDir>& posdir, vector<double> sommeScores);
    void loadResultsOnFile(const string& fileName, const vector<PosDir>& posdir);
    void loadResultsOnFile(const string& fileName);
    
/*!
* The loadMatrixOnFile method loads on a given file a given matrix.
*/
    
    void loadMatrixOnFile(const string& fileName, matrix matrice);

/*
* The find method takes a matrix, loads it, goes through sequences in fasta 
* and gets all motifs with sufficient scores with their positions and direction.

    
    void find(MatrixProtein&, const string&);*/
    
/*!
*  A getter method to get motifs4output to be exported into a .txt file.
* */
    
    vector<PosDir> getMotifs4Output() const;
    
/*!
* The loadSeq method reads a short .fasta and gets the appropriate sequences.
*/
    vector<string> loadSeq( string fileName);
       
/*!
 * The Convert_BedCtoC method converts the structBedCoordinate to Coordinate.
 * */    
 
    vector<Coordinate> Convert_BedCToC(vector<BedCoordinate> BedC);  //to delete ? no more BedCoordinate
    
/*!
 * 
 * */
 
    void Clean_Motif_Output(const string& fileName);
    
/*!
* The method interval_addition adds scores on an interval that is 50 bases upstream and downstream the position.
*/
    
    double interval_addition(int pos, vector<Coordinate> Coordinates);
    
/*!
 * 
 * */    
    
    void fillPosDir(MatrixProtein& MX, string chrn);

};

/*!
 * chromosomeNb is a local function used by the Sequence class to extract the chromosome number from headers in .fasta files
 * */
string chromosomeNb(const string&);

/*!
 * compare is a local function used by the Sequence class to compare vectors of characters; checks whether entries of both vectors are the same.
 * */
bool compare(const vector<char>&, const vector<char>&);

/*!
 * giveComplementaryBase is a local function used by the Sequence class for internal manipulations of DNA sequences, when necessary.
 * */
char giveComplementaryBase(const char&);

#endif
