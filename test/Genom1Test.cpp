// Genom1Test
#include <gtest/gtest.h>
#include "../src/matrixprotein.hpp"
#include "../src/sequence.hpp"
#include <iostream>

using namespace std;

//Const variables

const matrix pwmrel({{0, 1, 0},{0.25, 0.25, 0.5},{1, 0, 0},{0.5, 0.25, 0.25}});
const matrix pssmrel({{-100, 2, -100},{0, 0, 1},{2, -100, -100},{1, 0, 0}});
const matrix absolute_PWM({{0.25, 0.5, 0.1, 0.4},{0.25, 0.25, 0.5, 0.0}});
const matrix relative_PWM ({{1, 0.4, 0.4, 0.2},{1, 0.4, 0.4, 0.2}});
const matrix not_possible_matrix ({{1, 2, 3, 4, 5}});
const vector <bool> bool_absolute_PWM ({{1,0}});
const vector <bool> bool_relative_PWM ({{1,1}});
MatrixProtein mat;
const vector<char> seq1({'A','C','T','T','C','G','A','T','C'});
const vector<char> seq2({'A','T','T','T','G','A','A','C','C'});
const vector<char> seq1Complement({'T','G','A','A','G','C','T','A','G'});
const vector<char> reverseSeq1({'G','A','T','C','G','A','A','G','T'});
const vector<PosDir> goodVec({{6,1,"chrT",'+',"CCCTTTG",0.0},{38,1,"chrT",'-',"CCCTTTG",0.0}});
const vector<PosDir> fastaCheckPosDir({{7,1,"chrT",'+', "AAAACCC", 0.0}}); 
const vector<PosDir> fastaPlusMatrixCheck({{32,1,"chrT",'+',"AAAAAAC", 6.00794},{36,1,"chrT",'+',"AACAAAG", 6.01382}});

//Prototypes

bool isEqual(const vector<char>& v1, const vector<char>& v2);
bool isEqualPosDir(const vector<PosDir>& v1, const vector<PosDir>& v2);
bool isEqualM(const matrix& m1, const matrix& m2);
bool isEqualBool(const vector<bool>& v1, const vector<bool>& v2);

//Tests

TEST(SequenceTest, SequenceComparison)
{
	
	ASSERT_FALSE(isEqual(seq1, seq2));
	
}

TEST(SequenceTest, MotifRecognitionInFasta)
{
	Sequence seq;
	vector<PosDir> results = seq.motifRecognition("CCCTTTG", "fasta_test.fasta");
	ASSERT_TRUE(isEqualPosDir(goodVec, results));
	
}

TEST(SequenceTest, Complementary)
{
	vector<char> compl_;
	for(size_t i(0); i < seq1.size(); ++i)
	{
		compl_.push_back(giveComplementaryBase(seq1[i]));
	}
	ASSERT_TRUE(isEqual(seq1Complement, compl_));
}

TEST(SequenceTest, ReverseComplementary)
{
	Sequence seq;
	vector<char> rev_ = seq.giveReverseComplementarySeq(seq1);
	ASSERT_TRUE(isEqual(rev_, reverseSeq1));
}

TEST(SequenceTest, FastaCheckingSpaces)
{
	Sequence seq;
	try 
	{
		vector<PosDir> results = seq.motifRecognition("AAACCC", "SeqFail3.fasta");
	}
	catch (const runtime_error& err)
	{
		EXPECT_EQ(err.what(), string("Error: no spaces allowed!"));
	}
}

TEST(SequenceTest, FastaCheckingHeader)
{
	Sequence seq;
	try 
	{
		vector<PosDir> results = seq.motifRecognition("AAACCC", "SeqFail2.fasta");
	}
	catch (const runtime_error& err)
	{
		EXPECT_EQ(err.what(),string("Error: missing header in .fasta!"));
	}
}

TEST(SequenceTest, FastaCheckingNucl)
{
	Sequence seq;
	try 
	{
		vector<PosDir> results = seq.motifRecognition("AAACCC", "fasta_check1.fasta");
	}
	catch (const runtime_error& err)
	{
		EXPECT_EQ(err.what(), string("Error: Wrong nucleotide detected in the .fasta!"));
	}
}

TEST(SequenceTest, ScanFastaExtraction)
{
	Sequence seq;
	vector<Coordinate> coord = seq.readBed("BMAL1_sites.bed", "chr7");
	vector<string> regions = seq.scanFasta(coord,"chr7.fasta");
	EXPECT_GE(coord.size(),regions.size());
}

TEST(SequenceTest, ScanFastaAberrantMotif)
{
	Sequence seq;
	vector<Coordinate> coord = seq.readBed("BMAL1_sites.bed", "chr7");
	vector<string> regions = seq.scanFasta(coord,"chr7.fasta", 60);
	EXPECT_EQ(0,regions.size());
}

TEST(SequenceTest, FastaPlusMatrix)
{
	Sequence seq;
	MatrixProtein mat;
	mat.loadmatrix_fromfile("DBP_PPM.mat");
	matrix pssm = mat.getpssm_rel();
	seq.fastaPlusMatrix("fasta_test.fasta", pssm, 6.0);
	cout << seq.getMotifs4Output().size() << endl;
	cout << fastaPlusMatrixCheck.size() << endl;
	vector<PosDir> test = seq.getMotifs4Output();
	for(size_t i(0); i < seq.getMotifs4Output().size(); ++i)
	{
		cout << seq.getMotifs4Output()[i].pos << " " << seq.getMotifs4Output()[i].sequence;
		cout << " " << seq.getMotifs4Output()[i].seqNb << " " << seq.getMotifs4Output()[i].chrNb << " " << seq.getMotifs4Output()[i].dir << " " << seq.getMotifs4Output()[i].bindingscore << endl;
	}
	ASSERT_TRUE(isEqualPosDir(test,fastaPlusMatrixCheck));
}

TEST(MatrixProteinTest, SwapToPwm)
{
	MatrixProtein mp1;
	mp1.setpssm_rel(pssmrel);
	matrix pwmtest;
	mp1.swaptopwm(mp1.getpssm_rel());
	pwmtest = mp1.getmx();
	ASSERT_TRUE(isEqualM(pwmrel, pwmtest));	
}

TEST(MatrixProteinTest, SwapToPssm)
{
	MatrixProtein mp;
	mp.setpwm_rel(pwmrel);
	matrix pssmtest;
	mp.swaptopssm(mp.getpwm_rel());
	pssmtest = mp.getmx();
	ASSERT_TRUE(isEqualM(pssmrel, pssmtest));	
}

TEST(MatrixProteinTest, CheckStatusAbs)
{
	MatrixProtein mat;
	vector <bool> test;
	test = mat.matrix_status(absolute_PWM);
	ASSERT_TRUE(isEqualBool(test, bool_absolute_PWM));
}

TEST(MatrixProteinTest, CheckStatusRel)
{
	MatrixProtein mat;
	vector <bool> test;
	test = mat.matrix_status(relative_PWM);
	ASSERT_TRUE(isEqualBool(test, bool_relative_PWM));
	
}

TEST(MatrixProteinTest, CheckPossible)
{
	MatrixProtein mat;
	bool a(0);
	bool test;
	test = mat.possible(not_possible_matrix);
	EXPECT_EQ(a,test);
	
}
	
TEST(MatrixProteinTest, RandomlGeneratesDNA) 
{
	MatrixProtein mat;
	double length(100);
	vector <char> actg(mat.seq_generator(100));
	EXPECT_EQ(actg.size(), length);
}

//Definitions

bool isEqual(const vector<char>& v1, const vector<char>& v2)
{
	if(v1.size() != v2.size())   {return false;}
	else
	{
		for(size_t i(0); i < v1.size(); ++i)
		{
			if(v1[i] != v2[i]) {return false;}
		}
		return true;
	}
}


bool isEqualPosDir(const vector<PosDir>& v1, const vector<PosDir>& v2)
{
	if(v1.size() != v2.size())   {return false;}
	else 
	{
		for(size_t i(0); i < v1.size(); ++i)
		{
			if(v1[i].pos != v2[i].pos
				or v1[i].seqNb != v2[i].seqNb
				or v1[i].chrNb != v2[i].chrNb 
				or v1[i].dir != v2[i].dir
				or v1[i].sequence != v2[i].sequence
				or v1[i].bindingscore != v2[i].bindingscore)
				{
					return false;
				} 
		}
		return true;
	}
}

bool isEqualM(const matrix& m1, const matrix& m2)
{
	if(m1.size() != m2.size())    {return false;}
	else
	{
		for (unsigned int i(0); i < m1.size() ; ++i)
		{
			if(m1[i].size() != m2[i].size())	{return false;}
			else
			{
				for (unsigned int j(0); j < m1[i].size(); ++j)
				{
					//cout << m2[i][j] << " "; 
					if (m1[i][j] != m2[i][j])
					{
						return false;
					}
				}  
			} 
		}
	}
	return true; 
}

bool isEqualBool(const vector<bool>& v1, const vector<bool>& v2)
{
	bool res;
	if (v1.size() != v2.size())
	{
		res = false;
	}
	for (size_t i(0); i < v1.size(); ++i)
	{
		
		if (v1[i] == v2[i])
		{
			res = true;
		}
	}
	return res;
}


int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
