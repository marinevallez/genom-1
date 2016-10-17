#include <iostream>
#include <fstream> 						//provides file stream classes
#include <string>	
#include "sequence.hpp"	
using namespace std;


int main()
{
	//HERE MAIN ENUMERATES WHAT ARE THE POSSIBLE TASKS THAT THE USER CAN ASK
	//cout << "Hello,..."
	
	//MAIN INDICATES TO THE USER WHAT KIND OF FILES NEEDED TO PERFORM CHOSEN TASK(S) (one at a time?)
		
	//USER GIVE FILE
	
	
	string fileName; 
	cout << "Choose a file to work on :" << endl;
	cin >> fileName;
	
	ifstream currentFile(fileName);	 		// class instance from the file chosen by the user
	
	if(currentFile.fail()) 					// if it didnt open it will show an error 
	{
		cerr << "File could not be opened" << endl;
	}
		else  								//content will include functionalities already implemented by others
		{
			int n;
			cout << "You selected " << fileName << " as a file." << endl; 
			// HERE :redirect to appropriate functions according to file type
			cout << "Please choose a task with the corresponding number :" << endl;	 //then make user choose from a list of actions
			cout << "[list to be completed]" << endl;
			cout << "[1] Find a particular motif within a DNA sequence" << endl;
		
			cin >> n;
			
			Sequence seq;
			seq.loadFile("promoters.fasta");
			seq.display();
			
			cout << "Please select a base motif to search and its length: " << endl;
			string motif;
			int seqLength;
			cin >> motif;
			cin >> seqLength;
			
			for (auto p : seq.searchMotif(motif, seqLength))
			{
				cout << p << "  ";
			}
		}
	
	currentFile.close();					// make sure we close the strean
	return 0;
}
