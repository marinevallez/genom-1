#include <iostream>
#include <fstream> 						//provides file stream classes
#include <string>	
#include "sequence.hpp"	
using namespace std;


int main()
{
	//MAIN ENUMERATES WHAT ARE THE POSSIBLE TASKS THAT THE USER CAN ASK
	
	cout << "Hello, welcome to the Genom1 program ! Genom1 contains a package of functionalities that will help you analyse DNA binding sites." << endl;
	cout << endl;
	cout << "Here are the different tasks that this program performs :" << endl;
	cout << endl;
	cout << "[1] Display a list of binding sites along a genomic sequence for a protein" << endl;
	cout << "[2] Describe a consensus sequence through a motif logo" << endl;						//there could also be distribution graphs
	cout << "[3] Calculate an affinity score (along a genomic sequence or for a collection of sequences)" << endl;
	cout << "[4] Display a Position Weight Matrix as a motif description" << endl;
	cout << endl;
	
	int taskNumber;
	cout << "Please choose a task by typing the corresponding number (from 1 to 4) :" << endl;
	cin >> taskNumber;
	while ((taskNumber == 0) or (taskNumber >= 5))								//make sure that the user gives an int between. 1 and 4
	{
		cout << "Please choose a task by typing the corresponding number (from 1 to 4) :" << endl;
		cin >> taskNumber;	
	}
	
	//MAIN INDICATES TO THE USER WHAT KIND OF FILES NEEDED TO PERFORM CHOSEN TASK(S)
	
	
	switch(taskNumber)
	{
		case 1 : cout << "To execute this task, the user must provde a .fasta file and motif or list of motifs." << endl;
		break;
		
		case 2 : cout << "To execute this task, the user must provde a a genomic sequence or a list of sequences and an affinity score." << endl;
		break;
		
		case 3 : cout << "To execute this task, the user must provde a Position Weight Matrix " 
								"(Position-Probability Matrix or Position-Specific Scoring Matrix) and a genomic sequence." << endl;
		break;
		
		case 4 : cout << "To execute this task, the user must provde a .mat file or a list of sites and an affinity score." << endl;
		break;
	}
	
		
	//USER GIVE FILE
		
	return 0;
}
