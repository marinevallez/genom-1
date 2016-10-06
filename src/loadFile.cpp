#include <iostream>
#include <fstream> 						//provides file stream classes
#include <string>			
using namespace std;


int main()
{
	string fileName; 
	cout << "Choose a file to work on :" << endl;
	cin >> fileName;
	
	ifstream currentFile(fileName);	 	// class instance from the file chosen by the user
	if(file.fail())						// if it didnt open it will show an error 
	{
		cerr << "File could not be opened" << endl;
	}
	else  								//content will include functionalities already implemented by others
	{
		
		
	}
	
	
	currentFile.close();				// make sure we close the strean
	return 0;
}
