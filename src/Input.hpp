#ifndef INPUT_H
#define INPUT_H
#include <iostream>
using namespace std;

class Input		//a purely virtual class because we will never have JUST input, it will always be a specific sub-class
{
public:

    Input() {}	//constructor and destructor
    virtual ~Input() {}
	virtual void loadFile(string) = 0;	//a purely virtual method for loading files that requires a file name 
										//need to be redefined in every sub-class

};

#endif
