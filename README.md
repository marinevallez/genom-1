#genom-1
Project in genomics no 1

1 Introduction
============
We were given the objective to build a program with actual applications in current genomic problems.  Several resource files can be given by the user, such as a **genomic sequences**, a **Position Weight Matrix** (PWM) or **Position-Specific Scoring Matrix** (PSSM), **affinity scores** along a specific sequence...

The idea is to offer different outputs that can help the user to solve his/her problem.  The output can be a **list of binding sites** present in a chromosome, an **affinity score**/**PWM** (or a PSSM), or a **graphic representation of the consensus sequence** determined by the PMW (a logo showing the different nucleotides A, T, C and G, their sizes representing their respective probability).

2 A bit of biology
=============

3 Download the program
======================
To download the program, go to : 
https://github.com/EPFL-SV-cpp-projects/genom-1. 
Open a terminal and go to the directory where you want the program, 
then type :
`git clone https://github.com/EPFL-SV-cpp-projects/genom-1`
The genom-1 folder and all of its content will be copied into your directory, you can then compile and execute the program.

4 Compilation and execution
===========================
To compile and execute the program, a fews steps are required. Make sure you are in the genom-1 folder and then :

 - `rm –rf build` and `mkdir build` to make sure an empty, clean build
   folder is created
 - `cmake ..`
 - `make`

And then you have several options :

 - `make doc` to generate the API documentation
 - `make test` to compile the tests
 - `./Main` to execute the program
 - `./Genom1Test` to run the tests

5 The functionalities of the program
=====================================
When executing the program, a menu will appear offering a list of results. 
You can choose the task that the program will perform by hitting ’1’, ’2’ , ’3’ or ’4’. 
Then you will be asked to provide the resource files (if the program needs a file you can write its name like "example.fasta"). Make sure to also write the name of the extension, like .fasta or .mat. 

 Here are the main functionalities of the program :

 1. Read a DNA sequence and a PWM (or/and its logarithmic equivalent)and give as results the list of binding sites along the
    genome.
 2. Read a DNA sequence, a list of sites and their respective binding scores and output a PWM (or/and its logarithmic equivalent).
 3. Based on a matrix or on the binding scores and list of sites, produce the consensus sequence logo.
 4. Using a matrix and a sequence, give the list of all possibles sites that are above a threshold given by the user.




