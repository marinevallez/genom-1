
#include "Protein.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>


using namespace std;

double get_afinity_score_from_matrix(PWM matrix,vector<char> sequence) // calculates the binding score (double) of a certain sequence based on a Matrix of type PWM
// the function throws an exception which should be catched !
{
    if (matrix.size() < sequence.size()) // checks that the sequence is entirly contained in the matrix
    {
        throw std::runtime_error("Error : The Matrix doesn't fit the sequence"); //throw exeption
    }
    if (matrix.size()==0) { // checks if the Matrix is configurated
        throw std::runtime_error("Error : The Matrix isn't configurated "); //throw exeption
    }
    
    double score(1);
    for (unsigned int i(0); i <= sequence.size(); ++i) {
        switch (sequence[i]) {
            case 'A':
                score *= matrix[i][0];
                break;
            case 'C':
                score *= matrix[i][1];
                break;
            case 'G':
                score *= matrix[i][2];
                break;
            case 'T':
                score *= matrix[i][3];
                break;
                
                
            default:
                break;
        }
    }
    return score;
}


double get_afinity_score_from_matrix(PSSM matrix,vector<char> sequence) // surcharge of the function above :  calculates the binding score (double) of a certain sequence based on a Matrix of type PSSM (log base 2)
// the function throws an exception which should be catched !
{
    if (matrix.size() < sequence.size()) // checks that the sequence is entirly contained in the matrix
    {
        throw std::runtime_error("Error : The Matrix doesn't fit the sequence"); //throw exeption
    }
    if (matrix.size()==0) { // checks if the Matrix is configurated
        throw std::runtime_error("Error : The Matrix isn't configurated "); //throw exeption
    }
    
    double score(0);
    for (unsigned int i(0); i <= sequence.size(); ++i) { // assuming that the base are in colomn and in the ACGT order
        switch (sequence[i]) {
            case 'A':
                score += matrix[i][0]; //addition due to the log form of the probability
                break;
            case 'C':
                score += matrix[i][1];
                break;
            case 'G':
                score += matrix[i][2];
                break;
            case 'T':
                score += matrix[i][3];
                break;
                
                
            default:
                break;
        }
    }
    return score;
}
