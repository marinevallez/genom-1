

#ifndef Logo_hpp
#define Logo_hpp

#include <stdio.h>
#include <array>
#include <SFML/Graphics.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Letters.hpp"
#include "Axes.hpp"
#include "matrixprotein.hpp"

using namespace std;

class Logo {

/*!
 * contain the texture of each letter filled with the emplacment of the image by the constructory.
 * */
 
    array<sf::Texture, 4> lettre_texture; 
   
    
public:

/*!
 * Fills lettre_texture.
 * */
    Logo();
   
/*!
 * Display the logo.
 * */
    int afficher_logo(vector<Pattern> sequences, matrix); 
    
   
/*!
 * Calculates the size that must be displayed, based on its probability. This function calculates a arbitrary size and might result in a loss of information !
 * */
    double size( double size_b);


/*!
 * For a given letter, gives the appropriate size, texture and emplacement.
 * */
    sf::RectangleShape set_text( double taille, size_t lettre, double emplacement, int largeur); 
    
    

/*!
 * Determine the order of probability for a given place in the sequence.
 * */    
    vector<size_t> Letters_order(Letters Letters_, size_t j); 
    
};

#endif /* Logo_hpp */
 
