

#include "Logo.hpp"


Logo::Logo()
{
    lettre_texture[0].loadFromFile("../res/A.png");
    lettre_texture[1].loadFromFile("../res/C.png");
    lettre_texture[2].loadFromFile("../res/G.png");
    lettre_texture[3].loadFromFile("../res/T.png");
    
}

double Logo::size( double size_b)
{
    return pow(6,size_b)*10; //! arbitrary choice of size chosen by random testing
}

sf::RectangleShape Logo::set_text( double taille, size_t lettre, double emplacement, int largeur)
{
    
    sf::Texture* ptexture;
    
    
    switch (lettre) {
        case 0 :
        {
            
            
            ptexture = &lettre_texture[0];
            break;
        }
        case 1 :
        {
            
            
            ptexture = &lettre_texture[1];
            break;
        }
        case 2 :
        {
            
            
            ptexture = &lettre_texture[2];
            break;
            
        }
        case 3 :
        {
            
            
            ptexture = &lettre_texture[3];
            break;
        }
            
            
        default:
            break;
    }
    sf::RectangleShape rectangle;
    rectangle.setSize(sf::Vector2f(100,pow(6,taille)*10));
    rectangle.setFillColor(sf::Color(250,250,250,250));
    rectangle.setPosition(largeur, emplacement);
    rectangle.setTexture(ptexture);
    
    
    return rectangle;
}

vector<size_t> Logo::Letters_order(Letters Letters_, size_t j)
{
    
    double nbr1 (0.0);
    size_t lt1 (0);
    double nbr2(0.0);
    size_t lt2 (0);
    double nbr3(0.0);
    size_t lt3 (0);
    double nbr4(0.0);
    size_t lt4 (0);
    for (size_t i(0); i < 4; ++i) {
        if (Letters_.get_lettre_size()[j][i] > nbr1) {
            nbr4 = nbr3;
            lt4 = lt3;
            nbr3 = nbr2;
            lt3 = lt2;
            nbr2 = nbr1;
            lt2 = lt1;
            nbr1 = Letters_.get_lettre_size()[j][i];
            lt1 = i;
        } else if (Letters_.get_lettre_size()[j][i] > nbr2) {
            nbr4 = nbr3;
            lt4 = lt3;
            nbr3 = nbr2;
            lt3 = lt2;
            nbr2 = Letters_.get_lettre_size()[j][i];
            lt2 = i;
        } else if (Letters_.get_lettre_size()[j][i] > nbr3)
        {
            nbr4 = nbr3;
            lt4 = lt3;
            nbr3 = Letters_.get_lettre_size()[j][i];
            lt3 = i;
        } else  {
            nbr4 = Letters_.get_lettre_size()[j][i];
            lt4 = i;
        }
    }
    return {lt1,lt2,lt3,lt4};
}


int Logo::afficher_logo(vector<Motif> sequences)
{
    sf::RenderWindow window(sf::VideoMode(1400, 500), "PWM logo");
    
    //Creation axes
    sf::Font font;
    if (!font.loadFromFile("../res/arial.ttf")) return EXIT_FAILURE;
    window.clear(sf::Color(250,250,250));
    Axes Axes_;
    Axes_.draw(window);
    
    // Creation letters
    
    Letters Letters_(sequences);
    int largeur(95);
    
    
    int position_axe_x(348);
    for (size_t j(0); j < Letters_.get_lettre_size().size(); ++j) {
        
        vector<size_t> position(Letters_order(Letters_, j)); //! gets the order of probability for the letter
        
        //! Calculate the position of the letter from the 0 of the Bits axes, taking in account the size of the letter above it, the origine is on the left corner of the window
        
        sf::RectangleShape probabilite1(set_text(Letters_.get_lettre_size()[j][position[0]], position[0], position_axe_x - size(Letters_.get_lettre_size()[j][position[0]]) - size(Letters_.get_lettre_size()[j][position[1]]) - size(Letters_.get_lettre_size()[j][position[2]]) - size(Letters_.get_lettre_size()[j][position[3]]), largeur));
        
        window.draw(probabilite1);
        
        sf::RectangleShape probabilite2(set_text(Letters_.get_lettre_size()[j][position[1]], position[1], position_axe_x - size(Letters_.get_lettre_size()[j][position[3]]) - size(Letters_.get_lettre_size()[j][position[2]]) - size(Letters_.get_lettre_size()[j][position[1]]), largeur));
        
        window.draw(probabilite2);
        
        sf::RectangleShape probabilite3(set_text(Letters_.get_lettre_size()[j][position[2]], position[2], position_axe_x - size(Letters_.get_lettre_size()[j][position[3]]) - size(Letters_.get_lettre_size()[j][position[2]]), largeur));
        
        window.draw(probabilite3);
        
        sf::RectangleShape probabilite4(set_text(Letters_.get_lettre_size()[j][position[3]], position[3],position_axe_x - size(Letters_.get_lettre_size()[j][position[3]]) , largeur));
        
        window.draw(probabilite4);
        
        largeur+=75; //! spacing of 75 between each column on the axes (arbitary) 
    }
    
    //Display
    window.display();
    
    //Gestion Closing
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }
        
    }
    return EXIT_SUCCESS;
    
}
