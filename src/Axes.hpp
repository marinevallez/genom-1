

#ifndef Axes_hpp
#define Axes_hpp

#include <SFML/Graphics.hpp>
#include <vector>
#include <stdio.h>
#include <array>

using namespace std;
class Axes
{
public:
    void draw(sf::RenderTarget& target);
    Axes();
    private :
    sf::Texture texture;
    sf::Sprite sprite;
};


#endif /* Axes_hpp */
