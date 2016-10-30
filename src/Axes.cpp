
#include "Axes.hpp"

Axes::Axes()
{
    texture.loadFromFile("../res/graphPWM.png"); // loads the image of the axes
    sprite.setTexture(texture);
}
void Axes::draw(sf::RenderTarget& target)
{
    sprite.setOrigin(0, -50); //sets the axes on the middle of the window (this should be changed if the size of the window is changed)
    
    target.draw(sprite); 
}
