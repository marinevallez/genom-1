
#include "Letters.hpp"
#include <cmath>




void test(double& L) // assure that the size of a lettre isn't equal to 2 if the probability is 0 (for aesthtic consideration)
{
    if (L == 0) {
        L = pow(10, -10);
    }
}

vector<vector<double>>  Letters::get_lettre_size() // getter for the size of each lettre
{
    return Letters_size;
}

Letters::Letters(matrix m)
{
    
 
    for (size_t j(0) ; j < m.size() ; ++j)
    {
        double max (0);
        for (int i(0); i < 4; ++i) {
            if (m[j][i] > max) {
                max = m[j][i];
            }
        };
        
        double proportion (max * 300);
     
        Letters_size.push_back({m[j][0] * proportion ,m[j][1] * proportion, m[j][2] * proportion ,m[j][3] * proportion});
    
        
    }
    
}





