/*    This file is part of simplex
      Copyright (C) 2019  Julien Thevenon ( julien_thevenon at yahoo.fr )

      This program is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#ifndef _MY_EQUATION_SYSTEM_H_
#define _MY_EQUATION_SYSTEM_H_

#include "my_square_matrix.h"

class SystemEquation
{
  private:
	my_square_matrix matrice;
	my_matrix coef;

  public:


    SystemEquation(my_square_matrix p_matrix, my_matrix p_coef);

    double * solve();

    std::string to_string();
};

//-----------------------------------------------------------------------------
SystemEquation::SystemEquation(my_square_matrix p_matrix, my_matrix p_coef)
:matrice(p_matrix)
,coef(p_coef)
{
    if(p_matrix.get_height() != p_coef.get_height())
    {
        throw quicky_exception::quicky_logic_exception("Equation p_matrix and coefficient p_matrix have incompatible sizes", __LINE__, __FILE__);
    }
}

//-----------------------------------------------------------------------------
double * SystemEquation::solve()
{
    if(matrice.getDeterm()==0)
    {
        return(NULL);
    }

    int largeur = matrice.get_width();
    int i,i2,j;
    bool pivot;
    double * max = new double[2];
    double * result= new double[largeur];

    // Matrix triangularisation
    for (i=0;i<largeur ;i++ )
    {
        pivot=false;

        if(matrice.get_data(i,i) == 0)
        {

            pivot=true;
            max = matrice.max_abs_sub_column(i + 1,i);
            matrice.swap_line(i, (int)max[1]);
            coef.swap_line(i, (int)max[1]);
        }

        if(pivot == false || (pivot == true && matrice.get_data(i,i) != 0))
        {
            for(i2 = i+1; i2 < largeur; ++i2)
            {
                double m = matrice.get_data(i2, i) / matrice.get_data(i,i);

                coef.set_data(i2, 0, coef.get_data(i2, 0) - m * coef.get_data(i,0));
                for(j = 0; j < largeur; ++j)
                {
                    if(j == i)
                    {
                        matrice.set_data(i2, j, 0);
                    }
                    else
                    {
                        matrice.set_data(i2, j, matrice.get_data(i2, j) - m * matrice.get_data(i,j));
                    }
                }
            }
        }
    }

    // Triangular system resolution
    for(j = largeur - 1; j >= 0; --j)
    {
        result[j] = coef.get_data(j, 0) / matrice.get_data(j,j);
        for(i = j + 1; i < largeur; ++i)
        {
            result[j] = result[j] - (result[i] * matrice.get_data(j,i)) / matrice.get_data(j,j);
        }
    }
    return(result);
}

//-----------------------------------------------------------------------------
std::string SystemEquation::to_string()
{
    std::string l_string("dimension");
    l_string += std::to_string(matrice.get_width()) +"\n";
    for(unsigned int i = 0; i  < matrice.get_height(); ++i)
    {
        for(unsigned int j = 0; j <matrice.get_width(); ++j)
        {
            l_string += std::to_string(matrice.get_data(i,j)) +"\t";
        }
        l_string += "=" + std::to_string(coef.get_data(i,0)) +"\n";
    }
    return(l_string);
}

#ifdef SIMPLEX_SELF_TEST
void test_equation_system()
{
    my_square_matrix l_matrix(3);

    l_matrix.set_data(0,0,10);
    l_matrix.set_data(0,1,5);
    l_matrix.set_data(0,2,0);
    l_matrix.set_data(1,0,100);
    l_matrix.set_data(1,1,10);
    l_matrix.set_data(1,2,0);
    l_matrix.set_data(2,0,-20);
    l_matrix.set_data(2,1,100);
    l_matrix.set_data(2,2,0);

    std::cout << "Equation Matrix :" << std::endl;
    std::cout << l_matrix.to_string();


    my_matrix l_coef(3,1);

    l_coef.set_data(0,0,10);
    l_coef.set_data(1,0,10);
    l_coef.set_data(2,0,10);

    std::cout << "Coef matrix:" << std::endl;
    std::cout << l_coef.to_string();

    SystemEquation system(l_matrix,l_coef);

    std::cout << "Equation system:" << std::endl;
    std::cout << system.to_string();

    double * result = system.solve();

    if(result != NULL)
    {
        for(unsigned int i = 0; i < l_coef.get_height();++i)
        {
            std::cout << "Result[" << std::to_string(i) << "] : " << result[i] << std::endl;
        }
    }
}
#endif // SIMPLEX_SELF_TEST

#endif // _MY_EQUATION_SYSTEM_H_
// EOF
